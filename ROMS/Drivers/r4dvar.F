#include "cppdefs.h"
      MODULE r4dvar_mod

#ifdef R4DVAR
!
!svn $Id: r4dvar.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module splits the R4D-Var data assimilation algorithm into     !
!  its logical components routines:                                    !
!                                                                      !
!    background        Xb, nonlinear trajectory used to linearize the  !
!                          tangent linear and adjoint models           !
!                                                                      !
!    increment         dXa, 4D-Var data assimilation increment: inner  !
!                       loops and minimization solver                  !
!                                                                      !
!    analysis          Xa, analyzed state: Xa = Xb + dXa               !
!                                                                      !
!    prior_error       Process background and model prior error        !
!                      covariance and its normalization coefficients   !
!                                                                      !
!    posterior_error   Posterior error analysis                        !
!                                                                      !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Moore, A.M., H.G. Arango, G. Broquet, B.S. Powell, A.T. Weaver,   !
!      and J. Zavala-Garay, 2011: The Regional Ocean Modeling System   !
!      (ROMS)  4-dimensional variational data assimilations systems,   !
!      Part I - System overview and formulation, Prog. Oceanogr., 91,  !
!      34-49, doi:10.1016/j.pocean.2011.05.004.                        !
!                                                                      !
!    Moore, A.M., H.G. Arango, G. Broquet, C. Edward, M. Veneziani,    !
!      B. Powell, D. Foley, J.D. Doyle, D. Costa, and P. Robinson,     !
!      2011: The Regional Ocean Modeling System (ROMS) 4-dimensional   !
!      variational data assimilations systems, Part II - Performance   !
!      and application to the California Current System, Prog.         !
!      Oceanogr., 91, 50-73, doi:10.1016/j.pocean.2011.05.003.         !
!                                                                      !
!    Gurol, S., A.T. Weaver, A.M. Moore, A. Piacentini, H.G. Arango,   !
!      S. Gratton, 2014: B-preconditioned minimization algorithms for  !
!      data assimilation with the dual formulation, QJRMS, 140,        !
!      539-556.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC  :: background
      PUBLIC  :: increment
      PUBLIC  :: analysis
      PUBLIC  :: prior_error

# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F || \
     defined POSTERIOR_EOFS
      PUBLIC  :: posterior_error
# endif
!
!  Set module internal parameters.
!
# ifdef SPLIT_4DVAR
      logical :: Ldone
!
# endif
      integer :: Lini = 1       ! NLM initial conditions record in INI
      integer :: Lbck = 2       ! background record in INI
      integer :: Rec1 = 1
      integer :: Rec2 = 2
!
      CONTAINS
!
      SUBROUTINE background (my_outer, RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the backgound state trajectory, Xb_n-1(t),    !
!  used to linearize the tangent linear and adjoint models in the      !
!  inner loops. It interpolates the background at the observations     !
!  locations, and computes the accept/reject quality control flag,     !
!  ObsScale.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!     RunInterval     NLM kernel time stepping window (seconds)        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
      USE wrt_ini_mod, ONLY : wrt_ini
!
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: lstr, ng
      integer :: Fcount
# ifdef PROFILE
      integer :: thread
# endif
# if defined MODEL_COUPLING && !defined MCT_LIB
      integer :: NstrStep, NendStep, extra
!
      real(dp) :: ENDtime, NEXTtime
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", background"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize and set nonlinear model initial conditions.
!-----------------------------------------------------------------------

# ifdef PROFILE
!
!  Start profile clock.
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 86, __LINE__, MyFile)
        END DO
      END DO
# endif
!
!  Initialize the switch to gather weak constraint forcing.
!
      DO ng=1,Ngrids
        WRTforce(ng)=.FALSE.
      END DO
!
!  Initialize and set nonlinear model initial conditions.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lini" of INI(ng)%name NetCDF file. The
!  record "Lbck" becomes the background state record and the record
!  "Lini" becomes current nonlinear initial conditions.
!
      DO ng=1,Ngrids
        INI(ng)%Rindex=1
        Fcount=INI(ng)%load
        INI(ng)%Nrec(Fcount)=1
        CALL wrt_ini (ng, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Set nonlinear output history file as the initial basic state
!  trajectory.
!
      DO ng=1,Ngrids
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), my_outer
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
      END DO
!
!  Define output 4DVAR NetCDF file containing all processed data
!  at observation locations.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Run nonlinear model and compute basic state trajectory. It processes
!  and writes the observations accept/reject flag (ObsScale) once to
!  allow background quality control, if any.
# if defined MODEL_COUPLING && !defined MCT_LIB
!  Since the ROMS kernel has a delayed output and line diagnostics by
!  one timestep, subtact an extra value to the report of starting and
!  ending timestep for clarity. Usually, the model coupling interval
!  is of the same size as ROMS timestep.
# endif
!-----------------------------------------------------------------------
!
      MyRunInterval=RunInterval
      IF (Master) WRITE (stdout,'(1x)')
!
      DO ng=1,Ngrids
# ifdef RP_AVERAGES
        idefAVG(ng)=-1
        LdefAVG(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), my_outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
        wrtObsScale(ng)=.TRUE.
# if defined MODEL_COUPLING && !defined MCT_LIB
!
        NEXTtime=time(ng)+RunInterval
        ENDtime=INItime(ng)+(ntimes(ng)-1)*dt(ng)
        IF ((NEXTtime.eq.ENDtime).and.(ng.eq.1)) THEN
          extra=0                                   ! last time interval
        ELSE
          extra=1
        END IF
        step_counter(ng)=0
        NstrStep=iic(ng)
        NendStep=NstrStep+INT((MyRunInterval)/dt(ng))-extra
        IF (Master) WRITE (stdout,20) 'NL', ng, my_outer, 0,            &
     &                                NstrStep, NendStep
# else
        IF (Master) WRITE (stdout,20) 'NL', ng, my_outer, 0,            &
     &                                ntstart(ng), ntend(ng)
# endif
      END DO
!
# ifdef SOLVE3D
      CALL main3d (RunInterval)
# else
      CALL main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.FALSE.
        wrtObsScale(ng)=.FALSE.
      END DO

# ifdef PROFILE
!
!  Stop profile clock
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 86, __LINE__, MyFile)
        END DO
      END DO
# endif
!
 10   FORMAT (a,'_outer',i0,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS: started time-stepping:',                 &
     &        ' (Grid: ',i0,', Outer: ',i2.2,', Inner: ',i3.3,          &
              ', TimeSteps: ',i0,' - ',i0,')',/)
!
      RETURN
      END SUBROUTINE background
!
      SUBROUTINE increment (my_outer, RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the 4D-Var data assimilation state increment, !
!  dXa, by iterating the inner loops and minimizing the cost function. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!     RunInterval     TLM/ADM kernels time stepping window (seconds)   !
!                                                                      !
!                                                                      !
!  It solves the system:                                               !
!                                                                      !
!     (R_n + Cobs) * Beta_n = h_n                                      !
!                                                                      !
!     h_n = Xo - H * X_n                                               !
!                                                                      !
!  where R_n is the representer matrix, Cobs is the observation-error  !
!  covariance, Beta_n are the representer coefficients, h_n is the     !
!  misfit between observations (Xo) and model (H * X_n), and H is      !
!  the linearized observation operator. Here, _n denotes a sequence    !
!  of estimates.                                                       !
!                                                                      !
!  The system does not need to be solved explicitly by inverting the   !
!  symmetric stabilized representer matrix, P_n:                       !
!                                                                      !
!     P_n = R_n + Cobs                                                 !
!                                                                      !
!  but by computing the action of P_n on any vector PSI, such that     !
!                                                                      !
!     P_n * PSI = R_n * PSI + Cobs * PSI                               !
!                                                                      !
!  The representer matrix is not explicitly computed but evaluated by  !
!  one integration backward of the adjoint model and one integration   !
!  forward of the tangent linear model for any forcing vector PSI.     !
!                                                                      !
!  A preconditioned conjugate gradient algorithm is used to compute    !
!  an approximation PSI for Beta_n.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE congrad_mod,      ONLY : congrad
# ifdef SPLIT_4DVAR
      USE congrad_mod,      ONLY : cg_read_congrad
# endif
      USE convolve_mod,     ONLY : error_covariance
# ifdef ADJUST_BOUNDARY
      USE mod_boundary,     ONLY : initialize_boundary
# endif
      USE mod_forces,       ONLY : initialize_forces
      USE mod_ocean,        ONLY : initialize_ocean
      USE strings_mod,      ONLY : FoundError
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY : balance_ref, biconj
# endif
!
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: Lcgini, Linner, Lposterior, add
!
      integer :: i, ifile, lstr, my_inner, ng, tile
      integer :: Fcount, InpRec, Tindex
# ifdef PROFILE
      integer :: thread
# endif
# ifdef SPLIT_4DVAR
!
      real(dp) :: stime
# endif
!
      character (len=6  ) :: driver = 'r4dvar'
      character (len=10 ) :: suffix
      character (len=20 ) :: string
# ifdef SPLIT_4DVAR
      character (len=256) :: ncname
# endif
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", increment"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Run representer model and compute a "prior estimate" state
!  trajectory, X_n(t). Use linearized state trajectory (X_n-1) as
!  basic state.
!-----------------------------------------------------------------------

# ifdef PROFILE
!
!  Start profile clock.
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iTLM, 87, __LINE__, MyFile)
        END DO
      END DO
# endif

# ifdef SPLIT_4DVAR
!
!-----------------------------------------------------------------------
!  If split 4D-Var algorithm, set several variables that are computed
!  or assigned in other 4D-Var phase executable.
!-----------------------------------------------------------------------
!
!  Reset Nrun counter to a value greater than one.
!
      IF (my_outer.gt.1) THEN
        Nrun=1+(my_outer-1)*Ninner
      END IF
!
!  Open Nonlinear model initial conditions NetCDF file (INI struc) and
!  inquire about its variables IDs.
!
      DO ng=1,Ngrids
        LdefINI(ng)=.FALSE.
        CALL def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  If outer>1, open tangent linear initial conditions NetCDF file
!  (ITL struc) and inquire about its variables IDs.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefITL(ng)=.FALSE.
          CALL tl_def_ini (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!  If outer>1, open representer model initial conditions NetCDF file
!  (IRP struc) and inquire about its variables IDs.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefIRP(ng)=.FALSE.
          CALL rp_def_ini (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!  If outer>1, open impulse forcing NetCDF file (TLF struc) and inquire
!  about its variables IDs.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefTLF(ng)=.FALSE.
          CALL def_impulse (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!  If outer>1, open adjoint history NetCDF file (ADM struc) and inquire
!  about its variables IDs and set "FrcRec".
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefADJ(ng)=.FALSE.
          CALL ad_def_his (ng, LdefADJ(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Fcount=ADM(ng)%Fcount
          FrcRec(ng)=ADM(ng)%Nrec(Fcount)
        END DO
      END IF
!
!  Open 4D-Var NetCDF file (DAV struc) and inquire about its variables.
!  Deactivate "haveNLmod" since its values are read below. Otherwise,
!  the TLM will read zero values when calling "obs_read" for outer>1.
!  The DAV file does not have values for NLmodel_value(:,outer) yet.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.FALSE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        haveNLmod(ng)=.FALSE.        ! because is activated in "def_mod"
      END DO
!
!  If split 4D-Var and outer>1, read several variables from 4D-VAR
!  NetCDF (DAV struc) file needed in the conjugate gradient algorithm.
!  In the unsplit case, such values are available in memory.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          CALL cg_read_congrad (ng, iTLM, my_outer)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!  In the split 4D-Var algorithm, we need to read in initial NLmodVal
!  ("NLmodel_initial") and ObsScale ("obs_scale"), computed in the
!  background phase, from the DAV NetCDF file. Such data is in memory
!  in the unsplit algorithm.
!
!  HGA: What to do in 4D-Var with nested grids?
!
      IF (my_outer.eq.1) THEN
        ng=1                         ! initial values from "background"
        CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                   &
     &                        Vname(1,idNLmi), NLmodVal,                &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      ng=1
      CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                     &
     &                      Vname(1,idObsS), ObsScale,                  &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  In the split 4D-Var alorithm, we also need to read in ObsVal
!  ("obs_value") and ObsErr ("obs_error") from the observation file.
!  They are used in the initialization of the conjugate gradient
!  algorithm.
!
      ng=1
      CALL netcdf_get_fvar (ng, iTLM, OBS(ng)%name,                     &
     &                      Vname(1,idOval), ObsVal)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL netcdf_get_fvar (ng, iTLM, OBS(ng)%name,                     &
     &                      Vname(1,idOerr), ObsErr)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set nonlinear output history file to be used as the basic state
!  trajectory. The nonlinear model is run outside of the outer loops,
!  so we need outer=0 when updating HIS structure.  Notice that this
!  structure is needed to set the BLK structure below.
!
      DO ng=1,Ngrids
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), 0
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
        IF (HIS(ng)%Nfiles.gt.1) THEN
          DO ifile=1,HIS(ng)%Nfiles
            WRITE (suffix,"('_',i4.4,'.nc')") ifile
            HIS(ng)%files(ifile)=TRIM(HIS(ng)%base)//TRIM(suffix)
          END DO
          HIS(ng)%name=TRIM(HIS(ng)%files(1))
        ELSE
          HIS(ng)%files(1)=TRIM(HIS(ng)%name)
        END IF
      END DO
!
!  Read in 4D-Var starting time (sec) from forward trajectory NetCDF
!  file and set INItime. The starting time is need for boundary and
!  surface forcing adjustments, if any.
!
      InpRec=1
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iTLM, HIS(ng)%name,                   &
     &                        Vname(1,idtime), stime,                   &
     &                        start = (/InpRec/),                       &
     &                        total = (/1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        INItime(ng)=stime

#  ifdef ADJUST_BOUNDARY
!
!  Set time (sec) for the open boundary adjustment.
!
        OBC_time(1,ng)=stime
        DO i=2,Nbrec(ng)
          OBC_time(i,ng)=OBC_time(i-1,ng)+nOBC(ng)*dt(ng)
        END DO
#  endif

#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Set time (sec) for the surface forcing adjustment.
!
        SF_time(1,ng)=stime
        DO i=2,Nfrec(ng)
          SF_time(i,ng)=SF_time(i-1,ng)+nSFF(ng)*dt(ng)
        END DO
#  endif
      END DO

#  if defined ADJUST_BOUNDARY || \
      defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
!
!  If outer>1, read in adjust forcing arrays from IRP file (Rec2),
!  which are computed in the previous outer loop when calling the
!  "error_covarinace" after the inner loops. Such values are in
!  memory in the unsplit algorithm. Use "iTLM" descriptor so the
!  switch get_adjust=.TRUE.  Notice that all the TL arrays are
!  cleared in the next call to "rp_initial" except the "tl_*_obc"
!  arrays used when adjusting OBCs.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          InpRec=Rec2
          Tindex=1
          CALL get_state (ng, iTLM, 1, IRP(ng)%name, InpRec, Tindex)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#  endif
!
!  If outer>1, set previous outer loop representer model file.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          WRITE (TLM(ng)%name,10) TRIM(FWD(ng)%head), my_outer-1
          lstr=LEN_TRIM(TLM(ng)%name)
          TLM(ng)%base=TLM(ng)%name(1:lstr-3)
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  On first pass (outer=1), create NetCDF files and initialize all the
!  records in the ITL file to zero.
!-----------------------------------------------------------------------
!
      CHECK_OUTER1 : IF (my_outer.eq.1) THEN
!
!  If outer=1, define tangent linear initial conditions file.
!
        DO ng=1,Ngrids
          LdefITL(ng)=.TRUE.
          CALL tl_def_ini (ng)
          LdefITL(ng)=.FALSE.
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  If outer=1, define TLM/RPM impulse forcing NetCDF file.
!
        DO ng=1,Ngrids
          LdefTLF(ng)=.TRUE.
          CALL def_impulse (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

# if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
     defined POSTERIOR_ERROR_F
!
!  Define output Hessian NetCDF file that will eventually contain
!  the intermediate posterior analysis error covariance matrix
!  fields or its EOFs.
!
        DO ng=1,Ngrids
          LdefHSS(ng)=.TRUE.
          CALL def_hessian (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif

# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
!
!  Define output initial or final full posterior error covariance
!  (diagonal) matrix NetCDF.
!
        DO ng=1,Ngrids
          LdefERR (ng)=.TRUE.
          CALL def_error (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif
      END IF CHECK_OUTER1
!
!-----------------------------------------------------------------------
!  Run finite amplitude tangent linear (representers) model.
!-----------------------------------------------------------------------
!
!  Clear the nonlinear state arrays to make sure that the RHS terms
!  rzeta, rubar, rvbar, ru, and rv are zero.  Otherwise, those arrays
!  will have the last computed values when running the nonlinear model
!  if not processing the forward trajectory RHS terms. This needs to
!  be done to get identical solutions with the split schemes.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
        END DO
      END DO
!
!  Set structure for the nonlinear forward trajectory to be processed
!  by the tangent linear and adjoint models. Also, set switches to
!  process the FWD structure in routine "check_multifile". Notice that
!  it is possible to split solution into multiple NetCDF files to reduce
!  their size.
!
      IF (my_outer.eq.1) THEN
        CALL edit_multifile ('HIS2FWD')
      ELSE
        CALL edit_multifile ('TLM2FWD')
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      DO ng=1,Ngrids
        LreadFWD(ng)=.TRUE.
      END DO

# ifdef FORWARD_FLUXES
!
!  Set the BLK structure to contain the nonlinear model surface fluxes
!  needed by the tangent linear and adjoint models. Also, set switches
!  to process that structure in routine "check_multifile". Notice that
!  it is possible to split the solution into multiple NetCDF files to
!  reduce their size.
!
!  The switch LreadFRC is deactivated because all the atmospheric
!  forcing, including shortwave radiation, is read from the NLM
!  surface fluxes or is assigned during ESM coupling.  Such fluxes
!  are available from the QCK structure. There is no need for reading
!  and processing from the FRC structure input forcing-files.
!
      CALL edit_multifile ('QCK2BLK')
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      DO ng=1,Ngrids
        LreadBLK(ng)=.TRUE.
        LreadFRC(ng)=.FALSE.
      END DO
# endif
!
!  Set representer model output file name.  The strategy is to write
!  the representer solution at the beginning of each outer loop.
!
      DO ng=1,Ngrids
        idefTLM(ng)=-1
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        WRITE (TLM(ng)%name,10) TRIM(FWD(ng)%head), my_outer
        lstr=LEN_TRIM(TLM(ng)%name)
        TLM(ng)%base=TLM(ng)%name(1:lstr-3)
      END DO
!
!  Activate switch to write the representer model at observation points.
!  Turn off writing into history file and turn off impulse forcing.
!
      DO ng=1,Ngrids
        wrtRPmod(ng)=.TRUE.
        SporadicImpulse(ng)=.FALSE.
        FrequentImpulse(ng)=.FALSE.
      END DO

# ifndef DATALESS_LOOPS
!
!  As in the nonlinear model, initialize always the representer model
!  here with the background or reference state (IRP(ng)%name, record
!  Rec1).
!
      DO ng=1,Ngrids
        IRP(ng)%Rindex=Rec1
        CALL rp_initial (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Run representer model using the nonlinear trajectory as a basic
!  state.  Compute model solution at observation points, H * X_n.
!
      DO ng=1,Ngrids
#  ifdef RP_AVERAGES
        LdefAVG(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
#  endif
        IF (Master) THEN
          WRITE (stdout,20) 'RP', ng, my_outer, 0,                      &
     &                      ntstart(ng), ntend(ng)
        END IF
      END DO
!
#  ifdef SOLVE3D
      CALL rp_main3d (RunInterval)
#  else
      CALL rp_main2d (RunInterval)
#  endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Report data penalty function.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          DO i=0,NstateVar(ng)
            IF (i.eq.0) THEN
              string='Total'
            ELSE
              string=Vname(1,idSvar(i))
            END IF
            IF (FOURDVAR(ng)%DataPenalty(i).ne.0.0_r8) THEN
              WRITE (stdout,30) my_outer, inner, 'RPM',                 &
     &                          FOURDVAR(ng)%DataPenalty(i),            &
     &                          TRIM(string)
            END IF
          END DO
        END IF
!
!  Write out initial data penalty function to NetCDF file.
!
        SourceFile=MyFile
        CALL netcdf_put_fvar (ng, iRPM, DAV(ng)%name,                   &
     &                        'RP_iDataPenalty',                        &
     &                        FOURDVAR(ng)%DataPenalty(0:),             &
     &                        (/1,my_outer/), (/NstateVar(ng)+1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Clean penalty array before next run of RP model.
!
        FOURDVAR(ng)%DataPenalty=0.0_r8
      END DO
!
!  Turn off IO switches.
!
      DO ng=1,Ngrids
        LdefTLM(ng)=.FALSE.
        LwrtTLM(ng)=.FALSE.
        wrtRPmod(ng)=.FALSE.
      END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero after
!  running the representer model and must be zero when running the
!  tangent linear model.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_forces (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iTLM)
#  endif
        END DO
      END DO

#  if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
      IF (balance(isFsur)) THEN
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL balance_ref (ng, tile, Lini)
            CALL biconj (ng, tile, iNLM, Lini)
          END DO
          wrtZetaRef(ng)=.TRUE.
        END DO
      END IF
#  endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  4D-Var inner loops.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      INNER_LOOP : DO my_inner=0,Ninner
        inner=my_inner
!
! Initialize conjugate gradient algorithm depending on hot start or
! outer loop index.
!
        IF (my_inner.eq.0) THEN
          Lcgini=.TRUE.
          DO ng=1,Ngrids
            CALL congrad (ng, iRPM, my_outer, my_inner, Ninner, Lcgini)
          END DO
        END IF
!
!  If initialization step, skip the inner-loop computations.
!
        Linner=.FALSE.
        IF ((my_inner.ne.0).or.(Nrun.ne.1)) THEN
          IF (((my_inner.eq.0).and.LhotStart).or.(my_inner.ne.0)) THEN
            Linner=.TRUE.
          END IF
        END IF
!
!  Start inner loop computations.
!
        INNER_COMPUTE : IF (Linner) THEN
!
!-----------------------------------------------------------------------
!  Integrate adjoint model forced with any vector PSI at the observation
!  locations and generate adjoint trajectory, Lambda_n(t).
!-----------------------------------------------------------------------
!
!  Initialize the adjoint model from rest.
!
          DO ng=1,Ngrids
            CALL ad_initial (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            wrtMisfit(ng)=.FALSE.
          END DO

#  ifdef RPM_RELAXATION
!
!  Adjoint of representer relaxation is not applied during the
!  inner-loops.
!
          DO ng=1,Ngrids
            LweakRelax(ng)=.FALSE.
          END DO
#  endif
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file only once to avoid opening too many files.
!
          DO ng=1,Ngrids
            WRTforce(ng)=.TRUE.
            IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
            Fcount=ADM(ng)%load
            ADM(ng)%Nrec(Fcount)=0
            ADM(ng)%Rindex=0
          END DO
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'AD', ng, my_outer, my_inner,           &
     &                          ntstart(ng), ntend(ng)
            END IF
          END DO
!
#  ifdef SOLVE3D
          CALL ad_main3d (RunInterval)
#  else
          CALL ad_main2d (RunInterval)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Activate "LwrtState2d" switch to write the correct ad_ubar, ad_vbar,
!  ad_zeta here instead of ad_ubar_sol, ad_vbar_sol, and ad_zeta_sol in
!  the calls to "ad_wrt_his" below (HGA).
!  Here, ad_zeta(:,:,kout)=ad_zeta_sol. However, ad_ubar and ad_vbar
!  not equal to ad_ubar_sol and ad_vbar_sol, respectively. It is
!  irrelevant because ubar and vbar are not part of the state in
!  3D application.  Notice that "LwrtState2d" will be turned off at
!  the bottom of "error_covariance".
!
          DO ng=1,Ngrids
            LwrtState2d(ng)=.TRUE.
          END DO
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.
!
          DO ng=1,Ngrids
            CALL ad_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
          DO ng=1,Ngrids
            WRTforce(ng)=.FALSE.
            LwrtState2d(ng)=.TRUE.
            CALL ad_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!  Convolve adjoint trajectory with error covariances.  Notice that we
!  use iTLM to write the convolved solution into the ITL(ng)%name file
!  (record Rec1) as the TLM initial conditions for the next inner loop.
!
#  ifdef POSTERIOR_ERROR_I
          Lposterior=.TRUE.
#  else
          Lposterior=.FALSE.
#  endif
          CALL error_covariance (iTLM, driver, my_outer, my_inner,      &
     &                           Lbck, Lini, Lold, Lnew,                &
     &                           Rec1, Rec2, Lposterior)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
          IF (Master) THEN
            WRITE (stdout,40) my_outer, my_inner
          END IF
          DO ng=1,Ngrids
            TLF(ng)%Rindex=0
#  ifdef DISTRIBUTE
            tile=MyRank
#  else
            tile=-1
#  endif
            CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!-----------------------------------------------------------------------
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!-----------------------------------------------------------------------
!
          DO ng=1,Ngrids
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.TRUE.
          END DO
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
          DO ng=1,Ngrids
            IF (FrcRec(ng).gt.3) THEN
              FrequentImpulse(ng)=.TRUE.
            END IF
          END DO
!
!  Initialize tangent linear model from ITL(ng)%name, record Rec1.
!
          DO ng=1,Ngrids
            ITL(ng)%Rindex=Rec1
            CALL tl_initial (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!  Activate switch to write out initial misfit between model and
!  observations.
!
          IF ((my_outer.eq.1).and.(my_inner.eq.1)) THEN
            DO ng=1,Ngrids
              wrtMisfit(ng)=.TRUE.
            END DO
          END IF
!
!  Set tangent linear history NetCDF parameters.  Define tangent linear
!  history file at the beggining of each inner loop  to avoid opening
!  too many NetCDF files.
!
          DO ng=1,Ngrids
            IF (my_inner.gt.1) LdefTLM(ng)=.FALSE.
            Fcount=TLM(ng)%load
            TLM(ng)%Nrec(Fcount)=0
            TLM(ng)%Rindex=0
          END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute R_n * PSI at observation points which
!  are used in the conjugate gradient algorithm.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ng, my_outer, my_inner,           &
     &                          ntstart(ng), ntend(ng)
            END IF
          END DO
!
#  ifdef SOLVE3D
          CALL tl_main3d (RunInterval)
#  else
          CALL tl_main2d (RunInterval)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO ng=1,Ngrids
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.FALSE.
          END DO

#  ifdef POSTERIOR_ERROR_F
!
!  Copy the final time tl_var(Lold) into ad_var(Lold) so that it can be
!  written to the Hessian NetCDF file.
!
          add=.FALSE.
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL load_TLtoAD (ng, tile, Lold(ng), Lold(ng), add)
            END DO
          END DO
!
!  Write evolved tangent solution into hessian netcdf file for use
!  later.
!
          IF (my_inner.ne.0) THEN
            DO ng=1,Ngrids
              CALL wrt_hessian (ng, Lold(ng), Lold(ng))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
          END IF
#  endif
!
!-----------------------------------------------------------------------
!  Use conjugate gradient algorithm to find a better approximation
!  PSI to representer coefficients Beta_n.
!-----------------------------------------------------------------------
!
          Nrun=Nrun+1
          DO ng=1,Ngrids
            Lcgini=.FALSE.
            CALL congrad (ng, iRPM, my_outer, my_inner, Ninner, Lcgini)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO

        END IF INNER_COMPUTE

      END DO INNER_LOOP
!
!  Close tangent linear NetCDF file.
!
      SourceFile=MyFile
      DO ng=1,Ngrids
        CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Once that the representer coefficients, Beta_n, have been
!  approximated with sufficient accuracy, compute estimates of
!  Lambda_n and Xhat_n by carrying out one backward intergration
!  of the adjoint model and one forward itegration of the representer
!  model.
!-----------------------------------------------------------------------
!
!  Initialize the adjoint model always from rest.
!
      DO ng=1,Ngrids
        CALL ad_initial (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

#  ifdef RPM_RELAXATION
!
!  Adjoint of representer relaxation is applied during the
!  outer-loops.
!
      DO ng=1,Ngrids
        LweakRelax(ng)=.TRUE.
      END DO
#  endif
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
      DO ng=1,Ngrids
        WRTforce(ng)=.TRUE.
        IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
        Fcount=ADM(ng)%load
        ADM(ng)%Nrec(Fcount)=0
        ADM(ng)%Rindex=0
      END DO
!
!  Time-step adjoint model backwards forced with estimated representer
!  coefficients, Beta_n.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,20) 'AD', ng, my_outer, Ninner,                 &
     &                      ntstart(ng), ntend(ng)
        END IF
      END DO
!
#  ifdef SOLVE3D
      CALL ad_main3d (RunInterval)
#  else
      CALL ad_main2d (RunInterval)
#  endif

      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
      DO ng=1,Ngrids
        CALL ad_wrt_his (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
      DO ng=1,Ngrids
        WRTforce(ng)=.FALSE.
        CALL ad_wrt_his (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Convolve adjoint trajectory with error covariances.  Notice that we
!  use iRPM to write the convolved solution into the IRP(ng)%name file
!  (record Rec2) as the RPM initial conditions for the next outer loop.
!
      Lposterior=.FALSE.
      CALL error_covariance (iRPM, driver, my_outer, inner,             &
     &                       Lbck, Lini, Lold, Lnew,                    &
     &                       Rec1, Rec2, Lposterior)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
      IF (Master) THEN
        WRITE (stdout,40) my_outer, inner
      END IF
      DO ng=1,Ngrids
        TLF(ng)%Rindex=0
#  ifdef DISTRIBUTE
        tile=MyRank
#  else
        tile=-1
#  endif
        CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

# endif /* !DATALESS_LOOPS */

# ifdef PROFILE
!
!  Stop profile clock
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 87, __LINE__, MyFile)
        END DO
      END DO
# endif
!
 10   FORMAT (a,'_outer',i0,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS: started time-stepping:',                 &
     &        ' (Grid: ',i0,', Outer: ',i2.2,', Inner: ',i3.3,          &
              ', TimeSteps: ',i0,' - ',i0,')',/)
 30   FORMAT (' (',i3.3,',',i3.3,'): ',a,' data penalty, Jdata = ',     &
     &        1p,e17.10,0p,t68,a)
 40   FORMAT (/,' Converting Convolved Adjoint Trajectory to',          &
     &          ' Impulses: Outer = ',i3.3,' Inner = ',i3.3,/)
!
      RETURN
      END SUBROUTINE increment
!
      SUBROUTINE analysis (my_outer, RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes 4D-Var data assimilation analysis, Xa. The    !
!  nonlinear model initial conditions are computed by adding the       !
!  4D-Var increments to the current background:  Xa = Xb + dXa.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!     RunInterval     NLM kernel time stepping window (seconds)        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE strings_mod,  ONLY : FoundError
!
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: i, ifile, lstr, ng
      integer :: InpRec
# ifdef PROFILE
      integer :: thread
# endif
# ifdef SPLIT_4DVAR
!
      real(dp) :: stime
# endif
!
      character (len=10) :: suffix
      character (len=20) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", analysis"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Run representer model and compute a "new estimate" of the state
!  trajectory, X_n(t).
!-----------------------------------------------------------------------

# ifdef PROFILE
!
!  Start profile clock.
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 88, __LINE__, MyFile)
        END DO
      END DO
# endif

# ifdef SPLIT_4DVAR
!
!-----------------------------------------------------------------------
!  If split 4D-Var algorithm, set several variables that computed or
!  assigned in other 4D-Var phase executable.
!-----------------------------------------------------------------------
!
!  Set Nrun>1, to read in surface forcing and open boundary conditions
!  increments in "initial", if appropriate.
!
      Nrun=Ninner+1
!
!  Set ERstr=Nrun, to set the open boundary condition (OBC_time) and
!  surface forcing (SF_time) adjustment times, if needed.
!
      ERstr=Nrun
!
!  Open 4D-Var NetCDF file (DAV struc) and inquire about its variables.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.FALSE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  In the split 4D-Var algorithm, we need to read the values of ObsScale
!  computed in the background phase from the DAV NetCDF file. Such data
!  is in memory in the unsplit algorithm.
!
!  HGA: What to do in 4D-Var with nested grids?
!
      ng=1
      CALL netcdf_get_fvar (ng, iRPM, DAV(ng)%name,                     &
     &                      Vname(1,idObsS), ObsScale,                  &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set nonlinear output history file to be used as the basic state
!  trajectory. The nonlinear model is run outside of the outer loops,
!  so we need outer=0 when updating HIS structure.
!
      DO ng=1,Ngrids
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), 0
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
        IF (HIS(ng)%Nfiles.gt.1) THEN
          DO ifile=1,HIS(ng)%Nfiles
            WRITE (suffix,"('_',i4.4,'.nc')") ifile
            HIS(ng)%files(ifile)=TRIM(HIS(ng)%base)//TRIM(suffix)
          END DO
          HIS(ng)%name=TRIM(HIS(ng)%files(1))
        ELSE
          HIS(ng)%files(1)=TRIM(HIS(ng)%name)
        END IF
      END DO
!
!  Read in 4D-Var starting time (sec) from forward trajectory NetCDF
!  file and set INItime. The starting time is need for boundary and
!  surface forcing adjustments, if any.
!
      InpRec=1
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iRPM, HIS(ng)%name,                   &
     &                        Vname(1,idtime), stime,                   &
     &                        start = (/InpRec/),                       &
     &                        total = (/1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        INItime(ng)=stime

#  ifdef ADJUST_BOUNDARY
!
!  Set time (sec) for the open boundary adjustment.
!
        OBC_time(1,ng)=stime
        DO i=2,Nbrec(ng)
          OBC_time(i,ng)=OBC_time(i-1,ng)+nOBC(ng)*dt(ng)
        END DO
#  endif

#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Set time (sec) for the surface forcing adjustment.
!
        SF_time(1,ng)=stime
        DO i=2,Nfrec(ng)
          SF_time(i,ng)=SF_time(i-1,ng)+nSFF(ng)*dt(ng)
        END DO
#  endif
      END DO
!
!  If outer>1, set previous outer loop representer model file.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          WRITE (TLM(ng)%name,10) TRIM(FWD(ng)%head), my_outer-1
          lstr=LEN_TRIM(TLM(ng)%name)
          TLM(ng)%base=TLM(ng)%name(1:lstr-3)
        END DO
      END IF
!
!  Set basic state trajectory to either Nonlinear model (outer=1) or
!  previous outer loop representer model (outer>1). Also, set switches
!  to process the FWD structure in routine "check_multifile". Notice
!  that it is possible to split solution into multiple NetCDF files to
!  reduce their size.
!
      IF (my_outer.eq.1) THEN
        CALL edit_multifile ('HIS2FWD')
      ELSE
        CALL edit_multifile ('TLM2FWD')
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      DO ng=1,Ngrids
        LreadFWD(ng)=.TRUE.
      END DO

#  ifdef FORWARD_FLUXES
!
!  Set the BLK structure to contain the nonlinear model surface fluxes
!  needed by the tangent linear and adjoint models. Also, set switches
!  to process that structure in routine "check_multifile". Notice that
!  it is possible to split the solution into multiple NetCDF files to
!  reduce their size.
!
!  The switch LreadFRC is deactivated because all the atmospheric
!  forcing, including shortwave radiation, is read from the NLM
!  surface fluxes or is assigned during ESM coupling.  Such fluxes
!  are available from the QCK structure. There is no need for reading
!  and processing from the FRC structure input forcing-files.
!
      CALL edit_multifile ('QCK2BLK')
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      DO ng=1,Ngrids
        LreadBLK(ng)=.TRUE.
        LreadFRC(ng)=.FALSE.
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Run representer model with the new analysis initial conditions
!  cpmputed in the "increment" phase.
!-----------------------------------------------------------------------
!
!  Set new basic state trajectory for next outer loop.
!
      DO ng=1,Ngrids
        idefTLM(ng)=-1
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.TRUE.
        wrtRPmod(ng)=.TRUE.
        WRITE (TLM(ng)%name,10) TRIM(FWD(ng)%head), my_outer
        lstr=LEN_TRIM(TLM(ng)%name)
        TLM(ng)%base=TLM(ng)%name(1:lstr-3)
      END DO
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
      DO ng=1,Ngrids
        IF (FrcRec(ng).gt.3) THEN
          FrequentImpulse(ng)=.TRUE.
        END IF
      END DO
!
!  Initialize representer model IRP(ng)%name file, record Rec2.
!
      DO ng=1,Ngrids
# ifdef DATALESS_LOOPS
        IRP(ng)%Rindex=Rec1
# else
        IRP(ng)%Rindex=Rec2
# endif
        CALL rp_initial (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Activate switch to write out final misfit between model and
!  observations.
!
      IF (my_outer.eq.Nouter) THEN
        DO ng=1,Ngrids
          wrtMisfit(ng)=.TRUE.
        END DO
      END IF
!
!  Run representer model using previous linearized trajectory, X_n-1, as
!  basic state and forced with convolved adjoint trajectory impulses.
!
      DO ng=1,Ngrids
# ifdef RP_AVERAGES
        idefAVG(ng)=-1
        LdefAVG(ng)=.TRUE.
        LwrtAVG(ng)=.TRUE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), my_outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
        IF (Master) THEN
          WRITE (stdout,20) 'RP', ng, my_outer, Ninner,                 &
     &                      ntstart(ng), ntend(ng)
        END IF
      END DO
!
# ifdef SOLVE3D
      CALL rp_main3d (RunInterval)
# else
      CALL rp_main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
# ifdef RP_AVERAGES
      DO ng=1,Ngrids
        LdefAVG(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
      END DO
# endif
!
!  Report data penalty function.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          DO i=0,NstateVar(ng)
            IF (i.eq.0) THEN
              string='Total'
            ELSE
              string=Vname(1,idSvar(i))
            END IF
            IF (FOURDVAR(ng)%DataPenalty(i).ne.0.0_r8) THEN
              WRITE (stdout,30) my_outer, Ninner, 'RPM',                &
     &                          FOURDVAR(ng)%DataPenalty(i),            &
     &                          TRIM(string)
# ifdef DATALESS_LOOPS
              WRITE (stdout,30) my_outer, Ninner, 'NLM',                &
     &                          FOURDVAR(ng)%NLPenalty(i),              &
     &                          TRIM(string)
# endif
            END IF
          END DO
        END IF
      END DO
!
!  Write out final data penalty function to NetCDF file.
!
      SourceFile=MyFile
      DO ng=1,Ngrids
        CALL netcdf_put_fvar (ng, iRPM, DAV(ng)%name,                   &
     &                        'RP_fDataPenalty',                        &
     &                        FOURDVAR(ng)%DataPenalty(0:),             &
     &                        (/1,my_outer/), (/NstateVar(ng)+1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Clean array before next run of RP model.
!
      DO ng=1,Ngrids
        FOURDVAR(ng)%DataPenalty=0.0_r8
# ifdef DATALESS_LOOPS
        FOURDVAR(ng)%NLPenalty=0.0_r8
# endif
        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
      END DO
!
!  Close current forward NetCDF file.  If last outer loop, set history
!  file ID to closed state since we manipulated its indices with the
!  forward file ID, which it was closed.
!
      DO ng=1,Ngrids
        CALL netcdf_close (ng, iRPM, FWD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        IF (my_outer.eq.Nouter) HIS(ng)%ncid=-1
      END DO

# ifdef PROFILE
!
!  Stop profile clock
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 88, __LINE__, MyFile)
        END DO
      END DO
# endif
!
 10   FORMAT (a,'_outer',i0,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS: started time-stepping:',                 &
     &        ' (Grid: ',i0,', Outer: ',i2.2,', Inner: ',i3.3,          &
              ', TimeSteps: ',i0,' - ',i0,')',/)
 30   FORMAT (' (',i3.3,',',i3.3,'): ',a,' data penalty, Jdata = ',     &
     &        1p,e17.10,0p,t68,a)
!
      RETURN
      END SUBROUTINE analysis
!
      SUBROUTINE prior_error (ng)
!
!=======================================================================
!                                                                      !
!  This routine processes background prior error covariance standard   !
!  deviations and normalization coefficients.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng              Nested grid number                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE normalization_mod, ONLY : normalization
      USE strings_mod,       ONLY : FoundError
!
!  Imported variable declarations
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: tile
      integer :: NRMrec, STDrec, Tindex
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", prior_error"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set application grid, metrics, and associated variables and
!  parameters.
!-----------------------------------------------------------------------
!
!  The ROMS application grid configuration is done once. It is usually
!  done in the "initial" kernel routine. However, since we are calling
!  the "normalization" routine here, we need several grid variables and
!  parameter.  Also, if reading only water points, we need to know the
!  land/sea mask arrays to unpack.
!
      IF (SetGridConfig(ng)) THEN
        CALL set_grid (ng, iNLM)
      END IF
!
!-----------------------------------------------------------------------
!  Read in standard deviation factors for error covariance.
!-----------------------------------------------------------------------
!
!  Initial conditions standard deviation. They are loaded in Tindex=1
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 10, 10, STD(1,ng)%name, STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      IF (NSA.eq.2) THEN
        CALL get_state (ng, 11, 11, STD(2,ng)%name, STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 12, 12, STD(3,ng)%name, STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 13, 13, STD(4,ng)%name, STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Error covariance normalization coefficients.
!-----------------------------------------------------------------------
!
!  Compute or read in the error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid and decorrelation scales.
!
      IF (ANY(LwrtNRM(:,ng))) THEN
        CALL def_norm (ng, iNLM, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        IF (NSA.eq.2) THEN
          CALL def_norm (ng, iNLM, 2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef ADJUST_BOUNDARY
        CALL def_norm (ng, iNLM, 3)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        CALL def_norm (ng, iNLM, 4)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL normalization (ng, tile, 2)
        END DO
        LdefNRM(1:4,ng)=.FALSE.
        LwrtNRM(1:4,ng)=.FALSE.

      ELSE
        NRMrec=1
        CALL get_state (ng, 14, 14, NRM(1,ng)%name, NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        IF (NSA.eq.2) THEN
          CALL get_state (ng, 15, 15, NRM(2,ng)%name, NRMrec, 2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef ADJUST_BOUNDARY
        CALL get_state (ng, 16, 16, NRM(3,ng)%name, NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        CALL get_state (ng, 17, 17, NRM(4,ng)%name, NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
      END IF
!
      RETURN
      END SUBROUTINE prior_error

# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F || \
     defined POSTERIOR_EOFS
!
      SUBROUTINE posterior_error (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes posterior analysis error covariance matrix    !
!  using the Lanczos vectors.                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     RunInterval     Kernel time stepping window (seconds)            !
!                                                                      !
!  Warning:                                                            !
!                                                                      !
!  Currently, this code only works for a single outer-loop.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE convolve_mod,      ONLY : convolve
#  ifdef ADJUST_BOUNDARY
      USE mod_boundary,      ONLY : initialize_boundary
#  endif
      USE mod_forces,        ONLY : initialize_forces
      USE mod_ocean,         ONLY : initialize_ocean
      USE posterior_mod,     ONLY : posterior
#  if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      USE posterior_var_mod, ONLY : posterior_var
#  endif
      USE random_ic_mod,     ONLY : random_ic
      USE strings_mod,       ONLY : FoundError
!
!  Imported variable declarations
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: Ltrace
!
      integer :: my_inner, my_outer, ng, tile
      integer :: Fcount, Rec
!
      character (len=6) :: driver = 'w4dvar'

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", posterior_error"
!
      SourceFile=MyFile

#  if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
!
!-----------------------------------------------------------------------
!  Compute full (diagonal) posterior analysis error covariance matrix.
!-----------------------------------------------------------------------
!
!  Clear tangent and adjoint arrays because they are used as
!  work arrays below.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
          CALL initialize_ocean (ng, tile, iTLM)
          CALL initialize_forces (ng, tile, iADM)
          CALL initialize_forces (ng, tile, iTLM)
#   ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iADM)
          CALL initialize_boundary (ng, tile, iTLM)
#   endif
        END DO
      END DO
!
!   Compute the diagonal of the posterior/analysis error covariance
!   matrix. The result is written to record 2 of the ITL netcdf file.
!
      VAR_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL posterior_var (ng, tile, iTLM, my_outer)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END DO VAR_OLOOP
!
!  Write out the diagonal of the posterior/analysis covariance matrix
!  which is in tl_var(Rec1) to 4DVar error NetCDF file.
!
      DO ng=1,Ngrids
        CALL wrt_error (ng, Rec1, Rec1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Clear tangent and adjoint arrays because they are used as
!  work arrays below.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
          CALL initialize_ocean (ng, tile, iTLM)
          CALL initialize_forces (ng, tile, iADM)
          CALL initialize_forces (ng, tile, iTLM)
#   ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iADM)
          CALL initialize_boundary (ng, tile, iTLM)
#   endif
        END DO
      END DO
#  endif

#  ifdef POSTERIOR_EOFS
!
!-----------------------------------------------------------------------
!  Compute the posterior analysis error covariance matrix EOFs using a
!  Lanczos algorithm.
!
!  NOTE: Currently, this code only works for a single outer-loop.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,10)
!
!  Estimate first the trace of the posterior analysis error
!  covariance matrix since the evolved and convolved Lanczos
!  vectors stored in the Hessian NetCDF file will be destroyed
!  later.
!
      Ltrace=.TRUE.

      TRACE_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0

        TRACE_ILOOP : DO my_inner=1,NpostI
          inner=my_inner
!
!  Initialize the tangent linear variables with a random vector
!  comprised of +1 and -1 elements randomly chosen.
!
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL random_ic (ng, tile, iTLM, inner, outer,             &
     &                        Lold(ng), Ltrace)
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!  Apply horizontal convolution.
!
          CALL convolve (driver, Lini, Lold, Lnew)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Compute Lanczos vector and eigenvectors of the posterior analysis
!  error covariance matrix.
!
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL posterior (ng, tile, iTLM, inner, outer, Ltrace)
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO

        END DO TRACE_ILOOP

      END DO TRACE_OLOOP
!
!  Estimate posterior analysis error covariance matrix.
!
      Ltrace=.FALSE.

      POST_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0
!
!  The Lanczos algorithm requires to save all the Lanczos vectors.
!  They are used to compute the posterior EOFs.
!
        DO ng=1,Ngrids
          ADM(ng)%Rindex=0
          Fcount=ADM(ng)%load
          ADM(ng)%Nrec(Fcount)=0
        END DO

        POST_ILOOP : DO my_inner=0,NpostI
          inner=my_inner
!
!  Read first record of ITL file and apply convolutions.
!
!  NOTE: If inner=0, we would like to use a random starting vector.
!        For now we can use what ever is in record 1.
!
          IF (my_inner.ne.0) THEN
            DO ng=1,Ngrids
              Rec=1
              CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec, Lold(ng))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
          ELSE
            DO ng=1,Ngrids
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL random_ic (ng, tile, iTLM, inner, outer,           &
     &                          Lold(ng), Ltrace)
              END DO
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
          END IF
!
!  Apply horizontal convolution.
!
          CALL convolve (driver, Lini, Lold, Lnew)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Compute Lanczos vector and eigenvectors of the posterior analysis
!  error covariance matrix.
!
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL posterior (ng, tile, iTLM, inner, outer, Ltrace)
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!   Write the Lanczos vectors of the posterior error covariance
!   to the adjoint NetCDF file.
!
          DO ng=1,Ngrids
#   if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            Lfout(ng)=Lnew(ng)
#   endif
#   ifdef ADJUST_BOUNDARY
            Lbout(ng)=Lnew(ng)
#   endif
            kstp(ng)=Lnew(ng)
#   ifdef SOLVE3D
            nstp(ng)=Lnew(ng)
#   endif
            LwrtState2d(ng)=.TRUE.
            CALL ad_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtState2d(ng)=.FALSE.
          END DO
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing adjustments for next inner loop into
!  ITL(ng)%name (record Rec1).
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, Lnew(ng), Rec1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO

        END DO POST_ILOOP

      END DO POST_OLOOP

#  endif /* POSTERIOR_EOFS */
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        HIS(ng)%ncid=-1
      END DO

#  ifdef POSTERIOR_EOFS
!
 10   FORMAT (/,' <<<< Posterior Analysis Error Covariance Matrix',     &
     &          ' Estimation >>>>',/)
#  endif
!
      RETURN
      END SUBROUTINE posterior_error
# endif

#endif
      END MODULE r4dvar_mod
