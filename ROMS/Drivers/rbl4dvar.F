#include "cppdefs.h"
      MODULE rbl4dvar_mod

#ifdef RBL4DVAR
!
!git $Id$
!svn $Id: rbl4dvar.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module splits the RBL4D-Var data assimilation algorithm into   !
!  its logical components routines.                                    !
!                                                                      !
!    background_initialize:                                            !
!                                                                      !
!      Initializes the nonlinear model kernel. It is separated from    !
!      the 'background' phase to allow ESM coupling.                   !
!                                                                      !
!    background:                                                       !
!                                                                      !
!      Timesteps the nonlinear model to compute the basic state Xb(t)  !
!      used to linearize the tangent linear and adjoint models. It     !
!      interpolates the nonlinear model trajectory to the observations !
!      locations.                                                      !
!                                                                      !
!    increment:                                                        !
!                                                                      !
!      Minimizes of the 4D-Var cost function over Ninner inner loops   !
!      iterations to compute the data assimilation increment, dXa. It  !
!      also computes the nonlinear model initial conditions by adding  !
!      the increment to the background, Xa = Xb(t=0) + dXa. The new    !
!      NLM initial conditions is written inot INI(ng)%name NetCDF file.!
!                                                                      !
!    analysis_initialize:                                              !
!                                                                      !
!      Initializes the nonlinear model with the 4D-Var analysis, Xa.   !
!      It is separated from 'analysis' phase to allow ESM coupling.    !
!                                                                      !
!    analysis:                                                         !
!                                                                      !
!      Timesteps the nonlinear model to compute new state trajectory   !
!      that becomes the basic state for the next outer loop. It also   !
!      interpolates the solution at the observation locations.         !
!                                                                      !
!    prior_error:                                                      !
!                                                                      !
!      Processes the background and model prior error covariance and   !
!      its normalization coefficients.                                 !
!                                                                      !
!    posterior_error:                                                  !
!                                                                      !
!      Computes the posterior error analysis.                          !
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
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
      USE mod_scalars
      USE mod_stepping
!
# ifdef ADJUST_BOUNDARY
      USE mod_boundary,      ONLY : initialize_boundary
# endif
      USE mod_forces,        ONLY : initialize_forces
      USE mod_mixing,        ONLY : initialize_mixing
      USE mod_ocean,         ONLY : initialize_ocean
!
# ifdef SPLIT_4DVAR
      USE ad_def_his_mod,    ONLY : ad_def_his
# endif
      USE ad_wrt_his_mod,    ONLY : ad_wrt_his
# ifndef RPCG
      USE congrad_mod,       ONLY : congrad
#  ifdef SPLIT_4DVAR
      USE congrad_mod,       ONLY : cg_read_congrad
#  endif
# endif
# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F || \
     defined POSTERIOR_EOFS
      USE convolve_mod,      ONLY : convolve
# endif
      USE convolve_mod,      ONLY : error_covariance
# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      USE def_error_mod,     ONLY : def_error
# endif
# if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
     defined POSTERIOR_ERROR_F
      USE def_hessian_mod,   ONLY : def_hessian
# endif
# ifdef SPLIT_4DVAR
      USE def_ini_mod,       ONLY : def_ini
# endif
      USE def_impulse_mod,   ONLY : def_impulse
      USE def_mod_mod,       ONLY : def_mod
      USE def_norm_mod,      ONLY : def_norm
      USE get_state_mod,     ONLY : get_state
# ifdef POSTERIOR_ERROR_F
      USE ini_adjust_mod,    ONLY : load_TLtoAD
# endif
      USE normalization_mod, ONLY : normalization
# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F || \
     defined POSTERIOR_EOFS
      USE posterior_mod,     ONLY : posterior
# endif
#  if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      USE posterior_var_mod, ONLY : posterior_var
#  endif
# ifdef RPCG
      USE rpcg_lanczos_mod,  ONLY : rpcg_lanczos
#  ifdef SPLIT_4DVAR
      USE rpcg_lanczos_mod,  ONLY : cg_read_rpcg
#  endif
# endif
      USE strings_mod,       ONLY : FoundError
      USE tl_def_ini_mod,    ONLY : tl_def_ini
      USE tl_wrt_ini_mod,    ONLY : tl_wrt_ini
# ifdef RPCG
      USE wrt_aug_imp_mod,   ONLY : wrt_aug_imp
# endif
# if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      USE wrt_error_mod,     ONLY : wrt_error
# endif
# ifdef POSTERIOR_ERROR_F
      USE wrt_hessian_mod,   ONLY : wrt_hessian
# endif
      USE wrt_impulse_mod,   ONLY : wrt_impulse
      USE wrt_ini_mod,       ONLY : wrt_ini
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
      USE wrt_ini_mod,       ONLY : wrt_frc_AD
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod,  ONLY : balance_ref, biconj
# endif
!
      implicit none
!
      PUBLIC  :: background_initialize
      PUBLIC  :: background
      PUBLIC  :: increment
      PUBLIC  :: analysis_initialize
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
      integer :: Lbck = 2
      integer :: Lini = 1
      integer :: LTLM1 = 1
      integer :: LTLM2 = 2
      integer :: Rec1 = 1
      integer :: Rec2 = 2
      integer :: Rec3 = 3
      integer :: Rec4 = 4
      integer :: Rec5 = 5
!
      CONTAINS
!
      SUBROUTINE background_initialize (my_outer)
!
!=======================================================================
!                                                                      !
!  This routine initializes ROMS nonlinear model trajectory Xb_n-1(0)  !
!  for each (n-1) outer loops. It is separated from the 'background'   !
!  4D-Var phase in ESM coupling applications that use generic methods  !
!  for 'initialize', 'run', and 'finalize'.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
!  Local variable declarations.
!
      integer :: lstr, ng, tile
      integer :: Fcount, Tindex
# ifdef PROFILE
      integer :: thread
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", background_initialize"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize nonlinear model kernel.
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
!  Clear TLM mixing arrays.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_mixing (ng, tile, iTLM)
        END DO
      END DO
!
!  Get nonlinear model initial conditions.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
# ifdef FORWARD_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
# ifdef FRC_FILES
        LreadFRC(ng)=.TRUE.
# endif
        LreadFWD(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%load
        RST(ng)%Nrec(Fcount)=0
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lbck" of INI(ng)%name NetCDF file. The
!  record "Lbck" becomes the background state and "Lini" becomes
!  the current nonlinear initial conditions.
!
      Tindex=1
      DO ng=1,Ngrids
        INI(ng)%Rindex=1
        Fcount=INI(ng)%load
        INI(ng)%Nrec(Fcount)=1
# ifdef DISTRIBUTE
        CALL wrt_ini (ng, MyRank, Tindex, Lbck)
# else
        CALL wrt_ini (ng, -1, Tindex, Lbck)
# endif
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
!  Set the nonlinear model to output the quicksave history file as a
!  function of the outer loop. It may be used as the basic state
!  trajectory for the surface fluxes (wind stress, shortwave, heat
!  flux, and E-P) because they can be saved frequently to resolve
!  the daily cycle while avoiding large files.
!
      DO ng=1,Ngrids
        LdefQCK(ng)=.TRUE.
        LwrtQCK(ng)=.TRUE.
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), my_outer
        lstr=LEN_TRIM(QCK(ng)%name)
        QCK(ng)%base=QCK(ng)%name(1:lstr-3)
      END DO
!
!  Define output 4D-Var NetCDF file (DAV struc) containing all
!  processed data at observation locations.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!
      RETURN
      END SUBROUTINE background_initialize
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
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoneStepping
!
      integer :: i, lstr, ng, tile
      integer :: Fcount, Tindex
# ifdef PROFILE
      integer :: thread
# endif
# if defined MODEL_COUPLING && !defined MCT_LIB
      integer :: NstrStep, NendStep, extra
!
      real(dp) :: ENDtime, NEXTtime
# endif
!
      character (len=20) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", background"
!
      SourceFile=MyFile

# if !(defined MODEL_COUPLING && defined ESMF_LIB)
!
!-----------------------------------------------------------------------
!  Set nonlinear model initial conditions.
!-----------------------------------------------------------------------
!
      CALL background_initialize (my_outer)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Run nonlinear model and compute background state trajectory,
!  Xb_n-1(t), and the background values at the observation points
!  and times. It processes and writes the observations accept/reject
!  flag (ObsScale) once to allow background quality control, if any.
# if defined MODEL_COUPLING && !defined MCT_LIB
!  Since the ROMS kernel has a delayed output and line diagnostics by
!  one timestep, subtact an extra value to the report of starting and
!  ending timestep for clarity. Usually, the model coupling interval
!  is of the same size as ROMS timestep.
# endif
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
!  Initialize various parameters and switches.
!
      MyRunInterval=RunInterval
!
      DO ng=1,Ngrids
# ifdef AVERAGES
        LdefAVG(ng)=.TRUE.
        LwrtAVG(ng)=.TRUE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), my_outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
# ifdef DIAGNOSTICS
        LdefDIA(ng)=.TRUE.
        LwrtDIA(ng)=.TRUE.
        WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%head), my_outer
        lstr=LEN_TRIM(DIA(ng)%name)
        DIA(ng)%base=DIA(ng)%name(1:lstr-3)
# endif
        wrtMisfit(ng)=.TRUE.
        wrtObsScale(ng)=.TRUE.
        SporadicImpulse(ng)=.FALSE.
        FrequentImpulse(ng)=.FALSE.
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
!  Run nonlinear model kernel.
!
# ifdef SOLVE3D
      CALL main3d (RunInterval)
# else
      CALL main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  If completed stepping, reset switches and write out cost function.
# if defined MODEL_COUPLING && !defined MCT_LIB
!  In coupled applications, RunInterval is much less than ntimes*dt,
!  so we need to wait until the last coupling interval is finished.
!  Otherwise, the control switches will be turned off prematurely.
# endif
!-----------------------------------------------------------------------
!
# if defined MODEL_COUPLING && !defined MCT_LIB
      IF (NendStep.eq.ntend(1)) THEN
        DoneStepping=.TRUE.
      ELSE
        DoneStepping=.FALSE.
      END IF
# else
      DoneStepping=.TRUE.
# endif
!
      IF (DoneStepping) THEN
        DO ng=1,Ngrids
          LdefQCK(ng)=.FALSE.
          LwrtQCK(ng)=.FALSE.
# ifdef AVERAGES
          LdefAVG(ng)=.FALSE.
          LwrtAVG(ng)=.FALSE.
# endif
# ifdef DIAGNOSTICS
          LdefDIA(ng)=.FALSE.
          LwrtDIA(ng)=.FALSE.
# endif
          wrtNLmod(ng)=.FALSE.
          wrtObsScale(ng)=.FALSE.
        END DO
!
!  Report data penalty function.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            DO i=0,NobsVar(ng)
              IF (i.eq.0) THEN
                string='Total'
              ELSE
                string=ObsName(i)
              END IF
              IF (FOURDVAR(ng)%NLPenalty(i).ne.0.0_r8) THEN
                WRITE (stdout,30) my_outer, 0, 'NLM',                   &
     &                            FOURDVAR(ng)%NLPenalty(i),            &
     &                            TRIM(string)
              END IF
            END DO
          END IF
!
!  Write out initial data penalty function to NetCDF file.
!
          SourceFile=MyFile
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,             &
     &                              'NL_iDataPenalty',                  &
     &                              FOURDVAR(ng)%NLPenalty(0:),         &
     &                              (/1/), (/NobsVar(ng)+1/),           &
     &                              ncid = DAV(ng)%ncid)

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iNLM, DAV(ng)%name,         &
     &                                  'NL_iDataPenalty',              &
     &                                  FOURDVAR(ng)%NLPenalty(0:),     &
     &                                  (/1/), (/NobsVar(ng)+1/),       &
     &                                  pioFile = DAV(ng)%pioFile)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Clean penalty array before next run of NL model.
!
          FOURDVAR(ng)%NLPenalty=0.0_r8
        END DO
      END IF

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
 30   FORMAT (' (',i3.3,',',i3.3,'): ',a,' data penalty, Jdata = ',     &
     &        1p,e17.10,0p,t68,a)

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
!=======================================================================

# ifdef RPCG
!
      USE comp_Jb0_mod,   ONLY : aug_oper, comp_Jb0
      USE ini_adjust_mod, ONLY : ini_adjust
      USE sum_grad_mod,   ONLY : sum_grad
      USE sum_imp_mod,    ONLY : sum_imp
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
      logical :: Lcgini, Linner, Lposterior
!
      integer :: i, ifile, lstr, my_inner, ng, tile
      integer :: Fcount, InpRec
# ifdef RPCG
      integer :: ADrec, nLAST
      integer :: irec, jrec, jrec1, jrec2
      integer :: LiNL
# endif
# ifdef PROFILE
      integer :: thread
# endif
!
      integer, dimension(Ngrids) :: Nrec
# ifdef RPCG
      integer, dimension(Ngrids) :: nADrec
# endif
# ifdef SPLIT_4DVAR
!
      real(dp) :: stime
# endif
!
      character (len=8 ) :: driver = 'rbl4dvar'
      character (len=10) :: suffix

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", increment"
!
      SourceFile=MyFile
!
!=======================================================================
!  Compute 4D-Var increment.
!=======================================================================

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
#  ifdef RPCG
          CALL cg_read_rpcg (ng, iTLM, my_outer)
#  else
          CALL cg_read_congrad (ng, iTLM, my_outer)
#  endif
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
      IF (my_outer.eq.1) THEN        ! initial values from "background"
        ng=1
        SELECT CASE (DAV(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,               &
     &                            Vname(1,idNLmi), NLmodVal,            &
     &                            ncid = DAV(ng)%ncid)

#  if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,           &
     &                                Vname(1,idNLmi), NLmodVal,        &
     &                                pioFile = DAV(ng)%pioFile)

#  endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      ng=1
      SELECT CASE (DAV(ng)%IOtype)
        CASE (io_nf90)
          CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                 &
     &                          Vname(1,idObsS), ObsScale,              &
     &                          ncid = DAV(ng)%ncid)

#  if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              Vname(1,idObsS), ObsScale,          &
     &                              pioFile = DAV(ng)%pioFile)
#  endif
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  In the split 4D-Var alorithm, we also need to read in ObsVal
!  ("obs_value") and ObsErr ("obs_error") from the observation file.
!  They are used in the initialization of the conjugate gradient
!  algorithm.
!
      ng=1
      SELECT CASE (OBS(ng)%IOtype)
        CASE (io_nf90)
          CALL netcdf_get_fvar (ng, iTLM, OBS(ng)%name,                 &
     &                          Vname(1,idOval), ObsVal)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_get_fvar (ng, iTLM, OBS(ng)%name,                 &
     &                          Vname(1,idOerr), ObsErr)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_fvar (ng, iTLM, OBS(ng)%name,             &
     &                              Vname(1,idOval), ObsVal)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL pio_netcdf_get_fvar (ng, iTLM, OBS(ng)%name,             &
     &                              Vname(1,idOerr), ObsErr)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
      END SELECT
!
!  Set nonlinear output history file to be used as the basic state
!  trajectory. The 4D-Var increment phase is computed by a different
!  executable and needs to know some of the HIS structure information.
!
      DO ng=1,Ngrids
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), my_outer-1
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
!  Set the nonlinear model quicksave-history file as the basic state for
!  the surface fluxes computed in "bulk_flux", which may be available at
!  more frequent intervals while avoiding large files. Since the 4D-Var
!  increment phase is calculated by a different executable and needs to
!  know some of the QCK structure information.
!
      DO ng=1,Ngrids
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), my_outer-1
        lstr=LEN_TRIM(QCK(ng)%name)
        QCK(ng)%base=QCK(ng)%name(1:lstr-3)
        IF (QCK(ng)%Nfiles.gt.1) THEN
          DO ifile=1,QCK(ng)%Nfiles
            WRITE (suffix,"('_',i4.4,'.nc')") ifile
            QCK(ng)%files(ifile)=TRIM(QCK(ng)%base)//TRIM(suffix)
          END DO
          QCK(ng)%name=TRIM(QCK(ng)%files(1))
        ELSE
          QCK(ng)%files(1)=TRIM(QCK(ng)%name)
        END IF
      END DO
!
!  Read in 4D-Var starting time (sec) from nonlinear trajectory.
!  Initialize "tday" which are needed to write the correct time in
!  the ITL NetCDF file.  It is alse need for boundary and surface
!  forcing adjustments, if any.
!
      InpRec=1
      DO ng=1,Ngrids
        SELECT CASE (HIS(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, iTLM, HIS(ng)%name,               &
     &                            Vname(1,idtime), stime,               &
     &                            start = (/InpRec/),                   &
     &                            total = (/1/))

#  if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_fvar (ng, iTLM, HIS(ng)%name,           &
     &                                Vname(1,idtime), stime,           &
     &                                start = (/InpRec/),               &
     &                                total = (/1/))
#  endif
        END SELECT
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
# endif
!
!-----------------------------------------------------------------------
!  Set few variables needed in the "increment" phase.
!-----------------------------------------------------------------------
!
!  Set structure for the nonlinear forward trajectory to be processed
!  by the tangent linear and adjoint models. Also, set switches to
!  process the FWD structure in routine "check_multifile". Notice that
!  it is possible to split solution into multiple NetCDF files to reduce
!  their size.
!
      CALL edit_multifile ('HIS2FWD')
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
        LreadQCK(ng)=.FALSE.
      END DO
# endif

# ifdef RPCG
!
!  If weak constraint, set the number of TLF records used in the
!  augmented solver.
!
      DO ng=1,Ngrids
        IF (nADJ(ng).lt.ntimes(ng)) THEN
          nADrec(ng)=2*(1+ntimes(ng)/nADJ(ng))
        ELSE
          nADrec(ng)=0
        END IF
      END DO
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

# ifdef RPCG
!
!  If outer=1, Initialize all records of the ITL file to zero. The TLM
!  state variables are zero when outer=1.
!
        DO ng=1,Ngrids
          tdays(ng)=INItime(ng)*sec2day
#  ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rec1, Rec1)
#  else
          CALL tl_wrt_ini (ng, -1, Rec1, Rec1)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rec1, Rec2)
#  else
          CALL tl_wrt_ini (ng, -1, Rec1, Rec2)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rec1, Rec3)
#  else
          CALL tl_wrt_ini (ng, -1, Rec1, Rec3)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rec1, Rec4)
#  else
          CALL tl_wrt_ini (ng, -1, Rec1, Rec4)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rec1, Rec5)
#  else
          CALL tl_wrt_ini (ng, -1, Rec1, Rec5)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (nADJ(ng).lt.ntimes(ng)) THEN
            nLAST=Rec5
            DO irec=1,nADrec(ng)
#  ifdef DISTRIBUTE
              CALL tl_wrt_ini (ng, MyRank, Rec1, nLAST+irec)
#  else
              CALL tl_wrt_ini (ng, -1, Rec1, nLAST+irec)
#  endif
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
          END IF
        END DO
# endif
!
!  If outer=1, define impulse forcing NetCDF file.
!
        DO ng=1,Ngrids
          LdefTLF(ng)=.TRUE.
          CALL def_impulse (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

# if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
     defined POSTERIOR_ERROR_F
!
!  If outer=1, define output Hessian NetCDF file that will eventually
!  contain the intermediate posterior analysis error covariance matrix
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
!  If outer=1, define output initial or final full posterior error
!  covariance (diagonal) matrix NetCDF.
!
        DO ng=1,Ngrids
          LdefERR (ng)=.TRUE.
          CALL def_error (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif

      END IF CHECK_OUTER1

# ifdef RPCG
!
!-----------------------------------------------------------------------
!  If outer>1, integrate the tangent linear model to evaluate
!  G[x(k)-x(k-1)], which is written to the DAV NetCDF in the
!  variable "TLmodel_value". It is needed in routine "rpcg_lanczos".
!  In the old unsplit algorithm, the TLM was run in the "analysis"
!  phase. It is done here to allow lower resolution in the inner
!  loops in the future.
!-----------------------------------------------------------------------
!
      CHECK_OUTER2 : IF (my_outer.gt.1) THEN
!
!  Clear tangent linear forcing arrays. This is very important since
!  these arrays are non-zero and must be zero when running the tangent
!  linear model. We also need to clear the nonlinear state arrays to
!  make sure that the RHS terms rzeta, rubar, rvbar, ru, and rv are
!  zero.  Otherwise, those array will have the last computed values
!  when running the nonlinear model if not processing the forward
!  trajectory RHS terms.  This needs to be done to get identical
!  solutions with the split schemes.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iNLM)
            CALL initialize_forces (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
        END DO
!
!  Initialize tangent linear model from initial impulse which is now
!  stored in file ITL(ng)%name.
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
!  Initialize tangent linear model from ITLname, record Rec3.
!  The sum of the initial condition increments is in record 3.
!
        DO ng=1,Ngrids
          ITL(ng)%Rindex=Rec3
          CALL tl_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Run tangent linear model.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'TL', ng, my_outer, 0,                    &
     &                        ntstart(ng), ntend(ng)
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
      END IF CHECK_OUTER2
# endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  4D-Var inner loops.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero and must
!  be zero when running the tangent linear model.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iNLM)
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
        END DO

# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
        IF (balance(isFsur)) THEN
          DO ng=1,Ngrids
            CALL get_state (ng, iNLM, 2, INI(ng), Lini, Lini)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL balance_ref (ng, tile, Lini)
              CALL biconj (ng, tile, iNLM, Lini)
            END DO
            wrtZetaRef(ng)=.TRUE.
          END DO
        END IF
# endif
!
!  Inner loops iteration.
!
      INNER_LOOP : DO my_inner=0,Ninner
        inner=my_inner

# ifdef RPCG
        IF (inner.ne.Ninner) THEN
          Linner=.TRUE.
        ELSE
          Linner=.FALSE.
        END IF
!
!  Retrieve TLmodVal and NLmodVal when inner=0 and outer>1 for use as
!  Hbk and BCKmodVal, respectively.
!
        IF ((inner.eq.0).and.(my_outer.gt.1)) THEN
          DO ng=1,Ngrids
            SELECT CASE (DAV(ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,           &
     &                                'TLmodel_value', TLmodVal)
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) RETURN

                CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,           &
     &                                'NLmodel_value', NLmodVal,        &
     &                                ncid = DAV(ng)%ncid,              &
     &                                start = (/1,my_outer-1/),         &
     &                                total = (/Ndatum(ng),1/))
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) RETURN

#  if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,       &
     &                                    'TLmodel_value', TLmodVal)
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) RETURN

                CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,       &
     &                                    'NLmodel_value', NLmodVal,    &
     &                                    pioFile = DAV(ng)%pioFile,    &
     &                                    start = (/1,my_outer-1/),     &
     &                                    total = (/Ndatum(ng),1/))
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) RETURN
#  endif
            END SELECT
          END DO
        END IF
!
        IF (inner.eq.0) THEN
          Lcgini=.TRUE.
        END IF
        DO ng=1,Ngrids
          CALL rpcg_lanczos (ng, iTLM, my_outer, inner, Ninner, Lcgini)
        END DO
# else
!
!  Initialize conjugate gradient algorithm depending on hot start or
!  outer loop index.
!
        IF (inner.eq.0) THEN
          Lcgini=.TRUE.
          DO ng=1,Ngrids
            CALL congrad (ng, iRPM, my_outer, inner, Ninner, Lcgini)
          END DO
        END IF
!
!  If initialization step, skip the inner-loop computations.
!
        Linner=.FALSE.
        IF ((inner.ne.0).or.(Nrun.ne.1)) THEN
          IF (((inner.eq.0).and.LhotStart).or.(inner.ne.0)) THEN
            Linner=.TRUE.
          END IF
        END IF
# endif
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
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file only once to avoid opening too many files.
!
          DO ng=1,Ngrids
# ifdef WEAK_NOINTERP
            WRTforce(ng)=.FALSE.
# else
            WRTforce(ng)=.TRUE.
# endif
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
# ifdef SOLVE3D
          CALL ad_main3d (RunInterval)
# else
          CALL ad_main2d (RunInterval)
# endif
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

# ifndef WEAK_NOINTERP
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
          DO ng=1,Ngrids
#  ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
#  else
            CALL ad_wrt_his (ng, -1)
#  endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
# endif
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
          DO ng=1,Ngrids
            WRTforce(ng)=.FALSE.
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!  Convolve adjoint trajectory with error covariances.
!
# ifdef POSTERIOR_ERROR_I
          Lposterior=.TRUE.
# else
          Lposterior=.FALSE.
# endif
# ifdef RPCG
!
!  Set the flag that controls the augmentation of the model error
!  forcing terms. This is ONLY done in the outer-loop so
!  LaugWeak=.FALSE. here.
!
          LaugWeak=.FALSE.
# endif
          CALL error_covariance (iTLM, driver, my_outer, inner,         &
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
            WRITE (stdout,30) my_outer, inner
          END IF
          DO ng=1,Ngrids
            TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
            CALL wrt_impulse (ng, MyRank, iADM, ADM(ng)%name)
# else
            CALL wrt_impulse (ng, -1, iADM, ADM(ng)%name)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!-----------------------------------------------------------------------
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!-----------------------------------------------------------------------
!
!  Initialize tangent linear model from initial impulse which is now
!  stored in file ITL(ng)%name.
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
          IF ((my_outer.eq.1).and.(inner.eq.1)) THEN
            DO ng=1,Ngrids
              wrtMisfit(ng)=.TRUE.
            END DO
          END IF
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (H M B M' H')_n * PSI at observation
!  points which are used in the conjugate gradient algorithm.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ng, my_outer, my_inner,           &
     &                          ntstart(ng), ntend(ng)
            END IF
          END DO
!
# ifdef SOLVE3D
          CALL tl_main3d (RunInterval)
# else
          CALL tl_main2d (RunInterval)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO ng=1,Ngrids
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.FALSE.
          END DO

# ifdef POSTERIOR_ERROR_F
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
          IF (inner.ne.0) THEN
            DO ng=1,Ngrids
#  ifdef DISTRIBUTE
              CALL wrt_hessian (ng, tile, Lold(ng), Lold(ng))
#  else
              CALL wrt_hessian (ng, -1, Lold(ng), Lold(ng))
#  endif
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
          END IF
# endif
!
          Nrun=Nrun+1

# ifndef RPCG
!
!  Use conjugate gradient algorithm to find a better approximation
!  PSI to coefficients Beta_n.
!
          DO ng=1,Ngrids
            Lcgini=.FALSE.
            CALL congrad (ng, iTLM, my_outer, inner, Ninner, Lcgini)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
# else
          Lcgini=.FALSE.
# endif

        END IF INNER_COMPUTE

      END DO INNER_LOOP
!
!-----------------------------------------------------------------------
!  Once the w_n, have been approximated with sufficient accuracy,
!  compute estimates of Lambda_n and Xhat_n by carrying out one
!  backward intergration of the adjoint model and one forward
!  itegration of the nonlinear model.
!-----------------------------------------------------------------------
!
!  Initialize the adjoint model always from rest.
!
      DO ng=1,Ngrids
        CALL ad_initial (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
      DO ng=1,Ngrids
# ifdef WEAK_NOINTERP
        WRTforce(ng)=.FALSE.
# else
        WRTforce(ng)=.TRUE.
# endif
        IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
        Fcount=ADM(ng)%load
        ADM(ng)%Nrec(Fcount)=0
        ADM(ng)%Rindex=0
      END DO
!
!  Time-step adjoint model backwards forced with estimated coefficients,
!  Beta_n.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,20) 'AD', ng, my_outer, Ninner,                 &
     &                      ntstart(ng), ntend(ng)
        END IF
      END DO
!
# ifdef SOLVE3D
      CALL ad_main3d (RunInterval)
# else
      CALL ad_main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifndef WEAK_NOINTERP
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL ad_wrt_his (ng, MyRank)
#  else
        CALL ad_wrt_his (ng, -1)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
# endif
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
      DO ng=1,Ngrids
        WRTforce(ng)=.FALSE.
# ifdef DISTRIBUTE
        CALL ad_wrt_his (ng, MyRank)
# else
        CALL ad_wrt_his (ng, -1)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

# ifdef RPCG
!
!  Get number of records in adjoint NetCDF.
!  We need to do this here since ADM(ng)%Nrec is reset to zero in
!  error_covariance.
!
      DO ng=1,Ngrids
        Fcount=ADM(ng)%load
        Nrec(ng)=ADM(ng)%Nrec(Fcount)
      END DO
# endif
!
!  Convolve adjoint trajectory with error covariances.
!
      Lposterior=.FALSE.

# ifdef RPCG
!
!  Set the flag that controls the augmentation of the model error
!  forcing terms. This is ONLY done in the outer-loop so
!  LaugWeak=.TRUE. here.
!
      LaugWeak=.TRUE.
      LiNL=my_outer+1
      CALL error_covariance (iNLM, driver, my_outer, inner,             &
     &                       LiNL, Lini, Lold, Lnew,                    &
     &                       Rec1, Rec2, Lposterior)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# else
      CALL error_covariance (iNLM, driver, my_outer, inner,             &
     &                       Lbck, Lini, Lold, Lnew,                    &
     &                       Rec1, Rec2, Lposterior)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef RPCG
!
!  Augmented solver:
!
!  NOTES: The ITL file contains 5 records -
!         Rec2 = the new TL initial condition
!         Rec3 = the sum of the TL initial conditions
!         Rec4 = B^-1(xb-xk)=sum_j^k-1 G_j lambda_j
!         Rec5 = the augmented correction to Rec2.
!         Rec5+1 to Rec5+nADrec/2 = sum of the TLF forcing increments
!         Rec5+nADrec/2+1 to Rec5+nADrec = B^-1 of the sum of the TLF
!                                         forcing increments.
!  Reset the flag LaugWeak flag.
!
      LaugWeak=.FALSE.

# else
!
!  Write out nonlinear model initial conditions into INIname, record
!  INI(ng)%Rindex.
!
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL wrt_ini (ng, MyRank, Lnew(ng))
#  else
        CALL wrt_ini (ng, -1, Lnew(ng))
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
      defined ADJUST_BOUNDARY
#   ifdef DISTRIBUTE
        CALL wrt_frc_AD (ng, MyRank, Lold(ng), INI(ng)%Rindex)
#   else
        CALL wrt_frc_AD (ng, -1, Lold(ng), INI(ng)%Rindex)
#   endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
      END DO
# endif

# ifdef RPCG
!
! Compute the augmented correction to the adjoint propagator.
! We need to use sum (x(k)-x(k-1)) before it is updated. This is
! in record 3 of the ITL file.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iTLM, 4, ITL(ng), Rec3, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL aug_oper (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
! Save this augmented piece in record 5 of the ITL file.
!
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL tl_wrt_ini (ng, MyRank, LTLM2, Rec5)
#  else
        CALL tl_wrt_ini (ng, -1, LTLM2, Rec5)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
! Complete the computation of the TL initial condition by adding the
! contribution from the augmented adjoint propagator.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iTLM, 4, ITL(ng), Rec2, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        CALL get_state (ng, iTLM, 4, ITL(ng), Rec5, LTLM2)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL sum_grad (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
! Write the final TL increment to Rec2 of the ITL file.
!
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL tl_wrt_ini (ng, MyRank, LTLM2, Rec2)
#  else
        CALL tl_wrt_ini (ng, -1, LTLM2, Rec2)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
! Now update the non-linear model initial conditions.
!
      LiNL=my_outer+1
      DO ng=1,Ngrids
        CALL get_state (ng, iNLM, 9, INI(ng), LiNL, Lnew(ng))
        CALL get_state (ng, iADM, 4, ITL(ng), Rec2, LTLM1)
      END DO
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_adjust (ng, tile, LTLM1, Lnew(ng))
        END DO
      END DO
!
!  Write out nonlinear model initial conditions into INIname, record
!  INI(ng)%Rindex.
!
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL wrt_ini (ng, MyRank, Lnew(ng))
#  else
        CALL wrt_ini (ng, -1, Lnew(ng))
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
      defined ADJUST_BOUNDARY
#   ifdef DISTRIBUTE
        CALL wrt_frc_AD (ng, MyRank, LTLM1, INI(ng)%Rindex)
#   else
        CALL wrt_frc_AD (ng, -1, LTLM1, INI(ng)%Rindex)
#   endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
      END DO
!
! Compute the new B^-1(x(k)-x(k-1)) term.
! Gather the final adjoint solutions increments, sum and
! save in record 4 of the ITL file. Use the tl arrays as temporary
! storage.
!
! First add the augmented piece which is computed from the previous
! sum.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iTLM, 4, ITL(ng), Rec4, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL aug_oper (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL sum_grad (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
      DO ng=1,Ngrids
        ADrec=Nrec(ng)
        CALL get_state (ng, iTLM, 4, ADM(ng), ADrec, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL sum_grad (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
!  Write the current sum of adjoint solutions into record 4 of the ITL
!  file.
!
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL tl_wrt_ini (ng, MyRank, LTLM2, Rec4)
#  else
        CALL tl_wrt_ini (ng, -1, LTLM2, Rec4)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
# endif
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
      IF (Master) THEN
        WRITE (stdout,30) my_outer, inner
      END IF
      DO ng=1,Ngrids
        TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
        CALL wrt_impulse (ng, MyRank, iADM, ADM(ng)%name)
# else
        CALL wrt_impulse (ng, -1, iADM, ADM(ng)%name)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

# ifdef RPCG
!
!  Now Compute the augmented corrections to the weak constraint
!  forcing terms. The sums so far are in records 6 to
!  6+nADrec/2.
!
      DO ng=1,Ngrids
        DO i=1,nADrec(ng)/2
          irec=i
          jrec=Rec5+i
          CALL get_state (ng, iTLM, 4, ITL(ng), jrec, LTLM1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL aug_oper (ng, tile, LTLM1, LTLM2)
          END DO
!
!  Complete the computation of the TLF forcing term by adding the
!  contribution from the augmented adjoint propagator. Specify
!  the value of 7 for the model variable since this the special
!  case in get_state for reading the impulse forcing.
!
!!        CALL get_state (ng, iTLM, 4, TLF(ng), irec, LTLM1)      ! TEST
          CALL get_state (ng, 7, 4, TLF(ng), irec, LTLM1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO tile=first_tile(ng),last_tile(ng),+1
!!          CALL sum_grad (ng, tile, LTLM1, LTLM2)                ! TEST
            CALL sum_imp (ng, tile, LTLM2)
          END DO
!
!  Write the final forcing increment to he TLF file.
!  Note the original TLF file is overwritten at this point.
!
          TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
          CALL wrt_aug_imp (ng, MyRank, iTLM, LTLM2, i)
# else
          CALL wrt_aug_imp (ng, -1, iTLM, LTLM2, i)
# endif

        END DO
      END DO
!
!  Gather the increments from the final inner-loop and
!  save in record 3 of the ITL file. The current increment
!  is in record 2 and the sum so far is in record 3.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iTLM, 4, ITL(ng), Rec2, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        CALL get_state (ng, iTLM, 4, ITL(ng), Rec3, LTLM2)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL sum_grad (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
!  Write the current sum into record 3 of the ITL file.
!
      DO ng=1,Ngrids
# ifdef DISTRIBUTE
        CALL tl_wrt_ini (ng, MyRank, LTLM2, Rec3)
# else
        CALL tl_wrt_ini (ng, -1, LTLM2, Rec3)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
!
!  Write the current sum of the forcing and boundary increments
!  into the INI file into index INI(ng)%Rindex-1 (i.e. we are
!  overwriting the previous fields. Read the current sum into the
!  adjoint arrays since this is what wrt_frc_AD uses.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iADM, 4, ITL(ng), Rec3, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
      DO ng=1,Ngrids
#  ifdef DISTRIBUTE
        CALL wrt_frc_AD (ng, MyRank, LTLM1, INI(ng)%Rindex)
#  else
        CALL wrt_frc_AD (ng, -1, LTLM1, INI(ng)%Rindex)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
# endif
!
!  Gather the model error forcing increments and update the
!  sum in records 6 to 6+nADrec/2.
!
      DO ng=1,Ngrids
        DO i=1,nADrec(ng)/2
          irec=i
          jrec=Rec5+i
          CALL get_state (ng, iTLM, 4, ITL(ng), jrec, LTLM1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!!        CALL get_state (ng, iTLM, 4, TLF(ng)%, irec, LTLM2)     ! TEST
          CALL get_state (ng, 7, 4, TLF(ng), irec, LTLM2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO tile=first_tile(ng),last_tile(ng),+1
!!          CALL sum_grad (ng, tile, LTLM1, LTLM2)                ! TEST
            CALL sum_imp (ng, tile, LTLM1)
          END DO
!
!  Write the current sum into record jrec of the ITL file.
!
# ifdef DISTRIBUTE
!!        CALL tl_wrt_ini (ng, MyRank, LTLM2, jrec)               ! TEST
          CALL tl_wrt_ini (ng, MyRank, LTLM1, jrec)
# else
!!        CALL tl_wrt_ini (ng, -1, LTLM2, jrec)                   ! TEST
          CALL tl_wrt_ini (ng, -1, LTLM1, jrec)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END DO
!
!  Now compute the background cost function Jb0.
!  First compute the contribution from the increments in the
!  initial conditions, surface forcing, and boundary conditions.
!
      Jb0(my_outer)=0.0_r8
      DO ng=1,Ngrids
        CALL get_state (ng, iADM, 8, ITL(ng), Rec4, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        CALL get_state (ng, iTLM, 8, ITL(ng), Rec3, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL comp_Jb0 (ng, tile, iTLM, my_outer, LTLM1, LTLM1)
        END DO
      END DO
!
!  Now compute the contribution of model error terms to the
!  background cost function Jb0.
!
!  NOTE: I THINK WE NEED A BETTER WAY OF COMPUTING THE OBS ERROR
!  CONTRIBUTION TO Jb. THE FOLLOWING CALCULATION IS DANGEROUS BECAUSE
!  IT RELIES ON THE FACT THAT THE SURFACING FORCING AND OBC TL ARRAYS
!  READ FROM jrec2 ARE ZERO. THIS MEANS THAT ONLY THE MODEL STATE
!  CONTRIBUTES TO Jb FOR THE MODEL ERROR TERM, WHICH OF COURSE IS
!  WHAT SHOULD BE THE CASE. THE SURFACE FORCING AND OBC CONTRIBUTIONS
!  ARE COMPUTED IN THE PREVIOUS CALL TO COMP_JB0.
!
      DO ng=1,Ngrids
        DO irec=1,nADrec(ng)/2
          jrec1=Rec5+irec
          jrec2=Rec5+nADrec(ng)/2+irec
          CALL get_state (ng, iADM, 8, ITL(ng), jrec1, LTLM1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          CALL get_state (ng, iTLM, 8, ITL(ng), jrec2, LTLM1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL comp_Jb0 (ng, tile, iTLM, my_outer, LTLM1, LTLM1)
          END DO
        END DO
      END DO
!
!  Overwrite the TFL netcdf file with the sum of the model error
!  forcing increments - required for the next run of the NLM and
!  TLM.
!
      DO ng=1,Ngrids
        DO irec=1,nADrec(ng)/2
          jrec=Rec5+irec
          CALL get_state (ng, iTLM, 8, ITL(ng), jrec, LTLM1)
!
          TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
          CALL wrt_aug_imp (ng, MyRank, iTLM, LTLM1, irec)
# else
          CALL wrt_aug_imp (ng, -1, iTLM, LTLM1, irec)
# endif
        END DO
      END DO
# endif

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
 30   FORMAT (/,' Converting Convolved Adjoint Trajectory to',          &
     &          ' Impulses: Outer = ',i3.3,' Inner = ',i3.3,/)
!
      RETURN
      END SUBROUTINE increment
!
      SUBROUTINE analysis_initialize (my_outer)
!
!=======================================================================
!                                                                      !
!  This routine initializes the nonlinear kernel with the 4D-Var new   !
!  state estimate Xa = Xb + dXa.  It is separated from the 'analysis'  !
!  4D-Var phase in ESM coupling applications that use generic methods  !
!  for 'initialize', 'run', and 'finalize'.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
!  Local variable declarations.
!
      integer :: ifile, lstr, ng, tile
      integer :: Fcount
# ifdef PROFILE
      integer :: thread
# endif
!
      integer, dimension(Ngrids) :: indxSave
!
      character (len=10) :: suffix

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", analysis_initialize"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize nonlinear model kernel with new 4D-Var state estimate.
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
!
!  Clear tangent arrays and the nonlinear model mixing arrays
!  before running nonlinear model (important).
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iTLM)
          CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iTLM)
# endif
          CALL initialize_mixing (ng, tile, iNLM)
        END DO
      END DO

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
#  if defined MODEL_COUPLING && defined ESMF_LIB
      ERend=Nouter
#  endif
!
!  Open Nonlinear model initial conditions NetCDF file (INI struc) and
!  inquire about its variable IDs. In particular, set the current value
!  of INI(ng)%Rindex, which is needed in "initial".
!
      DO ng=1,Ngrids
        LdefINI(ng)=.FALSE.
        CALL def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
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
      SELECT CASE (DAV(ng)%IOtype)
        CASE (io_nf90)
          CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                 &
     &                          Vname(1,idObsS), ObsScale,              &
     &                          ncid = DAV(ng)%ncid)

#  if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              Vname(1,idObsS), ObsScale,          &
     &                              pioFile = DAV(ng)%pioFile)
#  endif
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  If outer=Nouter, read initial NLmodVal ("NLmodel_initial") and load
!  it in the scratch "NLincrement" array.
!
      IF (my_outer.eq.Nouter) THEN
        ng=1
        SELECT CASE (DAV(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,               &
     &                            Vname(1,idNLmi), NLincrement,         &
     &                            ncid = DAV(ng)%ncid)

#  if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,           &
     &                                Vname(1,idNLmi), NLincrement,     &
     &                                pioFile = DAV(ng)%pioFile)
#  endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef FORWARD_FLUXES
!
!  If not first NLM run, set BLK structure to the previous outer loop
!  quicksave trajectory.  There is a logic in "get_data" that reads
!  from BLK.
!  (HGA: This probably legacy code that it is no longer needed)
!
      DO ng=1,Ngrids
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), my_outer-1
        lstr=LEN_TRIM(QCK(ng)%name)
        QCK(ng)%base=QCK(ng)%name(1:lstr-3)
      END DO
!
      CALL edit_multifile ('QCK2BLK')
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Initialize nonlinear model.
!-----------------------------------------------------------------------
!
!  Set new basic state trajectory for next outer loop.  Notice that the
!  LreadFWD switch is turned off to suppress processing of the structure
!  when "check_multifile" during nonlinear model execution.
!
      DO ng=1,Ngrids
        idefHIS(ng)=-1
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
# ifdef FRC_FILES
        LreadFRC(ng)=.TRUE.
# endif
        LreadFWD(ng)=.FALSE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), my_outer
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
!  Set nonlinear output quicksave history file for the next outer loop.
!  Notice that the LreadBLK switch is turned off to suppress processing
!  of the structure when "check_multifile" during nonlinear model
!  execution.
!
      DO ng=1,Ngrids
        idefQCK(ng)=-1
        LdefQCK(ng)=.TRUE.
        LwrtQCK(ng)=.TRUE.
# ifdef FORWARD_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), my_outer
        lstr=LEN_TRIM(QCK(ng)%name)
        QCK(ng)%base=QCK(ng)%name(1:lstr-3)
        IF (QCK(ng)%Nfiles.gt.1) THEN
          DO ifile=1,QCK(ng)%Nfiles
            WRITE (suffix,"('_',i4.4,'.nc')") ifile
            QCK(ng)%files(ifile)=TRIM(QCK(ng)%base)//TRIM(suffix)
          END DO
          QCK(ng)%name=TRIM(QCK(ng)%files(1))
        ELSE
          QCK(ng)%files(1)=TRIM(QCK(ng)%name)
        END IF
      END DO
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
      DO ng=1,Ngrids
# ifdef SPLIT_4DVAR
        IF (nADJ(ng).lt.ntimes(ng)) THEN
          SELECT CASE (ADM(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_get_dim (ng, iNLM, ADM(ng)%name,              &
     &                             DimName = 'ocean_time',              &
     &                             DimSize = FrcRec(ng))
#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_get_dim (ng, iNLM, ADM(ng)%name,          &
     &                                 DimName = 'ocean_time',          &
     &                                 DimSize = FrcRec(ng))
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
        IF (FrcRec(ng).gt.3) THEN
          FrequentImpulse(ng)=.TRUE.
        END IF
      END DO
!
!  Initialize nonlinear model INI(ng)%name file, record outer+2.
!  Notice that NetCDF record index counter is saved because this
!  counter is used to write initial conditions.
!
      DO ng=1,Ngrids
        indxSave(ng)=INI(ng)%Rindex
        INI(ng)%Rindex=my_outer+2
        RST(ng)%Rindex=0
        Fcount=RST(ng)%load
        RST(ng)%Nrec(Fcount)=0
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      DO ng=1,Ngrids
        INI(ng)%Rindex=indxSave(ng)
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
 10   FORMAT (a,'_outer',i0,'.nc')
!
      RETURN
      END SUBROUTINE analysis_initialize
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
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoneStepping
!
      integer :: i, lstr, ng
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
      character (len=20) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", analysis"
!
      SourceFile=MyFile
!
# if !(defined MODEL_COUPLING && defined ESMF_LIB)
!
!-----------------------------------------------------------------------
!  Set nonlinear model initial conditions.
!-----------------------------------------------------------------------
!
      CALL analysis_initialize (my_outer)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Run nonlinear model and compute a "new estimate" of the state
!  trajectory, Xb(t)|n.
# if defined MODEL_COUPLING && !defined MCT_LIB
!  Since the ROMS kernel has a delayed output and line diagnostics by
!  one timestep, subtact an extra value to the report of starting and
!  ending timestep for clarity. Usually, the model coupling interval
!  is of the same size as ROMS timestep.
# endif
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
!
!  Initialize various parameters and switches.
!
      MyRunInterval=RunInterval
!
      DO ng=1,Ngrids
# ifdef AVERAGES
        idefAVG(ng)=-1
        LdefAVG(ng)=.TRUE.
        LwrtAVG(ng)=.TRUE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), my_outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
# ifdef DIAGNOSTICS
        idefDIA(ng)=-1
        LdefDIA(ng)=.TRUE.
        LwrtDIA(ng)=.TRUE.
        WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%head), my_outer
        lstr=LEN_TRIM(DIA(ng)%name)
        DIA(ng)%base=DIA(ng)%name(1:lstr-3)
# endif
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
        IF (Master) WRITE (stdout,20) 'NL', ng, my_outer, Ninner,       &
     &                                NstrStep, NendStep
# else
        IF (Master) WRITE (stdout,20) 'NL', ng, my_outer, Ninner,       &
     &                                ntstart(ng), ntend(ng)
# endif
      END DO
!
!  Run nonlinear model kernel.
!
# ifdef SOLVE3D
      CALL main3d (RunInterval)
# else
      CALL main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  If completed stepping, reset switches and write out cost function.
# if defined MODEL_COUPLING && !defined MCT_LIB
!  In coupled applications, RunInterval is much less than ntimes*dt,
!  so we need to wait until the last coupling interval is finished.
!  Otherwise, the control switches will be turned off prematurely.
# endif
!-----------------------------------------------------------------------
!
# if defined MODEL_COUPLING && !defined MCT_LIB
      IF (NendStep.eq.ntend(1)) THEN
        DoneStepping=.TRUE.
      ELSE
        DoneStepping=.FALSE.
      END IF
# else
      DoneStepping=.TRUE.
# endif
!
      IF (DoneStepping) THEN
        DO ng=1,Ngrids
# ifdef AVERAGES
          LdefAVG(ng)=.FALSE.
          LwrtAVG(ng)=.FALSE.
# endif
# ifdef DIAGNOSTICS
          LdefDIA(ng)=.FALSE.
          LwrtDIA(ng)=.FALSE.
# endif
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
        END DO
!
!  Set structure for the nonlinear forward trajectory to be processed
!  by the tangent linear and adjoint models. Also, set switches to
!  process the FWD structure in routine "check_multifile".  Notice that
!  it is possible to split solution into multiple NetCDF files to reduce
!  their size.
!
        CALL edit_multifile ('HIS2FWD')
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
        CALL edit_multifile ('QCK2BLK')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        DO ng=1,Ngrids
          LreadBLK(ng)=.TRUE.
        END DO
# endif
!
!  Report data penalty function.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            DO i=0,NobsVar(ng)
              IF (i.eq.0) THEN
                string='Total'
              ELSE
                string=ObsName(i)
              END IF
              IF (FOURDVAR(ng)%NLPenalty(i).ne.0.0_r8) THEN
                WRITE (stdout,30) my_outer, Ninner, 'NLM',              &
     &                            FOURDVAR(ng)%NLPenalty(i),            &
     &                            TRIM(string)
              END IF
            END DO
          END IF
!
!  Write out final data penalty function to NetCDF file.
!
          SourceFile=MyFile
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,             &
     &                              'NL_fDataPenalty',                  &
     &                              FOURDVAR(ng)%NLPenalty(0:),         &
     &                              (/1,my_outer/),                     &
     &                              (/NobsVar(ng)+1,1/),                &
     &                              ncid = DAV(ng)%ncid)

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iNLM, DAV(ng)%name,         &
     &                                  'NL_fDataPenalty',              &
     &                                  FOURDVAR(ng)%NLPenalty(0:),     &
     &                                  (/1,my_outer/),                 &
     &                                  (/NobsVar(ng)+1,1/),            &
     &                                  pioFile = DAV(ng)%pioFile)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Clean penalty array before next run of NL model.
!
          FOURDVAR(ng)%NLPenalty=0.0_r8
        END DO
      END IF

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
      USE def_norm_mod,      ONLY : def_norm
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
      CALL get_state (ng, 10, 10, STD(1,ng), STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      IF (NSA.eq.2) THEN
        CALL get_state (ng, 11, 11, STD(2,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 12, 12, STD(3,ng), STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 13, 13, STD(4,ng), STDrec, Tindex)
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        LdefNRM(1:4,ng)=.FALSE.
        LwrtNRM(1:4,ng)=.FALSE.

      ELSE

        NRMrec=1
        CALL get_state (ng, 14, 14, NRM(1,ng), NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        IF (NSA.eq.2) THEN
          CALL get_state (ng, 15, 15, NRM(2,ng), NRMrec, 2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef ADJUST_BOUNDARY
        CALL get_state (ng, 16, 16, NRM(3,ng), NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        CALL get_state (ng, 17, 17, NRM(4,ng), NRMrec, 1)
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
      USE random_ic_mod,     ONLY : random_ic
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
      character (len=8) :: driver = 'rbl4dvar'

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
!  Compute the diagonal of the posterior/analysis error covariance
!  matrix. The result is written to record 2 of the ITL netcdf file.
!
      VAR_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL posterior_var (ng, tile, iTLM, outer)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END DO VAR_OLOOP
!
!  Write out the diagonal of the posterior/analysis covariance matrix
!  which is in tl_var(Rec1) to 4DVar error NetCDF file.
!
      DO ng=1,Ngrids
#   ifdef DISTRIBUTE
        CALL wrt_error (ng, MyRank, Rec1, Rec1)
#   else
        CALL wrt_error (ng, -1, Rec1, Rec1)
#   endif
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
          IF (inner.ne.0) THEN
            DO ng=1,Ngrids
              Rec=1
              CALL get_state (ng, iTLM, 1, ITL(ng), Rec, Lold(ng))
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
#   ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
#   else
            CALL ad_wrt_his (ng, -1)
#   endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtState2d(ng)=.FALSE.
          END DO
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing adjustments for next inner
!  loop into ITL(ng)%name (record Rec1).
!
          DO ng=1,Ngrids
#   ifdef DISTRIBUTE
            CALL tl_wrt_ini (ng, MyRank, Lnew(ng), Rec1)
#   else
            CALL tl_wrt_ini (ng, -1, Lnew(ng), Rec1)
#   endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO

        END DO POST_ILOOP

      END DO POST_OLOOP

#  endif
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        IF (HIS(ng)%IOtype.eq.io_nf90) THEN
          HIS(ng)%ncid=-1
#  if defined PIO_LIB && defined DISTRIBUTE
        ELSE IF (HIS(ng)%IOtype.eq.io_pio)
          HIS(ng)%pioFile%fh=-1
#  endif
        END IF
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
      END MODULE rbl4dvar_mod
