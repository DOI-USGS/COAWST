#include "cppdefs.h"
      MODULE i4dvar_mod

#ifdef I4DVAR
!
!svn $Id: i4dvar.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module splits the I4D-Var data assimilation algorithm into     !
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
!      iterations to compute the data assimilation increment, dXa.     !
!                                                                      !
!    analysis:                                                         !
!                                                                      !
!      Computes the nonlinear model new initial conditions by adding   !
!      the 4D-Var data assimilation increment to the background        !
!      initial state, Xa = Xb(t=0) + dXa. Then, writes it into NLM     !
!      initial conditions NetCDF INI(ng)%name (record Lini).           !
!                                                                      !
!    posterior_analysis_initialize:                                    !
!                                                                      !
!      Initializes the nonlinear model with the 4D-Var analysis, Xa.   !
!      It is separated from 'posterior_analysis' phase to allow ESM    !
!      coupling.                                                       !
!                                                                      !
!    posterior_analysis:                                               !
!                                                                      !
!      Timesteps the nonlinear model to compute the posterior analysis !
!      state. It interpolates solution at observation locations        !
!                                                                      !
!    prior_error:                                                      !
!                                                                      !
!       Processes the prior background error covariance and its        !
!       normalization coefficients.                                    !
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
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC  :: background_initialize
      PUBLIC  :: background
      PUBLIC  :: increment
      PUBLIC  :: analysis
      PUBLIC  :: posterior_analysis_initialize
      PUBLIC  :: posterior_analysis
      PUBLIC  :: prior_error
!
!  Set module internal parameters.
!
      logical :: LgetNRM        ! process covariance normalization
      logical :: LgetSTD        ! process covariance standard deviation
# ifdef SPLIT_4DVAR
      logical :: Ldone          ! 4D-Var cycle finish switch
# endif
!
      integer :: LTLM1 = 1      ! trial x-space TLM IC record in ITL
      integer :: LTLM2 = 2      ! previous v-space TLM IC record in ITL
      integer :: LTLM3 = 3      ! trial v-space TLM IC record in ITL
      integer :: LADJ1 = 1      ! initial cost gradient
      integer :: LADJ2 = 2      ! new cost gradient (not normalized)
      integer :: Lini  = 1      ! NLM initial conditions record in INI
      integer :: Lbck  = 2      ! background record in INI
      integer :: Rec1  = 1
      integer :: Rec2  = 2
      integer :: Rec3  = 3
      integer :: Rec4  = 4
# ifdef SPLIT_4DVAR
      integer :: Rec5  = 5
# endif
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
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE mod_mixing,       ONLY : initialize_mixing
      USE strings_mod,      ONLY : FoundError
      USE wrt_ini_mod,      ONLY : wrt_ini
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY : balance_ref, biconj
# endif
!
!  Imported variable declarations
!
      integer, intent(in) :: my_outer
!
!  Local variable declarations.
!
      integer :: lstr, ng, tile
      integer :: Fcount, InpRec
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
!  Clear nonlinear mixing arrays.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
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
!  If outer>1, set ERstr=Nrun to compute the open boundary conditions
!  (OBC_time) and surface forcing (SF_time) adjustment times, if needed.
!
      IF (my_outer.gt.1) THEN
        Nrun=1+(my_outer-1)*(Ninner+1)
        ERstr=Nrun
      END IF

#  if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
      defined ADJUST_WSTRESS
!
!  If outer>1, reset reset surface forcing and lateral boundary output
!  time indices to LADJ2 since they are changed in the inner loops.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          Lfout(ng)=LADJ2
#   ifdef ADJUST_BOUNDARY
          Lbout(ng)=LADJ2
#   endif
        END DO
      END IF
#  endif
!
!  If outer>1, open Nonlinear model initial conditions NetCDF file
!  (INI struc) and inquire about its variable IDs. Reset the value of
!  INI(ng)%Rindex, which is needed in "initial" to the 4D-Var analysis
!  record.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefINI(ng)=.FALSE.
          CALL def_ini (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          INI(ng)%Rindex=Lini
        END DO
      END IF
!
!  If outer>1, open 4D-Var NetCDF file (DAV struc) and inquire about its
!  variables.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefMOD(ng)=.FALSE.
          CALL def_mod (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          haveNLmod(ng)=.FALSE.                   ! do not read NLmodVal
        END DO
      END IF
!
!  If outer>1, read in the global observation screening and quality
!  control scale ObsScaleGlobal ("obs_scale") that it is written in the
!  first outer loop. Such data is in memory in the unsplit algorithm.
!  Recall that in I4D-Var, observation variables are read and load into
!  the arrays elements 1:Nobs(ng) for each survey time.  That is, only
!  the values needed are read.
!
!  HGA: What to do in 4D-Var with nested grids?
!
      IF (my_outer.gt.1) THEN
        ng=1
        CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                   &
     &                        Vname(1,idObsS), ObsScaleGlobal,          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  If outer>1, read in initial (outer=1) nonlinear cost function values
!  from DAV file and load it into "CostNorm", its values are used for
!  reporting normalized values.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'NLcost_function',                      &
     &                          FOURDVAR(ng)%CostNorm(0:),              &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/1,1/),                        &
     &                          total = (/NobsVar(ng)+1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF

# if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
!
!  If outer>1, read in surface forcing and or lateral boundary
!  increments (X-space) from ITL file Rec5.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefITL(ng)=.FALSE.
          CALL tl_def_ini (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          InpRec=Rec5
          CALL get_state (ng, iTLM, 1, ITL(ng)%name, InpRec, LTLM1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#  endif

#  ifdef FORWARD_FLUXES
!
!  If outer>1, set structure to process nonlinear surface fluxes from
!  the initial (outer=1) NLM quicksave trajectory. Notice that 2 is
!  substracted since the counter starts at 'outer0'.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), my_outer-2
          lstr=LEN_TRIM(QCK(ng)%name)
          QCK(ng)%base=QCK(ng)%name(1:lstr-3)
        END DO
!
        CALL edit_multifile ('QCK2BLK')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Initialize nonlinear model.
!-----------------------------------------------------------------------
!
!  Activate switch to write out initial misfit between model and
!  observations.
!
      DO ng=1,Ngrids
        IF (my_outer.eq.1) THEN
          wrtMisfit(ng)=.TRUE.
        ELSE
          wrtMisfit(ng)=.FALSE.
        END IF
      END DO
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop.
!
      DO ng=1,Ngrids
        idefHIS(ng)=-1
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        LreadFWD(ng)=.FALSE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), my_outer-1
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
# ifdef FORWARD_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), my_outer-1
        lstr=LEN_TRIM(QCK(ng)%name)
        QCK(ng)%base=QCK(ng)%name(1:lstr-3)
      END DO
!
!  Initialize nonlinear model. If outer=1, the model is initialized
!  with the background or reference state. Otherwise, the model is
!  initialized with the estimated initial conditions from previous
!  iteration, X(0) = X(0) + deltaX(0).
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%load
        RST(ng)%Nrec(Fcount)=0
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  If first pass, save nonlinear initial conditions (currently in time
!  index 1, background) into next record (Lbck) of INI(ng)%name NetCDF
!  file. The record "Lbck" becomes the background state record and the
!  record "Lini" becomes current nonlinear initial conditions.  Both
!  records are used in the algorithm below.
!
      IF (my_outer.eq.1) THEN
        DO ng=1,Ngrids
          INI(ng)%Rindex=1
          Fcount=INI(ng)%load
          INI(ng)%Nrec(Fcount)=1
          CALL wrt_ini (ng, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF

# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
      IF (balance(isFsur)) THEN
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL balance_ref (ng, tile, Lini)
            CALL biconj (ng, tile, iNLM, Lini)
          END DO
          wrtZetaRef(ng)=.TRUE.
        END DO
      END IF
# endif
!
!  If first pass, define output 4DVAR NetCDF file containing all
!  processed data at observation locations.
!
      IF (my_outer.eq.1) THEN
        DO ng=1,Ngrids
          LdefMOD(ng)=.TRUE.
          CALL def_mod (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
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
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE strings_mod, ONLY : FoundError
# if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
     defined ADJUST_WSTRESS
      USE wrt_ini_mod, ONLY : wrt_frc
# endif
!
!  Imported variable declarations
!
      integer, intent(in)  :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoneStepping
!
      integer :: i, lstr, ng
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
!  Run nonlinear model. Save nonlinear tracjectory needed by the
!  adjoint and tangent linear models. Interpolate nonlinear model
!  to observation locations (compute and save H x). It processes
!  and writes the observations accept/reject flag (ObsScale) once
!  to allow background quality control, if any.
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
!-----------------------------------------------------------------------
!  If completed stepping, write out into NetCDF files.
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

# if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
     defined ADJUST_WSTRESS
!
!  Write out initial and background surface forcing into initial
!  INI(ng)%name NetCDF file for latter use.
!
        DO ng=1,Ngrids
          CALL wrt_frc (ng, Lfout(ng), Lini)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (my_outer.eq.1) THEN
            CALL wrt_frc (ng, Lfout(ng), Lbck)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
!
!  Write out nonlinear model misfit cost function into DAV(ng)%name
!  NetCDF file.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'NLcost_function',                      &
     &                          FOURDVAR(ng)%NLobsCost(0:),             &
     &                          (/1,my_outer/), (/NobsVar(ng)+1,1/),    &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
# ifdef BALANCE_OPERATOR
      USE ad_balance_mod,     ONLY : ad_balance
# endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE back_cost_mod,      ONLY : back_cost
      USE cgradient_mod,      ONLY : cgradient
# ifdef SPLIT_4DVAR
      USE cgradient_mod,      ONLY : cg_read_cgradient
# endif
      USE cost_grad_mod,      ONLY : cost_grad
# ifdef ADJUST_BOUNDARY
      USE mod_boundary,       ONLY : initialize_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces,         ONLY : initialize_forces
# endif
      USE mod_ocean,          ONLY : initialize_ocean
      USE strings_mod,        ONLY : FoundError
# ifdef BALANCE_OPERATOR
      USE tl_balance_mod,     ONLY : tl_balance
# endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Imported variable declarations
!
      logical :: Lweak = .FALSE.
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: i, ifile, lstr, my_inner, ng, tile
      integer :: Fcount, InpRec, Lcon, Lsav
# ifdef PROFILE
      integer :: thread
# endif
!
      real(r8) :: rate
# ifdef SPLIT_4DVAR
      real(dp) :: stime
# endif
!
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
        Nrun=1+(my_outer-1)*(Ninner+1)
      END IF
!
!  Open 4D-Var NetCDF file (DAV struc) and inquire about its variables.
!  Activate "haveNLmod" to read its values in calls to "obs_read". Its
!  values were written in the "background" phase.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.FALSE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        haveNLmod(ng)=.TRUE.
      END DO
!
!  In the split 4D-Var algorithm, we need to read the global observation
!  screening and quality control scale ObsScaleGlobal ("obs_scale") that
!  it is written in the "background" phase. Such data is in memory in
!  the unsplit algorithm.  Recall that in I4D-Var, observation variables
!  are read and load into the arrays elements 1:Nobs(ng) for each survey
!  time.  That is, only the values needed are read.
!
!  HGA: What to do in 4D-Var with nested grids?
!
      ng=1
      CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                     &
     &                      Vname(1,idObsS), ObsScaleGlobal,            &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in first value of NLM cost function computed at outer=1, and
!  load it into "CostNorm" that is used in the reporting of normalized
!  values.
!
      SourceFile=MyFile
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name,                   &
     &                        'NLcost_function',                        &
     &                        FOURDVAR(ng)%CostNorm(0:),                &
     &                        ncid = DAV(ng)%ncid,                      &
     &                        start = (/1,1/),                          &
     &                        total = (/NobsVar(ng)+1,1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Read in current outer loop NLM cost function.
!
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name,                   &
     &                        'NLcost_function',                        &
     &                        FOURDVAR(ng)%NLobsCost(0:),               &
     &                        ncid = DAV(ng)%ncid,                      &
     &                        start = (/1,my_outer/),                   &
     &                        total = (/NobsVar(ng)+1,1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  If outer>1, read in previous background and TLM cost function to
!  compute FOURDVAR(ng)%CostFun(0).  Its value is in memory in the
!  unsplit driver.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'TLcost_function',                      &
     &                          FOURDVAR(ng)%CostFunOld(0:),            &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/Nrun-1/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'back_function',                        &
     &                          FOURDVAR(ng)%BackCost(0:),              &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          start = (/Nrun-1/),                     &
     &                          total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          FOURDVAR(ng)%CostFun(0)=FOURDVAR(ng)%CostFunOld(0)+           &
     &                            FOURDVAR(ng)%BackCost(0)
        END DO
      END IF
!
!  If outer>1, read several conjugate gradient variables from 4D-VAR
!  NetCDF (DAV struc) file.  In the unsplit case, such values are
!  available in memory.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          CALL cg_read_cgradient (ng, iTLM, my_outer)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
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
!  If outer>1, open tangent linear history NetCDF file (TLM struc) and
!  inquire about its variables IDs.
!
      IF (my_outer.gt.1) THEN
        DO ng=1,Ngrids
          LdefTLM(ng)=.FALSE.
          CALL tl_def_his (ng, LdefTLM(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!  Set nonlinear output history file to be used as the basis state
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
# endif
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
!  Set various switches.
!
      DO ng=1,Ngrids
        LwrtCost(ng)=.TRUE.
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
        wrtTLmod(ng)=.TRUE.
      END DO
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
!
!  Clear the nonlinear state arrays to make sure that the RHS terms
!  rzeta, rubar, rvbar, ru, and rv are zero. Otherwise, those arrays
!  will have the last computed values when running the nonlinear model
!  if not processing the forward trajectory RHS terms. It needs to be
!  done to get identical solutions with the split schemes.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
        END DO
      END DO
!
!  The minimization algorithm requires to save all the gradient
!  solutions for each inner loop iteration.  They are used for
!  orthogonalization in the conjugate gradient algorithm.  Thus,
!  we need to reset adjoint file record indices.
!
      DO ng=1,Ngrids
        ADM(ng)%Rindex=0
        Fcount=ADM(ng)%load
        ADM(ng)%Nrec(Fcount)=0
      END DO
!
!  An adjoint NetCDF is created for each outer loop.
!
      DO ng=1,Ngrids
        idefADJ(ng)=-1
        LdefADJ(ng)=.TRUE.
        WRITE (ADM(ng)%name,10) TRIM(ADM(ng)%head), my_outer
        lstr=LEN_TRIM(ADM(ng)%name)
        ADM(ng)%base=ADM(ng)%name(1:lstr-3)
      END DO
!
!  Define output Hessian NetCDF file containing the eigenvectors
!  approximation to the Hessian matrix computed from the Lanczos
!  algorithm. Notice that the file name is a function of the
!  outer loop. That is, a file is created for each outer loop.
!
      DO ng=1,Ngrids
        LdefHSS(ng)=.TRUE.
        WRITE (HSS(ng)%name,10) TRIM(HSS(ng)%head), my_outer
        lstr=LEN_TRIM(HSS(ng)%name)
        HSS(ng)%base=HSS(ng)%name(1:lstr-3)
        CALL def_hessian (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Notice that inner loop iteration start from zero. This is needed to
!  compute the minimization initial increment deltaX(0), its associated
!  gradient G(0), and descent direction d(0) used in the conjugate
!  gradient algorithm.
!
      INNER_LOOP : DO my_inner=0,Ninner
        inner=my_inner
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time-step tangent linear model: compute cost function.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  If first pass inner=0, initialize tangent linear state (increments,
!  deltaX) from rest. Otherwise, use trial initial conditions estimated
!  by the conjugate gradient algorithm in previous inner loop. The TLM
!  initial conditions are read from ITL(ng)%name, record 1.
!
        DO ng=1,Ngrids
          ITL(ng)%Rindex=1
          CALL tl_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  On first pass, initialize records 2, 3 and 4 of the ITL file to zero.
!
        IF ((my_inner.eq.0).and.(my_outer.eq.1)) THEN
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, LTLM1, Rec2)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            CALL tl_wrt_ini (ng, LTLM1, Rec3)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            CALL tl_wrt_ini (ng, LTLM1, Rec4)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
        END IF

# ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, activate writing and determine
!  output filename. The multiple file option is used to perturb initial
!  state and create ensembles.  The TLM final trajectory is written for
!  each inner loop on separated NetCDF files.
!
        DO ng=1,Ngrids
          idefTLM(ng)=-1
          LdefTLM(ng)=.TRUE.
          LwrtTLM(ng)=.TRUE.
          WRITE (TLM(ng)%name,20) TRIM(TLM(ng)%head), Nrun
          lstr=LEN_TRIM(TLM(ng)%name)
          TLM(ng)%base=TLM(ng)%name(1:lstr-3)
        END DO
# endif
!
!  Run tangent linear model. Compute misfit observation cost function,
!  Jo.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,30) 'TL', ng, my_outer, my_inner,             &
     &                        ntstart(ng), ntend(ng)
          END IF
        END DO
!
# ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
# else
        CALL tl_main2d (RunInterval)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef EVOLVED_LCZ
!
!  Write evolved tangent Lanczos vector into hessian netcdf file for use
!  later.
!
!  NOTE: When using this option, it is important to set LhessianEV and
!  Lprecond to FALSE in s4dvar.in, otherwise the evolved Lanczos vectors
!  with be overwritten by the Hessian eigenvectors. The fix to this
!  is to define a new NetCDF file that contains the evolved Lanczos
!  vectors.
!
        IF (my_inner.ne.0) THEN
          DO ng=1,Ngrids
            CALL wrt_evolved (ng, kstp(ng), nrhs(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
        END IF
# endif

# ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, close current NetCDF file.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          IF (TLM(ng)%ncid.ne.-1) THEN
            CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step adjoint model backwards: compute cost function gradient.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model always from rest.
!
        DO ng=1,Ngrids
          CALL ad_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Time-step adjoint model backwards. The adjoint model is forced with
!  the adjoint of the observation misfit (Jo) term.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,30) 'AD', ng, my_outer, my_inner,             &
     &                        ntstart(ng), ntend(ng)
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
!  Clear adjoint arrays.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iADM)
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            CALL initialize_forces (ng, tile, iADM)
# endif
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iADM)
# endif
          END DO
        END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Conjugate gradient algorithm.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read TLM v-space initial conditions, record 3 in ITL(ng)%name, and
!  load it into time index LTLM1. This is needed to compute background
!  cost function. Also read in new (x-space) gradient vector, GRADx(Jo),
!  from adjoint history file ADM(ng)%name.  Read in the sum of all the
!  previous outer-loop increments which are always in record 4 of
!  the ITL file.
!
        DO ng=1,Ngrids
          IF (my_inner.eq.0) THEN
            InpRec=Rec1
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, InpRec, LTLM1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          ELSE
            InpRec=Rec3
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, InpRec, LTLM1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

          InpRec=Rec4
          CALL get_state (ng, iTLM, 8, ITL(ng)%name, InpRec, LTLM2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          InpRec=ADM(ng)%Rindex
          CALL get_state (ng, iADM, 4, ADM(ng)%name, InpRec, LADJ2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef BALANCE_OPERATOR
          InpRec=Lini
          CALL get_state (ng, iNLM, 2, INI(ng)%name, InpRec, Lini)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          nrhs(ng)=Lini
# endif
        END DO
!
!  Convert observation cost function gradient, GRADx(Jo), from model
!  space (x-space) to minimization space (v-space):
!
!     GRADv(Jo) = B^(T/2) GRADx(Jo),  operator: S G L^(T/2) W^(-1/2)
!
!  First, multiply the adjoint solution, GRADx(Jo), by the background-
!  error standard deviations, S.  Second, convolve result with the
!  adjoint diffusion operator, G L^(T/2) W^(-1/2). Then, backgound
!  cost function contribution (BackCost) and cost function gradient
!  (v-space) by adding background and observation contributions:
!
!     GRADv(J) = GRADv(Jb) + GRADv(Jo) = deltaV + GRADv(Jo)
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
# ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, tile, Lini, LADJ2)
# endif
            CALL ad_variability (ng, tile, LADJ2, Lweak)
            CALL ad_convolution (ng, tile, LADJ2, Lweak, 2)
            CALL cost_grad (ng, tile, LTLM1, LTLM2, LADJ2)
          END DO
        END DO
!
!  Compute current total cost function.
!
        DO ng=1,Ngrids
          IF (Nrun.eq.1) THEN
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
              FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
            END DO
          ELSE
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
            END DO
          END IF
        END DO
!
!  Prepare for background cost function (Jb) calculation:
!
!  Read the convolved gradient from inner=0 (which is permanently
!  saved in record 1 of the adjoint file)  ALWAYS into record 1.
!
        IF (my_inner.gt.0) THEN
          DO ng=1,Ngrids
            InpRec=LADJ1
            CALL get_state (ng, iADM, 3, ADM(ng)%name, InpRec, LADJ1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
        END IF
!
!  Compute background cost function (Jb) for inner=0:
!
!  If first pass of inner loop, read in the sum of previous v-space
!  gradients from record 4 of ITL file using the TLM model variables
!  as temporary storage. Also add background cost function to Cost0.
!
        IF (my_inner.eq.0) THEN
          DO ng=1,Ngrids
            InpRec=Rec4
            CALL get_state (ng, iTLM, 2, ITL(ng)%name, InpRec, LTLM2)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL back_cost (ng, tile, LTLM2)
            END DO
!
            FOURDVAR(ng)%Cost0(my_outer)=FOURDVAR(ng)%Cost0(my_outer)+  &
     &                                   FOURDVAR(ng)%BackCost(0)
          END DO
        END IF
!
!  Compute current total cost function.
!
        DO ng=1,Ngrids
          IF (Nrun.eq.1) THEN
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostNorm(i)=FOURDVAR(ng)%CostNorm(i)+        &
     &                                 FOURDVAR(ng)%BackCost(i)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
              FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
            END DO
          ELSE
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
            END DO
          END IF
        END DO
!
!  Determine the descent direction in which the quadractic total cost
!  function decreases. Then, determine the TLM initial conditions,
!  deltaV(LTLM1), and its gradient, GRADv{J(Lnew)} at the new
!  direction.  Also, Compute TLM v-space trial initial conditions for
!  next inner loop, deltaV(LTLM2). The new gradient minimize the
!  quadratic cost function spanned by current and previous inner loop
!  iterations.  This is achieved by orthogonalizing (Gramm-Schmidt
!  algorithm) against all previous inner loop gradients.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL cgradient (ng, tile, iTLM, my_inner, my_outer)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Report background (Jb) and observations (Jo) cost function values
!  normalized by their first minimization value. It also reports the
!  percentage change on total cost function value with respect to
!  previous iteration. Compute the optimality of the minimization to
!  check the statistical hypotheses between the background and
!  observations errors: the cost function value at the minimum, Jmin,
!  is ideally equal to half the number of observations assimilated
!  (Optimality=1=2*Jmin/Nobs), for a linear system.
!
        IF (Master) THEN
          DO ng=1,Ngrids
            IF (Nrun.gt.1) THEN
              rate=100.0_r8*ABS(FOURDVAR(ng)%CostFun(0)-                &
     &                          FOURDVAR(ng)%CostFunOld(0))/            &
     &                      FOURDVAR(ng)%CostFunOld(0)
            ELSE
              rate=0.0_r8
            END IF
            Optimality(ng)=2.0_r8*FOURDVAR(ng)%CostFun(0)/              &
     &                     (FOURDVAR(ng)%ObsCount(0)-                   &
     &                      FOURDVAR(ng)%ObsReject(0))
            WRITE (stdout,40) my_outer, my_inner,                       &
     &                        FOURDVAR(ng)%BackCost(0)/                 &
     &                        FOURDVAR(ng)%CostNorm(0),                 &
     &                        FOURDVAR(ng)%ObsCost(0)/                  &
     &                        FOURDVAR(ng)%CostNorm(0),                 &
     &                        rate
            IF (my_inner.eq.0) THEN
              DO i=0,NobsVar(ng)
                IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
                  IF (i.eq.0) THEN
                    WRITE (stdout,50) my_outer, my_inner,               &
     &                                FOURDVAR(ng)%NLobsCost(i)/        &
     &                                FOURDVAR(ng)%CostNorm(i)
                  ELSE
                    WRITE (stdout,60) my_outer, my_inner,               &
     &                                FOURDVAR(ng)%NLobsCost(i)/        &
     &                                FOURDVAR(ng)%CostNorm(i),         &
     &                                TRIM(ObsName(i))
                  END IF
                END IF
                FOURDVAR(ng)%NLobsCost(i)=0.0
              END DO
            END IF
            WRITE (stdout,70) my_outer, my_inner, Optimality(ng)
          END DO
        END IF
!
!  Save total v-space cost function gradient, GRADv{J(Lnew)}, into
!  ADM(ng)%name history NetCDF file. Noticed that the lastest adjoint
!  solution record is over-written in the NetCDF file for future use.
!  The switch "LwrtState2d" is activated to write out state arrays
!  instead ad_*_sol arrays (HGA).
!
        DO ng=1,Ngrids
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          Lfout(ng)=LADJ2
# endif
# ifdef ADJUST_BOUNDARY
          Lbout(ng)=LADJ2
# endif
          kstp(ng)=LADJ2
# ifdef SOLVE3D
          nstp(ng)=LADJ2
# endif
          ADM(ng)%Rindex=ADM(ng)%Rindex-1
          LwrtState2d(ng)=.TRUE.
          CALL ad_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          LwrtState2d(ng)=.FALSE.
        END DO
!
!  Write out trial v-space TLM initial conditions, currently in time
!  index LTM2, into record 3 of ITL(ng)%name NetCDF file.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, LTLM2, Rec3)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Read current outer loop nonlinear model initial conditions and
!  background state vectors.
!
        DO ng=1,Ngrids
          InpRec=Lini
          CALL get_state (ng, iNLM, 2, INI(ng)%name, InpRec, Lini)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          InpRec=Lbck
          CALL get_state (ng, iNLM, 9, INI(ng)%name, InpRec, Lbck)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Convert increment vector, deltaV, from minimization space (v-space)
!  to model space (x-space):
!
!     deltaX = B^(1/2) deltaV
!  or
!     deltaX = W^(1/2) L^(1/2) G S
!
!  First, convolve estimated increment vector (v-space) by with the
!  tangent linear diffusion operator, W^(1/2) L^(1/2) G.  Second,
!  multiply result by the background-error standard deviation, S.
!
        Lcon=LTLM2
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_convolution (ng, tile, Lcon, Lweak, 2)
            CALL tl_variability (ng, tile, Lcon, Lweak)
# ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, tile, Lini, Lcon)
# endif
          END DO
        END DO
!
!  Write out trial x-space (convolved) TLM initial conditions, currently
!  in time index Lcon, into record 1 of ITL(ng)%name NetCDF file.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, Lcon, Rec1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Update counters.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          Lsav=Lnew(ng)
          Lnew(ng)=Lold(ng)
          Lold(ng)=Lsav
          Nrun=Nrun+1
        END DO
!
      END DO INNER_LOOP
!
!  Turn of switch to write cost functions in the DAV NetCDF file. They
!  are written in calls to "tl_wrt_ini" only inside the inner loops.
!
      DO ng=1,Ngrids
        LwrtCost(ng)=.FALSE.
      END DO
!
!  Close adjoint NetCDF file.
!
      SourceFile=MyFile
      DO ng=1,Ngrids
        IF (ADM(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, ADM(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Close Hessian NetCDF file.
!
      DO ng=1,Ngrids
        IF (HSS(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, HSS(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO

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
 20   FORMAT (a,'_member',i3.3,'.nc')
 30   FORMAT (/,1x,a,1x,'ROMS: started time-stepping:',                 &
     &        ' (Grid: ',i0,', Outer: ',i2.2,', Inner: ',i3.3,          &
              ', TimeSteps: ',i0,' - ',i0,')',/)
 40   FORMAT (/,' (',i3.3,',',i3.3,'): TLM Cost Jb, J  = ',             &
     &        1p,e17.10,0p,1x,1p,e17.10,0p,t68,1p,e11.4,' %')
 50   FORMAT (/,'>(',i3.3,',',i3.3,'): NLM Cost     J  = ',             &
     &        18x,1p,e17.10,0p)
 60   FORMAT (' (',i3.3,',',i3.3,'): NLM Cost     J  = ',               &
     &        18x,1p,e17.10,0p,t69,a)
 70   FORMAT (/,1x,'(',i3.3,',',i3.3,'): Optimality (2*J/Nobs) = ',     &
     &        1p,e17.10,/)
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
      USE mod_scalars
      USE mod_stepping
!
      USE ini_adjust_mod,     ONLY : ini_adjust
!!    USE ini_fields_mod,     ONLY : ini_fields
      USE mod_ocean,          ONLY : initialize_ocean
#if defined MASKING && defined SPLIT_4DVAR
      USE set_masks_mod,      ONLY : set_masks
#endif
      USE strings_mod,        ONLY : FoundError
      USE sum_grad_mod,       ONLY : sum_grad
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
      USE wrt_ini_mod,        ONLY : wrt_ini
!
!  Imported variable declarations
!
      logical :: Lweak = .FALSE.
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: Fcount, Lcon, Nrec
      integer :: InpRec, OutRec
# ifdef PROFILE
      integer :: thread
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", analysis"
!
      SourceFile=MyFile
!
!=======================================================================
!  Compute new nonlinear initial conditions by adding minimization
!  increment to previous outer loop initial conditions:
!
!         Xi(outer+1) = Xi(outer) + deltaX(Lcon)
!
!=======================================================================

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
!  Clear nonlinear state variables.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
        END DO
      END DO

# ifdef SPLIT_4DVAR
!
!-----------------------------------------------------------------------
!  If split 4D-Var algorithm, set several variables that computed or
!  assigned in other 4D-Var phase executable.
!-----------------------------------------------------------------------
!
!  Open nonlinear model initial conditions NetCDF file (INI struc) and
!  inquire about its variable IDs.
!
      DO ng=1,Ngrids
        LdefINI(ng)=.FALSE.
        CALL def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Open tangent linear model initial conditions NetCDF file (ITL struc)
!  and inquire about its variable IDs.
!
      DO ng=1,Ngrids
        LdefITL(ng)=.FALSE.
        CALL tl_def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Open 4D-Var NetCDF file (DAV struc) and inquire about its variables.
!  It is used in "tl_wrt_ini" to write out observation and background
!  cost functions.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.FALSE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        Nrec=DAV(ng)%Rindex
      END DO
!
!  Set Nrun>1, set Nrun couter. The counter is only allow to increase
!  in the "increment" phase.
!
      Nrun=1+(my_outer-1)*(Ninner+1)
      ERstr=Nrun

#  ifdef MASKING
!
!  Set internal mask arrays to process I/O NetCDF files.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_masks (ng, tile, iNLM)
        END DO
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute nonlinear model initial from the 4D-Var analysis.
!-----------------------------------------------------------------------
!
!  In order to use the correct fields, the model time indices are set
!  to Lini.
!
!  The appropriate TLM correction for the NLM model resides in record 1
!  of the ITL file.
!
      DO ng=1,Ngrids
        kstp(ng)=Lini
# ifdef SOLVE3D
        nstp(ng)=Lini
# endif
        InpRec=Lini
        CALL get_state (ng, iNLM, 1, INI(ng)%name, InpRec, Lini)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        InpRec=LTLM1
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, InpRec, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_adjust (ng, tile, LTLM1, Lini)
!!        CALL ini_fields (ng, tile, iNLM)
        END DO
      END DO
!
!  Write out new nonlinear model initial conditions into record Lini
!  of INI(ng)%name.
!
      DO ng=1,Ngrids
        Fcount=INI(ng)%load
        INI(ng)%Nrec(Fcount)=1
        OutRec=Lini
        CALL wrt_ini (ng, Lini, OutRec)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Gather the v-space increments from the final inner-loop and
!  save in record 4 of the ITL file. The current v-space increment
!  is in record 3 and the sum so far is in record 4.
!
      DO ng=1,Ngrids
        InpRec=Rec3
        CALL get_state (ng, iTLM, 8, ITL(ng)%name, InpRec, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        InpRec=Rec4
        CALL get_state (ng, iTLM, 8, ITL(ng)%name, InpRec, LTLM2)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL sum_grad (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
! Write the current sum into record 4 of the ITL file.
!
      DO ng=1,Ngrids
        CALL tl_wrt_ini (ng, LTLM2, Rec4)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

# if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
!
!  Set index containing the surface forcing increments used the run
!  the nonlinear model in the outer loop and read the forcing
!  increments. For bulk fluxes, we read Rec1 because the stress
!  fluxes change by virtue of the changing initial conditions.
!  When not using bulk fluxes, we read Rec4 because the background
!  stress and flux is prescribed by input files which are not
!  overwritten so we need to correct the background using the
!  sum of the increments from all previous outer-loops.
!  If using Rec4 we need to convert from v-space to x-space
!  by applying the convolution.
!  Note that Lfinp=Lbinp so the the forcing and boundary
!  adjustments are both processsed correctly.
#  ifdef BALANCE_OPERATOR
!  Currently, We do not need the call to tl_balance below, but we
!  might later if we impose a balance constraint on the wind stress
!  corrections.
#  endif
!
!  AMM: CHECK WHAT HAPPENS WITH SECONDARY PRECONDITIONING.
!
      DO ng=1,Ngrids
        Lfinp(ng)=LTLM1
#  if defined BULK_FLUXES && !defined FORWARD_FLUXES
        InpRec=Rec1
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, InpRec, Lfinp(ng))
#  endif
#  if defined FORWARD_FLUXES || !defined BULK_FLUXES
        InpRec=Rec4
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec4, Lfinp(ng))
        Lcon=Lfinp(ng)
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL tl_convolution (ng, tile, Lcon, Lweak, 2)
          CALL tl_variability (ng, tile, Lcon, Lweak)
#   ifdef BALANCE_OPERATOR
!!        CALL tl_balance (ng, tile, Lini, Lcon)
#   endif
        END DO
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

#  ifdef SPLIT_4DVAR
!
!  If split 4D-Var, write out surface forcing and lateral boundary
!  condition increments (X-space) needed for the nonlinear model into
!  Rec5.
!
      DO ng=1,Ngrids
        CALL tl_wrt_ini (ng, Lfinp(ng), Rec5)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#  endif
# endif
!
!  Clear tangent linear state variables.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iTLM)
        END DO
      END DO
!
!  Close current forward NetCDF file.
!
      SourceFile=MyFile
      DO ng=1,Ngrids
!       CALL netcdf_close (ng, iNLM, FWD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      RETURN
      END SUBROUTINE analysis
!
      SUBROUTINE posterior_analysis_initialize
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
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE mod_mixing,  ONLY : initialize_mixing
      USE strings_mod, ONLY : FoundError
!
!  Local variable declarations.
!
      integer :: i, ifile, lstr, ng, tile
      integer :: Fcount, InpRec
# ifdef PROFILE
      integer :: thread
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", posterior_analysis_initialize"
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
!  Clear nonlinear mixing arrays.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
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
      outer=Nouter
      inner=Ninner
!
!  Set Nrun>1, to read in surface forcing and open boundary conditions
!  increments in "initial", if appropriate.
!
      Nrun=1+Nouter*(Ninner+1)
!
!  Set ERstr=Nrun, to set the open boundary condition (OBC_time) and
!  surface forcing (SF_time) adjustment times, if needed.
!
      ERstr=Nrun
!
!  Open Nonlinear model initial conditions NetCDF file (INI struc) and
!  inquire about its variable IDs.
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
        haveNLmod(ng)=.FALSE.                     ! do not read NLmodVal
      END DO
!
!  In the split 4D-Var algorithm, we need to read the global observation
!  screening and quality control scale ObsScaleGlobal ("obs_scale") that
!  it is written in the "background" phase. Such data is in memory in
!  the unsplit algorithm.  Recall that in I4D-Var, observation variables
!  are read and load into the arrays elements 1:Nobs(ng) for each survey
!  time.  That is, only the values needed are read.
!
!  HGA: What to do in 4D-Var with nested grids?
!
      ng=1
      CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,                     &
     &                      Vname(1,idObsS), ObsScaleGlobal,            &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read initial nonlinear cost function values from DAV file and load it
!  into CostNorm, its values are used for reporting normalized values.
!
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name,                   &
     &                        'NLcost_function',                        &
     &                        FOURDVAR(ng)%CostNorm(0:),                &
     &                        ncid = DAV(ng)%ncid,                      &
     &                        start = (/1,1/),                          &
     &                        total = (/NobsVar(ng)+1,1/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

# if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
!
!  If split 4D-Var, read in surface forcing and or lateral boundary
!  increments (X-space) from ITL file Rec5.
!
      DO ng=1,Ngrids
        LdefITL(ng)=.FALSE.
        CALL tl_def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        InpRec=Rec5
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, InpRec, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  If split 4D-Var, reset reset surface forcing and lateral boundary
!  output time indices to LADJ2 since they are changed in the inner
!  loops.
!
      DO ng=1,Ngrids
        Lfout(ng)=LADJ2
#   ifdef ADJUST_BOUNDARY
        Lbout(ng)=LADJ2
#   endif
      END DO
#  endif

#  ifdef FORWARD_FLUXES
!
!  If not first NLM run, set BLK structure to the previous quicksave
!  trajectory.  There is a logic in "get_data" that reads from BLK.
!  (HGA: This probably legacy code that it is no longer needed)
!
      DO ng=1,Ngrids
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), Nouter-1
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
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop. Notice that the LreadBLK and LreadFWD switches
!  are turned off to suppress processing of the structures when
!  "check_multifile" during nonlinear model execution.
!
      DO ng=1,Ngrids
        idefHIS(ng)=-1
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
# ifdef FORWARD_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
        LreadFWD(ng)=.FALSE.
        HIS(ng)%Rindex=0
        Fcount=HIS(ng)%load
        HIS(ng)%Nrec(Fcount)=0
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), Nouter
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
      END DO
!
!  Set the nonlinear model output quicksave-history file for the
!  posterior analysis trajectory. Notice that the LreadBLK switch
!  is turned off to suppress processing of the structure in
!  "check_multifile" during nonlinear model execution.
!
      DO ng=1,Ngrids
        idefQCK(ng)=-1
        LdefQCK(ng)=.TRUE.
        LwrtQCK(ng)=.TRUE.
# ifdef FORWARD_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
        WRITE (QCK(ng)%name,10) TRIM(QCK(ng)%head), Nouter
        lstr=LEN_TRIM(QCK(ng)%name)
        QCK(ng)%base=QCK(ng)%name(1:lstr-3)
      END DO
!
!  Initialize nonlinear model with estimated initial conditions.  Reset
!  the value of INI(ng)%Rindex, which is needed in "initial" to the
!  4D-Var analysis record.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        wrtMisfit(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%load
        RST(ng)%Nrec(Fcount)=0
        INI(ng)%Rindex=Lini
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Clear NLobsCost. Activate switch to write out final misfit between
!  model and observations.
!
      DO ng=1,Ngrids
        DO i=0,NobsVar(ng)
          FOURDVAR(ng)%NLobsCost(i)=0.0_r8
        END DO
        wrtMisfit(ng)=.TRUE.
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
!
      RETURN
      END SUBROUTINE posterior_analysis_initialize
!
      SUBROUTINE posterior_analysis (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine initialize the NLM with estimated 4D-Var state and     !
!  interpolates solution at observation locations for posterior        !
!  analysis.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
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
      USE mod_mixing,  ONLY : initialize_mixing
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoneStepping
!
      integer :: i, lstr, ng
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
     &  __FILE__//", posterior_analysis"
!
      SourceFile=MyFile

# if !(defined MODEL_COUPLING && defined ESMF_LIB)
!
!-----------------------------------------------------------------------
!  Set nonlinear model initial conditions.
!-----------------------------------------------------------------------
!
      CALL posterior_analysis_initialize
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Run nonlinear model. Interpolate nonlinear model to observation
!  locations.
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
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
# ifdef DIAGNOSTICS
        idefDIA(ng)=-1
        LdefDIA(ng)=.TRUE.
        LwrtDIA(ng)=.TRUE.
        WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%head), outer
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
        IF (Master) WRITE (stdout,20) 'NL', ng, NstrStep, NendStep
# else
        IF (Master) WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
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
!-----------------------------------------------------------------------
!  If completed stepping, write out and report cost function.
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
!
!  Write out nonlinear model final misfit cost function into
!  DAV(ng)%name NetCDF file. Notice that it is written in the
!  Nouter+1 record.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'NLcost_function',                      &
     &                          FOURDVAR(ng)%NLobsCost(0:),             &
     &                          (/1,Nouter+1/), (/NobsVar(ng)+1,1/),    &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Report the final value of the nonlinear model misfit cost function.
!
        IF (Master) THEN
          DO ng=1,Ngrids
            DO i=0,NobsVar(ng)
              IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
                IF (i.eq.0) THEN
                  WRITE (stdout,30) outer, inner,                       &
     &                              FOURDVAR(ng)%NLobsCost(i)/          &
     &                              FOURDVAR(ng)%CostNorm(i)
                ELSE
                  WRITE (stdout,40) outer, inner,                       &
     &                              FOURDVAR(ng)%NLobsCost(i)/          &
     &                              FOURDVAR(ng)%CostNorm(i),           &
     &                              TRIM(ObsName(i))
                END IF
              END IF
            END DO
          END DO
        END IF
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
        DO ng=1,Ngrids
          HIS(ng)%ncid=-1
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
     &        ' (Grid: ',i0,', TimeSteps: ',i0,' - ',i0,')',/)
 30   FORMAT (/,'>(',i3.3,',',i3.3,'): NLM Cost     J  = ',             &
     &        18x,1p,e17.10,0p)
 40   FORMAT (' (',i3.3,',',i3.3,'): NLM Cost     J  = ',               &
     &        18x,1p,e17.10,0p,t69,a)
!
      RETURN
      END SUBROUTINE posterior_analysis
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
      IF (LgetSTD) THEN
!
!  Initial conditions standard deviation. They are loaded in Tindex=1
!  of the e_var(...,Tindex) state variables.
!
        STDrec=1
        Tindex=1
        CALL get_state (ng, 10, 10, STD(1,ng)%name, STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
      END IF
!
!-----------------------------------------------------------------------
!  Error covariance normalization coefficients.
!-----------------------------------------------------------------------
!
      IF (LgetNRM) THEN
!
!  Compute or read in background-error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid.
!
        IF (ANY(LwrtNRM(:,ng))) THEN
          CALL def_norm (ng, iNLM, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
          CALL get_state (ng, 14, 14, NRM(1,ng)%name, NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
          CALL get_state (ng, 16, 16, NRM(3,ng)%name, NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL get_state (ng, 17, 17, NRM(4,ng)%name, NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        END IF
      END IF
!
      RETURN
      END SUBROUTINE prior_error
#endif
      END MODULE i4dvar_mod
