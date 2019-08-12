      MODULE ocean_control_mod
!
!svn $Id: obs_sen_w4dpsas_forecast.h$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Strong/Weak Constraint 4-Dimensional Variational Data     !
!    Assimilation and Observation Sensitivity Driver: Physical-space   !
!    Statistical Analysis System (4D-PSAS).                            !
!    Dual formulation in observarion space.                            !
!                                                                      !
!  This driver is used for weak constraint 4D-Var where errors are     !
!  considered in both model and observations. It also computes the     !
!  the sensitivity of the assimilation system to each observation.     !
!  It measures the degree to which each observation contributes to     !
!  the uncertainty in the estimate.  This analysis  can be used to     !
!  determine the type of measurements that need to be made,  where     !
!  to observe, and when.                                               !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
!  conventions:                                                        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Moore, A.M., H.G. Arango, G. Broquet, B.S. Powell, A.T. Weaver,     !
!    and J. Zavala-Garay, 2011: The Regional Ocean Modeling System     !
!    (ROMS)  4-dimensional variational data assimilations systems,     !
!    Part I - System overview and formulation, Prog. Oceanogr., 91,    !
!    34-49, doi:10.1016/j.pocean.2011.05.004.                          !
!                                                                      !
!  Moore, A.M., H.G. Arango, G. Broquet, C. Edward, M. Veneziani,      !
!    B. Powell, D. Foley, J.D. Doyle, D. Costa, and P. Robinson,       !
!    2011: The Regional Ocean Modeling System (ROMS) 4-dimensional     !
!    variational data assimilations systems, Part II - Performance     !
!    and application to the California Current System, Prog.           !
!    Oceanogr., 91, 50-73, doi:10.1016/j.pocean.2011.05.003.           !
!                                                                      !
!  Moore, A.M., H.G. Arango, G. Broquet, C. Edward, M. Veneziani,      !
!    B. Powell, D. Foley, J.D. Doyle, D. Costa, and P. Robinson,       !
!    2011: The Regional Ocean Modeling System (ROMS) 4-dimensional     !
!    variational data assimilations systems, Part III - Observation    !
!    impact and observation sensitivity in the California Current      !
!    System, Prog. Oceanogr., 91, 74-94,                               !
!    doi:10.1016/j.pocean.2011.05.005.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
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
#ifdef MCT_LIB
# ifdef ATM_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAV_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
      USE strings_mod,       ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: STDrec, Tindex
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(mpiCOMM)) THEN
        OCN_COMM_WORLD=mpiCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel control switches. These scalars switches are
!  independent from standard input parameters.
!
        CALL initialize_parallel
!
!  Read in model tunable parameters from standard input. Allocate and
!  initialize variables in several modules after the number of nested
!  grids and dimension parameters are known.
!
        CALL inp_par (iNLM)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
!$OMP PARALLEL
#if defined _OPENMP
      MyThread=my_threadnum()
#elif defined DISTRIBUTE
      MyThread=MyRank
#else
      MyThread=0
#endif
      DO ng=1,Ngrids
        chunk_size=(NtileX(ng)*NtileE(ng)+numthreads-1)/numthreads
        first_tile(ng)=MyThread*chunk_size
        last_tile (ng)=first_tile(ng)+chunk_size-1
      END DO
!$OMP END PARALLEL
!
!  Initialize internal wall clocks. Notice that the timings does not
!  includes processing standard input because several parameters are
!  needed to allocate clock variables.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Process Information:',/)
        END IF
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iNLM, 0, __LINE__, __FILE__)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL mod_arrays (allocate_vars)
!$OMP END PARALLEL
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar

      END IF

#if defined MCT_LIB && (defined ATM_COUPLING || defined WAV_COUPLING)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef ATM_COUPLING
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAV_COUPLING
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif

#ifndef OBS_SPACE
!
!-----------------------------------------------------------------------
!  If the required vectors and arrays from congrad from a previous
!  run of the assimilation cycle are available, read them here from
!  LCZ(ng)%name NetCDF file.
!-----------------------------------------------------------------------
!
      SourceFile=__FILE__ //                                            &
     &           ", ROMS_initialize"
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_beta',        &
     &                        cg_beta)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_delta',       &
     &                        cg_delta)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_Gnorm_v',     &
     &                        cg_Gnorm_v)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_dla',         &
     &                        cg_dla)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_QG',          &
     &                        cg_QG)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'zgrad0',         &
     &                        zgrad0)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'zcglwk',         &
     &                        zcglwk)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'TLmodVal_S',     &
     &                        TLmodVal_S,                               &
     &                        broadcast = .FALSE.)   ! Master use only
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# if defined RPCG && !defined OBS_IMPACT
        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'Hbk',            &
     &                        Hbk)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'Jb0',            &
     &                        Jb0)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'vcglwk',         &
     &                        vcglwk)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif
      END DO
!
!-----------------------------------------------------------------------
!  If skiping runing nonlinear model, read in observation screening and
!  quality control flag.
!-----------------------------------------------------------------------
!
      SourceFile=__FILE__ //                                            &
     &           ", ROMS_initialize"
      wrtObsScale(1:Ngrids)=.FALSE.
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, Vname(1,idObsS),  &
     &                        ObsScale)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
#endif
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
      DO ng=1,Ngrids
        CALL get_state (ng, 10, 10, STD(1,ng)%name, STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      IF (NSA.eq.2) THEN
        DO ng=1,Ngrids
          CALL get_state (ng, 11, 11, STD(2,ng)%name, STDrec, Tindex)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 12, 12, STD(3,ng)%name, STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 13, 13, STD(4,ng)%name, STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
#endif

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS nonlinear, tangent linear and     !
!  adjoint models.                                                     !
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
      USE convolve_mod,      ONLY : error_covariance
#ifdef ADJUST_BOUNDARY
      USE mod_boundary,      ONLY : initialize_boundary
#endif
      USE mod_forces,        ONLY : initialize_forces
      USE mod_ocean,         ONLY : initialize_ocean
      USE normalization_mod, ONLY : normalization
      USE strings_mod,       ONLY : FoundError, uppercase
#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod,  ONLY : balance_ref, biconj
#endif
      USE ini_adjust_mod,     ONLY : load_ADtoTL
      USE ini_adjust_mod,     ONLY : load_TLtoAD!
!
!  Imported variable declarations
!
      real(dp), intent(inout) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: Lcgini, Linner, Lposterior, add

      integer :: my_inner, my_outer
      integer :: Lbck, Lini, Rec1, Rec2, ImpOrd
      integer :: i, ng, status, tile
      integer :: Fcount, NRMrec

      integer, dimension(Ngrids) :: indxSave
      integer, dimension(Ngrids) :: Nrec

      real(r8) :: str_day, end_day, dstartS, rtime

      character (len=1) :: charA, charB
#ifdef OBS_SPACE
      character (len=1) :: charC
#endif
      character (len=15) :: driver
      character (len=20) :: string
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
!  Initialize relevant parameters.
!
      DO ng=1,Ngrids
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Lfinp(ng)=1         ! forcing index for input
        Lfout(ng)=1         ! forcing index for output history files
#endif
#ifdef ADJUST_BOUNDARY
        Lbinp(ng)=1         ! boundary index for input
        Lbout(ng)=1         ! boundary index for output history files
#endif
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
      END DO
      Lini=1                ! NLM initial conditions record in INI
      Lbck=2                ! background record in INI
      Rec1=1
      Rec2=2
      Nrun=1
      outer=0
      inner=0
      ERstr=1
      ERend=Nouter
      driver='obs_sen_w4dpsas_forecast'
      charA='A'
      charB='B'
#ifdef OBS_SPACE
      charC='C'
#endif
!
!  Choose the desired order approximation of the forecast impact
!  calculations:
!
!          ImpOrd: 1 (1st order)
!                  2 (2nd order) or
!                  3 (3rd order)
!
!  ImpOrd=1 is the least accurate and can be in error by a factor of 2.
!
      ImpOrd=3
!
!  For this driver DSTART is assumed to be the start time of the
!  analysis cycle that yields the analysis for the subsequent forecast
!  cycle. The number of time steps in the analysis cycle is NTIMES_ANA,
!  and the number of time steps in the forecast cycle is NTIMES_FCT.
!
!  Initially we will be running the adjoint model over the forecast
!  time interval so we temporarily reset dstart and RunInterval to
!  reflect the forecast interval start and run time.
!
      dstartS=dstart
      rtime=0.0_r8
      DO ng=1,Ngrids
        rtime=MAX(rtime, dt(ng)*ntimes_ana(ng))
      END DO
      dstart=dstart+rtime*sec2day
!
      rtime=0.0_r8
      DO ng=1,Ngrids
        rtime=MAX(rtime, dt(ng)*ntimes_fct(ng))
        ntimes(ng)=ntimes_fct(ng)
      END DO
      RunInterval=rtime
!
!-----------------------------------------------------------------------
!  Configure weak constraint 4DVAR algorithm: PSAS Approach.
!-----------------------------------------------------------------------
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
        LsenPSAS(ng)=.FALSE.
#ifdef OBS_SPACE
        Lobspace(ng)=.FALSE.
# ifndef OBS_IMPACT
        LadjVAR(ng)=.FALSE.
# endif
#endif
      END DO

!$OMP PARALLEL
      CALL initial
!$OMP END PARALLEL
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lbck" of INI(ng)%name NetCDF file. The
!  record "Lbck" becomes the background state record and the record
!  "Lini" becomes current nonlinear initial conditions.
!
      DO ng=1,Ngrids
        INI(ng)%Rindex=1
        Fcount=INI(ng)%Fcount
        INI(ng)%Nrec(Fcount)=1
        CALL wrt_ini (ng, 1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Set nonlinear output history file as the initial basic state
!  trajectory for the specified outer loop (Nimpact) used in the
!  observation impact/sensitivity analysis.
!
      DO ng=1,Ngrids
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%base), Nimpact-1
      END DO

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
      DO ng=1,Ngrids
        BLK(ng)%name=HIS(ng)%name
      END DO
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Model-error covariance normalization and stardard deviation factors.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Compute or read in the error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid and decorrelation scales.
!
      DO ng=1,Ngrids
        IF (ANY(LwrtNRM(:,ng))) THEN
          CALL def_norm (ng, iNLM, 1)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          IF (NSA.eq.2) THEN
            CALL def_norm (ng, iNLM, 2)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
          CALL def_norm (ng, iNLM, 3)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL def_norm (ng, iNLM, 4)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
#endif

!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL normalization (ng, tile, 2)
          END DO
!$OMP END PARALLEL
          LdefNRM(1:4,ng)=.FALSE.
          LwrtNRM(1:4,ng)=.FALSE.
        ELSE
          NRMrec=1
          CALL get_state (ng, 14, 14, NRM(1,ng)%name, NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          IF (NSA.eq.2) THEN
            CALL get_state (ng, 15, 15, NRM(2,ng)%name, NRMrec, 2)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
          CALL get_state (ng, 16, 16, NRM(3,ng)%name, NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL get_state (ng, 17, 17, NRM(4,ng)%name, NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
#endif
        END IF
      END DO

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
      IF (balance(isFsur)) THEN
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
!
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL balance_ref (ng, tile, Lini)
            CALL biconj (ng, tile, iNLM, Lini)
          END DO
!$OMP END PARALLEL
          wrtZetaRef(ng)=.TRUE.
        END DO
      END IF
#endif
!
!  Define tangent linear initial conditions file.
!
      DO ng=1,Ngrids
        LdefITL(ng)=.TRUE.
        CALL tl_def_ini (ng)
        LdefITL(ng)=.FALSE.
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Define impulse forcing NetCDF file.
!
      DO ng=1,Ngrids
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO

#ifndef OBS_SPACE
!
!  Define output 4DVAR NetCDF file containing all processed data
!  at observation locations.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Write out outer loop beeing processed.
!
      SourceFile=__FILE__ // ", ROMS_run"
      DO ng=1,Ngrids
        CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, 'Nimpact',        &
     &                        Nimpact, (/0/), (/0/),                    &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
#endif
!
!  Set forward basic state NetCDF ID to nonlinear model trajectory to
!  avoid the inquiring stage.
!
      DO ng=1,Ngrids
        FWD(ng)%ncid=HIS(ng)%ncid
      END DO
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        HIS(ng)%ncid=-1
        FWD(ng)%ncid=-1
#ifdef OBS_SPACE
        OBS(ng)%ncid=-1
#endif
      END DO
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!  Run the adjoint model forced by forecast minus analysis difference
!  using appropriate solutioh trajectories.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Reset the start and end times for the adjoint forcing.
!
      DO ng=1,Ngrids
        str_day=tdays(ng)+ntimes(ng)*dt(ng)*sec2day
        end_day=tdays(ng)
        IF ((DstrS(ng).eq.0.0_r8).and.(DendS(ng).eq.0.0_r8)) THEN
          DstrS(ng)=end_day
          DendS(ng)=str_day
        END IF
        IF (Master) THEN
          WRITE (stdout,70) 'AD', DendS(ng), DstrS(ng)
        END IF
      END DO
!
!  Set basic state trajectory and adjoint forcing file.
!
!     ads_xxxx_a.nc    if obs_space is off
!     ads_xxxx_b.nc
!
!     obs_xxxx_a.nc    if obs_space is on
!     obs_xxxx_b.nc
!
      DO ng=1,Ngrids
        WRITE (FWD(ng)%name,90) TRIM(FCTB(ng)%base)
        BLK(ng)%name=FWD(ng)%name
#ifdef OBS_SPACE
        IF (ImpOrd.ne.2) THEN
           WRITE (OBS(ng)%name,90) TRIM(OIFB(ng)%base)
        ELSE
           WRITE (OBS(ng)%name,90) TRIM(OIFA(ng)%base)
        END IF
#else
        IF (ImpOrd.ne.2) THEN
           WRITE (ADS(ng)%name,90) TRIM(FOIB(ng)%base)
        ELSE
           WRITE (ADS(ng)%name,90) TRIM(FOIA(ng)%base)
        END IF
#endif
      END DO
      IF (Master) THEN
       WRITE (stdout,50)
      END IF
!
!  Initialize the adjoint model: initialize using dI/dxf is appropriate.
!
      DO ng=1,Ngrids
        Lstiffness=.FALSE.
#ifdef OBS_SPACE
        LsenPSAS(ng)=.FALSE.
        LsenFCT(ng)=.TRUE.
        Lobspace(ng)=.TRUE.
# ifndef OBS_IMPACT
        LadjVAR(ng)=.FALSE.
# endif
#else
        LsenPSAS(ng)=.TRUE.
        LsenFCT(ng)=.FALSE.
#endif
!$OMP PARALLEL
        CALL ad_initial (ng)
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                 __FILE__)) RETURN
        END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
      DO ng=1,Ngrids
        WRTforce=.TRUE.
        IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
        Fcount=ADM(ng)%Fcount
        ADM(ng)%Nrec(Fcount)=0
        ADM(ng)%Rindex=0
      END DO

!  Time-step adjoint model backwards.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
        END IF
      END DO

!$OMP PARALLEL
#ifdef SOLVE3D
      CALL ad_main3d (RunInterval)
#else
      CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
      IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &               __FILE__)) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
      DO ng=1,Ngrids
        CALL ad_wrt_his (ng)
        IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                 __FILE__)) RETURN
      END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
      DO ng=1,Ngrids
        WRTforce(ng)=.FALSE.
        LwrtState2d(ng)=.FALSE.
        CALL ad_wrt_his (ng)
        IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                 __FILE__)) RETURN
      END DO
!
!   Retrieve adjoint solution into tl_var(1).
!
      DO ng=1,Ngrids
         CALL get_state (ng, iTLM, 4, ADM(ng)%name, ADM(ng)%Rindex,   &
     &                      Rec1)
         IF (FoundError(exit_flag, NoError, __LINE__,                 &
     &                  __FILE__)) RETURN

      END DO
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        HIS(ng)%ncid=-1
        FWD(ng)%ncid=-1
#ifdef OBS_SPACE
        OBS(ng)%ncid=-1
#else
        ADS(ng)%ncid=-1
#endif
      END DO
!
!  Set basic state trajectory and adjoint forcing file.
!
      IF (ImpOrd.gt.1) THEN
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,90) TRIM(FCTA(ng)%base)
          BLK(ng)%name=FWD(ng)%name
#ifdef OBS_SPACE
          IF (ImpOrd.eq.2) THEN
             WRITE (OBS(ng)%name,90) TRIM(OIFB(ng)%base)
          ELSE
             WRITE (OBS(ng)%name,90) TRIM(OIFA(ng)%base)
          END IF
#else
          IF (ImpOrd.eq.2) THEN
             WRITE (ADS(ng)%name,90) TRIM(FOIB(ng)%base)
          ELSE
             WRITE (ADS(ng)%name,90) TRIM(FOIA(ng)%base)
          END IF
#endif
        END DO
        IF (Master) THEN
         WRITE (stdout,50)
        END IF
!
!  Initialize the adjoint model: initialize using dI/dxf is appropriate.
!
        DO ng=1,Ngrids
          Lstiffness=.FALSE.
#ifdef OBS_SPACE
          LsenPSAS(ng)=.FALSE.
          LsenFCT(ng)=.TRUE.
          Lobspace(ng)=.TRUE.
# ifndef OBS_IMPACT
          LadjVAR(ng)=.FALSE.
# endif
#else
          LsenPSAS(ng)=.TRUE.
          LsenFCT(ng)=.FALSE.
#endif
!$OMP PARALLEL
          CALL ad_initial (ng)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
        DO ng=1,Ngrids
          WRTforce=.TRUE.
          IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
          Fcount=ADM(ng)%Fcount
          ADM(ng)%Nrec(Fcount)=0
          ADM(ng)%Rindex=0
        END DO
!
!  Time-step adjoint model backwards.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL ad_main3d (RunInterval)
#else
        CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
        DO ng=1,Ngrids
          CALL ad_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
        DO ng=1,Ngrids
          WRTforce(ng)=.FALSE.
          LwrtState2d(ng)=.FALSE.
          CALL ad_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!   Retrieve adjoint solution.
!
        DO ng=1,Ngrids
           CALL get_state (ng, iADM, 4, ADM(ng)%name, ADM(ng)%Rindex,   &
     &                      Rec1)
           IF (FoundError(exit_flag, NoError, __LINE__,                 &
     &                    __FILE__)) RETURN
        END DO
!
!   Add the retrieved adjoint solution to the previous solution saved
!   in tl_var(1).
!
        add=.TRUE.
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_ADtoTL (ng, tile, Rec1, Rec1, add)
          END DO
!$OMP END PARALLEL
        END DO
!
      END IF
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        HIS(ng)%ncid=-1
        FWD(ng)%ncid=-1
#ifdef OBS_SPACE
        OBS(ng)%ncid=-1
#else
        ADS(ng)%ncid=-1
#endif
      END DO
!!
!! Compute and report model-observation comparison statistics.
!!
!!    DO ng=1,Ngrids
!!      CALL stats_modobs (ng)
!!      IF (FoundError(exit_flag, NoError, __LINE__,                    &
!!   &                 __FILE__)) RETURN
!!    END DO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Adjoint of 4D-PSAS to compute the observation sensitivity.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  Reset the start and end times for the adjoint forcing.
!
      DO ng=1,Ngrids
        str_day=tdays(ng)+ntimes(ng)*dt(ng)*sec2day
        end_day=tdays(ng)
        IF ((DstrS(ng).eq.0.0_r8).and.(DendS(ng).eq.0.0_r8)) THEN
          DstrS(ng)=end_day
          DendS(ng)=str_day
        END IF
        IF (Master) THEN
          WRITE (stdout,70) 'AD', DendS(ng), DstrS(ng)
        END IF
      END DO
!
!  WARNING: ONLY one outer loop can be used for this application.
!  =======  For more than 1 outer-loop, we require the second
!  derivative of each model operator (i.e. the tangent linear
!  of the tangent linear operator).
!
      AD_OUTER_LOOP : DO my_outer=Nimpact,Nimpact
        outer=my_outer
        inner=0
!
!-----------------------------------------------------------------------
!  Run the adjoint model initialized and forced by dI/dx where I is the
!  chosen function of the analysis/forecast state x.
!-----------------------------------------------------------------------
!
!  Reset DSTART and RunIterval for the analysis cycle interval using
!  NTIMES_ANA.
!
        dstart=dstartS
!
        rtime=0.0_r8
        DO ng=1,Ngrids
          rtime=MAX(rtime, dt(ng)*ntimes_ana(ng))
          ntimes(ng)=ntimes_ana(ng)
        END DO
        RunInterval=rtime
!
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX || \
    defined ADJUST_BOUNDARY
!
!  Call initial again to reset the forcing and obc adjustment time
!  counters.
!
!$OMP PARALLEL
        CALL initial
!$OMP END PARALLEL
#endif
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
          BLK(ng)%name=FWD(ng)%name
#ifdef OBS_SPACE
          OBS(ng)%name=TRIM(OBS(ng)%base)
#endif
        END DO
        IF (Master) THEN
          WRITE (stdout,50)
        END IF
!
!  Initialize the adjoint model: initialize using dI/dxf is
!  appropriate.
!
        DO ng=1,Ngrids
          Lstiffness=.FALSE.
!!AMM
#ifdef OBS_SPACE
          LsenPSAS(ng)=.TRUE.
          LsenFCT(ng)=.FALSE.
          Lobspace(ng)=.FALSE.
# ifndef OBS_IMPACT
          LadjVAR(ng)=.FALSE.
# endif
#else
          LsenPSAS(ng)=.FALSE.
          LsenFCT(ng)=.FALSE.
#endif
!!AMM
!$OMP PARALLEL
          CALL ad_initial (ng)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!    Copy tl_var(1) solution to ad_var(1).
!    It is important to clear the tl surface forcing and boundary
!    arrays first since only the initial fields must be used to
!    initialize the next run of the adjoint.
!
        add=.FALSE.
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
            CALL initialize_forces (ng, tile, iTLM)
#endif
#ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#endif
            CALL load_TLtoAD (ng, tile, Rec1, Rec1, add)
          END DO
!$OMP END PARALLEL
        END DO
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX || \
    defined ADJUST_BOUNDARY
!
!  Clear the adjoint forcing and open boundary arrays
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iADM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iADM)
# endif
          END DO
!$OMP END PARALLEL
        END DO
#endif
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
        DO ng=1,Ngrids
          WRTforce=.TRUE.
          IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
          Fcount=ADM(ng)%Fcount
          ADM(ng)%Nrec(Fcount)=0
          ADM(ng)%Rindex=0
        END DO
!
!  Time-step adjoint model backwards.
!  ??? What do we do in the case of model error? Save forcing for TLM?
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL ad_main3d (RunInterval)
#else
        CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
        DO ng=1,Ngrids
          CALL ad_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
        DO ng=1,Ngrids
          WRTforce(ng)=.FALSE.
          LwrtState2d(ng)=.FALSE.
          CALL ad_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Convolve adjoint trajectory with error covariances.
!
        Lposterior=.FALSE.
        CALL error_covariance (iTLM, driver, outer, inner,              &
     &                         Lbck, Lini, Lold, Lnew,                  &
     &                         Rec1, Rec2, Lposterior)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
!  AMM: Do not know what to do in the weak constraint case yet.
!
!!      IF (Master) THEN
!!        WRITE (stdout,40) outer, inner
!!      END IF
!!      DO ng=1,Ngrids
!!        TLF(ng)%Rindex=0
#ifdef DISTRIBUTE
!!        tile=MyRank
#else
!!        tile=-1
#endif
!!        CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
!!        IF (FoundError(exit_flag, NoError, __LINE__,                  &
!!   &                   __FILE__)) RETURN
!!      END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
#ifdef OBS_SPACE
!
!  Reinitialize fourdvar arrays using the analysis interval observation
!  file.
!
        DO ng=1,Ngrids
!!AMM     WRITE (OBS(ng)%name,90) TRIM(OBS(ng)%base)
          WRITE (OBS(ng)%name,80) TRIM(OBS(ng)%base), charC
        END DO
!
        CALL deallocate_fourdvar
!
        CALL initialize_fourdvar
!
!-----------------------------------------------------------------------
!  If the required vectors and arrays from congrad from a previous
!  run of the assimilation cycle are available, read them here from
!  LCZ(ng)%name NetCDF file.
!-----------------------------------------------------------------------
!
      SourceFile=__FILE__ // ", ROMS_initialize"
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_beta',        &
     &                        cg_beta)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_delta',       &
     &                        cg_delta)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_Gnorm_v',     &
     &                        cg_Gnorm_v)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_dla',         &
     &                        cg_dla)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'cg_QG',          &
     &                        cg_QG)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'zgrad0',         &
     &                        zgrad0)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'zcglwk',         &
     &                        zcglwk)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'TLmodVal_S',     &
     &                        TLmodVal_S,                               &
     &                        broadcast = .FALSE.)   ! Master use only
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# ifdef RPCG
        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'Hbk',            &
     &                        Hbk)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'Jb0',            &
     &                        Jb0)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, 'vcglwk',         &
     &                        vcglwk)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif
      END DO
!
!-----------------------------------------------------------------------
!  If skiping runing nonlinear model, read in observation screening and
!  quality control flag.
!-----------------------------------------------------------------------
!
      SourceFile=__FILE__ // ", ROMS_initialize"
      wrtObsScale(1:Ngrids)=.FALSE.
      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name, Vname(1,idObsS),  &
     &                        ObsScale)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Define output 4DVAR NetCDF file containing all processed data
!  at observation locations.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO

#endif
!
        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.TRUE.
          LwrtTLM(ng)=.FALSE.
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
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
!$OMP PARALLEL
          CALL tl_initial (ng)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (HMBM'H')_n * PSI at observation points
!  which are used in the conjugate gradient algorithm.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
#else
        CALL tl_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
        END DO

#ifdef OBS_IMPACT
!
!  Compute observation impact to the data assimilation system.
!
        DO ng=1,Ngrids
# ifdef RPCG
          CALL rep_matrix (ng, iTLM, outer, Ninner-1)
# else
          CALL rep_matrix (ng, iTLM, outer, Ninner)
# endif
        END DO
#else
!
!  Set basic state trajectory for adjoint inner-loops.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Set obs file flag for next call to ad_initial.
!
        DO ng=1,Ngrids
          OBS(ng)%ncid=-1
        END DO
!
# ifdef RPCG
        AD_INNER_LOOP : DO my_inner=Ninner,0,-1
# else
        AD_INNER_LOOP : DO my_inner=Ninner,1,-1
# endif
          inner=my_inner

# ifdef RPCG
!
!  Retrieve NLmodVal when inner=0 for use as BCKmodVal.
!
          IF (inner.eq.0) THEN
             DO ng=1,Ngrids
              CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,             &
     &                               'NLmodel_value', NLmodVal)
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                 __FILE__)) RETURN
            END DO
          END IF
          IF (inner.ne.Ninner) THEN
            Linner=.TRUE.
          ELSE
            Linner=.FALSE.
          END IF
# endif

          IF (Master) THEN
            WRITE (stdout,60) 'Adjoint of', uppercase('w4dpsas'),       &
     &                        outer, inner
          END IF
# ifdef RPCG
!
          INNER_COMPUTE : IF (Linner) THEN
!
# else
!
!  Call adjoint conjugate gradient algorithm.
!
          Lcgini=.FALSE.
          DO ng=1,Ngrids
            CALL ad_congrad (ng, iTLM, outer, inner, Ninner, Lcgini)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
# endif
!
!  Initialize the adjoint model from rest.
!
          DO ng=1,Ngrids
            LsenPSAS(ng)=.FALSE.
            LsenFCT(ng)=.TRUE.
#ifdef OBS_SPACE
            Lobspace(ng)=.TRUE.
# ifndef OBS_IMPACT
            LadjVAR(ng)=.TRUE.
# endif
#endif
!$OMP PARALLEL
            CALL ad_initial (ng)
!$OMP END PARALLEL
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
            wrtMisfit(ng)=.FALSE.
          END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file only once to avoid opening too many files.
!
          WRTforce=.TRUE.
          DO ng=1,Ngrids
            IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
            Fcount=ADM(ng)%Fcount
            ADM(ng)%Nrec(Fcount)=0
            ADM(ng)%Rindex=0
          END DO
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
            END IF
          END DO

!$OMP PARALLEL
# ifdef SOLVE3D
          CALL ad_main3d (RunInterval)
# else
          CALL ad_main2d (RunInterval)
# endif
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
          DO ng=1,Ngrids
            CALL ad_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
          DO ng=1,Ngrids
            WRTforce(ng)=.FALSE.
            LwrtState2d(ng)=.FALSE.
            CALL ad_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
!
!  Convolve adjoint trajectory with error covariances.
!
          Lposterior=.FALSE.
          CALL error_covariance (iTLM, driver, outer, inner,            &
     &                           Lbck, Lini, Lold, Lnew,                &
     &                           Rec1, Rec2, Lposterior)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
!!
!! AMM: Do not know what to do in the weak constraint case yet.
!!
!!        IF (Master) THEN
!!          WRITE (stdout,40) outer, inner
!!        END IF
!!        DO ng=1,Ngrids
!!          TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
!!          tile=MyRank
# else
!!          tile=-1
# endif
!!          CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
!!          IF (FoundError(exit_flag, NoError, __LINE__,                &
!!   &                     __FILE__)) RETURN
!!        END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
!  Initialize tangent linear model from ITL(ng)%name, record 1.
!
          DO ng=1,Ngrids
            ITL(ng)%Rindex=Rec1
!$OMP PARALLEL
            CALL tl_initial (ng)
!$OMP END PARALLEL
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
!
!  Set tangent linear history NetCDF parameters.  Define tangent linear
!  history file at the beggining of each inner loop  to avoid opening
!  too many NetCDF files.
!
          DO ng=1,Ngrids
            IF (inner.gt.Ninner) LdefTLM(ng)=.FALSE.
            Fcount=TLM(ng)%Fcount
            TLM(ng)%Nrec(Fcount)=0
            TLM(ng)%Rindex=0
          END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
            END IF
          END DO

!$OMP PARALLEL
# ifdef SOLVE3D
          CALL tl_main3d (RunInterval)
# else
          CALL tl_main2d (RunInterval)
# endif
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          DO ng=1,Ngrids
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.FALSE.
          END DO
# ifdef RPCG
          END IF INNER_COMPUTE
!
          DO ng=1,Ngrids
            CALL ad_rpcg_lanczos (ng, iRPM, outer, inner, Ninner,       &
     &                            Lcgini)
          END DO
# endif

        END DO AD_INNER_LOOP
# ifndef RPCG
!
!  Call adjoint conjugate gradient algorithm.
!
        inner=0
        Lcgini=.TRUE.
        DO ng=1,Ngrids
          CALL ad_congrad (ng, iTLM, outer, inner, Ninner, Lcgini)
        END DO
# endif

#endif /* !OBS_IMPACT */

#ifdef OBS_IMPACT
!
!  Write out total observation impact.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'ObsImpact_total', ad_ObsVal,           &
# ifdef IMPACT_INNER
     &                          (/1,1/), (/Mobs,Ninner/),               &
# else
     &                          (/1/), (/Mobs/),                        &
# endif
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
#else
!
!  Write out observation sensitivity.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,                 &
     &                          'ObsSens_total', ad_ObsVal,             &
# ifdef IMPACT_INNER
     &                          (/1,1/), (/Mobs,Ninner/),               &
# else
     &                          (/1/), (/Mobs/),                        &
# endif
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
#endif
!
!  Close tangent linear NetCDF file.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
          TLM(ng)%ncid=-1
        END DO

#if defined OBS_IMPACT && defined OBS_IMPACT_SPLIT
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with initial condition increments
!  only to compute the observation impact associated with the initial
!  conditions.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.TRUE.
          LwrtTLM(ng)=.FALSE.
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero after
!  running the representer model and must be zero when running the
!  tangent linear model.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
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
!  Initialize tangent linear model from ITL(ng)%name, record 1.
!
        DO ng=1,Ngrids
          ITL(ng)%Rindex=Rec1
!$OMP PARALLEL
          CALL tl_initial (ng)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Clear tangent linear forcing arrays and boundary arrays
!  before the obs impact initial condition calculation.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (HMBM'H')_n * PSI at observation points
!  which are used in the conjugate gradient algorithm.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
# ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
# else
        CALL tl_main2d (RunInterval)
# endif
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
        END DO
!
!  Compute observation impact to the data assimilation system.
!
        DO ng=1,Ngrids
# ifdef RPCG
          CALL rep_matrix (ng, iTLM, outer, Ninner-1)
# else
          CALL rep_matrix (ng, iTLM, outer, Ninner)
# endif
        END DO
!
!  Write out observation sentivity.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,                 &
     &                          'ObsImpact_IC', ad_ObsVal,              &
# ifdef IMPACT_INNER
     &                          (/1,1/), (/Mobs,Ninner/),               &
# else
     &                          (/1/), (/Mobs/),                        &
# endif
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with surface forcing increments
!  only to compute the observation impact associated with the surface
!  forcing.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.TRUE.
          LwrtTLM(ng)=.FALSE.
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero after
!  running the representer model and must be zero when running the
!  tangent linear model.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
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
!  Initialize tangent linear model from ITL(ng)%name, record 1.
!
        DO ng=1,Ngrids
          ITL(ng)%Rindex=Rec1
!$OMP PARALLEL
          CALL tl_initial (ng)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Clear tangent initial condition arrays and boundary arrays
!  before the obs impact initial condition calculation.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (HMBM'H')_n * PSI at observation points
!  which are used in the conjugate gradient algorithm.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#  ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
#  else
        CALL tl_main2d (RunInterval)
#  endif
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
        END DO
!
!  Compute observation impact to the data assimilation system.
!
        DO ng=1,Ngrids
# ifdef RPCG
          CALL rep_matrix (ng, iTLM, outer, Ninner-1)
# else
          CALL rep_matrix (ng, iTLM, outer, Ninner)
# endif
        END DO
!
!  Write out observation sentivity.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,                 &
     &                          'ObsImpact_FC', ad_ObsVal,              &
# ifdef IMPACT_INNER
     &                          (/1,1/), (/Mobs,Ninner/),               &
# else
     &                          (/1/), (/Mobs/),                        &
# endif
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
# endif

# if defined ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with boundary condition increments
!  only to compute the obs impact associated with the boundary
!  conditions.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.TRUE.
          LwrtTLM(ng)=.FALSE.
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero after
!  running the representer model and must be zero when running the
!  tangent linear model.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
!$OMP END PARALLEL
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
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
!$OMP PARALLEL
          CALL tl_initial (ng)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Clear tangent linear initial condition and forcing arrays
!  before the obs impact initial condition calculation.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
            CALL initialize_forces (ng, tile, iTLM)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (HMBM'H')_n * PSI at observation points
!  which are used in the conjugate gradient algorithm.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#  ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
#  else
        CALL tl_main2d (RunInterval)
#  endif
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        DO ng=1,Ngrids
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
        END DO
!
!  Compute observation impact to the data assimilation system.
!
        DO ng=1,Ngrids
#  ifdef RPCG
          CALL rep_matrix (ng, iTLM, outer, Ninner-1)
#  else
          CALL rep_matrix (ng, iTLM, outer, Ninner)
#  endif
        END DO
!
!  Write out observation sentivity.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,                 &
     &                          'ObsImpact_BC', ad_ObsVal,              &
#  ifdef IMPACT_INNER
     &                          (/1,1/), (/Mobs,Ninner/),               &
#  else
     &                          (/1/), (/Mobs/),                        &
#  endif
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
# endif
#endif /* OBS_IMPACT_SPLIT */
!
!  Close current forward NetCDF file.
!
        SourceFile=__FILE__ // ", ROMS_run"
        DO ng=1,Ngrids
          CALL netcdf_close (ng, iNLM, FWD(ng)%ncid)
          FWD(ng)%ncid=-1
        END DO

      END DO AD_OUTER_LOOP
!
 10   FORMAT (a,'_',i3.3,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (' (',i3.3,',',i3.3,'): ',a,' data penalty, Jdata = ',     &
     &        1p,e17.10,0p,t68,a)
 40   FORMAT (/,' Converting Convolved Adjoint Trajectory to',          &
     &          ' Impulses: Outer = ',i3.3,' Inner = ',i3.3,/)
 50   FORMAT (/,'ROMS/TOMS: Started adjoint Sensitivity calculation',   &
     &          ' ...',/)
 60   FORMAT (/,'ROMS/TOMS: ',a,1x,a,', Outer = ',i3.3,                 &
     &          ' Inner = ',i3.3,/)
 70   FORMAT (/,1x,a,1x,'ROMS/TOMS: adjoint forcing time range: ',      &
     &        f12.4,' - ',f12.4 ,/)
 80   FORMAT (a,'_',a,'.nc')
 90   FORMAT (a,'.nc')

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear, tangent linear, and    !
!  adjoint models execution.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: Fcount, ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into record 3.
!
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
            Fcount=RST(ng)%Fcount
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
            CALL wrt_rst (ng)
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed CPU time (seconds):',/)
      END IF
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0, __LINE__, __FILE__)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Report dynamic memory and automatic memory requirements.
!
!$OMP PARALLEL
      CALL memory
!$OMP END PARALLEL
!
!  Close IO files.
!
      CALL close_out

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
