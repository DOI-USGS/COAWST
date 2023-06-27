      MODULE roms_kernel_mod
!
!git $Id$
!svn $Id: obs_sen_rbl4dvar_analysis.h 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Strong/Weak Constraint 4-Dimensional Variational Data     !
!    Assimilation and Observation Sensitivity Driver: Restricted       !
!    B-preconditioned Lanczos (RBL4D-Var).                             !
!                                                                      !
!  This driver is used for the dual formulation (observation space),   !
!  strong or weak constraint 4D-Var where errors may be considered     !
!  in both model and observations.                                     !
!                                                                      !
!  It also computes the sensitivity of the assimilation system to      !
!  each observation. It measures the degree to which each observation  !
!  contributes to the uncertainty in the estimate. This analysis can   !
!  be used to determine the type of measurements that need to be made, !
!  where to observe, and when.                                         !
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
      USE mod_param
      USE mod_parallel
      USE mod_arrays
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
      USE mod_stepping
!
#ifdef ADJUST_BOUNDARY
      USE mod_boundary,      ONLY : initialize_boundary
#endif
      USE mod_forces,        ONLY : initialize_forces
      USE mod_ocean,         ONLY : initialize_ocean
!
      USE ad_wrt_his_mod,    ONLY : ad_wrt_his
      USE close_io_mod,      ONLY : close_file, close_inp, close_out
#ifndef RPCG
      USE congrad_mod,       ONLY : congrad
#endif
      USE convolve_mod,      ONLY : error_covariance
      USE def_impulse_mod,   ONLY : def_impulse
      USE def_mod_mod,       ONLY : def_mod
      USE def_norm_mod,      ONLY : def_norm
      USE get_state_mod,     ONLY : get_state
      USE inp_par_mod,       ONLY : inp_par
      USE normalization_mod, ONLY : normalization
#ifdef MCT_LIB
# ifdef ATM_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAV_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
#ifdef RPCG
      USE rpcg_lanczos_mod,  ONLY : rpcg_lanczos
#endif
      USE stats_modobs_mod,  ONLY : stats_modobs
      USE strings_mod,       ONLY : FoundError, uppercase
      USE tl_def_ini_mod,    ONLY : tl_def_ini
      USE wrt_ini_mod,       ONLY : wrt_ini
      USE wrt_rst_mod,       ONLY : wrt_rst
#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod,  ONLY : balance_ref, biconj
#endif
!
      implicit none
!
      PUBLIC :: ROMS_initialize
      PUBLIC :: ROMS_run
      PUBLIC :: ROMS_finalize
!
      CONTAINS
!
      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first
!
      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.
!
#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: STDrec, Tindex
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_initialize"

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (mpi) world communictor.
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
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
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iNLM, 0, __LINE__, MyFile)
          END DO
        END DO
!
!  Allocate and initialize modules variables.
!
        CALL ROMS_allocate_arrays (allocate_vars)
        CALL ROMS_initialize_arrays
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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

#if !defined RECOMPUTE_4DVAR
!
!-----------------------------------------------------------------------
!  If the required vectors and arrays from congrad from a previous
!  run of the assimilation cycle are available, read them here from
!  LCZ(ng)%name NetCDF file.
!-----------------------------------------------------------------------
!
      SourceFile=MyFile
      DO ng=1,Ngrids
        SELECT CASE (LCZ(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'cg_beta', cg_beta)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'cg_delta', cg_delta)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'cg_Gnorm_v', cg_Gnorm_v)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'cg_dla', cg_dla)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'cg_QG', cg_QG)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'zgrad0', zgrad0)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'zcglwk', zcglwk)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'TLmodVal_S', TLmodVal_S,             &
     &                            broadcast = .FALSE.) ! Master use only
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef RPCG
            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'Hbk', Hbk)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'Jb0', Jb0)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            'vcglwk', vcglwk)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'cg_beta', cg_beta)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'cg_delta', cg_delta)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'cg_Gnorm_v', cg_Gnorm_v)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'cg_dla', cg_dla)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'cg_QG', cg_QG)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'zgrad0', zgrad0)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'zcglwk', zcglwk)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'TLmodVal_S', TLmodVal_S)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef RPCG
            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'Hbk', Hbk)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'Jb0', Jb0)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                'vcglwk', vcglwk)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
        END SELECT
      END DO
#endif

#ifdef SKIP_NLM
!
!-----------------------------------------------------------------------
!  If skiping runing nonlinear model, read in observation screening and
!  quality control flag.
!-----------------------------------------------------------------------
!
      SourceFile=MyFile
      wrtObsScale(1:Ngrids)=.FALSE.
      DO ng=1,Ngrids
        SELECT CASE (LCZ(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,               &
     &                            Vname(1,idObsS), ObsScale)

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_fvar (ng, iTLM, LCZ(ng)%name,           &
     &                                Vname(1,idObsS), ObsScale)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        CALL get_state (ng, 10, 10, STD(1,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      IF (NSA.eq.2) THEN
        DO ng=1,Ngrids
          CALL get_state (ng, 11, 11, STD(2,ng), STDrec, Tindex)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 12, 12, STD(3,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 13, 13, STD(4,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#endif
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS nonlinear, tangent linear and     !
!  adjoint models.                                                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      real(dp), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: Lcgini, Linner, Lposterior
!
      integer :: my_inner, my_outer
      integer :: Lbck, Lini, Rec1, Rec2
      integer :: i, lstr, ng, status, tile
      integer :: Fcount, NRMrec

      integer, dimension(Ngrids) :: indxSave
      integer, dimension(Ngrids) :: Nrec
!
      real(r8) :: str_day, end_day
!
      character (len=25) :: driver
      character (len=20) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
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
      driver='obs_sen_rbl4dvar_analysis'
!
!-----------------------------------------------------------------------
!  Configure weak constraint RBL4D-Var algorithm.
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
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lbck" of INI(ng)%name NetCDF file. The
!  record "Lbck" becomes the background state record and the record
!  "Lini" becomes current nonlinear initial conditions.
!
      DO ng=1,Ngrids
        INI(ng)%Rindex=1
        Fcount=INI(ng)%load
        INI(ng)%Nrec(Fcount)=1
#ifdef DISTRIBUTE
        CALL wrt_ini (ng, MyRank, 1)
#else
        CALL wrt_ini (ng, -1, 1)
#endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Set nonlinear output history file as the initial basic state
!  trajectory for the specified outer loop (Nimpact) used in the
!  observation impact/sensitivity analysis.
!
      DO ng=1,Ngrids
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), Nimpact-1
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
      END DO
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          IF (NSA.eq.2) THEN
            CALL def_norm (ng, iNLM, 2)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
          CALL def_norm (ng, iNLM, 3)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL def_norm (ng, iNLM, 4)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

          DO tile=first_tile(ng),last_tile(ng),+1
            CALL normalization (ng, tile, 2)
          END DO
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

#ifdef ADJUST_BOUNDARY
          CALL get_state (ng, 16, 16, NRM(3,ng), NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL get_state (ng, 17, 17, NRM(4,ng), NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
        END IF
      END DO

#if !defined RECOMPUTE_4DVAR && defined BALANCE_OPERATOR && \
     defined ZETA_ELLIPTIC
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
#endif
!
!  Define tangent linear initial conditions file.
!
      DO ng=1,Ngrids
        LdefITL(ng)=.TRUE.
        CALL tl_def_ini (ng)
        LdefITL(ng)=.FALSE.
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Define impulse forcing NetCDF file.
!
      DO ng=1,Ngrids
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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

#ifdef SKIP_NLM
!
!  If skipping NLM, write out screening flag to 4DVAR NetCDF file.
!
      DO ng=1,Ngrids
        SELECT CASE (DAV(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,               &
     &                            Vname(1,idObsS), ObsScale,            &
     &                            (/1/), (/Ndatum(ng)/),                &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idObsS))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, DAV(ng)%name,           &
     &                                Vname(1,idObsS), ObsScale,        &
     &                                (/1/), (/Ndatum(ng)/),            &
     &                               pioFile = DAV(ng)%pioFile,         &
     &                               pioVar = DAV(ng)%pioVar(idObsS)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#endif
!
!  Write out outer loop beeing processed.
!
      SourceFile=MyFile
      DO ng=1,Ngrids
        SELECT CASE (DAV(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name,               &
     &                            'Nimpact', Nimpact,                   &
     &                            (/0/), (/0/),                         &
     &                            ncid = DAV(ng)%ncid)

#if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_ivar (ng, iNLM, DAV(ng)%name,           &
     &                                'Nimpact', Nimpact,               &
     &                                (/0/), (/0/),                     &
     &                                pioFile = DAV(ng)%pioFile)
#endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO

#ifndef SKIP_NLM
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute background state trajectory, X_n-1(t)
!  and the background values at the observation points and times. It
!  processes and writes the observations accept/reject flag (ObsScale)
!  once to allow background quality control, if any.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      DO ng=1,Ngrids
        wrtObsScale(ng)=.TRUE.
        SporadicImpulse(ng)=.FALSE.
        FrequentImpulse(ng)=.FALSE.
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
        END IF
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
              WRITE (stdout,30) outer, inner, 'NLM',                    &
     &                          FOURDVAR(ng)%NLPenalty(i),              &
     &                          TRIM(string)
            END IF
          END DO
        END IF
!
!  Write out initial data penalty function to NetCDF file.
!
        SourceFile=MyFile
        SELECT CASE (DAV(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,               &
     &                            'NL_iDataPenalty',                    &
     &                            FOURDVAR(ng)%NLPenalty(0:),           &
     &                            (/1/), (/NobsVar(ng)+1/),             &
     &                            ncid = DAV(ng)%ncid)

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iNLM, DAV(ng)%name,           &
     &                                'NL_iDataPenalty',                &
     &                                FOURDVAR(ng)%NLPenalty(0:),       &
     &                                (/1/), (/NobsVar(ng)+1/),         &
     &                                pioFile = DAV(ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Clean penalty array before next run of NL model.
!
        FOURDVAR(ng)%NLPenalty=0.0_r8
      END DO

#endif /* !SKIP_NLM */
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

#ifdef FORWARD_FLUXES
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
#endif

#ifdef RECOMPUTE_4DVAR
!
!-----------------------------------------------------------------------
!  Solve the system (following Courtier, 1997):
!
!              (H M_n B (M_n)' H' + Cobs) * w_n = d_n
!
!              d_n = yo - H * Xb_n
!
!  where M_n is the tangent linear model matrix, Cobs is the
!  observation-error covariance, B is the background error covariance
!  and dx_n=B M' H' w_n is the analysis increment so that Xa=Xb+dx_n.
!  d_n is the misfit between observations (yo) and model (H * Xb_n),
!  and H is the linearized observation operator.
!
!  Here, _n denotes a sequence of outer-loop estimates.
!
!  The system does not need to be solved explicitly by inverting the
!  symmetric matrix, P_n:
!
!              P_n = H M_n B (M_n)' H' + Cobs
!
!  but by computing the action of P_n on any vector PSI, such that
!
!              P_n * PSI =  H M_n B (M_n)' H' * PSI + Cobs * PSI
!
!  The (H M_n B (M_n)' H') matrix is not explicitly computed but
!  evaluated by one integration backward of the adjoint model and
!  one integration forward of the tangent linear model for any
!  forcing vector PSI.
!
!  A preconditioned conjugate gradient algorithm is used to compute
!  an approximation PSI for w_n.
!
!-----------------------------------------------------------------------
!
!  If the required vectors and arrays from congrad from a previous run
!  of the assimilation cycle are not available, rerun the 4D-Var cycle.
!
      OUTER_LOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0
!
!  Set basic state trajectory (X_n-1) file to previous outer loop file
!  (outer-1).
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero and must
!  be zero when running the tangent linear model.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
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
        INNER_LOOP : DO my_inner=0,Ninner
          inner=my_inner
!
!  Initialize conjugate gradient algorithm depending on hot start or
!  outer loop index.
!
          IF (inner.eq.0) THEN
            Lcgini=.TRUE.
            DO ng=1,Ngrids
              CALL congrad (ng, iRPM, outer, inner, Ninner, Lcgini)
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
!
!  Start inner loop computations.
!
          INNER_COMPUTE : IF (Linner) THEN
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate adjoint model forced with any vector PSI at the observation
!  locations and generate adjoint trajectory, Lambda_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model from rest.
!
            DO ng=1,Ngrids
              Lsen4DVAR(ng)=.FALSE.
              CALL ad_initial (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              wrtMisfit(ng)=.FALSE.
            END DO
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
                WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
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
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
            DO ng=1,Ngrids
# ifdef DISTRIBUTE
              CALL ad_wrt_his (ng, MyRank)
# else
              CALL ad_wrt_his (ng, -1)
# endif
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
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
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!  Convolve adjoint trajectory with error covariances.
!
            Lposterior=.FALSE.
            CALL error_covariance (iTLM, driver, outer, inner,          &
     &                             Lbck, Lini, Lold, Lnew,              &
     &                             Rec1, Rec2, Lposterior)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
            IF (Master) THEN
              WRITE (stdout,40) outer, inner
            END IF
            DO ng=1,Ngrids
              TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
              CALL wrt_impulse (ng, MyRank, iADM, ADM(ng)%name)
# else
              CALL wrt_impulse (ng, -1, iADM, ADM(ng)%name)
# endif
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
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
!  Initialize tangent linear model from ITL(ng)%name, record Rec1.
!
            DO ng=1,Ngrids
              ITL(ng)%Rindex=Rec1
              CALL tl_initial (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!  Activate switch to write out initial misfit between model and
!  observations.
!
            IF ((outer.eq.1).and.(inner.eq.1)) THEN
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
                WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
              END IF
            END DO
!
# ifdef SOLVE3D
            CALL tl_main3d (RunInterval)
# else
            CALL tl_main2d (RunInterval)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            DO ng=1,Ngrids
              wrtNLmod(ng)=.FALSE.
              wrtTLmod(ng)=.FALSE.
            END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Use conjugate gradient algorithm to find a better approximation
!  PSI to coefficients Beta_n.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            Nrun=Nrun+1
            DO ng=1,Ngrids
              Lcgini=.FALSE.
              CALL congrad (ng, iTLM, outer, inner, Ninner, Lcgini)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO

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
          Lsen4DVAR(ng)=.FALSE.
          CALL ad_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
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
!  Time-step adjoint model backwards forced with estimated coefficients,
!  Beta_n.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
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
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          CALL ad_wrt_his (ng, MyRank)
# else
          CALL ad_wrt_his (ng, -1)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
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
        Lposterior=.FALSE.
        CALL error_covariance (iNLM, driver, outer, inner,              &
     &                         Lbck, Lini, Lold, Lnew,                  &
     &                         Rec1, Rec2, Lposterior)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
        IF (Master) THEN
          WRITE (stdout,40) outer, inner
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
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute a "new estimate" of the state
!  trajectory, X_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set new basic state trajectory for next outer loop.
!
        DO ng=1,Ngrids
          LdefHIS(ng)=.TRUE.
          LwrtHIS(ng)=.TRUE.
          wrtNLmod(ng)=.TRUE.
          wrtTLmod(ng)=.FALSE.
          WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), outer
          lstr=LEN_TRIM(HIS(ng)%name)
          HIS(ng)%base=HIS(ng)%name(1:lstr-3)
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
!  Clear tangent arrays before running nonlinear model (important).
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
        END DO
!
!  Initialize nonlinear model INI(ng)%name file, record outer+2.
!  Notice that NetCDF record index counter is saved because this
!  counter is used to write initial conditions.
!
        DO ng=1,Ngrids
          indxSave(ng)=INI(ng)%Rindex
          INI(ng)%Rindex=outer+2
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
        IF (outer.eq.Nouter) THEN
          DO ng=1,Ngrids
            wrtMisfit(ng)=.TRUE.
          END DO
        END IF
!
!  Run nonlinear forced by convolved adjoint trajectory impulses and
!  compute new basic state trajectory X_n.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
          END IF
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
          wrtTLmod(ng)=.FALSE.
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
                WRITE (stdout,30) outer, inner, 'NLM',                  &
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
     &                              (/1,outer/),                        &
     &                              (/NobsVar(ng)+1,1/),                &
     &                              ncid = DAV(ng)%ncid)

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iNLM, DAV(ng)%name,         &
     &                                  'NL_fDataPenalty',              &
     &                                  FOURDVAR(ng)%NLPenalty(0:),     &
     &                                  (/1,outer/),                    &
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
!
!  Close current forward NetCDF file.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          CALL close_file (ng, iNLM, FWD(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          HIS(ng)%ncid=-1
        END DO

      END DO OUTER_LOOP

#endif /* RECOMPUTE_4DVAR */
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
!       HIS(ng)%ncid=-1
      END DO
!!
!! Compute and report model-observation comparison statistics.
!!
!!    DO ng=1,Ngrids
!!      CALL stats_modobs (ng)
!!      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!!    END DO
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!  Adjoint of RBL4D-Var to compute the observation sensitivity.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Reset the start and end times for the adjoint forcing.
!
      DO ng=1,Ngrids
#ifdef SKIP_NLM
        str_day=tdays(ng)+ntimes(ng)*dt(ng)*sec2day
        end_day=tdays(ng)
#else
        str_day=tdays(ng)
        end_day=str_day-ntimes(ng)*dt(ng)*sec2day
#endif
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
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
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
          Lsen4DVAR(ng)=.TRUE.
          CALL ad_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
        DO ng=1,Ngrids
          WRTforce=.TRUE.
          IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
          Fcount=ADM(ng)%load
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
!
#ifdef SOLVE3D
        CALL ad_main3d (RunInterval)
#else
        CALL ad_main2d (RunInterval)
#endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
        DO ng=1,Ngrids
#ifdef DISTRIBUTE
          CALL ad_wrt_his (ng, MyRank)
#else
          CALL ad_wrt_his (ng, -1)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
        DO ng=1,Ngrids
          WRTforce(ng)=.FALSE.
#ifdef DISTRIBUTE
          CALL ad_wrt_his (ng, MyRank)
#else
          CALL ad_wrt_his (ng, -1)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Convolve adjoint trajectory with error covariances.
!
        Lposterior=.FALSE.
        CALL error_covariance (iTLM, driver, outer, inner,              &
     &                         Lbck, Lini, Lold, Lnew,                  &
     &                         Rec1, Rec2, Lposterior)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!!        CALL wrt_impulse (ng, MyRank, iADM, ADM(ng)%name)
#else
!!        CALL wrt_impulse (ng, -1, iADM, ADM(ng)%name)
#endif
!!        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!!      END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#endif
          END DO
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
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
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (HMBM'H')_n * PSI at observation points
!  which are used in the conjugate gradient algorithm.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO
!
#ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
#else
        CALL tl_main2d (RunInterval)
#endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
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
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
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
              SELECT CASE (DAV(ng)%IOtype)
                CASE (io_nf90)
                  CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,         &
     &                                  'NLmodel_value', NLmodVal)

#  if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  CALL pio_netcdf_get_fvar (ng, iTLM, DAV(ng)%name,     &
     &                                      'NLmodel_value', NLmodVal)
#  endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
          END IF
          IF (inner.ne.Ninner) THEN
            Linner=.TRUE.
          ELSE
            Linner=.FALSE.
          END IF
# endif

          IF (Master) THEN
            WRITE (stdout,60) 'Adjoint of', uppercase('rbl4dvar'),      &
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
            IF (FoundError(exit_flag, NoError,                          &
     &                     __LINE__, MyFile)) RETURN
          END DO
# endif
!
!  Initialize the adjoint model from rest.
!
          DO ng=1,Ngrids
            Lsen4DVAR(ng)=.FALSE.
            CALL ad_initial (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            wrtMisfit(ng)=.FALSE.
          END DO
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
              WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
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
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
          DO ng=1,Ngrids
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
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
          Lposterior=.FALSE.
          CALL error_covariance (iTLM, driver, outer, inner,            &
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
!!
!! AMM: Do not know what to do in the weak constraint case yet.
!!
!!        IF (Master) THEN
!!          WRITE (stdout,40) outer, inner
!!        END IF
!!        DO ng=1,Ngrids
!!          TLF(ng)%Rindex=0
# ifdef DISTRIBUTE
!!          CALL wrt_impulse (ng, MyRank, iADM, ADM(ng)%name)
# else
!!          CALL wrt_impulse (ng, -1, iADM, ADM(ng)%name)
# endif
!!          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
            CALL tl_initial (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
!
!  Set tangent linear history NetCDF parameters.  Define tangent linear
!  history file at the beggining of each inner loop  to avoid opening
!  too many NetCDF files.
!
          DO ng=1,Ngrids
            IF (inner.gt.Ninner) LdefTLM(ng)=.FALSE.
            Fcount=TLM(ng)%load
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
        SourceFile=MyFile
        DO ng=1,Ngrids
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,             &
     &                              'ObsImpact_total', ad_ObsVal,       &
# ifdef IMPACT_INNER
     &                              (/1,1/), (/Mobs,Ninner/),           &
# else
     &                              (/1/), (/Mobs/),                    &
# endif
     &                              ncid = DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iNLM, DAV(ng)%name,         &
     &                                  'ObsImpact_total', ad_ObsVal,   &
#  ifdef IMPACT_INNER
     &                                  (/1,1/), (/Mobs,Ninner/),       &
#  else
     &                                  (/1/), (/Mobs/),                &
#  endif
     &                                  pioFile = DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL pio_netcdf_sync (ng, iNLM, DAV(ng)%name,             &
     &                              DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
          END SELECT
        END DO
#else
!
!  Write out observation sensitivity.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              'ObsSens_total', ad_ObsVal,         &
     &                              (/1/), (/Mobs/),                    &
     &                              ncid = DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, DAV(ng)%name,         &
     &                                  'ObsSens_total', ad_ObsVal,     &
     &                                  (/1/), (/Mobs/),                &
     &                                  pioFile = DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL pio_netcdf_sync (ng, iNLM, DAV(ng)%name,             &
     &                              DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
          END SELECT
        END DO
#endif
!
!  Close tangent linear NetCDF file.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          CALL close_file (ng, iTLM, TLM(ng))
          IF (FoundError(exit_flag, NoError, __LINE__,  MyFile)) RETURN
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
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
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
          CALL tl_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Clear tangent linear forcing arrays and boundary arrays
!  before the obs impact initial condition calculation.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
# endif
          END DO
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
        SourceFile=MyFile
        DO ng=1,Ngrids
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              'ObsImpact_IC', ad_ObsVal,          &
# ifdef IMPACT_INNER
     &                              (/1,1/), (/Mobs,Ninner/),           &
# else
     &                              (/1/), (/Mobs/),                    &
# endif
     &                              ncid = DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, DAV(ng)%name,         &
     &                                  'ObsImpact_IC', ad_ObsVal,      &
#  ifdef IMPACT_INNER
     &                                  (/1,1/), (/Mobs,Ninner/),       &
#  else
     &                                  (/1/), (/Mobs/),                &
#  endif
     &                                  pioFile = DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL pio_netcdf_sync (ng, iNLM, DAV(ng)%name,             &
     &                              DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
          END SELECT
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
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
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
          CALL tl_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  Clear tangent initial condition arrays and boundary arrays
!  before the obs impact initial condition calculation.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
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
        SourceFile=MyFile
        DO ng=1,Ngrids
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              'ObsImpact_FC', ad_ObsVal,          &
# ifdef IMPACT_INNER
     &                              (/1,1/), (/Mobs,Ninner/),           &
# else
     &                              (/1/), (/Mobs/),                    &
# endif
     &                              ncid = DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, DAV(ng)%name,         &
     &                                  'ObsImpact_FC', ad_ObsVal,      &
#  ifdef IMPACT_INNER
     &                                  (/1,1/), (/Mobs,Ninner/),       &
#  else
     &                                  (/1/), (/Mobs/),                &
#  endif
     &                                  pioFile = DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL pio_netcdf_sync (ng, iNLM, DAV(ng)%name,             &
     &                              DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
          END SELECT
        END DO
# endif

# if defined ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with boundary condition increments
!  only to compute the observation impact associated with the boundary
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
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#  endif
          END DO
        END DO
!
!  Set basic state trajectory.
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%head), outer-1
          lstr=LEN_TRIM(FWD(ng)%name)
          FWD(ng)%base=FWD(ng)%name(1:lstr-3)
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
!  Clear tangent linear initial condition and forcing arrays
!  before the obs impact initial condition calculation.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
            CALL initialize_forces (ng, tile, iTLM)
          END DO
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
        SourceFile=MyFile
        DO ng=1,Ngrids
          SELECT CASE (DAV(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              'ObsImpact_BC', ad_ObsVal,          &
#  ifdef IMPACT_INNER
     &                              (/1,1/), (/Mobs,Ninner/),           &
#  else
     &                              (/1/), (/Mobs/),                    &
#  endif
     &                              ncid = DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, DAV(ng)%name,         &
     &                                  'ObsImpact_BC', ad_ObsVal,      &
#   ifdef IMPACT_INNER
     &                                  (/1,1/), (/Mobs,Ninner/),       &
#   else
     &                                  (/1/), (/Mobs/),                &
#   endif
     &                                  pioFile = DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              CALL pio_netcdf_sync (ng, iNLM, DAV(ng)%name,             &
     &                              DAV(ng)%pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
          END SELECT
        END DO
# endif
#endif /* OBS_IMPACT_SPLIT */
!
!  Close current forward NetCDF file.
!
        SourceFile=MyFile
        DO ng=1,Ngrids
          CALL close_file (ng, iNLM, FWD(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (HIS(ng)%IOtype.eq.io_nf90) THEN
            HIS(ng)%ncid=-1
#if defined PIO_LIB && defined DISTRIBUTE
          ELSE IF (HIS(ng)%IOtype.eq.io_pio) THEN
            HIS(ng)%pioFile%fh=-1
#endif
          END IF
        END DO

      END DO AD_OUTER_LOOP
!
 10   FORMAT (a,'_outer',i0,'.nc')
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
!
      RETURN
      END SUBROUTINE ROMS_run
!
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear, tangent linear, and    !
!  adjoint models execution.                                           !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: Fcount, ng, thread
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_finalize"
!
!-----------------------------------------------------------------------
!  Read and write observation variables for completeness.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#ifdef DISTRIBUTE
        CALL stats_modobs (ng, MyRank)
#else
        CALL stats_modobs (ng, -1)
#endif
      END DO
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
            Fcount=RST(ng)%load
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
#ifdef DISTRIBUTE
            CALL wrt_rst (ng, MyRank)
#else
            CALL wrt_rst (ng, -1)
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks, report memory requirements, and
!  close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed wall CPU time for each process (seconds):',/)
      END IF
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0, __LINE__, MyFile)
        END DO
      END DO
!
!  Report dynamic memory and automatic memory requirements.
!
      CALL memory
!
!  Close IO files.
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iNLM)
      END DO
      CALL close_out
!
      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE roms_kernel_mod
