      MODULE roms_kernel_mod
!
!git $Id$
!git $Id$
!svn $Id: jedi_roms.h 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMSJEDI Interface: Joint Effort for Data-assimilation Integration  !
!                                                                      !
!  The routines in this driver control the initialization, time-       !
!  stepping, and finalization of ROMSJEDI interface following          !
!  ESMF/NUOPC conventions:                                             !
!                                                                      !
!     ROMS_initialize           Phase 1 ROMSJEDI Initialization        !
!     ROMS_initializeP2         Phase 2 ROMSJEDI Initialization        !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_arrays
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
!
      USE close_io_mod,      ONLY : close_inp, close_out
      USE dateclock_mod,     ONLY : time_string
#ifdef DISTRIBUTE
      USE distribute_mod,    ONLY : mp_bcasti
#endif
#ifdef WET_DRY
      USE get_wetdry_mod,    ONLY : get_wetdry
#endif
      USE ini_hmixcoef_mod,  ONLY : ini_hmixcoef
      USE inp_par_mod,       ONLY : inp_par
#ifdef NESTING
      USE nesting_mod,       ONLY : nesting
#endif
#ifdef SOLVE3D
      USE set_depth_mod,     ONLY : set_depth0, set_depth
      USE omega_mod,         ONLY : omega
      USE rho_eos_mod,       ONLY : rho_eos
      USE set_massflux_mod,  ONLY : set_massflux
#endif
#ifdef MASKING
      USE set_masks_mod,     ONLY : set_masks
#endif
      USE stiffness_mod,     ONLY : stiffness
      USE strings_mod,       ONLY : FoundError
#ifdef TANGENT
      USE tl_set_depth_mod,  ONLY : tl_bath
#endif
#ifdef WET_DRY
      USE wetdry_mod,        ONLY : wetdry
#endif
      USE wrt_rst_mod,       ONLY : wrt_rst
!
      implicit none
!
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_initializeP2
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize
!
      PRIVATE :: nlm_initial
#ifdef TANGENT
      PRIVATE :: tlm_initial
#endif
#ifdef ADJOINT
      PRIVATE :: adm_initial
#endif
!
      CONTAINS
!
      SUBROUTINE ROMS_initialize (first, mpiCOMM, kernel)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS state variables and     !
!  internal parameters. It reads standard input parameters.            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first
!
      integer, intent(in), optional :: mpiCOMM
      integer, intent(in), optional :: kernel
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: Phase, chunk_size, ng, thread, tile
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
!  Initialize counters.
!
        Nrun=1                ! run counter
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
#if defined DISTRIBUTE
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
!
!-----------------------------------------------------------------------
!  Initialize ROMS, phase 1.
!-----------------------------------------------------------------------
!
      Phase=1
!
      SELECT CASE (kernel)
        CASE (iNLM)
          CALL nlm_initial (Phase)
#ifdef TANGENT
        CASE (iTLM)
          CALL tlm_initial (Phase)
#endif
#ifdef ADJOINT
        CASE (iADM)
          CALL adm_initial (Phase)
#endif
      END SELECT
!
!-----------------------------------------------------------------------
!  Set ROMS application grid configuration. It is done once.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (SetGridConfig(ng)) THEN
          CALL set_grid (ng, iNLM)
          SetGridConfig(ng)=.FALSE.
        END IF
      END DO
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_initializeP2 (kernel)
!
!=======================================================================
!                                                                      !
!  ROMSJEDI phase 2 initialization. It requires the initial state      !
!  (set elsewhere) to complete the full initialization of the kernel.  !
!  It computes depths, density, and horizontal mass fluxes.            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in) :: kernel
!
!  Local variable declarations.
!
      integer :: Phase
!
!-----------------------------------------------------------------------
!  Initialize ROMSJEDI, Phase 2.
!-----------------------------------------------------------------------
!
      Phase=2
!
      SELECT CASE (kernel)
        CASE (iNLM)
          CALL nlm_initial (Phase)
#ifdef TANGENT
        CASE (iTLM)
          CALL tlm_initial (Phase)
#endif
#ifdef ADJOINT
        CASE (iADM)
          CALL adm_initial (Phase)
#endif
      END SELECT
!
      RETURN
      END SUBROUTINE ROMS_initializeP2
!
      SUBROUTINE ROMS_run (RunInterval, kernel)
!
!=======================================================================
!                                                                      !
!  This routine advance ROMS kernel for the specified time interval.   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     RunInterval     Execution time stepping window (seconds)         !
!     kernel          Dynamical/numerical kernel (integer)             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in), optional :: kernel
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng
      integer :: NstrStep, NendStep, extra
!
      real(dp) :: ENDtime, NEXTtime
!
      character (len=2), dimension(4) :: MID = (/'NL','TL','RP','AD'/)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
!
      SourceFile=MyFile
!
!=======================================================================
!  Advance ROMS dynamical/numerical kernel: NLM, TLM, or ADM
!=======================================================================
!
!  OOPS will advance ROMS kernels by smaller intervals, usually a single
!  timestep. The strategy here is different to that used for coupling
!  since ROMS delayed output. The delayed ouput last half timestep will
!  affect the OOPS trajectory logic needed to save the NLM background
!  fields needed to linearize the TLM and ADM kernels.
!
      MyRunInterval=RunInterval
      IF (Master) WRITE (stdout,'(1x)')
!
      SELECT CASE (kernel)
        CASE (iNLM)
          DO ng=1,Ngrids
            NEXTtime=time(ng)+RunInterval
            ENDtime=INItime(ng)+ntimes(ng)*dt(ng)
            extra=1
            step_counter(ng)=0
            NstrStep=iic(ng)-1
            NendStep=NstrStep+INT((MyRunInterval)/dt(ng))-extra
            IF (Master) WRITE (stdout,10) MID(kernel), ng,              &
     &                                    NstrStep, MAX(0,NendStep)
          END DO
          IF (Master) WRITE (stdout,'(1x)')
        CASE (iTLM)
          DO ng=1,Ngrids
            NEXTtime=time(ng)+RunInterval
            ENDtime=INItime(ng)+ntimes(ng)*dt(ng)
            IF (NEXTtime.eq.ENDtime) THEN
              extra=0                              ! last time interval
            ELSE
              extra=1
            END IF
            step_counter(ng)=0
            NstrStep=iic(ng)-1
            NendStep=NstrStep+INT((MyRunInterval)/dt(ng))-extra
            IF (Master) WRITE (stdout,10) MID(kernel), ng,              &
     &                                    NstrStep, MAX(0,NendStep)
          END DO
          IF (Master) WRITE (stdout,'(1x)')
        CASE (iADM)
          DO ng=1,Ngrids
            NEXTtime=time(ng)-RunInterval
            ENDtime=INItime(ng)+ntimes(ng)*dt(ng)
            IF (time(ng).eq.ENDtime) THEN
              extra=0                              ! first time interval
            ELSE
              extra=1
            END IF
            step_counter(ng)=0
            NstrStep=iic(ng)-1
            NendStep=NstrStep-INT((MyRunInterval)/dt(ng))+extra
            IF (Master) WRITE (stdout,10) MID(kernel), ng,              &
     &                                    NstrStep, MAX(0,NendStep)
          END DO
          IF (Master) WRITE (stdout,'(1x)')
      END SELECT
!
!  Time-step ROMS kernel.
!
      SELECT CASE (kernel)
        CASE (iNLM)
          CALL main3d (RunInterval)
#ifdef TANGENT
        CASE (iTLM)
          CALL tl_main3d (RunInterval)
#endif
#ifdef ADJOINT
        CASE (iADM)
          CALL ad_main3d (RunInterval)
#endif
      END SELECT
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
 10   FORMAT (1x,a,1x,'ROMS: started time-stepping:',                   &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i0,' - ',i0,')')
!
      RETURN
      END SUBROUTINE ROMS_run
!
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS W4D-RBLanczos execution.               !
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
!  Stop model and time profiling clocks, report memory requirements,
!  and close output NetCDF files.
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
!
      SUBROUTINE nlm_initial (Phase)
!
!=======================================================================
!                                                                      !
!  ROMSJEDI nonlinear model (NLM) kernel initialization Phases:        !
!                                                                      !
!     Phase 1: Set time-stepping parameters                            !
!                                                                      !
!     Phase 2: Computes initial depths, density, horizontal mass       !
!              fluxes, and other configuration arrays. It reads        !
!              forcing snapshots.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in) :: Phase
!
!  Local variable declarations.
!
      integer :: Fcount
      integer :: ng, thread, tile
#ifdef NESTING
      integer :: ig, nl
      integer :: cr, i, m
#endif
      integer, dimension(Ngrids) :: IniRec, Tindex
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", nlm_initial"

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 2, __LINE__, MyFile)
        END DO
      END DO
#endif
!
!=======================================================================
!  ROMSJEDI NLM kernel, Phase 1 Initialization.
!=======================================================================
!
      PHASE1 : IF (Phase.eq.1) THEN
!
        IF (Master) THEN
          WRITE (stdout,10) 'NLM_INITIAL: Phase 1 Initialization: ',    &
     &                      'Configuring ROMS nonlinear kernel ...'
        END IF
!
!  Initialize time stepping indices and counters.
!
        DO ng=1,Ngrids
          iif(ng)=1
          indx1(ng)=1
          kstp(ng)=1
          krhs(ng)=1
          knew(ng)=1
          PREDICTOR_2D_STEP(ng)=.FALSE.
!
          iic(ng)=0
# ifdef JEDI
          jic(ng)=0
# endif
          nstp(ng)=1
          nrhs(ng)=1
          nnew(ng)=1
#ifdef FLOATS
          nf(ng)=0
          nfp1(ng)=1
          nfm1(ng)=4
          nfm2(ng)=3
          nfm3(ng)=2
#endif
!
          IniRec(ng)=nrrec(ng)
          Tindex(ng)=1
!
          synchro_flag(ng)=.TRUE.
          first_time(ng)=0
          tdays(ng)=dstart
          time(ng)=tdays(ng)*day2sec
#ifdef JEDI
          time4jedi(ng)=time(ng)
#endif
          ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
          ntfirst(ng)=ntstart(ng)
          step_counter(ng)=0
        END DO
!
!  Initialize global diagnostics variables.
!
        avgke=0.0_dp
        avgpe=0.0_dp
        avgkp=0.0_dp
        volume=0.0_dp
!
!  Reset output history files time record counters. These counters are
!  reset on every iteration pass. This file is created on the first
!  iteration pass.
!
        DO ng=1,Ngrids
          LdefHIS(ng)=.TRUE.
          LwrtHIS(ng)=.TRUE.
          HIS(ng)%Rindex=0
          Fcount=HIS(ng)%Fcount
          HIS(ng)%Nrec(Fcount)=0

          LdefQCK(ng)=.TRUE.
          LwrtQCK(ng)=.TRUE.
          QCK(ng)%Rindex=0
          Fcount=QCK(ng)%Fcount
          QCK(ng)%Nrec(Fcount)=0

          LdefRST(ng)=.TRUE.
          LwrtRST(ng)=.TRUE.
          RST(ng)%Rindex=0
          Fcount=RST(ng)%Fcount
          RST(ng)%Nrec(Fcount)=0
        END DO

      END IF PHASE1
!
!=======================================================================
!  ROMSJEDI NLM kernel, Phase 2 Initialization.
!=======================================================================
!
      PHASE2: IF (Phase.eq.2) THEN
!
        IF (Master) THEN
          WRITE (stdout,10) 'NLM_INITIAL: Phase 2 Initialization: ',    &
     &                      'Get/Set required applications fields ...'
        END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
#ifndef ANA_SPONGE
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
#endif
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_hmixcoef (ng, tile, iNLM)
          END DO
        END DO

#ifdef ANA_SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in sponge areas using
!  analytical functions.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lsponge(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_sponge (ng, tile, iNLM)
            END DO
          END IF
        END DO
#endif

#ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
        END DO
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time independent (Zt_avg1=0) anf initial time dependent
!  depths and level thicknesses.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_depth0 (ng, tile, iNLM)
            CALL set_depth  (ng, tile, iNLM)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_massflux (ng, tile, iNLM)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL omega (ng, tile, iNLM)
            CALL rho_eos (ng, tile, iNLM)
          END DO
        END DO
#endif

#ifdef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile with analytcal expressions.
!  Point sources are at U- and V-points. We need to get their positions
!  to process internal Land/Sea masking arrays during initialization.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_psource (ng, tile, iNLM)
            END DO
          END IF
        END DO
#endif
!
!-----------------------------------------------------------------------
!  If applicable, close all input boundary, climatology, and forcing
!  NetCDF files and set associated parameters to the closed state. This
!  step is essential in iterative algorithms that run the full TLM
!  repetitively. Then, Initialize several parameters in their file
!  structure, so the appropriate input single or multi-file is selected
!  during initialization/restart.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL close_inp (ng, iNLM)
          CALL check_multifile (ng, iNLM)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!  If applicable, read in input data.
!
        DO ng=1,Ngrids
          CALL get_idata (ng)
          CALL get_data (ng)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_masks (ng, tile, iNLM)
          END DO
        END DO
#endif

#ifdef NESTING
# if defined MASKING || defined WET_DRY
!
!-----------------------------------------------------------------------
!  If nesting and Land/Sea masking, scale horizontal interpolation
!  weights to account for land contact points.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL nesting (ng, iNLM, nmask)
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If nesting, process state fields initial conditions in the contact
!  regions.
!-----------------------------------------------------------------------
!
!  Free-surface and 2D-momentum.
!
        DO nl=1,NestLayers
          DO ig=1,GridsInLayer(nl)
            ng=GridNumber(ig,nl)
            IF (ANY(CompositeGrid(:,ng))) THEN
              CALL nesting (ng, iNLM, nFSIC)      ! free-surface
# ifndef SOLVE3D
              CALL nesting (ng, iNLM, n2dIC)      ! 2d momentum
# endif
            END IF
          END DO
        END DO

# ifdef SOLVE3D
!
!  Determine vertical indices and vertical interpolation weights in
!  the contact zone using initial unperturbed depth arrays.
!
        DO ng=1,Ngrids
          CALL nesting (ng, iNLM, nzwgt)
        END DO
!
!  3D-momentum and tracers.
!
        DO nl=1,NestLayers
          DO ig=1,GridsInLayer(nl)
            ng=GridNumber(ig,nl)
            IF (ANY(CompositeGrid(:,ng))) THEN
              CALL nesting (ng, iNLM, n3dIC)      ! 3D momentum
              CALL nesting (ng, iNLM, nTVIC)      ! Tracer variables
            END IF
          END DO
        END DO
# endif
#endif

#if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
        IF (Nrun.eq.ERstr) THEN
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_drag (ng, tile, iNLM)
            END DO
          END DO
        END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
        IF (Lstiffness) THEN
          Lstiffness=.FALSE.
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL stiffness (ng, tile, iNLM)
            END DO
          END DO
        END IF

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL grid_coords (ng, iNLM)
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and time/date string. Save NLM
!  initial conditions time.
!
!  Notice that it is allowed to modify the "simulation length" in the
!  roms-jedi YAML file, which will affect the values of "ntimes" and
!  "ntend".
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          INItime(ng)=time(ng)
          iic(ng)=ntstart(ng)
          ntend(ng)=ntstart(ng)+ntimes(ng)-1
#ifdef JEDI
          jic(ng)=ntstart(ng)-1
          time4jedi(ng)=time(ng)-dt(ng)
#endif
          CALL time_string (time(ng), time_code(ng))
        END DO

      END IF PHASE2

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 2, __LINE__, MyFile)
        END DO
      END DO
#endif
!
  10  FORMAT (/,1x,a,a,/,1x,'***********',/)
!
      RETURN
      END SUBROUTINE nlm_initial

#ifdef TANGENT
!
      SUBROUTINE tlm_initial (Phase)
!
!=======================================================================
!                                                                      !
!  ROMSJEDI tangent linear model (TLM) kernel initialization Phases:   !
!                                                                      !
!     Phase 1: Set time-stepping parameters                            !
!                                                                      !
!     Phase 2: Computes masks and other configuration arrays. It reads !
!              initial forcing snapshots.                              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in) :: Phase
!
!  Local variable declarations.
!
      integer :: Fcount
      integer :: ng, thread, tile
      integer :: lstr, ifile

      integer, dimension(Ngrids) :: IniRec, Tindex

# ifdef GENERIC_DSTART
!
      real(dp) :: my_dstart
# endif
!
      character (len=10) :: suffix

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tlm_initial"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iTLM, 2, __LINE__, MyFile)
        END DO
      END DO
# endif
!
!=======================================================================
!  ROMSJEDI TLM kernel, Phase 1 Initialization.
!=======================================================================
!
      PHASE1 : IF (Phase.eq.1) THEN
!
        IF (Master) THEN
          WRITE (stdout,10) 'TLM_INITIAL: Phase 1 Initialization: ',    &
     &                      'Configuring ROMS tangent linear kernel ...'
        END IF
!
!  Initialize time stepping indices and counters.
!
        DO ng=1,Ngrids
          iif(ng)=1
          indx1(ng)=1
          kstp(ng)=1
          krhs(ng)=1
          knew(ng)=1
          PREDICTOR_2D_STEP(ng)=.FALSE.
!
          iic(ng)=0
# ifdef JEDI
          jic(ng)=0
# endif
          nstp(ng)=1
          nrhs(ng)=1
          nnew(ng)=1
# ifdef FLOATS_NOT_YET
          nf(ng)=0
          nfp1(ng)=1
          nfm1(ng)=4
          nfm2(ng)=3
          nfm3(ng)=2
# endif
!
          synchro_flag(ng)=.TRUE.
          first_time(ng)=0
          IF (ANY(tl_VolCons(:,ng))) THEN
            tl_ubar_xs=0.0_r8
          END IF

# if defined GENERIC_DSTART
!
!  Rarely, the tangent linear model is initialized from a NetCDF file,
!  so we do not know its actual initialization time. Usually, it is
!  computed from DSTART, implying that its value is correct in the ROMS
!  input script. Therefore, the user needs to check and update its value
!  to every time that ROMS is executed. Alternatively, if available, we
!  can use the initialization time from the nonlinear model, INItime.
!  This variable is assigned when computing or processing the basic
!  state trajectory needed to linearize the adjoint model.
!
          IF (INItime(ng).lt.0.0_dp) THEN
            my_dstart=dstart                  ! ROMS input script
          ELSE
            my_dstart=INItime(ng)/86400.0_dp  ! NLM IC time is known
          END IF
          tdays(ng)=my_dstart
# else
          tdays(ng)=dstart
# endif
          time(ng)=tdays(ng)*day2sec
# ifdef JEDI
          time4jedi(ng)=time(ng)
# endif
          ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
          ntend(ng)=ntstart(ng)+ntimes(ng)-1
          ntfirst(ng)=ntstart(ng)
!
          IniRec(ng)=nrrec(ng)
          Tindex(ng)=1
        END DO
!
!  Initialize global diagnostics variables.
!
        avgke=0.0_dp
        avgpe=0.0_dp
        avgkp=0.0_dp
        volume=0.0_dp

      END IF PHASE1
!
!=======================================================================
!  ROMSJEDI TLM kernel, Phase 2 Initialization.
!=======================================================================
!
      PHASE2: IF (Phase.eq.2) THEN
!
        IF (Master) THEN
          WRITE (stdout,10) 'TLM_INITIAL: Phase 2 Initialization: ',   &
     &                      'Get/Set required applications fields ...'
        END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
# ifndef ANA_SPONGE
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
# endif
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_hmixcoef (ng, tile, iTLM)
          END DO
        END DO

# ifdef ANA_SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in sponge areas using
!  analytical functions.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lsponge(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_sponge (ng, tile, iTLM)
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry to zero.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_bath (ng, tile)
          END DO
        END DO

# ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
        END DO
# endif

# ifdef FORWARD_FLUXES
!
!-----------------------------------------------------------------------
!  Set the BLK structure to contain the nonlinear model surface fluxes
!  needed by the tangent linear and adjoint models. Also, set switches
!  to process that structure in routine "check_multifile". Notice that
!  it is possible to split the solution into multiple NetCDF files to
!  reduce their size.
!-----------------------------------------------------------------------
!
!  Set the nonlinear model quicksave-history file as the basic state for
!  the surface fluxes computed in "bulk_flux", which may be available at
!  more frequent intervals while avoiding large files. Since the 4D-Var
!  increment phase is calculated by a different executable and needs to
!  know some of the QCK structure information.
!
      DO ng=1,Ngrids
        WRITE (QCK(ng)%name,"(a,'.nc')") TRIM(QCK(ng)%head)
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
!-----------------------------------------------------------------------
!  If applicable, close all input boundary, climatology, and forcing
!  NetCDF files and set associated parameters to the closed state. This
!  step is essential in iterative algorithms that run the full TLM
!  repetitively. Then, Initialize several parameters in their file
!  structure, so the appropriate input single or multi-file is selected
!  during initialization/restart.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL close_inp (ng, iTLM)
          CALL check_multifile (ng, iTLM)
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL tl_get_idata (ng)
          CALL tl_get_data (ng)
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_masks (ng, tile, iTLM)
          END DO
        END DO
# endif

# if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
        IF (Nrun.eq.ERstr) THEN
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_drag (ng, tile, iTLM)
            END DO
          END DO
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
        IF (Lstiffness) THEN
          Lstiffness=.FALSE.
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL stiffness (ng, tile, iTLM)
            END DO
          END DO
        END IF

# if defined FLOATS_NOT_YET || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL grid_coords (ng, iTLM)
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and time/date string.
!
!  Notice that it is allowed to modify the "simulation length" in the
!  roms-jedi YAML file, which will affect the values of "ntimes" and
!  "ntend".
!-----------------------------------------------------------------------
!
!  Subsract one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
        DO ng=1,Ngrids
          iic(ng)=ntstart(ng)
          ntend(ng)=ntstart(ng)+ntimes(ng)-1
#ifdef JEDI
          jic(ng)=ntstart(ng)-1
          time4jedi(ng)=time(ng)-dt(ng)
#endif
          CALL time_string (time(ng), time_code(ng))
          LdefTLM(ng)=.TRUE.
          LwrtTLM(ng)=.TRUE.
        END DO

      END IF PHASE2

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 2, __LINE__, MyFile)
        END DO
      END DO
# endif
!
  10  FORMAT (/,1x,a,a,/,1x,'***********',/)
!
      RETURN
      END SUBROUTINE tlm_initial
#endif

#ifdef ADJOINT
!
      SUBROUTINE adm_initial (Phase)
!
!=======================================================================
!                                                                      !
!  ROMSJEDI adjoint model (ADM) kernel initialization Phases:          !
!                                                                      !
!     Phase 1: Set time-stepping parameters                            !
!                                                                      !
!     Phase 2: Computes masks and other configuration arrays. It reads !
!              initial forcing snapshots.                              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      integer, intent(in) :: Phase
!
!  Local variable declarations.
!
      integer :: Fcount
      integer :: ng, thread, tile

      integer, dimension(Ngrids) :: IniRec, Tindex

# ifdef GENERIC_DSTART
!
      real(dp) :: my_dstart
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", adm_initial"

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iADM, 2, __LINE__, MyFile)
        END DO
      END DO
# endif
!
!=======================================================================
!  ROMSJEDI ADM kernel, Phase 1 Initialization.
!=======================================================================
!
      PHASE1 : IF (Phase.eq.1) THEN
!
        IF (Master) THEN
          WRITE (stdout,10) 'ADM_INITIAL: Phase 1 Initialization: ',    &
     &                      'Configuring ROMS nonlinear kernel ...'
        END IF
!
!  Initialize time stepping indices and counters.
!
        DO ng=1,Ngrids
          iif(ng)=1
          indx1(ng)=1
          kstp(ng)=1
          krhs(ng)=3
          knew(ng)=2
          PREDICTOR_2D_STEP(ng)=.FALSE.
!
          iic(ng)=0
# ifdef JEDI
          jic(ng)=0
# endif
# ifdef SOLVE3D
          nstp(ng)=1
          nnew(ng)=2
          nrhs(ng)=nstp(ng)
# endif
# ifdef FLOATS_NOT_YET
          nf(ng)=0
          nfp1(ng)=1
          nfm1(ng)=4
          nfm2(ng)=3
          nfm3(ng)=2
# endif
!
          synchro_flag(ng)=.TRUE.
          first_time(ng)=0
          ad_ubar_xs=0.0_r8

# ifdef GENERIC_DSTART
!
!  Rarely, the adjoint model is initialized from a NetCDF file, so we do
!  not know its actual initialization time. Usually, it is computed from
!  DSTART, implying that its value is correct in the ROMS input script.
!  Therefore, the user needs to check and update its value to every time
!  that ROMS is executed. Alternatively, if available, we can use the
!  initialization time from the nonlinear model, INItime. This variable
!  is assigned when computing or processing the basic state trajectory
!  needed to linearize the adjoint model.
!
          IF (INItime(ng).lt.0.0_dp) THEN
            my_dstart=dstart                ! ROMS input script
          ELSE
            my_dstart=INItime(ng)/86400.0_dp    ! NLM IC time is known
          END IF
          tdays(ng)=my_dstart+dt(ng)*REAL(ntimes(ng),r8)*sec2day
          time(ng)=tdays(ng)*day2sec
          ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
          ntend(ng)=ntstart(ng)-ntimes(ng)
          ntfirst(ng)=ntend(ng)
# else
          tdays(ng)=dstart+                                             &
     &              dt(ng)*REAL(ntimes(ng)-ntfirst(ng)+1,r8)*sec2day
          time(ng)=tdays(ng)*day2sec
          ntstart(ng)=ntimes(ng)+1
          ntend(ng)=ntfirst(ng)
          ntfirst(ng)=ntend(ng)
# endif
# ifdef JEDI
          time4jedi(ng)=time(ng)
# endif
          IniRec(ng)=nrrec(ng)
          Tindex(ng)=1
        END DO
!
!  Initialize global diagnostics variables.
!
        avgke=0.0_dp
        avgpe=0.0_dp
        avgkp=0.0_dp
        volume=0.0_dp

      END IF PHASE1
!
!=======================================================================
!  ROMSJEDI ADM kernel, Phase 2 Initialization.
!=======================================================================
!
      PHASE2: IF (Phase.eq.2) THEN
!
        IF (Master) THEN
          WRITE (stdout,10) 'ADM_INITIAL: Phase 2 Initialization: ',    &
     &                      'Get/Set required applications fields ...'
        END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
# ifndef ANA_SPONGE
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
# endif
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_hmixcoef (ng, tile, iADM)
          END DO
        END DO

# ifdef ANA_SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in sponge areas using
!  analytical functions.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lsponge(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_sponge (ng, tile, iADM)
            END DO
          END IF
        END DO
# endif

# ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
        END DO
# endif

# ifdef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile with analytcal expressions.
!  Point sources are at U- and V-points. We need to get their positions
!  to process internal Land/Sea masking arrays during initialization.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_psource (ng, tile, iADM)
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If applicable, close all input boundary, climatology, and forcing
!  NetCDF files and set associated parameters to the closed state. This
!  step is essential in iterative algorithms that run the full TLM
!  repetitively. Then, Initialize several parameters in their file
!  structure, so the appropriate input single or multi-file is selected
!  during initialization/restart.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL close_inp (ng, iADM)
          CALL check_multifile (ng, iADM)
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iADM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL ad_get_idata (ng)
          CALL ad_get_data (ng)
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iADM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_masks (ng, tile, iADM)
          END DO
        END DO
# endif

# if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
        IF (Nrun.eq.ERstr) THEN
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_drag (ng, tile, iADM)
            END DO
          END DO
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
        IF (Lstiffness) THEN
          Lstiffness=.FALSE.
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL stiffness (ng, tile, iADM)
            END DO
          END DO
        END IF

# if defined FLOATS_NOT_YET || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL grid_coords (ng, iADM)
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!!        ntstart(ng)=ntstart(ng)-1
          iic(ng)=ntstart(ng)
#ifdef JEDI
          jic(ng)=ntstart(ng)+1
          time4jedi(ng)=time(ng)+dt(ng)
#endif
          CALL time_string (time(ng), time_code(ng))
          LdefADJ(ng)=.TRUE.
          LwrtADJ(ng)=.TRUE.
        END DO

      END IF PHASE2

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iADM, 2, __LINE__, MyFile)
        END DO
      END DO
# endif
!
  10  FORMAT (/,1x,a,a,/,1x,'***********',/)
!
      RETURN
      END SUBROUTINE adm_initial
#endif

      END MODULE roms_kernel_mod
