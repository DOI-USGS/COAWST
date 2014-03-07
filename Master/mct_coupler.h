      PROGRAM mct_coupler
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  Master program to couple ROMS/TOMS to other models using the Model  !
!  Coupling Toolkit (MCT) library.                                     !
!                                                                      !
!  The following models are coupled to ROMS/TOMS:                      !
!                                                                      !
#ifdef WRF_COUPLING
!  WRF, Weather Research and Forecasting model:                        !
!       http://www.wrf-model.org                                       !
!                                                                      !
#endif
#ifdef SWAN_COUPLING
!  SWAN, Simulating WAves Nearshore model:                             !
!        http://vlm089.citg.tudelft.nl/swan/index.htm                  !
!                                                                      !
#endif
!=======================================================================
!
#if defined ROMS_COUPLING
!     USE mod_param
      USE mod_iounits
      USE mod_scalars
#endif
#if defined SWAN_COUPLING
      USE swan_iounits
#endif
#ifdef WRF_COUPLING
      USE module_wrf_top
#endif
!     USE mod_parallel
      USE mct_coupler_params
      USE mod_coupler_iounits
!     USE read_couplepar_mod
!
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean

#if defined ROMS_COUPLING
      USE ocean_control_mod, ONLY : ROMS_initialize
      USE ocean_control_mod, ONLY : ROMS_run
      USE ocean_control_mod, ONLY : ROMS_finalize
#endif
#if defined SWAN_COUPLING
      USE waves_control_mod, ONLY : SWAN_driver_init
      USE waves_control_mod, ONLY : SWAN_driver_run
      USE waves_control_mod, ONLY : SWAN_driver_finalize
#endif
#if defined WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : finalize_ocn2wav_coupling
#endif
#if defined AIR_OCEAN
      USE ocean_coupler_mod, ONLY : finalize_ocn2atm_coupling
#endif
!
      implicit none
      include 'mpif.h'
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: MyColor, MyCOMM, MyError, MyKey, Nnodes
      integer :: MyRank
      integer :: ng, pelast
      integer :: Ocncolor, Wavcolor, Atmcolor

#if defined ROMS_COUPLING
      integer, dimension(Ngrids) :: Tstr   ! starting ROMS time-step
      integer, dimension(Ngrids) :: Tend   ! ending   ROMS time-step
#endif

      real(m4) :: CouplingTime             ! single precision
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration
!-----------------------------------------------------------------------
!
!  Initialize MPI execution environment.
! 
      CALL mpi_init (MyError)
!
!  Get rank of the local process in the group associated with the
!  comminicator.
!
      CALL mpi_comm_size (MPI_COMM_WORLD, Nnodes, MyError)
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!  Read in coupled model parameters from standard input.
!    
      CALL read_CouplePar
!
!  Allocate several coupling variables.
!
#ifdef REFINED_GRID
# ifdef ROMS_COUPLING
      allocate(ocnids(Ngrids))
# endif
# ifdef SWAN_COUPLING
      allocate(wavids(Ngridss))
# endif
#endif
!
      N_mctmodels=0
      OCNid=0
      WAVid=0
      ATMid=0
#ifdef ROMS_COUPLING
# ifdef REFINED_GRID
      DO ng=1,Ngrids
        N_mctmodels=N_mctmodels+1
        ocnids(ng)=N_mctmodels
      END DO
      OCNid=ocnids(1)
# else
      N_mctmodels=N_mctmodels+1
      OCNid=N_mctmodels
# endif
#endif
#ifdef SWAN_COUPLING
# ifdef REFINED_GRID
      DO ng=1,Ngridss
        N_mctmodels=N_mctmodels+1
        wavids(ng)=N_mctmodels
      END DO
      WAVid=wavids(1)
# else
      N_mctmodels=N_mctmodels+1
      WAVid=N_mctmodels
# endif
#endif
#ifdef WRF_COUPLING
!     DO ng=1,max_dom
      N_mctmodels=N_mctmodels+1
      ATMid=N_mctmodels
!     END DO
#endif
!
!  Assign processors to the models.
!
      pelast=-1
#ifdef ROMS_COUPLING
      peOCN_frst=pelast+1
      peOCN_last=peOCN_frst+NnodesOCN-1
      pelast=peOCN_last
#endif
#ifdef SWAN_COUPLING
      peWAV_frst=pelast+1
      peWAV_last=peWAV_frst+NnodesWAV-1
      pelast=peWAV_last
#endif
#ifdef WRF_COUPLING
      peATM_frst=pelast+1
      peATM_last=peATM_frst+NnodesATM-1
      pelast=peATM_last
#endif
      IF (pelast.ne.Nnodes-1) THEN
        IF (MyRank.eq.0) THEN
          WRITE (stdout,10) pelast+1, Nnodes
 10       FORMAT (/,' mct_coupler - Number assigned processors: '       &
     &            ,i3.3,/,15x,'not equal to spawned MPI nodes: ',i3.3)
        END IF
        STOP
      ELSE
        IF (MyRank.eq.0) THEN
          WRITE (stdout,19)
 19       FORMAT (/,' Model Coupling: ',/)
#ifdef ROMS_COUPLING
          WRITE (stdout,20) peOCN_frst, peOCN_last
 20       FORMAT (/,7x,'Ocean Model MPI nodes: ',i3.3,' - ', i3.3)
#endif
#ifdef SWAN_COUPLING
          WRITE (stdout,21) peWAV_frst, peWAV_last
 21       FORMAT (/,7x,'Waves Model MPI nodes: ',i3.3,' - ', i3.3)
#endif
#ifdef WRF_COUPLING
          WRITE (stdout,22) peATM_frst, peATM_last
 22       FORMAT (/,7x,'Atmos Model MPI nodes: ',i3.3,' - ', i3.3)
#endif
        END IF
      END IF
!
!  Split the communicator into SWAN, WRF, and ROMS subgroups based 
!  on color and key.
!
      Atmcolor=1
      Ocncolor=2
      Wavcolor=3
      MyKey=0
#ifdef ROMS_COUPLING
      IF ((peOCN_frst.le.MyRank).and.(MyRank.le.peOCN_last)) THEN
        MyColor=OCNcolor
      END IF
#endif
#ifdef SWAN_COUPLING
      IF ((peWAV_frst.le.MyRank).and.(MyRank.le.peWAV_last)) THEN
        MyColor=WAVcolor
      END IF
#endif
#ifdef WRF_COUPLING
      IF ((peATM_frst.le.MyRank).and.(MyRank.le.peATM_last)) THEN
        MyColor=ATMcolor
      END IF
#endif
      CALL mpi_comm_split (MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
     &                     MyError)
!
!-----------------------------------------------------------------------
!  Run coupled models according to the processor rank.
!-----------------------------------------------------------------------
!
#if defined SWAN_COUPLING
      IF (MyColor.eq.WAVcolor) THEN
        CALL SWAN_driver_init (MyCOMM, REAL(TI_WAV_OCN))
        CALL SWAN_driver_run (REAL(TI_WAV_OCN))
        CALL SWAN_driver_finalize
      END IF
#elif defined REFDIF_COUPLING
      IF (MyColor.eq.WAVcolor) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iwaves))
        CALL refdif_initialize (MyCOMM)
        CALL refdif_run (CouplingTime, INPname(Iwaves))
        CALL refdif_finalize
      END IF
#endif
#ifdef WRF_COUPLING
      IF (MyColor.eq.ATMcolor) THEN
        CALL wrf_init (MyCOMM, N_mctmodels,                             &
#  ifdef REFINED_GRID
     &                 ocnids,                                          &
#  endif
     &                 OCNid, ATMid, WAVid,                             &
     &                 WRF_CPL_GRID,                                    &
     &                 REAL(TI_ATM_OCN),REAL(TI_ATM_WAV))
        CALL wrf_run
        CALL wrf_finalize
      END IF
#endif
#ifdef ROMS_COUPLING
      IF (MyColor.eq.OCNcolor) THEN
        first=.TRUE.
        Nrun=1
        IF (exit_flag.eq.NoError) THEN
          CALL ROMS_initialize (first, MyCOMM)
        END IF
        DO ng=1,Ngrids
          Tstr(ng)=ntstart(ng)
          Tend(ng)=ntend(ng)+1
        END DO
        IF (exit_flag.eq.NoError) THEN
          CALL ROMS_run (Tstr, Tend)
        END IF
        CALL ROMS_finalize
# if defined SWAN_COUPLING || defined REFDIF_COUPLING
        CALL finalize_ocn2wav_coupling
# endif
# ifdef WRF_COUPLING
        CALL finalize_ocn2atm_coupling
# endif
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Terminates all the mpi-processing and coupling.
!-----------------------------------------------------------------------
!
      CALL mpi_barrier (MPI_COMM_WORLD, MyError)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyError)

      STOP

      END PROGRAM mct_coupler
