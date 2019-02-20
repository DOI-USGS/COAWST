      PROGRAM mct_driver
!
!svn $Id: mct_driver.h 830 2017-01-24 21:21:11Z arango $
!=======================================================================
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
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
#ifdef WW3_COUPLING
!  WW3, Wave Watch 3:                                                  !
!        http://polar.ncep.noaa.gov/waves/wavewatch/                   !
!                                                                      !
#endif
!=======================================================================
!
#if defined ROMS_COUPLING
      USE mod_iounits
      USE mod_scalars
#endif
#if defined SWAN_COUPLING
      USE swan_iounits
#endif
#ifdef WRF_COUPLING
      USE module_wrf_top, ONLY : wrf_init
      USE module_wrf_top, ONLY : wrf_run
      USE module_wrf_top, ONLY : wrf_finalize
#endif
      USE mct_coupler_params
      USE mod_coupler_iounits
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
      integer :: MyRank, pelast
      integer :: Ocncolor, Wavcolor, Atmcolor
      integer :: ng, iw, io, ia, icc
      real(m8) :: lcm, gcdlcm

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
      CALL read_coawst_par(1)
!
!  Now that we know the input file names and locations for each model,
!  for each model read in the number of grids and the grid time steps.
!
      CALL read_model_inputs
!
      CALL allocate_coupler_params
!
#if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT || \
     defined MCT_INTERP_WV2AT
!
!  Read coupled model sparse matrix file names from standard input.
!
        CALL allocate_coupler_iounits
        CALL read_coawst_par(2)
#endif
!
!  Compute the mct send and recv instances.
!
!  For each model grid, determine the number of steps it should
!  compute before it sends data out.
!  For example, nWAV2OCN(1,2) is the number of steps the wave model
!  grid 1 should take before it sends data to the ocn grid 2.
!
#ifdef WAVES_OCEAN
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_WAV2OCN,lcm).eq.0) THEN
            nWAV2OCN(iw,io)=INT(TI_WAV2OCN/dtwav(iw))
          ELSE
            lcm=gcdlcm(TI_WAV2OCN,lcm)
            nWAV2OCN(iw,io)=INT(lcm/dtwav(iw))
          END IF
        END DO
      END DO
!
      DO io=1,Nocn_grids
        DO iw=1,Nwav_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_OCN2WAV,lcm).eq.0) THEN
            nOCN2WAV(io,iw)=INT(TI_OCN2WAV/dtocn(io))
          ELSE
            lcm=gcdlcm(TI_OCN2WAV,lcm)
            nOCN2WAV(io,iw)=INT(lcm/dtocn(io))
          END IF
        END DO
      END DO
#endif
#ifdef AIR_OCEAN
      DO ia=1,Natm_grids
        DO io=1,Nocn_grids
          lcm=gcdlcm(dtatm(ia),dtocn(io))
          IF (MOD(TI_ATM2OCN,lcm).eq.0) THEN
            nATM2OCN(ia,io)=INT(TI_ATM2OCN/dtatm(ia))
          ELSE
            lcm=gcdlcm(TI_ATM2OCN,lcm)
            nATM2OCN(ia,io)=INT(lcm/dtatm(ia))
          END IF
        END DO
      END DO
!
      DO io=1,Nocn_grids
        DO ia=1,Natm_grids
          lcm=gcdlcm(dtatm(ia),dtocn(io))
          IF (MOD(TI_OCN2ATM,lcm).eq.0) THEN
            nOCN2ATM(io,ia)=INT(TI_OCN2ATM/dtocn(io))
          ELSE
            lcm=gcdlcm(TI_OCN2ATM,lcm)
            nOCN2ATM(io,ia)=INT(lcm/dtocn(io))
          END IF
        END DO
      END DO
#endif
#ifdef AIR_WAVES
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          lcm=gcdlcm(dtatm(ia),dtwav(iw))
          IF (MOD(TI_ATM2WAV,lcm).eq.0) THEN
            nATM2WAV(ia,iw)=INT(TI_ATM2WAV/dtatm(ia))
          ELSE
            lcm=gcdlcm(TI_ATM2WAV,lcm)
            nATM2WAV(ia,iw)=INT(lcm/dtatm(ia))
          END IF
        END DO
      END DO
!
      DO iw=1,Nwav_grids
        DO ia=1,Natm_grids
          lcm=gcdlcm(dtatm(ia),dtwav(iw))
          IF (MOD(TI_WAV2ATM,lcm).eq.0) THEN
            nWAV2ATM(iw,ia)=INT(TI_WAV2ATM/dtwav(iw))
          ELSE
            lcm=gcdlcm(TI_WAV2ATM,lcm)
            nWAV2ATM(iw,ia)=INT(lcm/dtwav(iw))
          END IF
        END DO
      END DO
#endif
!
!  Similarly, for each model grid, determine the number of steps
!  it should compute before it recvs data from somewhere.
!  For example, nWAVFOCN(1,2) is the number of steps the wave model
!  grid 1 should take before it gets data from ocn grid 2.
!
#ifdef WAVES_OCEAN
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_OCN2WAV,lcm).eq.0) THEN
            nWAVFOCN(iw,io)=INT(TI_OCN2WAV/dtwav(iw))
          ELSE
            lcm=gcdlcm(TI_OCN2WAV,lcm)
            nWAVFOCN(iw,io)=INT(lcm/dtwav(iw))
          END IF
        END DO
      END DO
!
      DO io=1,Nocn_grids
        DO iw=1,Nwav_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_WAV2OCN,lcm).eq.0) THEN
            nOCNFWAV(io,iw)=INT(TI_WAV2OCN/dtocn(io))
          ELSE
            lcm=gcdlcm(TI_WAV2OCN,lcm)
            nOCNFWAV(io,iw)=INT(lcm/dtocn(io))
          END IF
        END DO
      END DO
#endif
#ifdef AIR_OCEAN
      DO ia=1,Natm_grids
        DO io=1,Nocn_grids
          lcm=gcdlcm(dtatm(ia),dtocn(io))
          IF (MOD(TI_OCN2ATM,lcm).eq.0) THEN
            nATMFOCN(ia,io)=INT(TI_OCN2ATM/dtatm(ia))
          ELSE
            lcm=gcdlcm(TI_OCN2ATM,lcm)
            nATMFOCN(ia,io)=INT(lcm/dtatm(ia))
          END IF
        END DO
      END DO
!
      DO io=1,Nocn_grids
        DO ia=1,Natm_grids
          lcm=gcdlcm(dtatm(ia),dtocn(io))
          IF (MOD(TI_ATM2OCN,lcm).eq.0) THEN
            nOCNFATM(io,ia)=INT(TI_ATM2OCN/dtocn(io))
          ELSE
            lcm=gcdlcm(TI_ATM2OCN,lcm)
            nOCNFATM(io,ia)=INT(lcm/dtocn(io))
          END IF
        END DO
      END DO
#endif
#ifdef AIR_WAVES
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          lcm=gcdlcm(dtatm(ia),dtwav(iw))
          IF (MOD(TI_WAV2ATM,lcm).eq.0) THEN
            nATMFWAV(ia,iw)=INT(TI_WAV2ATM/dtatm(ia))
          ELSE
            lcm=gcdlcm(TI_WAV2ATM,lcm)
            nATMFWAV(ia,iw)=INT(lcm/dtatm(ia))
          END IF
        END DO
      END DO
!
      DO iw=1,Nwav_grids
        DO ia=1,Natm_grids
          lcm=gcdlcm(dtatm(ia),dtwav(iw))
          IF (MOD(TI_ATM2WAV,lcm).eq.0) THEN
            nWAVFATM(iw,ia)=INT(TI_ATM2WAV/dtwav(iw))
          ELSE
            lcm=gcdlcm(TI_ATM2WAV,lcm)
            nWAVFATM(iw,ia)=INT(lcm/dtwav(iw))
          END IF
        END DO
      END DO
#endif
!
!  Allocate several coupling variables.
!
#ifdef ROMS_COUPLING
      allocate(ocnids(Nocn_grids))
#endif
#if defined SWAN_COUPLING || defined WW3_COUPLING
      allocate(wavids(Nwav_grids))
#endif
#ifdef WRF_COUPLING
      allocate(atmids(Natm_grids))
#endif
!
      N_mctmodels=0
#ifdef ROMS_COUPLING
      DO ng=1,Nocn_grids
        N_mctmodels=N_mctmodels+1
        ocnids(ng)=N_mctmodels
      END DO
#endif
#if defined SWAN_COUPLING || defined WW3_COUPLING
      DO ng=1,Nwav_grids
        N_mctmodels=N_mctmodels+1
        wavids(ng)=N_mctmodels
      END DO
#endif
#ifdef WRF_COUPLING
      DO ng=1,Natm_grids
        N_mctmodels=N_mctmodels+1
        atmids(ng)=N_mctmodels
     END DO
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
#if defined SWAN_COUPLING || defined WW3_COUPLING
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
#if defined SWAN_COUPLING || defined WW3_COUPLING
          WRITE (stdout,21) peWAV_frst, peWAV_last
 21       FORMAT (/,7x,'Waves Model MPI nodes: ',i3.3,' - ', i3.3)
#endif
#ifdef WRF_COUPLING
          WRITE (stdout,22) peATM_frst, peATM_last
 22       FORMAT (/,7x,'Atmos Model MPI nodes: ',i3.3,' - ', i3.3)
#endif
!
!  Write out some coupled model info.
!

#ifdef WAVES_OCEAN
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          WRITE (stdout,25) iw, dtwav(iw),io, dtocn(io),                &
     &                      TI_WAV2OCN, nWAV2OCN(iw,io)
 25       FORMAT (/,7x,'WAVgrid ',i2.2,' dt= ',f5.1,' -to- OCNgrid ',   &
     &            i2.2,' dt= ',f5.1,', CplInt: ',f7.1,' Steps: ',i3.3)
          WRITE (stdout,26) io, dtocn(io),iw, dtwav(iw),                &
     &                      TI_OCN2WAV, nOCN2WAV(io,iw)
 26       FORMAT (/,7x,'OCNgrid ',i2.2,' dt= ',f5.1,' -to- WAVgrid ',   &
     &            i2.2,' dt= ',f5.1,', CplInt: ',f7.1,' Steps: ',i3.3)
        END DO
      END DO
#endif
#ifdef AIR_OCEAN
      DO ia=1,Natm_grids
        DO io=1,Nocn_grids
          WRITE (stdout,27) ia, dtatm(ia),io, dtocn(io),                &
     &                      TI_ATM2OCN, nATM2OCN(ia,io)
 27       FORMAT (/,7x,'ATMgrid ',i2.2,' dt= ',f5.1,' -to- OCNgrid ',   &
     &            i2.2,' dt= ',f5.1,', CplInt: ',f7.1,' Steps: ',i3.3)
          WRITE (stdout,28) io, dtocn(io),ia, dtatm(ia),                &
     &                      TI_OCN2ATM, nOCN2ATM(io,ia)
 28       FORMAT (/,7x,'OCNgrid ',i2.2,' dt= ',f5.1,' -to- ATMgrid ',   &
     &            i2.2,' dt= ',f5.1,', CplInt: ',f7.1,' Steps: ',i3.3)
        END DO
      END DO
#endif
#ifdef AIR_WAVES
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          WRITE (stdout,29) ia, dtatm(ia),iw, dtwav(iw),                &
     &                      TI_ATM2WAV, nATM2WAV(ia,iw)
 29       FORMAT (/,7x,'ATMgrid ',i2.2,' dt= ',f5.1,' -to- WAVgrid ',   &
     &            i2.2,' dt= ',f5.1,', CplInt: ',f7.1,' Steps: ',i3.3)
          WRITE (stdout,30) iw, dtwav(iw),ia, dtatm(ia),                &
     &                      TI_WAV2ATM, nWAV2ATM(iw,ia)
 30       FORMAT (/,7x,'WAVgrid ',i2.2,' dt= ',f5.1,' -to- ATMgrid ',   &
     &            i2.2,' dt= ',f5.1,', CplInt: ',f7.1,' Steps: ',i3.3)
        END DO
      END DO
#endif
        END IF
      END IF
!     CALL flush_coawst (stdout)
!
!  Split the communicator into SWAN or WW3, WRF, and ROMS subgroups based
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
#if defined SWAN_COUPLING || defined WW3_COUPLING
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
        CALL SWAN_driver_init (MyCOMM)
        CALL SWAN_driver_run
        CALL SWAN_driver_finalize
      END IF
#elif defined WW3_COUPLING
      IF (MyColor.eq.WAVcolor) THEN
        CALL WW3_init (MyCOMM)
!        CALL WW3_driver_run
!        CALL WW3_driver_finalize
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
        CALL wrf_init (MyCOMM)
        CALL wrf_run
        CALL wrf_finalize(.TRUE.)
      END IF
#endif
#ifdef ROMS_COUPLING
      IF (MyColor.eq.OCNcolor) THEN
        first=.TRUE.
        Nrun=1
        IF (exit_flag.eq.NoError) THEN
          CALL ROMS_initialize (first, MyCOMM)
        END IF
        IF (exit_flag.eq.NoError) THEN
          run_time=0.0_m8
          DO ng=1,Ngrids
            run_time=MAX(run_time, dt(ng)*ntimes(ng))
          END DO
          CALL ROMS_run (run_time)
        END IF
        CALL ROMS_finalize
# if defined SWAN_COUPLING || defined REFDIF_COUPLING || \
     defined WW3_COUPLING
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
      END PROGRAM mct_driver

      FUNCTION gcdlcm (dtAin, dtBin)
!
!=======================================================================
!                                                                      !
!  This function computes the greatest common denominator              !
!  and lowest common multiple.                                         !
!                                                                      !
!  On Input:                                                           !
!     dtA        time step of model A                                  !
!     dtB        time step of model B                                  !
!                                                                      !
!  On Output:                                                          !
!     lcm        least common multiple                                 !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      real(m8), intent(in) :: dtAin, dtBin
      real(m8) :: gcdlcm
!
!  Local variable declarations.
!
      logical :: stayin
      real(m8) :: r, m, n, p, gcd, dtA, dtB
!
!-----------------------------------------------------------------------
!  Compute greatest common denominator and least common multiplier.
!-----------------------------------------------------------------------
      dtA=dtAin
      dtB=dtBin
      m=dtA
      n=dtB
      IF (dtA.gt.dtB) THEN
        p=dtA
        dtA=dtB
        dtB=p
      END IF
      stayin=.true.
      DO WHILE (stayin)
        r=mod(dtB,dtA)
        IF (r.eq.0) THEN
          gcd=dtA
          stayin=.false.
        ELSE
          dtB=dtA
          dtA=r
        END IF
      END DO
      gcdlcm=m*n/dtA

      RETURN
      END FUNCTION gcdlcm

      !SUBROUTINE flush_coawst (unit)
	  !-----------------------------------------------------------------------
	  !
	  !     USE mod_kinds
	  !
	  !      implicit none
	  !
	  !  Imported variable declarations.
	  !
	  !      integer, intent(in) :: unit
	  !
	  !      RETURN
	  !      END SUBROUTINE flush_coawst

