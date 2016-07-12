      MODULE ocean_control_mod
!
!svn $Id: symmetry.h 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Representer Driver:                                       !
!                                                                      !
!  This driver executes ROMS/TOMS weak constraint 4DVAR inner loop and !
!  checks the symmetry of the H R R' H' operator,  where R' and H' are !
!  transpose of R and H, respectively. The  R' H'  term is computed by !
!  integrating the adjoint model backwards while the  H R  is computed !
!  integrating forward the tangent linear model.                       !
!                                                                      !
!  It controls the initialization,  time-stepping, and finalization of !
!  the model execution following ESMF conventions:                     !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
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
      USE mod_scalars
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
        IF (exit_flag.ne.NoError) RETURN
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
            CALL wclock_on (ng, iNLM, 0)
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
        CALL get_state (ng, 6, 6, STD(1,ng)%name, STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      IF (NSA.eq.2) THEN
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STD(2,ng)%name, STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 8, 8, STD(3,ng)%name, STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 9, 9, STD(4,ng)%name, STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS weak constraint 4DVAR inner loop. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE convolve_mod, ONLY : error_covariance
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
#endif
      USE normalization_mod, ONLY : normalization
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: BOUNDED_TL

      integer :: i, j, ng, tile
      integer :: Lbck, Lini, Lstate, NRMrec
      integer :: IperAD, JperAD, KperAD, ivarAD
      integer :: IoutTL, JoutTL, KoutTL, ivarTL
#ifdef DISTRIBUTE
      integer :: Istr, Iend, Jstr, Jend
#endif
      integer, dimension(Ngrids) :: Nrec

      integer, allocatable :: StateVar(:)

      real(r8), allocatable :: R(:,:,:), Rerr(:,:)

      character (len=8 ) :: driver
      character (len=20) :: frmt
!
!-----------------------------------------------------------------------
!  Run nonlinear model and compute basic state tracjectory.
!-----------------------------------------------------------------------
!
      Lini=1                ! NLM initial conditions record in INI
      Lbck=1                ! background record in INI
      Rec1=1
      Rec2=2
      driver='symmetry'
!
!  Initialize and set nonlinear model initial conditions.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.FALSE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
      END DO

!$OMP PARALLEL
      CALL initial
!$OMP END PARALLEL
      IF (exit_flag.ne.NoError) RETURN
!
!  Run nonlinear model and compute basic state trajectory.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,10) 'NL', ng, ntstart(ng), ntend(ng)
        END IF
      END DO

!$OMP PARALLEL
#ifdef SOLVE3D
      CALL main3d (RunInterval)
#else
      CALL main2d (RunInterval)
#endif
!$OMP END PARALLEL
      IF (exit_flag.ne.NoError) RETURN
!
!  Set forward file trajectory.
!
      DO ng=1,Ngrids
        FWD(ng)%name=HIS(ng)%name
        FWD(ng)%ncid=HIS(ng)%ncid
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
!-----------------------------------------------------------------------
!  Compute or read in background-error correlations normalization
!  factors.
!-----------------------------------------------------------------------
!
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for a particular
!  application grid.
!
      DO ng=1,Ngrids
        IF (ANY(LwrtNRM(:,ng))) THEN
          CALL def_norm (ng, iNLM, 1)
          IF (exit_flag.ne.NoError) RETURN

          IF (NSA.eq.2) THEN
            CALL def_norm (ng, iNLM, 2)
          IF (exit_flag.ne.NoError) RETURN
          END IF
#ifdef ADJUST_BOUNDARY
          CALL def_norm (ng, iNLM, 3)
          IF (exit_flag.ne.NoError) RETURN
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL def_norm (ng, iNLM, 4)
          IF (exit_flag.ne.NoError) RETURN
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
          CALL get_state (ng, 5, 5, NRM(1,ng)%name, NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN

          IF (NSA.eq.2) THEN
            CALL get_state (ng, 5, 5, NRM(2,ng)%name, NRMrec, 2)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#ifdef ADJUST_BOUNDARY
          CALL get_state (ng, 10, 10, NRM(3,ng)%name, NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL get_state (ng, 11, 11, NRM(4,ng)%name, NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN
#endif
        END IF
      END DO
!
!  Define TLM impulse forcing NetCDF file.
!
      DO ng=1,Ngrids
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!=======================================================================
!  Perturb the adjoint state variable, one at the time, with a delta
!  function at the specified perturbation point. Integrate the adjoint
!  model backwards and the force the tangent linear model with the
!  result. Then, store tangent linear solution at the requested point
!  for each state variable.
!=======================================================================
!
!  Determine variables to perturb.
!
#ifdef SOLVE3D
      allocate ( StateVar(MAXVAL(NstateVar)-2) )
      StateVar(1)=isFsur
      StateVar(2)=isUvel
      StateVar(3)=isVvel
      DO i=1,MT
        StateVar(3+i)=isTvar(i)
      END DO
      Lstate=MAXVAL(StateVar)
#else
      allocate ( StateVar(MAXVAL(NstateVar)) )
      StateVar(1)=isFsur
      StateVar(2)=isUbar
      StateVar(3)=isVbar
      Lstate=3
#endif
      allocate ( R(Lstate,Lstate,Ngrids) )
      allocate ( Rerr(Lstate,Lstate) )
!
!  Initialize sample representer matrix to report.
!
      R(1:Lstate,1:Lstate,1:Ngrids)=0.0_r8
!
!  Set point to perturb and point to report.
!
      ivarTL=INT(user(1))
      ivarAD=INT(user(2))
      IoutTL=INT(user(3))
      IperAD=INT(user(4))
      JoutTL=INT(user(5))
      JperAD=INT(user(6))
      KoutTL=INT(user(7))
      KperAD=INT(user(8))
!
!  Determine number of perturbation iterations.
!
      ERstr=1
#ifdef SOLVE3D
      ERend=MINVAL(NstateVar)-2
#else
      ERend=MINVAL(NstateVar)
#endif
!
      PERT_LOOP : DO Nrun=ERstr,ERend
        user(2)=StateVar(Nrun)
!
!-----------------------------------------------------------------------
!  Time-step the adjoint model.
!-----------------------------------------------------------------------
!
!  Perturb the adjoint model at specified point.
!
        ADmodel=.TRUE.
        TLmodel=.FALSE.
        DO ng=1,Ngrids
          Lold(ng)=1
!$OMP PARALLEL
          CALL ad_initial (ng)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
          Fcount=ADM(ng)%Fcount
          ADM(ng)%Nrec(Fcount)=0
          ADM(ng)%Rindex=0
        END DO
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,10) 'AD', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL ad_main3d (RunInterval)
#else
        CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Convolve adjoint trajectory with model-error covariance and convert
!  to impulse forcing.
!-----------------------------------------------------------------------
!
        Lposterior=.FALSE.
        CALL error_covariance (iTLM, driver, -1, -1,                    &
     &                         Lbck, Lini, Lold, Lnew,                  &
     &                         Rec1, Rec2, Lposterior)
        IF (exit_flag.ne.NoError) RETURN
!
!  Convert convolved adjoint solution to impulse forcing. Write out
!  impulse forcing into TLF(ng)%name NetCDF file. To facilitate the
!  forcing by the TLM and RPM, the forcing is process and written in
!  increasing time coordinates.
!
        DO ng=1,Ngrids
          TLF(ng)%Rindex=0
#ifdef DISTRIBUTE
          tile=MyRank
#else
          tile=-1
#endif
          CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!-----------------------------------------------------------------------
!
!  Initialize tangent linear model from rest. The initial contribution
!  from the adjoint model will be added as impulse forcing in
!  "tl_forcing".
!
        ADmodel=.FALSE.
        TLmodel=.FALSE.
        DO ng=1,Ngrids
          wrtTLmod(ng)=.FALSE.
!$OMP PARALLEL
          CALL tl_initial (ng)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute R_n * PSI at observation points which
!  are used in the conjugate gradient algorithm.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,10) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
#else
        CALL tl_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN
!
!  Extract solution at requested points.
!
        DO ng=1,Ngrids
#ifdef DISTRIBUTE
          Istr=BOUNDS(ng)%Istr(MyRank)
          Iend=BOUNDS(ng)%Iend(MyRank)
          Jstr=BOUNDS(ng)%Jstr(MyRank)
          Jend=BOUNDS(ng)%Jend(MyRank)
          BOUNDED_TL=((Istr.le.IoutTL).and.(IoutTL.le.Iend)).and.       &
     &               ((Jstr.le.JoutTL).and.(JoutTL.le.Jend))
#else
          BOUNDED_TL=.TRUE.
#endif
          R(1,Nrun,ng)=OCEAN(ng)%tl_zeta(IoutTL,JoutTL,knew(ng))
#ifndef SOLVE3D
          R(2,Nrun,ng)=OCEAN(ng)%tl_ubar(IoutTL,JoutTL,knew(ng))
          R(3,Nrun,ng)=OCEAN(ng)%tl_vbar(IoutTL,JoutTL,knew(ng))
#else
          R(2,Nrun,ng)=OCEAN(ng)%tl_u(IoutTL,JoutTL,KoutTL,nstp(ng))
          R(3,Nrun,ng)=OCEAN(ng)%tl_v(IoutTL,JoutTL,KoutTL,nstp(ng))
          DO i=1,NT(ng)
            R(i+3,Nrun,ng)=OCEAN(ng)%tl_t(IoutTL,JoutTL,KoutTL,         &
     &                                    nstp(ng),i)
          END DO
#endif
        END DO

      END DO PERT_LOOP
!
!  Report sampled representer matrix and report symmetry.
!
      DO ng=1,Ngrids
#ifdef DISTRIBUTE
        CALL mp_bcastf (ng, iTLM, R(:,:,ng))
#endif
        IF (Master) THEN
          WRITE (stdout,20) 'Representer Matrix Symmetry Test: ',       &
     &                      'Perturbing Point: ',                       &
     &                      IperAD, JperAD, KperAD,                     &
     &                      'Sampling   Point: ',                       &
     &                      IoutTL, JoutTL, KoutTL
          WRITE (stdout,30) 'Sampled Representer Matrix: '
          IF (Lstate.lt.10) THEN
            WRITE (frmt,'(i1,a)') Lstate, '(1x,1p,e14.7,0p)'
          ELSE
            WRITE (frmt,'(i2,a)') Lstate, '(1x,1p,e14.7,0p)'
          END IF
          DO i=1,Lstate
            DO j=1,Lstate
              Rerr(i,j)=R(i,j,ng)-R(j,i,ng)
            END DO
          END DO
          DO i=1,Lstate
            WRITE (stdout,frmt) (R(i,j,ng),j=1,Lstate)
          END DO
          WRITE (stdout,30) 'Representer Matrix Symmetry Error: '
          DO i=1,Lstate
            WRITE (stdout,frmt) (Rerr(i,j),j=1,Lstate)
          END DO
        END IF
      END DO
!
 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
 20   FORMAT (/,1x,a,/,/,3x,a,' i = ',i4.4,' j = ',i4.4,' k = ',i4.4,   &
     &                 /,3x,a,' i = ',i4.4,' j = ',i4.4,' k = ',i4.4)
 30   FORMAT (/,1x,a,/)

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear model execution.        !
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
!  Compute and report model-observation comparison statistics.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        DO ng=1,Ngrids
          CALL stats_modobs (ng)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
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
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Close IO files.
!
      CALL close_out

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
