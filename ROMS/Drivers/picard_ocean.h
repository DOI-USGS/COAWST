      MODULE ocean_control_mod
!
!svn $Id: picard_ocean.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Picard Iterations Driver:                                 !
!                                                                      !
!  This driver is used to perform the Picard iterations test for the   !
!  representers tangent linear model used in IOMs weak constraint 4D   !
!  variational data assimilation (W4DVAR).  Recall that all  tangent   !
!  linear variables are in term of the full fields and the model can   !
!  expressed symbolically as:                                          !
!                                                                      !
!         d(S')/d(t) = N(So) + A(S' - So)                              !
!                                                                      !
!  where S' is the tangent linear state and So is the "basic state".   !
!  The "basic state" here is the solution of previous tangent linear   !
!  model iteration.                                                    !
!                                                                      !
!  This driver uses ESMF conventions for the  initialization,  time-   !
!  stepping,  and  finalization  of the representer  tangent  linear   !
!  model via:                                                          !
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

      SUBROUTINE ROMS_initialize (first, MyCOMM)
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
      USE mod_iounits
      USE mod_scalars
!
#ifdef AIR_OCEAN 
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

      integer :: ng, thread

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MyCOMM)) THEN
        OCN_COMM_WORLD=MyCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
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
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)
      END IF

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS representer tangent linear model. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      integer :: lstr, ng, my_iic
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------
!
      NEST_LOOP : DO ng=1,Ngrids
!
!  Use ensemble parameters for Picard itereations.
!
        ITER_LOOP : DO Nrun=ERstr,ERend
!
!  Cycle history and forward file names in such a way that the history
!  from the previous iteration becomes the basic state for the next.
!
          lstr=LEN_TRIM(TLMbase(ng))
          WRITE (TLMname(ng),10) TLMbase(ng)(1:lstr-3), Nrun
          WRITE (FWDname(ng),10) TLMbase(ng)(1:lstr-3), Nrun-1

          IF (Master) THEN
            WRITE (stdout,20) 'ROMS/TOMS Picard Iteration: ', Nrun,     &
     &                        TRIM(TLMname(ng)),                        &
     &                        TRIM(FWDname(ng))
          END IF
!
!  Activate defining history an restart files on each iteration. The
!  restart file is used to the store the solution of each iteration.
!
          iic(ng)=0
          LdefTLM(ng)=.TRUE.
          LwrtTLM(ng)=.TRUE.
          LdefRST(ng)=.FALSE.
!
!  Initialize representer tangent linear model.
!
          CALL rp_initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Time-step representers tangent linear model
!
          IF (Master) THEN
            WRITE (stdout,30) ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          RP_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL rp_main3d (ng)
#else
            CALL rp_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO RP_LOOP
!
!  Close IO and re-initialize NetCDF switches.
!
          CALL close_io

        END DO ITER_LOOP

      END DO NEST_LOOP
!
 10   FORMAT (a,'_',i2.2,'.nc')
 20   FORMAT (/,a,i3,/,/,5x,'  History file: ',a,                       &
     &                 /,5x,'  Forward file: ',a,/)
 30   FORMAT (/,'RP ROMS/TOMS: started time-stepping:',                 &
     &            '( TimeSteps: ',i8.8,' - ',i8.8,')',/)

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
      integer :: ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      DO ng=1,Ngrids
        IF (LwrtRST(ng).and.(exit_flag.eq.1)) THEN
          IF (Master) WRITE (stdout,10)
 10       FORMAT (/,' Blowing-up: Saving latest model state into ',     & 
     &              ' RESTART file',/)
          IF (LcycleRST(ng).and.(NrecRST(ng).ge.2)) THEN
            tRSTindx(ng)=2
            LcycleRST(ng)=.FALSE.
          END IF
          blowup=exit_flag
          exit_flag=NoError
          CALL wrt_rst (ng)
        END IF
      END DO
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
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
