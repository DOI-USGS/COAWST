      PROGRAM waves
!
!svn $Id: waves.h 594 2008-04-01 18:11:31Z jcwarner $
!================================================== John C. Warner   ===
!                                                                      !
!  This is a higher level driver to run SWAN or WaveWatch III          !
!                                                                      !
!=======================================================================

#if defined SWAN_MODEL
      USE mod_swan_kinds
      USE swan_iounits
      USE M_MPI
      USE waves_control_mod, ONLY : SWAN_driver_init
      USE waves_control_mod, ONLY : SWAN_driver_run
      USE waves_control_mod, ONLY : SWAN_driver_finalize
#endif
      implicit none
#if defined WW3_MODEL
      include 'mpif.h'
#endif
!
      integer :: MyError, MyRank, exit_flag
      integer, parameter :: stdout = 6
      integer, parameter :: NoError = 0
!
      exit_flag=0
#ifdef DISTRIBUTE
# ifdef MPI
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory MPI configuration.
!-----------------------------------------------------------------------
!
      CALL mpi_init (MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' WAVES - Unable to initialize MPI.')
        exit_flag=6
      END IF
!
!  Get rank of the local process in the group associated with the
!  comunicator.
!
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,20)
  20    FORMAT (/,' WAVES - Unable to inquire rank of local',            &
     &              ' processor.')
        exit_flag=6
      END IF
# endif
#endif

#if defined SWAN_MODEL
      IF (exit_flag.eq.NoError) THEN
        CALL SWAN_driver_init (MPI_COMM_WORLD)
      END IF

      IF (exit_flag.eq.NoError) THEN
        CALL SWAN_driver_run
      END IF

      CALL SWAN_driver_finalize

#elif defined WW3_MODEL
      IF (exit_flag.eq.NoError) THEN
        CALL WW3_init (MPI_COMM_WORLD)
      END IF
!     CALL WW3_driver_run
!     CALL WW3_driver_finalize
#endif

#if defined DISTRIBUTE && defined MPI
      CALL mpi_finalize (MyError)
#endif

      END PROGRAM waves
