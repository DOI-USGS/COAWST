      PROGRAM inwave
!
!svn $Id: inwave.h 594 2008-04-01 18:11:31Z warner $
!=======================================================================
!                                                                      !
!  InFragravity Wave model                                             !
!                                                                      !
!  Dr. John C. Warner                                                  !
!    U.S. Geological Survey                                            !
!    Woods Hole, MA, USA                                               !
!    (jcwarner@usgs.gov)                                               !
!
!   Maitane Olabarrieta
!                                                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE ocean_control_mod, ONLY : ROMS_initialize
      USE ocean_control_mod, ONLY : ROMS_run
      USE ocean_control_mod, ONLY : ROMS_finalize
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: ng, MyError

      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend

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
  10    FORMAT (/,' ROMS/TOMS - Unable to initialize MPI.')
        exit_flag=6
      END IF
!
!  Get rank of the local process in the group associated with the
!  comunicator.
!
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,20)
  20    FORMAT (/,' ROMS/TOMS - Unable to inquire rank of local',       &
     &              ' processor.')
        exit_flag=6
      END IF
# endif
#endif
!
!-----------------------------------------------------------------------
!  Initialize ocean internal and external parameters and state
!  variables.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        first=.TRUE.
        CALL ROMS_initialize (first)
!
        CALL inwave_init
      END IF
!
!-----------------------------------------------------------------------
!  Time-step ocean model.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        Tstr(ng)=ntstart(ng)
        Tend(ng)=ntend(ng)+1
      END DO
      IF (exit_flag.eq.NoError) THEN
!       CALL ROMS_run (Tstr, Tend) 
!
        CALL inwave_run
      END IF
!
!-----------------------------------------------------------------------
!  Terminate ocean model execution: flush and close all IO files.
!-----------------------------------------------------------------------
!
      CALL inwave_finalize
      CALL ROMS_finalize
#if defined DISTRIBUTE && defined MPI
      CALL mpi_finalize (MyError)
#endif

      END PROGRAM inwave
