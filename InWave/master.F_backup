      PROGRAM InWave
!
!svn $Id: InWave 001 2008-04-01 18:11:31Z jcwarner $
!=======================================================================
!    Licensed under a MIT/X style license                              !
!    See License_InWave.txt                                            !
!======================================================================= 
!                                                                      !
!  Infragravity Wave propagation model (InWave).                       !
!                                                                      !
!  This wave model solves the infragravity wave envelope of the        !
!  free surface over variable topography and orthogonal curvilinear    !
!  coordinates in the horizontal.                                      !
!                                                                      !
!  Dr. Maitane Olabarrieta                                             !
!    University Cantabria, Spain                                       !
!                                                                      !
!  Dr. John C. Warner                                                  !
!    U.S. Geological Survey                                            !
!    Woods Hole, MA, USA                                               !
!    (jcwarner@usgs.gov)                                               !
!                                                                      !
!=======================================================================
!
      USE mod_iparam

#if defined MPI
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory MPI configuration.
!-----------------------------------------------------------------------
!
      CALL mpi_init (MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' InWave - Unable to initialize MPI.')
        exit_flag=6
      END IF
!
!  Get rank of the local process in the group associated with the
!  comunicator.
!
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,20)
  20    FORMAT (/,' InWave - Unable to inquire rank of local',          &
     &              ' processor.')
        exit_flag=6
      END IF
#endif

       call InWave_init

       call InWave_run

       call InWave_finalize

#if defined MPI
      CALL mpi_finalize (MyError)
#endif

       END PROGRAM InWave


