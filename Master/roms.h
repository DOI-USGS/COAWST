      PROGRAM MyROMS
!
!git $Id$
!svn $Id: roms.h 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Regional Ocean Model System (ROMS)                                  !
!  Terrain-following Ocean Model System (TOMS)                         !
!                                                                      !
!  Master program to execute  ROMS/TOMS  drivers in ocean mode only    !
!  without coupling (sequential or concurrent) to  any  atmospheric    !
!  model.                                                              !
!                                                                      !
!  This ocean model solves the free surface, hydrostatic, primitive    !
!  equations  over  variable  topography  using  stretched terrain-    !
!  following coordinates in the vertical and orthogonal curvilinear    !
!  coordinates in the horizontal.                                      !
!                                                                      !
!  Nonlinear Model Developers:                                         !
!                                                                      !
!  Dr. Hernan G. Arango                                                !
!    Institute of Marine and Coastal Sciences                          !
!    Rutgers University, New Brunswick, NJ, USA                        !
!    (arango@marine.rutgers.edu)                                       !
!                                                                      !
!  Dr. Alexander F. Shchepetkin                                        !
!    Institute of Geophysics and Planetary Physics                     !
!    UCLA, Los Angeles, CA, USA                                        !
!    (alex@atmos.ucla.edu)                                             !
!                                                                      !
!  Dr. John C. Warner                                                  !
!    U.S. Geological Survey                                            !
!    Woods Hole, MA, USA                                               !
!    (jcwarner@usgs.gov)                                               !
!                                                                      !
!  Tangent linear and Adjoint Models and Algorithms Developers:        !
!                                                                      !
!    Dr. Hernan G. Arango    (arango@marine.rutgers.edu)               !
!    Dr. Bruce Cornuelle     (bcornuelle@ucsd.edu)                     !
!    Dr. Emanuele Di Lorenzo (edl@eas.gatech.edu)                      !
!    Dr. Arthur J. Miller    (ajmiller@ucsd.edu)                       !
!    Dr. Andrew M. Moore     (ammoore@ucsc.edu)                        !
!    Dr. Brian Powell        (powellb@uscs.edu)                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE mod_arrays,      ONLY : ROMS_deallocate_arrays
!
      USE roms_kernel_mod, ONLY : ROMS_initialize
      USE roms_kernel_mod, ONLY : ROMS_run
      USE roms_kernel_mod, ONLY : ROMS_finalize
#if defined PIO_LIB && defined DISTRIBUTE
      USE set_pio_mod,     ONLY : finalize_pio
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: ng, MyError

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
!  variables for all nested grids, if applicable.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        first=.TRUE.
        CALL ROMS_initialize (first)
      END IF
!
!-----------------------------------------------------------------------
!  Time-step ocean model over all nested grids, if applicable, by the
!  specified time interval in seconds.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        run_time=0.0_r8
        DO ng=1,Ngrids
          run_time=MAX(run_time, dt(ng)*ntimes(ng))
        END DO
        CALL ROMS_run (run_time)
      END IF
!
!-----------------------------------------------------------------------
!  Terminate ocean model execution: flush and close all IO files.
!-----------------------------------------------------------------------
!
      CALL ROMS_finalize
#if defined PIO_LIB && defined DISTRIBUTE
      CALL finalize_pio
#endif
      CALL ROMS_deallocate_arrays
#if defined DISTRIBUTE && defined MPI
      CALL mpi_finalize (MyError)
#endif
!
!-----------------------------------------------------------------------
!  If error, issue abort signal. It is expected that the Unix or POSIX
!  environment should handle this as a return errorcode from the ROMS
!  main to abort. It facilitates error handling in scripting.
!-----------------------------------------------------------------------
!
!!    IF ((exit_flag.ne.NoError).or.(blowup.ne.0)) STOP 911   ! F-2003
      IF ((exit_flag.ne.NoError).or.(blowup.ne.0)) ERROR STOP ! F-2008
!
      END PROGRAM MyROMS
