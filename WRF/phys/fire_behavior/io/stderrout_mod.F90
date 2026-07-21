  module stderrout_mod

    use, intrinsic :: iso_fortran_env, only : ERROR_UNIT, OUTPUT_UNIT

    implicit none

    private

    public :: Print_message, Stop_simulation

  contains

    subroutine Print_message (msg)

      implicit none

      character (len = *), intent (in) :: msg


      write (OUTPUT_UNIT, *) trim (msg)

    end subroutine Print_message

    subroutine Stop_simulation (msg)

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      character (len = *), intent (in) :: msg
      integer :: ierr


      write (ERROR_UNIT, *) 'STOP:' // trim (msg)
#ifdef DM_PARALLEL
      call mpi_abort (MPI_COMM_WORLD, 1, ierr)
      call mpi_finalize (ierr)
#else
      stop
#endif

    end subroutine Stop_simulation

  end module stderrout_mod
