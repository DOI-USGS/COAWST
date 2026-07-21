program build_tables
  !! program to build 3 lookup tables for tempo microphysics
  !!
  !! example build `make build_tables -f Makefile.intel mpi=1`
  !!
  !> timing information (ursa) to build the largest lookup table (rain-graupel collection)
  !>
  !> | Description | Time (s) |
  !> |-----------|-------------|
  !> | gfortran (-O2)                               | 3000 |
  !> | mpif90 (-O2) 1 process                       | 3000 |
  !> | mpif90 (-O2) 13 processes                    | 250 |
  !> | ifort (-fp-model consistent)                 | 1300 |
  !> | mpiifort (-fp-model consistent) 1 process    | 1300 |
  !> | mpiifort (-fp-model consistent) 13 processes | 120 |

  use module_mp_tempo_params, only : table_dp
  use module_mp_tempo_tables, only : tempo_build_tables

#ifdef build_tables_with_mpi
  use mpi_f08 
#endif

  implicit none

  integer :: build_tables_rank, build_tables_num_proc

#ifdef build_tables_with_mpi
    integer :: ierror, mpi_dp_size  
#endif

  build_tables_rank = 0
  build_tables_num_proc = 1

#ifdef build_tables_with_mpi
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, build_tables_num_proc, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, build_tables_rank, ierror)
  call MPI_Type_size(MPI_DOUBLE_PRECISION, mpi_dp_size, ierror)
  if (mpi_dp_size /= table_dp) then
    write(*,'(A,I1,A,I1,A)') 'MPI double precision size (', mpi_dp_size, &
      ') does not match fortran double precision size (', table_dp, ')'
    error stop 'MPI and fortran double precision sizes do not match.'
  endif
#endif

  call tempo_build_tables(build_tables_rank, build_tables_num_proc)

#ifdef build_tables_with_mpi
  call MPI_Finalize(ierror)
#endif

end program build_tables