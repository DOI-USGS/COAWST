  program fire_behavior

#ifdef DM_PARALLEL
    use mpi
#endif
    use state_mod, only : state_fire_t
    use namelist_mod, only : namelist_t
    use initialize_mod, only : Init_fire_state, Init_atm_state
    use advance_mod, only : Advance_state
    use wrfdata_mod, only : wrfdata_t
    use, intrinsic :: iso_fortran_env, only : ERROR_UNIT, OUTPUT_UNIT

    implicit none

    integer :: ierr, rank, mpi_comm_cfbm
    type (state_fire_t) :: grid
    type (wrfdata_t) :: atm_state
    type (namelist_t) :: config_flags
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'Running fire_behavior...'

#ifdef DM_PARALLEL
    call mpi_init (ierr)
    if (ierr /= MPI_SUCCESS) then
      write (ERROR_UNIT, *) 'ERROR: mpi_init failed'
      stop
    end if

    call MPI_Comm_dup (MPI_COMM_WORLD, mpi_comm_cfbm, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (ERROR_UNIT, *) 'ERROR: mpi_comm_dup failed'
      stop
    end if
#endif

      ! Read namelist
#ifdef DM_PARALLEL
    call Mpi_comm_rank (mpi_comm_cfbm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (ERROR_UNIT, *) 'ERROR: mpi_comm_rank failed'
      stop
    end if
#else
    rank = 0
#endif

    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) '  Reading namelist...'
    if (rank == 0) call config_flags%Initialization (file_name = 'namelist.fire')

#ifdef DM_PARALLEL
    call config_flags%Broadcast_nml (mpi_comm_cfbm)
#endif

    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) '  Initialization fire state...'
#ifdef DM_PARALLEL
    call grid%Set_mpi_comm_cfbm (mpi_comm_cfbm)
#endif
    select case (config_flags%ideal_opt)
      case (0)
        call Init_atm_state (atm_state, config_flags)
        call Init_fire_state (grid, config_flags, atm_state)

      case (1)
        call Init_fire_state (grid, config_flags)

      case default
        write (ERROR_UNIT, *) 'ERROR: ideal_opt option not supported: ', config_flags%ideal_opt
        stop

    end select

    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) '  Saving fire state...'
    call grid%Save_state ()

    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) '  Starting temporal loop...'
    do while (grid%datetime_now < grid%datetime_end)
      call Advance_state (grid, config_flags)
      call grid%Handle_output (config_flags)
      if (config_flags%ideal_opt == 0) call grid%Handle_wrfdata_update (atm_state, config_flags)
    end do
    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) '  Completed temporal loop'

#ifdef DM_PARALLEL
    call mpi_finalize (ierr)
    if (ierr /= MPI_SUCCESS) then
      write (ERROR_UNIT, *) 'ERROR: mpi_finalize failed'
      stop
    end if
#endif

    if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'Completed running fire_behavior'

  end program fire_behavior
