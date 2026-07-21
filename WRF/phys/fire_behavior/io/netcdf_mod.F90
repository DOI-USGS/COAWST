  module netcdf_mod

    use stderrout_mod, only : Stop_simulation, Print_message
    use mpi_mod, only : Gather_var2d

    implicit none

    private

    character (len = 2), parameter :: NAME_DIM_X = 'nx', NAME_DIM_Y = 'ny'
    public :: Get_netcdf_var, Get_netcdf_att, Get_netcdf_dim, Create_netcdf_file, Add_netcdf_dim, Add_netcdf_var, &
       Is_netcdf_file_present, Is_netcdf_var_present, Add_netcdf_var_mpi, NAME_DIM_X, NAME_DIM_Y

    interface Add_netcdf_var
      module procedure Add_netcdf_var_real32_2d
      module procedure Add_netcdf_var_real32_3d
    end interface Add_netcdf_var

    interface Add_netcdf_var_mpi
      module procedure Add_netcdf_var_real32_2d_mpi
    end interface Add_netcdf_var_mpi

    interface Get_netcdf_var
      module procedure Get_netcdf_var_char_1d
      module procedure Get_netcdf_var_int32_3d
      module procedure Get_netcdf_var_real32_2d
      module procedure Get_netcdf_var_real32_3d
      module procedure Get_netcdf_var_real32_4d
    end interface Get_netcdf_var

    interface Get_netcdf_att
      module procedure Get_netcdf_att_int32
      module procedure Get_netcdf_att_real32
    end interface Get_netcdf_att

  contains

    subroutine Add_netcdf_dim (file_name, name_dim, val_dim)

      use netcdf

      implicit none

      character (len = *), intent(in) :: file_name, name_dim
      integer, intent (in) :: val_dim

      integer :: ncidout, status, dimid


      status = nf90_open (trim(file_name), NF90_WRITE, ncidout)
      call Check_status (status)

      status = nf90_redef (ncidout)
      call Check_status (status)

      status = nf90_def_dim (ncidout, name_dim, val_dim, dimid)
      call Check_status (status)

      status = nf90_enddef (ncidout)
      call Check_status (status)

      status = nf90_close (ncidout)
      call Check_status (status)

    end subroutine Add_netcdf_dim

    subroutine Add_netcdf_var_real32_2d (file_name, name_dims, varname, var)

      use netcdf

      implicit none

      character (len = *), intent (in) :: file_name, varname
      character (len = *), dimension(:), intent (in) :: name_dims
      real, dimension(:, :), intent (in) :: var

      integer, parameter :: N_DIMS = 2
      integer :: ncidout, status, varid, n
      integer :: dimids(N_DIMS)


      status = nf90_open (file_name, NF90_WRITE, ncidout)
      call Check_status (status)

      status = nf90_redef (ncidout)
      call Check_status (status)

      do n = 1, N_DIMS
        status = nf90_inq_dimid(ncidout, name_dims(n), dimids(n))
        call Check_status (status)
      end do

      status = nf90_inq_varid (ncidout, varname, varid)
      if (status == NF90_ENOTVAR) then
        status = nf90_def_var (ncidout, varname, NF90_FLOAT, dimids, varid)
        call Check_status (status)
      end if

      status = nf90_enddef (ncidout)
      call Check_status (status)

      status = nf90_put_var (ncidout, varid, var)
      call Check_status (status)

      status = nf90_close (ncidout)
      call Check_status (status)

    end subroutine Add_netcdf_var_real32_2d

    subroutine Add_netcdf_var_real32_2d_mpi (file_name, cfbm_comm, nx, ny, ifps, ifpe, jfps, jfpe, var_name, var2d_local)

#ifdef DM_PARALLEL
      use mpi
#endif
      implicit none

      integer, intent (in) :: cfbm_comm, nx, ny, ifps, ifpe, jfps, jfpe
      character (len = *), intent (in) :: file_name, var_name
      real, dimension(ifps:ifpe, jfps:jfpe), intent (in) :: var2d_local

      real, dimension(nx, ny) :: var2d
      integer :: rank, ierr
      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('Entering Add_netcdf_var_real32_2d_mpi...')

      if (DEBUG_LOCAL) call Print_message ('file name = ' // trim (file_name))
      if (DEBUG_LOCAL) call Print_message ('var name = ' // trim (var_name))
      if (DEBUG_LOCAL) call Print_message ('dim X name = ' // trim (NAME_DIM_X))
      if (DEBUG_LOCAL) call Print_message ('dim Y name = ' // trim (NAME_DIM_Y))

#ifdef DM_PARALLEL
      call Mpi_comm_rank (cfbm_comm, rank, ierr)
      if (ierr /= MPI_SUCCESS) call Stop_simulation ('Problems with Mpi_comm_rank ')

      call Gather_var2d (cfbm_comm, nx, ny, ifps, ifpe, jfps, jfpe, var2d_local(ifps:ifpe, jfps:jfpe), var2d)

      if (rank == 0) then
        call Add_netcdf_var (file_name, [NAME_DIM_X, NAME_DIM_Y], var_name, var2d(1:nx, 1:ny))
      end if
#else
      call Add_netcdf_var (file_name, [NAME_DIM_X, NAME_DIM_Y], var_name, var2d_local(1:nx, 1:ny))
#endif

      if (DEBUG_LOCAL) call Print_message ('Leaving Add_netcdf_var_real32_2d_mpi...')

    end subroutine Add_netcdf_var_real32_2d_mpi

    subroutine Add_netcdf_var_real32_3d (file_name, name_dims, varname, var)

      use netcdf

      implicit none

      character (len = *), intent (in) :: file_name, varname
      character (len = *), dimension(:), intent (in) :: name_dims
      real, dimension(:, :, :), intent (in) :: var

      integer, parameter :: N_DIMS = 3
      integer :: ncidout, status, varid, n
      integer :: dimids(N_DIMS)


      status = nf90_open (file_name, NF90_WRITE, ncidout)
      call Check_status (status)

      status = nf90_redef (ncidout)
      call Check_status (status)

      do n = 1, N_DIMS
        status = nf90_inq_dimid(ncidout, name_dims(n), dimids(n))
        call Check_status (status)
      end do

      status = nf90_inq_varid (ncidout, varname, varid)
      if (status == NF90_ENOTVAR) then
        status = nf90_def_var (ncidout, varname, NF90_FLOAT, dimids, varid)
        call Check_status (status)
      end if

      status = nf90_enddef (ncidout)
      call Check_status (status)

      status = nf90_put_var (ncidout, varid, var)
      call Check_status (status)

      status = nf90_close (ncidout)
      call Check_status (status)

    end subroutine Add_netcdf_var_real32_3d


    subroutine Check_status (status)

      use netcdf

      implicit none

      integer, intent (in) :: status


      if (status /= NF90_NOERR) call Stop_simulation (trim (nf90_strerror (status)))

    end subroutine Check_status

    subroutine Create_netcdf_file (file_name)

      use netcdf

      implicit none

      character (len = *), intent(in) :: file_name
      integer :: status, ncid

      status = nf90_create (trim(file_name), NF90_NETCDF4, ncid)
      call Check_status (status)

      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Create_netcdf_file

    subroutine Get_netcdf_att_int32 (file_name, var_name, att_name, att_value)

      use netcdf

      use, intrinsic :: iso_fortran_env, only :  INT32

      implicit none

      character (len = *), intent(in) :: file_name, var_name, att_name
      integer (kind = INT32), intent(out) :: att_value

      integer :: status, ncid, varid


        ! Open file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get att length
      if (var_name == 'global') then
        status = nf90_get_att (ncid, NF90_GLOBAL, att_name, att_value)
        call Check_status (status)
      else
        status = nf90_inq_varid (ncid, var_name, varid)
        call Check_status (status)

        status = nf90_get_att (ncid, varid, att_name, att_value)
        call Check_status (status)
      end if

        ! Close file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_att_int32

    subroutine Get_netcdf_att_real32 (file_name, var_name, att_name, att_value)

      use netcdf

      use, intrinsic :: iso_fortran_env, only :  REAL32

      implicit none

      character (len = *), intent(in) :: file_name, var_name, att_name
      real (kind = REAL32), intent(out) :: att_value

      integer :: status, ncid, varid


        ! Open file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get att length
      if (var_name == 'global') then
        status = nf90_get_att (ncid, NF90_GLOBAL, att_name, att_value)
        call Check_status (status)
      else
        status = nf90_inq_varid (ncid, var_name, varid)
        call Check_status (status)

        status = nf90_get_att (ncid, varid, att_name, att_value)
        call Check_status (status)
      end if

        ! Close file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_att_real32

    subroutine Get_netcdf_dim (file_name, dim_name, dim_value)

      use netcdf

      implicit none

      character (len = *), intent(in) :: file_name, dim_name
      integer, intent(out) :: dim_value

      integer :: status, ncid, dimid


        ! Opens file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get dim ID
      status = nf90_inq_dimid (ncid, dim_name, dimid)
      call Check_status (status)

        ! Get dim len
      status = nf90_inquire_dimension (ncid, dimid, len = dim_value)
      call Check_status (status)

        ! Closing file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_dim

    subroutine Get_netcdf_var_char_1d (file_name, name_var, var1d)

      use netcdf

      implicit none

      character (len = *), intent(in) :: file_name, name_var
      character (len = :), dimension (:), allocatable, intent (out) :: var1d

      integer :: status, ncid, ivar, io_stat, nf_type, nvdims, i, n1, n2
      integer, dimension (NF90_MAX_VAR_DIMS) :: dimids, idims
      character (len = :), allocatable :: io_errmsg, msg


      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

      status = nf90_inq_varid (ncid, name_var, ivar)
      call Check_status (status)

      status = nf90_inquire_variable (ncid, ivar, xtype = nf_type, ndims = nvdims, dimids = dimids)
      call Check_status (status)

      if (nf_type == NF90_CHAR) then
        if (nvdims == 2) then
          do i = 1, nvdims
            status = nf90_inquire_dimension (ncid, dimids(i), len = idims(i))
            call Check_status (status)
          end do

          n1 = idims(1)
          n2 = idims(2)

          allocate (character (len = n1) :: var1d(n2), stat = io_stat, errmsg = io_errmsg)
          if (io_stat /= 0) then
            call Print_message ('Problems allocating char var (1D)')
            call Stop_simulation (io_errmsg)
          end if

          status = nf90_get_var (ncid, ivar, var1d)
          call Check_status (status)
        else
          msg = 'FATAL ERROR: ' // trim (name_var) // ' is not a 2D array'
          call Stop_simulation (msg)
        end if
      else
        msg = 'FATAL ERROR: ' // name_var // ' is not a CHAR variable'
        call Stop_simulation (msg)
      end if

      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_var_char_1d

    subroutine Get_netcdf_var_int32_3d (file_name, var_name, output)

      use netcdf
      use, intrinsic :: iso_fortran_env, only : INT32

      implicit none

      character (len = *), intent(in) :: file_name, var_name
      integer (kind = INT32), dimension (:, :, :), allocatable :: output

      integer :: status, ncid, ivar, nf_type, nvdims
      integer :: i, n1, n2, n3, io_stat
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimids, idims
      character (len = :), allocatable :: msg


        ! Open file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get var
      status = nf90_inq_varid (ncid, trim (var_name), ivar)
      call Check_status (status)
      status = nf90_inquire_variable (ncid, ivar, xtype = nf_type, ndims = nvdims, dimids = dimids)
      call Check_status (status)

      if (nf_type == NF90_INT) then
        if (nvdims == 3) then
          do i = 1, nvdims
            status = nf90_inquire_dimension (ncid, dimids(i), len = idims(i))
            call Check_status (status)
          end do

          n1 = idims(1)
          n2 = idims(2)
          n3 = idims(3)
          allocate (output(n1, n2, n3), stat = io_stat)
          if (io_stat /= 0) call Stop_simulation ('Problems allocating output variable')

          status = nf90_get_var (ncid, ivar, output)
          call Check_status (status)
        else
          msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a 3D array'
          Call Stop_simulation (msg)
        end if
      else
        msg = 'FATAL ERROR: ' // trim (var_name) // ' is not an int32 variable'
        call Stop_simulation (msg)
      end if

        ! Closing file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_var_int32_3d

    subroutine Get_netcdf_var_real32_2d (file_name, var_name, output)

      use netcdf
      use, intrinsic :: iso_fortran_env, only : REAL32

      implicit none

      character (len = *), intent(in) :: file_name, var_name
      real (kind = REAL32), dimension (:, :), allocatable :: output

      integer :: status, ncid, ivar, nf_type, nvdims
      integer :: i, n1, n2, io_stat
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimids, idims
      character (len = :), allocatable :: msg


        ! Open file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get var
      status = nf90_inq_varid (ncid, trim (var_name), ivar)
      call Check_status (status)
      status = nf90_inquire_variable (ncid, ivar, xtype = nf_type, ndims = nvdims, dimids = dimids)
      call Check_status (status)

      if (nf_type == NF90_FLOAT) then
        if (nvdims == 2) then
          do i = 1, nvdims
            status = nf90_inquire_dimension (ncid, dimids(i), len = idims(i))
            call Check_status (status)
          end do

          n1 = idims(1)
          n2 = idims(2)
          allocate (output(n1, n2), stat = io_stat)
          if (io_stat /= 0) Call Stop_simulation ('Problems allocating output variable')

          status = nf90_get_var (ncid, ivar, output)
          call Check_status (status)
        else
          msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a 2D array'
          call Stop_simulation (msg)
        end if
      else
        msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a real32 variable'
        call Stop_simulation (msg)
      end if

        ! Closing file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_var_real32_2d

    subroutine Get_netcdf_var_real32_3d (file_name, var_name, output)

      use netcdf
      use, intrinsic :: iso_fortran_env, only : REAL32

      implicit none

      character (len = *), intent(in) :: file_name, var_name
      real (kind = REAL32), dimension (:, :, :), allocatable :: output

      integer :: status, ncid, ivar, nf_type, nvdims
      integer :: i, n1, n2, n3, io_stat
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimids, idims
      character (len = :), allocatable :: msg


        ! Open file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get var
      status = nf90_inq_varid (ncid, trim (var_name), ivar)
      call Check_status (status)
      status = nf90_inquire_variable (ncid, ivar, xtype = nf_type, ndims = nvdims, dimids = dimids)
      call Check_status (status)

      if (nf_type == NF90_FLOAT) then
        if (nvdims == 3) then
          do i = 1, nvdims
            status = nf90_inquire_dimension (ncid, dimids(i), len = idims(i))
            call Check_status (status)
          end do

          n1 = idims(1)
          n2 = idims(2)
          n3 = idims(3)
          allocate (output(n1, n2, n3), stat = io_stat)
          if (io_stat /= 0) call Stop_simulation ('Problems allocating output variable')

          status = nf90_get_var (ncid, ivar, output)
          call Check_status (status)
        else
          msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a 3D array'
          call Stop_simulation (msg)
        end if
      else
        msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a real32 variable'
        call Stop_simulation (msg)
      end if

        ! Closing file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_var_real32_3d

    subroutine Get_netcdf_var_real32_4d (file_name, var_name, output)

      use netcdf
      use, intrinsic :: iso_fortran_env, only : REAL32

      implicit none

      character (len = *), intent(in) :: file_name, var_name
      real (kind = REAL32), dimension (:, :, :, :), allocatable :: output

      integer :: status, ncid, ivar, nf_type, nvdims
      integer :: i, n1, n2, n3, n4, io_stat
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimids, idims
      character (len = :), allocatable :: msg


        ! Open file
      status = nf90_open (trim(file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

        ! Get var
      status = nf90_inq_varid (ncid, trim (var_name), ivar)
      call Check_status (status)
      status = nf90_inquire_variable (ncid, ivar, xtype = nf_type, ndims = nvdims, dimids = dimids)
      call Check_status (status)

      if (nf_type == NF90_FLOAT) then
        if (nvdims == 4) then
          do i = 1, nvdims
            status = nf90_inquire_dimension (ncid, dimids(i), len = idims(i))
            call Check_status (status)
          end do

          n1 = idims(1)
          n2 = idims(2)
          n3 = idims(3)
          n4 = idims(4)
          allocate (output(n1, n2, n3, n4), stat = io_stat)
          if (io_stat /= 0) Call Stop_simulation ('Problems allocating output variable')

          status = nf90_get_var (ncid, ivar, output)
          call Check_status (status)
        else
          msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a 4D array'
          call Stop_simulation (msg)
        end if
      else
        msg = 'FATAL ERROR: ' // trim (var_name) // ' is not a real32 variable'
        call Stop_simulation (msg)
      end if

        ! Closing file
      status = nf90_close (ncid)
      call Check_status (status)

    end subroutine Get_netcdf_var_real32_4d

    subroutine Is_netcdf_file_present (file_name)

      use netcdf

      implicit none

      character (len = *), intent (in) :: file_name

      integer :: ncidout, status
      character (len = :), allocatable :: msg


      status = nf90_open (trim(file_name), NF90_NOWRITE, ncidout)
      if (status /= NF90_NOERR) then
        msg = 'Problems opening file ' // trim (file_name)
        call Print_message (msg)
        call Stop_simulation (trim (nf90_strerror (status)))
      end if

      status = nf90_close (ncidout)
      call Check_status (status)

    end subroutine Is_netcdf_file_present

    function Is_netcdf_var_present (file_name, name_var) result (return_value)

      use netcdf

      implicit none

      character (len = *), intent (in) :: file_name, name_var
      logical :: return_value

      integer :: status, ncid, varid


      status = nf90_open (trim (file_name), NF90_NOWRITE, ncid)
      call Check_status (status)

      status = nf90_inq_varid (ncid, trim (name_var), varid)
      if (status == NF90_NOERR) then
        return_value = .true.
      else
        return_value = .false.
      end if

      status = nf90_close (ncid)
      call Check_status (status)

    end function Is_netcdf_var_present

  end module netcdf_mod
