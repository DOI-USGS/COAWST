module NoahmpReadLandMod

  use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE, C_PTR, C_CHAR
  use netcdf
  use Machine
  use NoahmpIOVarType

  implicit none

  public :: NoahmpReadLandHeader, NoahmpReadLandMain

  private :: FATAL, NOT_FATAL, get_2d_netcdf, error_handler, get_landuse_netcdf, &
             get_soilcat_netcdf, get_netcdf_soillevel, init_interp, get_2d_netcdf_cfloat, &
             get_2d_netcdf_ffloat

  logical, parameter :: FATAL = .TRUE.
  logical, parameter :: NOT_FATAL = .FALSE.

  interface get_2d_netcdf
    module procedure get_2d_netcdf_cfloat
    module procedure get_2d_netcdf_ffloat
  end interface get_2d_netcdf

contains

subroutine NoahmpReadLandHeader(NoahmpIO)

    implicit none
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

    integer :: ncid, dimid, varid, ierr
    real, allocatable, dimension(:,:) :: dum2d
    character(len=256) :: units
    integer :: i
    integer :: rank
    integer :: ilev, is, js, ratio, xoffset, yoffset

    if (NoahmpIO%rank == 0) write(*,'("Noah-MP reading ''", A, "'' headers")') trim(NoahmpIO%erf_setup_file_lev)

    ierr = nf90_open(NoahmpIO%erf_setup_file_lev, NF90_NOWRITE, ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_lev))

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION", NoahmpIO%xsglobal)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'WEST-EAST_GRID_DIMENSION'")
    NoahmpIO%xsglobal = NoahmpIO%xsglobal-1

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", NoahmpIO%ysglobal)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'SOUTH-NORTH_GRID_DIMENSION'")
    NoahmpIO%ysglobal = NoahmpIO%ysglobal-1

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "DX", NoahmpIO%dx)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'DX'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "DY", NoahmpIO%dy)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'DY'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", NoahmpIO%truelat1)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'TRUELAT1'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2", NoahmpIO%truelat2)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'TRUELAT2'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", NoahmpIO%cen_lon)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'STAND_LON'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ", NoahmpIO%mapproj)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'MAP_PROJ'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "GRID_ID", NoahmpIO%igrid)
    if (ierr /= 0) then
       ierr = nf90_get_att(ncid, NF90_GLOBAL, "grid_id", NoahmpIO%igrid)
       call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'GRID_ID' or 'grid_id'")
    endif

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISWATER", NoahmpIO%iswater)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISWATER'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISLAKE", NoahmpIO%islake)
    if(ierr /= 0) then
      if (NoahmpIO%rank == 0) write(*,*) "Problems finding global attribute: ISLAKE; setting to -1"
      NoahmpIO%islake = -1
    end if

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISURBAN", NoahmpIO%isurban)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISURBAN'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "ISICE", NoahmpIO%isice)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'ISICE'")

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "MMINLU", NoahmpIO%llanduse)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'MMINLU'")
    
    ! IBM XLF seems to need something like this:
    do i = 1, 256
       if (ichar(NoahmpIO%llanduse(i:i)) == 0) NoahmpIO%llanduse(i:i) = " "
    enddo

    allocate(dum2d(NoahmpIO%xstart:NoahmpIO%xstart,NoahmpIO%ystart:NoahmpIO%ystart))
    call get_2d_netcdf("XLAT", ncid, dum2d,  units, NoahmpIO%xstart, NoahmpIO%xstart, NoahmpIO%ystart, NoahmpIO%ystart, FATAL, ierr)
    NoahmpIO%lat1 = dum2d(NoahmpIO%xstart,NoahmpIO%ystart)

    call get_2d_netcdf("XLONG", ncid, dum2d,  units, NoahmpIO%xstart, NoahmpIO%xstart, NoahmpIO%ystart, NoahmpIO%ystart, FATAL, ierr)
    NoahmpIO%lon1 = dum2d(NoahmpIO%xstart,NoahmpIO%ystart)
    deallocate (dum2d)

    ierr = nf90_close(ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO:  Problems closing NetCDF file.")

    NoahmpIO%xoffset = 0
    NoahmpIO%yoffset = 0

    do ilev=0,NoahmpIO%level

      select case(ilev)
      case(0)
        ierr = nf90_open(NoahmpIO%erf_setup_file_01, NF90_NOWRITE, ncid)
        call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_01))
      case(1)
        ierr = nf90_open(NoahmpIO%erf_setup_file_02, NF90_NOWRITE, ncid)
        call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_02))
      case(2)
        ierr = nf90_open(NoahmpIO%erf_setup_file_03, NF90_NOWRITE, ncid)
        call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_03)) 
      case default
        print *, "Error: unsupported level: ", ilev
        stop
      end select

      ierr = nf90_get_att(ncid, NF90_GLOBAL, "I_PARENT_START", is)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'I_PARENT_START'")

      ierr = nf90_get_att(ncid, NF90_GLOBAL, "J_PARENT_START", js)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'J_PARENT_START'")

      ierr = nf90_get_att(ncid, NF90_GLOBAL, "PARENT_GRID_RATIO", ratio)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems finding global attribute 'PARENT_GRID_RATIO'")

      ierr = nf90_close(ncid)
      call error_handler(ierr, "READ_ERF_HDRINFO:  Problems closing NetCDF file.")        

      NoahmpIO%xoffset = ratio*(NoahmpIO%xoffset+is-1)
      NoahmpIO%yoffset = ratio*(NoahmpIO%yoffset+js-1)

    end do

end subroutine NoahmpReadLandHeader

subroutine NoahmpReadLandMain(NoahmpIO)
    implicit none
    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

    character(len=256) :: units
    integer :: ierr
    integer :: ncid
    real, dimension(NoahmpIO%xstart-NoahmpIO%xoffset:NoahmpIO%xend-NoahmpIO%xoffset, \
                    NoahmpIO%ystart-NoahmpIO%yoffset:NoahmpIO%yend-NoahmpIO%yoffset) :: xdum
    integer :: rank

    character(len=256) :: titlestr
    character(len=8)   :: name
    character(len=256) :: llanduse

    integer :: ierr_snodep, varid
    integer :: idx, isoil
    real, dimension(100) :: layer_bottom
    real, dimension(100) :: layer_top
    real, dimension(NoahmpIO%nsoil)   :: dzs

    real, dimension(NoahmpIO%xstart-NoahmpIO%xoffset:NoahmpIO%xend-NoahmpIO%xoffset, \
                    NoahmpIO%ystart-NoahmpIO%yoffset:NoahmpIO%yend-NoahmpIO%yoffset, NoahmpIO%nsoil) :: insoil

    real, dimension(NoahmpIO%xstart-NoahmpIO%xoffset:NoahmpIO%xend-NoahmpIO%xoffset, \
                    NoahmpIO%nsoil, \
                    NoahmpIO%ystart-NoahmpIO%yoffset:NoahmpIO%yend-NoahmpIO%yoffset) :: soildummy

    integer :: ierr_vegfra
    integer :: ierr_lai

    integer :: i, j
    integer :: iret
    integer :: xstart, ystart, xend, yend

    if (NoahmpIO%rank == 0) write(*,'("Noah-MP reading ''", A, "'' variables")') trim(NoahmpIO%erf_setup_file_lev)

    ierr = nf90_open(NoahmpIO%erf_setup_file_lev, NF90_NOWRITE, ncid)
    call error_handler(ierr, "READ_ERF_HDRINFO: Problem opening wrfinput file: "//trim(NoahmpIO%erf_setup_file_lev))

    xstart = NoahmpIO%xstart-NoahmpIO%xoffset
    ystart = NoahmpIO%ystart-NoahmpIO%yoffset

    xend = NoahmpIO%xend-NoahmpIO%xoffset
    yend = NoahmpIO%yend-NoahmpIO%yoffset

    ! Get Latitude (lat)
    call get_2d_netcdf("XLAT", ncid, NoahmpIO%xlat,  units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Get Longitude (lon)
    call get_2d_netcdf("XLONG", ncid, NoahmpIO%xlong, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Get land mask (xland)
    call get_2d_netcdf("XLAND", ncid, NoahmpIO%xland, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)

    ! Get seaice (seaice)
    call get_2d_netcdf("SEAICE", ncid, NoahmpIO%seaice, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)

    ! Get Terrain (avg)
    call get_2d_netcdf("HGT", ncid, NoahmpIO%terrain, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Get Deep layer temperature (TMN)
    call get_2d_netcdf("TMN", ncid, NoahmpIO%TMN, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Get Map Factors (MAPFAC_MX)
    call get_2d_netcdf("MAPFAC_MX", ncid, NoahmpIO%msftx, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)
    if (ierr /= 0) print*, 'Did not find MAPFAC_MX, only needed for iopt_run=5'

    ! Get Map Factors (MAPFAC_MY)
    call get_2d_netcdf("MAPFAC_MY", ncid, NoahmpIO%msfty, units, xstart, xend, ystart, yend, NOT_FATAL, ierr)
    if (ierr /= 0) print*, 'Did not find MAPFAC_MY, only needed for iopt_run=5'

    ! Get Dominant Land Use categories (use)
    call get_landuse_netcdf(ncid, xdum , units, xstart, xend, ystart, yend)
    NoahmpIO%ivgtyp = nint(xdum)

    ! Get Dominant Soil Type categories in the top layer (stl)
    call get_soilcat_netcdf(ncid, xdum , units, xstart, xend, ystart, yend)
    NoahmpIO%isltyp = nint(xdum)

    where (NoahmpIO%SEAICE > 0.0) NoahmpIO%XICE = 1.0
 
    NoahmpIO%CROPTYPE   = 0       ! make default 0% crops everywhere

    NoahmpIO%TD_FRACTION = 0.0

    NoahmpIO%SLOPETYP  =  1 ! it was 2 here and 1 in the noahmpdrv- pvk
    NoahmpIO%DZS       =  NoahmpIO%SOIL_THICK_INPUT(1:NoahmpIO%NSOIL)
    NoahmpIO%ITIMESTEP = 1
    NoahmpIO%restart_flag = .false.

    ierr = nf90_get_att(ncid, NF90_GLOBAL, "MMINLU", llanduse)
    if (ierr /= 0) then
       if (NoahmpIO%rank == 0) write(*,'("WARNING:  Noah-MP input file does not have MMINLU attribute.")')
       if (NoahmpIO%rank == 0) write(*,'("          This probably means that the file  is from an older release.")')
       if (NoahmpIO%rank == 0) write(*,'("          I assume you know what you are doing.")')
    else
       if (NoahmpIO%rank == 0) write(*,'("MMNINLU attribute: ", A)') llanduse
    endif

    call get_2d_netcdf("CANWAT", ncid, NoahmpIO%canwat, units, xstart, xend, ystart, yend, FATAL, ierr)
    call get_2d_netcdf("TSK",    ncid, NoahmpIO%tsk, units, xstart, xend, ystart, yend, FATAL, ierr)
    call get_2d_netcdf("SNOW",   ncid, NoahmpIO%snow, units, xstart, xend, ystart, yend, FATAL, ierr)
    call get_2d_netcdf("SNOWC",  ncid, NoahmpIO%snowc, units, xstart, xend, ystart, yend, FATAL, ierr)

    NoahmpIO%snowh = 0.0
    call get_2d_netcdf("SNOWH", ncid, NoahmpIO%snowh, units, xstart, xend, ystart, yend, NOT_FATAL, ierr_snodep)
    NoahmpIO%fndsnowh = .true.
    if (ierr_snodep /= 0) NoahmpIO%fndsnowh = .false.
   
    ierr = nf90_inq_varid(ncid,  "DZS",  varid)
    ierr = nf90_get_var(ncid, varid, values=dzs, start=(/1/), count=(/NoahmpIO%nsoil/))    

    layer_top(1) = 0.0
    layer_bottom(1) = dzs(1)
    do isoil = 2, NoahmpIO%nsoil
      layer_top(isoil) = layer_bottom(isoil-1)
      layer_bottom(isoil) = layer_top(isoil) + dzs(isoil)
    end do

    call get_netcdf_soillevel("TSLB", ncid, NoahmpIO%nsoil, soildummy, units,  xstart, xend, ystart, yend, FATAL, ierr)
    
    !if (NoahmpIO%rank == 0) write(*, '("layer_bottom(1:nsoil) = ", 4F9.4)') layer_bottom(1:NoahmpIO%nsoil)
    !if (NoahmpIO%rank == 0) write(*, '("layer_top(1:nsoil)    = ", 4F9.4)') layer_top(1:NoahmpIO%nsoil)
    !if (NoahmpIO%rank == 0) write(*, '("Soil depth = ", 10F12.6)') NoahmpIO%dzs
    
    call init_interp(NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, NoahmpIO%nsoil, &
                     NoahmpIO%dzs, NoahmpIO%tslb, NoahmpIO%nsoil, soildummy, layer_bottom(1:NoahmpIO%nsoil), layer_top(1:NoahmpIO%nsoil), NoahmpIO%rank)

    call get_netcdf_soillevel("SMOIS", ncid, NoahmpIO%nsoil, soildummy, units,  xstart, xend, ystart, yend, FATAL, ierr)

    call init_interp(NoahmpIO%xstart, NoahmpIO%xend, NoahmpIO%ystart, NoahmpIO%yend, NoahmpIO%nsoil, \
    NoahmpIO%dzs, NoahmpIO%smois, NoahmpIO%nsoil, soildummy, layer_bottom(1:NoahmpIO%nsoil), layer_top(1:NoahmpIO%nsoil), NoahmpIO%rank)

    NoahmpIO%VEGFRA =  0.0

    call get_2d_netcdf("VEGFRA", ncid, NoahmpIO%vegfra, units, xstart, xend, ystart, yend, NOT_FATAL, ierr_vegfra)
    call get_2d_netcdf("LAI", ncid, NoahmpIO%lai, units, xstart, xend, ystart, yend, NOT_FATAL, ierr_lai)

    ! Get Minimum Green Vegetation Fraction SHDMIN
    call get_2d_netcdf("SHDMIN", ncid, NoahmpIO%gvfmin, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Get Minimum Green Vegetation Fraction SHDMAX
    call get_2d_netcdf("SHDMAX", ncid, NoahmpIO%gvfmax, units, xstart, xend, ystart, yend, FATAL, ierr)

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= 0) stop "MODULE_NOAHLSM_ERF_INPUT:  READLAND_ERF:  NF90_CLOSE"

end subroutine NoahmpReadLandMain

subroutine get_2d_netcdf_cfloat(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real(c_double), dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    ! FATAL_IF_ERROR:  an input code value:
    !      .TRUE. if an error in reading the data should stop the program.
    !      Otherwise the, IERR error flag is set, but the program continues.
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr
 
    units = " "
    
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid = ', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file.")
       else
          ierr = iret
         return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=array, start=(/xstart+1,ystart+1/), count=(/xend-xstart+1,yend-ystart+1/))

    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid =', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem retrieving variable '"//trim(name)//"' from NetCDF file.")
       else
          ierr = iret
          return
       endif
    endif

    ierr = 0;

end subroutine get_2d_netcdf_cfloat

subroutine get_2d_netcdf_ffloat(name, ncid, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=*), intent(out) :: units
    integer :: iret, varid
    ! FATAL_IF_ERROR:  an input code value:
    !      .TRUE. if an error in reading the data should stop the program.
    !      Otherwise the, IERR error flag is set, but the program continues.
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr
 
    units = " "
    
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid = ', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file.")
       else
          ierr = iret
         return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=array, start=(/xstart+1,ystart+1/), count=(/xend-xstart+1,yend-ystart+1/))

    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'ncid =', ncid
          call error_handler(iret, "MODULE_ERF_NETCDF_IO:  Problem retrieving variable '"//trim(name)//"' from NetCDF file.")
       else
          ierr = iret
          return
       endif
    endif

    ierr = 0;

end subroutine get_2d_netcdf_ffloat


subroutine error_handler(status, failure, success)
    !
    ! Check the error flag from a NetCDF function call, and print appropriate
    ! error message.
    !
    implicit none
    integer,                    intent(in) :: status
    character(len=*), optional, intent(in) :: failure
    character(len=*), optional, intent(in) :: success

    if (status .ne. NF90_NOERR) then
       write(*,'(/,A)') nf90_strerror(status)
       if (present(failure)) then
          write(*,'(/," ***** ", A,/)') failure
       endif
       stop 'Stopped'
    endif

    if (present(success)) then
       write(*,'(A)') success
    endif
end subroutine error_handler

subroutine get_landuse_netcdf(ncid, array, units, xstart, xend, ystart, yend)
    implicit none
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    character(len=24), parameter :: name = "IVGTYP"

    units = " "

    iret = nf90_inq_varid(ncid,  trim(name),  varid)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       stop "MODULE_NOAHLSM_ERF_INPUT:  get_landuse_netcdf:  nf90_inq_varid"
    endif

    iret = nf90_get_var(ncid, varid, array, (/xstart+1, ystart+1/), (/xend-xstart+1, yend-ystart+1/))
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       stop "MODULE_NOAHLSM_ERF_INPUT:  get_landuse_netcdf:  nf90_get_var"
    endif
end subroutine get_landuse_netcdf

subroutine get_soilcat_netcdf(ncid, array, units, xstart, xend, ystart, yend)
    implicit none
    integer, intent(in) :: ncid
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,ystart:yend), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    character(len=24), parameter :: name = "ISLTYP"

    units = " "

    iret = nf90_inq_varid(ncid,  trim(name),  varid)
    call error_handler(iret, "Problem finding variable '"//trim(name)//"' in the wrfinput file.")

    iret = nf90_get_var(ncid, varid, array, (/xstart+1, ystart+1/), (/xend-xstart+1, yend-ystart+1/))
    call error_handler(iret, "Problem retrieving variable "//trim(name)//" from the wrfinput file.")

end subroutine get_soilcat_netcdf

subroutine get_netcdf_soillevel(name, ncid, nsoil, array, units, xstart, xend, ystart, yend, fatal_if_error, ierr)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: nsoil
    integer, intent(in) :: xstart, xend, ystart, yend
    real, dimension(xstart:xend,nsoil,ystart:yend), intent(out) :: array
    character(len=256), intent(out) :: units
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr

    integer :: iret, varid, isoil
    real:: insoil(xstart:xend,ystart:yend,nsoil)

    units = " "

    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'name = "', trim(name)//'"'
          stop "MODULE_NOAHLSM_HRLDAS_INPUT:  get_2d_netcdf:  nf90_inq_varid"
       else
          ierr = iret
          return
       endif
    endif

    iret = nf90_get_att(ncid, varid, "units", units)
    if (iret /= 0) units = "units unknown"

    iret = nf90_get_var(ncid, varid, values=insoil, start=(/xstart+1,ystart+1,1,1/), count=(/xend-xstart+1,yend-ystart+1,nsoil,1/))
    do isoil = 1,nsoil
      array(:,isoil,:) = insoil(:,:,isoil)
    end do

    if (iret /= 0) then
       if (FATAL_IF_ERROR) then
          print*, 'name = "', trim(name)//'"'
          print*, 'varid =', varid
          print*, trim(nf90_strerror(iret))
          stop "MODULE_NOAHLSM_HRLDAS_INPUT:  get_2d_netcdf:  nf90_get_var"
       else
          ierr = iret
          return
       endif
    endif

    ierr = 0;
end subroutine get_netcdf_soillevel

subroutine init_interp(xstart, xend, ystart, yend, nsoil, sldpth, var, nvar, src, layer_bottom, layer_top, rank)
    implicit none
    integer, intent(in)    :: xstart, xend, ystart, yend, nsoil, nvar
    real, dimension(nsoil) :: sldpth ! the thickness of each layer
    real, dimension(xstart:xend, nsoil, ystart:yend), intent(out) :: var
    real, dimension(xstart:xend, nvar, ystart:yend ), intent(in)  :: src
    real, dimension(nvar),                            intent(in)  :: layer_bottom ! The depth from the surface of each layer bottom.
    real, dimension(nvar),                            intent(in)  :: layer_top    ! The depth from the surface of each layer top.
    integer :: i, j, k, kk, ktop, kbottom
    real, dimension(nsoil) :: dst_centerpoint
    real, dimension(nvar)  :: src_centerpoint
    real :: fraction
    integer :: ierr
    integer, intent(in) :: rank

    do k = 1, nsoil
       if (k==1) then
          dst_centerpoint(k) = sldpth(k)/2.
       else
          dst_centerpoint(k) = sldpth(k)/2. + sum(sldpth(1:k-1))
       endif
       !if (rank == 0) print*, 'k, dst_centerpoint(k) = ', k, dst_centerpoint(k)
    enddo
    print*

    do k = 1, nvar
       src_centerpoint(k) = 0.5*(layer_bottom(k)+layer_top(k))
       !if (rank == 0) print*, 'k, src_centerpoint(k) = ', k, src_centerpoint(k)
    enddo

    KLOOP : do k = 1, nsoil

       if (dst_centerpoint(k) < src_centerpoint(1)) then
          ! If the center of the destination layer is closer to the surface than
          ! the center of the topmost source layer, then simply set the 
          ! value of the destination layer equal to the topmost source layer:
          !if (rank == 0) then
          !   print'("Shallow destination layer:  Taking destination layer at ",F7.4, " from source layer at ", F7.4)', &
          !        dst_centerpoint(k), src_centerpoint(1)
          !endif
          var(:,k,:) = src(:,1,:)
          cycle KLOOP
       endif

       if (dst_centerpoint(k) > src_centerpoint(nvar)) then
          ! If the center of the destination layer is deeper than
          ! the center of the deepest source layer, then simply set the 
          ! value of the destination layer equal to the deepest source layer:
          !if (rank == 0) then
          !   print'("Deep destination layer:  Taking destination layer at ",F7.4, " from source layer at ", F7.4)', &
          !        dst_centerpoint(k), src_centerpoint(nvar)
          !endif
          var(:,k,:) = src(:,nvar,:)
          cycle KLOOP
       endif

       ! Check if the center of the destination layer is "close" to the center
       ! of a source layer.  If so, simply set the value of the destination layer
       ! equal to the value of that close soil layer:
       do kk = 1, nvar
          if (abs(dst_centerpoint(k)-src_centerpoint(kk)) < 0.01) then
             !if (rank == 0) then
             !   print'("(Near) match for destination layer:  Taking destination layer at ",F7.4, " from source layer at ", F7.4)', &
             !        dst_centerpoint(k), src_centerpoint(kk)
             !endif
             var(:,k,:) = src(:,kk,:)
             cycle KLOOP
          endif
       enddo

       ! Otherwise, do a linear interpolation

       ! Get ktop, the index of the top bracketing layer from the source dataset.
       ! Which from the bottom up, will be the first source level that is closer 
       ! to the surface than the destination level
       ktop = -99999
       TOPLOOP : do kk = nvar,1,-1
          if (src_centerpoint(kk) < dst_centerpoint(k)) then
             ktop = kk
             exit TOPLOOP
          endif
       enddo TOPLOOP
       if (ktop < -99998) stop "ktop problem"



       ! Get kbottom, the index of the bottom bracketing layer from the source dataset.
       ! Which, from the top down, will be the first source level that is deeper than
       ! the destination level
       kbottom = -99999
       BOTTOMLOOP : do kk = 1, nvar
          if ( src_centerpoint(kk) > dst_centerpoint(k) ) then
             kbottom = kk
             exit BOTTOMLOOP
          endif
       enddo BOTTOMLOOP
       if (kbottom < -99998) stop "kbottom problem"

       fraction = (src_centerpoint(kbottom)-dst_centerpoint(k)) / (src_centerpoint(kbottom)-src_centerpoint(ktop))

       ! print '(I2, 1x, 3F7.3, F8.5)', k, src_centerpoint(ktop), dst_centerpoint(k), src_centerpoint(kbottom), fraction

       !if (rank == 0) then
       !   print '("dst(",I1,") = src(",I1,")*",F8.5," + src(",I1,")*",F8.5)', k, ktop, fraction, kbottom, (1.0-fraction)
       !endif

       var(:,k,:) = (src(:,ktop,:)*fraction) + (src(:,kbottom,:)*(1.0-fraction))

    enddo KLOOP     
end subroutine init_interp


end module NoahmpReadLandMod
