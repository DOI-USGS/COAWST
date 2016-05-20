!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines reading ROMS grids  
!
!-----------------------------------------------------------------------
!     This is a wrapper around SCRIP package developed at USGS, 
!     Woods Hole Science And Marine Center. This routine calls the 
!     SCRIP's package routine that is the driver for computing the 
!     addresses and weights for interpolating between two grids on 
!     a sphere.
!
!---- Written by John C. Warner-------------------------------------- 
!-----         Tarandeep S. Kalra -------------------------------------
!--------------Date: 08/04/2015----------------------------------------
!***********************************************************************

      module read_roms

!-----------------------------------------------------------------------
      use kinds_mod      ! defines common data types 
      use scripwrap_mod  ! defines input file variable names for SCRIP-COAWST wrapper
      use create_fullgrid! call the subroutine to create full grid of psi points 
      use netcdf_mod     !netCDF include file and a netcdf error handling routine.

      contains 

      subroutine load_roms_grid()

      implicit none

!      include 'netcdf.inc'

      integer(int_kind) :: i, j, iunit
      integer(int_kind) :: nx, ny, mo 
      integer(int_kind)  :: ncstat, nc_file_id, nc_grdsize_id, 
     &                      nc_grdlat_id, nc_grdlon_id, 
     &                      nc_grdcrnrlat_id, nc_grdcrnrlon_id,
     &                      nc_grdmsk_id

      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max, scale
      real(dbl_kind), allocatable :: lon_psi_o(:,:), lat_psi_o(:,:)

      integer(int_kind) :: checksphere_id
      integer(int_kind) :: my_type
      integer(int_kind) :: AI
      integer(int_kind) :: RHTYPE            ! variable type
      integer(int_kind) :: RHN               ! number of dimensions
      integer(int_kind) :: RHDIMS(1)         ! variable shape
      integer(int_kind) :: RHNATT            ! number of attributes
      logical :: spherical
      character(char_len) :: RHNAME          ! variable name
      character(len=1) :: Achr

      allocate(ngrd_rm(Ngrids_roms))

      do mo=1,Ngrids_roms

!     Open the file. 
        ncstat=nf_open(roms_grids(mo),nf_nowrite,nc_file_id)
        call netcdf_error_handler(ncstat)

!     Read dimension id 
        ncstat=nf_inq_dimid(nc_file_id,'xi_rho',nc_grdsize_id) 
        call netcdf_error_handler(ncstat)
!     Get the grid size in each direction 
        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,                  
     &                            ngrd_rm(mo)%xi_size)
        call netcdf_error_handler(ncstat)

        ncstat=nf_inq_dimid(nc_file_id,'eta_rho',nc_grdsize_id) 
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,
     &                           ngrd_rm(mo)%eta_size)
        call netcdf_error_handler(ncstat)

!     Read spherical variable (it can be True/False or 1/0)
        ncstat=nf_inq_varid(nc_file_id, 'spherical', checksphere_id)
        call netcdf_error_handler(ncstat)

        ncstat=nf_inq_var(nc_file_id, checksphere_id,  RHNAME, RHTYPE,
     &                    RHN, RHDIMS, RHNATT)
        call netcdf_error_handler(ncstat)

!       The valid netCDF external data types are NF_BYTE, NF_CHAR, 
!       NF_SHORT, NF_INT, NF_FLOAT, AND NF_DOUBLE
!       RHTYPE gets the type of data 
!       RHTYPE =2 refers to NF_CHAR and RHTYPE=4 refers to NF_INT 

        if (RHTYPE.eq.2) then
          my_type=nf_char 
        else if (RHTYPE.eq.4) then 
          my_type=nf_int
        endif

!    Use scale to convert m to degrees, if spherical=T/t/1
        scale=1.
        if (my_type.eq.nf_int) then
          ncstat=nf_get_var_int(nc_file_id, checksphere_id, AI)
          call netcdf_error_handler(ncstat)
            if (AI.eq.0) then
              spherical=.FALSE.
              scale=1./6371000.
            else
              spherical=.TRUE.
              scale=1.
            endif
        else if (my_type.eq.nf_char) then
          ncstat=nf_get_var_text(nc_file_id, checksphere_id, Achr)
          call netcdf_error_handler(ncstat)
            if ((Achr.eq."t").or.(Achr.eq."T")) then
              spherical=.TRUE.
              scale=1.
            else 
              spherical=.FALSE.
              scale=1./6371000.
            endif
        endif 
!         
!     Allocate arrays
!
        allocate(ngrd_rm(mo)%
     &           lon_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%
     &           lat_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%
     &           mask_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%x_full_grid(ngrd_rm(mo)%xi_size+1,
     &                                   ngrd_rm(mo)%eta_size+1))
        allocate(ngrd_rm(mo)%y_full_grid(ngrd_rm(mo)%xi_size+1
     &                                  ,ngrd_rm(mo)%eta_size+1))
!   
!       Making lon_psi and lat_psi as local arrays 
!       Local arrays need to deallocated, no_psi_pts=no_rho_pts-1
!
        allocate(lon_psi_o(ngrd_rm(mo)%xi_size-1,
     &                     ngrd_rm(mo)%eta_size-1))
        allocate(lat_psi_o(ngrd_rm(mo)%xi_size-1,
     &                     ngrd_rm(mo)%eta_size-1))


        if (.not.spherical) then 
          ncstat=nf_inq_varid(nc_file_id,'x_rho',nc_grdlon_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'y_rho',nc_grdlat_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'x_psi',nc_grdcrnrlon_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'y_psi',nc_grdcrnrlat_id)
          call netcdf_error_handler(ncstat)
        else
          ncstat=nf_inq_varid(nc_file_id,'lon_rho',nc_grdlon_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'lat_rho',nc_grdlat_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'lon_psi',nc_grdcrnrlon_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'lat_psi',nc_grdcrnrlat_id)
          call netcdf_error_handler(ncstat)
        endif 
        ncstat=nf_inq_varid(nc_file_id,'mask_rho',nc_grdmsk_id)
        call netcdf_error_handler(ncstat)

!     Get variables
        ncstat=nf_get_var_double(nc_file_id, nc_grdlon_id,
     &                                     ngrd_rm(mo)%lon_rho_o)
        call netcdf_error_handler(ncstat)
        ngrd_rm(mo)%lon_rho_o=ngrd_rm(mo)%lon_rho_o*scale

        ncstat=nf_get_var_double(nc_file_id, nc_grdlat_id,
     &                                     ngrd_rm(mo)%lat_rho_o)
        call netcdf_error_handler(ncstat)
        ngrd_rm(mo)%lat_rho_o=ngrd_rm(mo)%lat_rho_o*scale

        ncstat=nf_get_var_int(nc_file_id, nc_grdmsk_id,
     &                                     ngrd_rm(mo)%mask_rho_o)
        call netcdf_error_handler(ncstat)

        ncstat=nf_get_var_double(nc_file_id, nc_grdcrnrlon_id,
     &                                       lon_psi_o)
        call netcdf_error_handler(ncstat)
        lon_psi_o=lon_psi_o*scale

        ncstat=nf_get_var_double(nc_file_id, nc_grdcrnrlat_id,
     &                                       lat_psi_o)
        call netcdf_error_handler(ncstat)
        lat_psi_o=lat_psi_o*scale
 
        call create_extra_rho_grid(ngrd_rm(mo)%xi_size-1,
     &                             ngrd_rm(mo)%eta_size-1,
     &                             lon_psi_o, lat_psi_o,
     &                             ngrd_rm(mo)%x_full_grid,
     &                             ngrd_rm(mo)%y_full_grid)

        deallocate(lon_psi_o, lat_psi_o)
!
      end do

!   Find child grid with respect to parent grid
!
      do mo=1,Ngrids_roms-1
        allocate(ngrd_rm(mo)%istr_o,ngrd_rm(mo)%jstr_o,
     &           ngrd_rm(mo)%iend_o,ngrd_rm(mo)%jend_o)

        nx=ngrd_rm(mo)%xi_size
        ny=ngrd_rm(mo)%eta_size 
        dist_max=10e6
!       First the (2,2) point 
        xx2=ngrd_rm(mo+1)%lon_rho_o(2,2)
        yy2=ngrd_rm(mo+1)%lat_rho_o(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_rm(mo)%lon_rho_o(i,j)-xx2)**2+
     &                 (ngrd_rm(mo)%lat_rho_o(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_rm(mo)%istr_o=i
              ngrd_rm(mo)%jstr_o=j
            endif
          enddo
        enddo
        dist_max=10e6
!       The second last point 
        xxend_1=ngrd_rm(mo+1)%lon_rho_o(ngrd_rm(mo+1)%xi_size-1,        
     &                                ngrd_rm(mo+1)%eta_size-1 )
        yyend_1=ngrd_rm(mo+1)%lat_rho_o(ngrd_rm(mo+1)%xi_size-1,
     &                                ngrd_rm(mo+1)%eta_size-1 )
        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_rm(mo)%lon_rho_o(i,j)-xxend_1)**2+
     &                 (ngrd_rm(mo)%lat_rho_o(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_rm(mo)%iend_o=i
              ngrd_rm(mo)%jend_o=j
            endif
          enddo
        enddo
         print*,"ROMS starting i & j index of parent w.r.t child grid"
         print*,ngrd_rm(mo)%istr_o, ngrd_rm(mo)%jstr_o
         print*,"ending i & j index of parent w.r.t child grid--"
         print*,ngrd_rm(mo)%iend_o, ngrd_rm(mo)%jend_o
      enddo



      end subroutine load_roms_grid
!======================================================================

      end module read_roms 



