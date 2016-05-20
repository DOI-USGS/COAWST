!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines reading WRF grids  
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

      module read_wrf

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use scripwrap_mod ! defines input file variable names 
      use netcdf_mod    !netCDF include file and a netcdf error handling routine
      use create_fullgrid ! adds additional psi points on all four sides of mesh

      contains

      subroutine load_wrf_grid()

      implicit none

      ! local variables 
      integer(int_kind) :: i, j, iunit
      integer(int_kind) :: t, nx, ny, ma, time_size
      integer(int_kind) :: ncstat, nc_file_id, nc_grdsize_id,
     &         nc_grdlat_id, nc_grdlon_id, nc_grdmsk_id,
     &         nc_grdcrnrlat_id, nc_grdcrnrlon_id

      real(dbl_kind)    :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind)    :: dist1, dist_max
      real(dbl_kind), allocatable :: lon_3drho_a(:,:,:)  
      real(dbl_kind), allocatable :: lat_3drho_a(:,:,:) 
      real(dbl_kind), allocatable :: mask_3drho_a(:,:,:) 
      real(dbl_kind), allocatable :: lon_psi_a(:,:), lat_psi_a(:,:)


      allocate(ngrd_wr(Ngrids_wrf))

      do ma=1,Ngrids_wrf
!     Open the file. 
        if (wrf_grids(ma)=="moving") then
          call create_wrf_moving_grid(ma)
        else
          ncstat=nf_open(wrf_grids(ma),nf_nowrite,nc_file_id)
          call netcdf_error_handler(ncstat)

!     Read dimension id, these are the 'rho' point sizes
          ncstat=nf_inq_dimid(nc_file_id,'west_east',nc_grdsize_id)
          call netcdf_error_handler(ncstat)
!     Get the grid size in each direction
          ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,
     &                                    ngrd_wr(ma)%we_size)

          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_dimid(nc_file_id,'south_north',nc_grdsize_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,
     &                                    ngrd_wr(ma)%sn_size)
          call netcdf_error_handler(ncstat)
          ncstat = nf_inq_dimid(nc_file_id,'Time',nc_grdsize_id)
          call netcdf_error_handler(ncstat)
          ncstat = nf_inq_dimlen(nc_file_id,nc_grdsize_id,
     &                                      time_size)
        call netcdf_error_handler(ncstat)

!      Read variable for rho points and masking id 
          ncstat=nf_inq_varid(nc_file_id,'XLONG',nc_grdlon_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'XLAT',nc_grdlat_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'LANDMASK',nc_grdmsk_id)
          call netcdf_error_handler(ncstat)

          nx=ngrd_wr(ma)%we_size
          ny=ngrd_wr(ma)%sn_size
!      Allocate arrays
          allocate( lon_3drho_a(time_size,nx,ny) )
          allocate( lat_3drho_a(time_size,nx,ny) )
          allocate( mask_3drho_a(time_size,nx,ny) )

!     Get rho points and masking values
          ncstat=nf_get_var_double(nc_file_id,nc_grdlon_id,lon_3drho_a)
          call netcdf_error_handler(ncstat)
          ncstat=nf_get_var_double(nc_file_id,nc_grdlat_id,lat_3drho_a)
          call netcdf_error_handler(ncstat)
           ncstat=nf_get_var_int(nc_file_id,nc_grdmsk_id,mask_3drho_a)
          call netcdf_error_handler(ncstat)

!      Allocate 2d arrays
          allocate( ngrd_wr(ma)%lon_rho_a(nx,ny) )
          allocate( ngrd_wr(ma)%lat_rho_a(nx,ny) )
          allocate( ngrd_wr(ma)%mask_rho_a(nx,ny) )

!     For WRF the variables in NETCDF file are also required from
!     3D array to 2d array for creating masks and SCRIP package
!     Select to use t=1
          do j=1,ny
            do i=1,nx
              ngrd_wr(ma)%lon_rho_a(i,j)=lon_3drho_a(1,i,j)
              ngrd_wr(ma)%lat_rho_a(i,j)=lat_3drho_a(1,i,j)
!      Convert atm masking to be oceanographic (convert to sea=1, land=0)
              ngrd_wr(ma)%mask_rho_a(i,j)=1-mask_3drho_a(1,i,j)
            end do
          end do

          deallocate(lon_3drho_a,lat_3drho_a,mask_3drho_a)
        endif
      end do 
!  come up to here for the moving nest full parent grid.

      do ma=1,Ngrids_wrf
        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
!
!       Need to create lat_psi and lon_psi points
!       because WRF input does not carry these values.
!       nx and ny are the number of 'rho' points.
!       First make lon_psi and lat_psi as local arrays
!       with num points =num_rho_pts-1.  This is 2 less than needed
!       but first step is to just interpolate these interior psi points
!       from the rho points.
        allocate(lon_psi_a(nx-1,ny-1), lat_psi_a(nx-1,ny-1))
        allocate(ngrd_wr(ma)%x_full_grid(nx+1,ny+1))
        allocate(ngrd_wr(ma)%y_full_grid(nx+1,ny+1))

        call create_psimesh(nx, ny, ngrd_wr(ma)%lon_rho_a,
     &                              ngrd_wr(ma)%lat_rho_a,
     &                               lon_psi_a, lat_psi_a)

        call create_extra_rho_grid(nx-1, ny-1,
     &                             lon_psi_a, lat_psi_a,
     &                             ngrd_wr(ma)%x_full_grid,
     &                             ngrd_wr(ma)%y_full_grid)

        deallocate(lon_psi_a, lat_psi_a)
      end do

!     Find the parent grid indices w.r.t child grid
!
      do ma=1,Ngrids_wrf-1
        allocate(ngrd_wr(ma)%istr_a,ngrd_wr(ma)%jstr_a,
     &           ngrd_wr(ma)%iend_a,ngrd_wr(ma)%jend_a)

        nx=ngrd_wr(ma)%we_size
        ny=ngrd_wr(ma)%sn_size
        dist_max=10e6
        xx2=ngrd_wr(ma+1)%lon_rho_a(2,2)
        yy2=ngrd_wr(ma+1)%lat_rho_a(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_wr(ma)%lon_rho_a(i,j)-xx2)**2+
     &                 (ngrd_wr(ma)%lat_rho_a(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_wr(ma)%istr_a=i
              ngrd_wr(ma)%jstr_a=j
            endif
          enddo
        enddo

        dist_max=10e6
        xxend_1=ngrd_wr(ma+1)%lon_rho_a(ngrd_wr(ma+1)%we_size-1,
     &                                  ngrd_wr(ma+1)%sn_size-1)
        yyend_1=ngrd_wr(ma+1)%lat_rho_a(ngrd_wr(ma+1)%we_size-1,
     &                                    ngrd_wr(ma+1)%sn_size-1)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_wr(ma)%lon_rho_a(i,j)-xxend_1)**2+      
     &                 (ngrd_wr(ma)%lat_rho_a(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_wr(ma)%iend_a=i
              ngrd_wr(ma)%jend_a=j
            endif
          enddo
        enddo
         print*, "WRF starting i & j index of parent w.r.t child grid--"
         print*,ngrd_wr(ma)%istr_a, ngrd_wr(ma)%jstr_a
         print*, "ending i & j index of parent w.r.t child grid--"
         print*,ngrd_wr(ma)%iend_a, ngrd_wr(ma)%jend_a
      enddo

      end subroutine load_wrf_grid
!
!======================================================================
!
      subroutine create_psimesh(nx, ny, lon_rho,lat_rho,
     &                      lon_psi, lat_psi) 

      implicit none 
!     This subroutine refines the incoming mesh by two times
!     by using bilinear interpolation formula. 
!     Then take every other point of the resulting mesh so that
!     the resulting grid points are at corners of the mesh.  

      integer (int_kind), intent(in) :: nx, ny
      real (dbl_kind), intent(in) :: lon_rho(nx,ny), lat_rho(nx,ny)
      real (dbl_kind), intent(out) :: lon_psi(nx-1,ny-1), 
     &                                lat_psi(nx-1,ny-1)

      integer (int_kind) :: i, j
      real (dbl_kind) :: x1, y1, x2, y2
    
!       Save the longitudnal direction first
!       x1,y1 are mid points of (ln1,la1 and ln1,la2)
!       x2,y2 are mid points of (ln2,la1 and ln2,la2)
!      
!ln1,la2 __________ ln2,la2
!       |          |  
!       |          |
!x1,y1  |    *     |x2,y2
!       |          | 
!       |__________|
!ln1,la1           ln2,la1
!  *=mid point of x1,y1 and x2,y2

      do j=1,ny-1
        do i=1,nx-1
          x1=(lon_rho(i,j)+lon_rho(i,j+1))*half
          y1=(lat_rho(i,j)+lat_rho(i,j+1))*half
          x2=(lon_rho(i+1,j)+lon_rho(i+1,j+1))*half
          y2=(lat_rho(i+1,j)+lat_rho(i+1,j+1))*half

          lon_psi(i,j)=(x1+x2)*half
          lat_psi(i,j)=(y1+y2)*half
        end do
      end do
!
      end subroutine create_psimesh

!======================================================================
!
      subroutine create_wrf_moving_grid(ma)

      implicit none

      integer (int_kind), intent(in) :: ma

      integer (int_kind) :: i, j, nx, ny, mm, nn, pgr
      integer (int_kind) :: icount, jcount, pid
      real (dbl_kind) :: x1, y1, x2, y2, dx, dy
      real (dbl_kind), allocatable :: lon_rho_t(:,:), lat_rho_t(:,:)

!  Create a grid that is the size of the parent x refined ratio.
      pgr=parent_grid_ratio(ma)
      pid=parent_id(ma)
      nx=ngrd_wr(pid)%we_size*pgr
      ny=ngrd_wr(pid)%sn_size*pgr
      ngrd_wr(ma)%we_size=nx
      ngrd_wr(ma)%sn_size=ny
!     t is temporary arrays
      allocate( lon_rho_t(nx-2,ny-2) )
      allocate( lat_rho_t(nx-2,ny-2) )
      allocate( ngrd_wr(ma)%lon_rho_a(nx,ny) )
      allocate( ngrd_wr(ma)%lat_rho_a(nx,ny) )
      allocate( ngrd_wr(ma)%mask_rho_a(nx,ny) )

! need to revisit masking
      ngrd_wr(ma)%mask_rho_a=1

! this gets every pgr th point
      icount=0
      do i=2,nx-1,pgr
        icount=icount+1
        jcount=0
        do j=2,ny-1,pgr
          jcount=jcount+1
          ngrd_wr(ma)%lon_rho_a(i,j)=                                   &
     &                             ngrd_wr(pid)%lon_rho_a(icount,jcount)
          ngrd_wr(ma)%lat_rho_a(i,j)=                                   &
     &                             ngrd_wr(pid)%lat_rho_a(icount,jcount)
        end do
      end do

! now get the refined points.
      icount=0
      do i=2,nx-1-pgr,pgr
        icount=icount+1
        jcount=0
        do j=2,ny-1-pgr,pgr
          jcount=jcount+1

!  bottom
          dx=ngrd_wr(pid)%lon_rho_a(icount+1,jcount)-                   &
     &       ngrd_wr(pid)%lon_rho_a(icount,jcount)
          dy=ngrd_wr(pid)%lat_rho_a(icount+1,jcount)-                   &
     &       ngrd_wr(pid)%lat_rho_a(icount,jcount)
          do mm=1,pgr-1
            ngrd_wr(ma)%lon_rho_a(i+mm,j)=                              &
     &                        ngrd_wr(pid)%lon_rho_a(icount,jcount)+    &
     &                        dx/pgr*mm
            ngrd_wr(ma)%lat_rho_a(i+mm,j)=                              &
     &                        ngrd_wr(pid)%lat_rho_a(icount,jcount)+    &
     &                        dy/pgr*mm
          enddo

!  top
          dx=ngrd_wr(pid)%lon_rho_a(icount+1,jcount+1)-                 &
     &       ngrd_wr(pid)%lon_rho_a(icount,jcount+1)
          dy=ngrd_wr(pid)%lat_rho_a(icount+1,jcount+1)-                 &
     &       ngrd_wr(pid)%lat_rho_a(icount,jcount+1)
          do mm=1,pgr-1
            ngrd_wr(ma)%lon_rho_a(i+mm,j+pgr)=                          &
     &               ngrd_wr(pid)%lon_rho_a(icount,jcount+1)+dx/pgr*mm
            ngrd_wr(ma)%lat_rho_a(i+mm,j+pgr)=                          &
     &               ngrd_wr(pid)%lat_rho_a(icount,jcount+1)+dy/pgr*mm
          enddo

!  left
          dx=ngrd_wr(pid)%lon_rho_a(icount,jcount+1)-                   &
     &       ngrd_wr(pid)%lon_rho_a(icount,jcount)
          dy=ngrd_wr(pid)%lat_rho_a(icount,jcount+1)-                   &
     &       ngrd_wr(pid)%lat_rho_a(icount,jcount)
          do mm=1,pgr-1
            ngrd_wr(ma)%lon_rho_a(i,j+mm)=                              &
     &                 ngrd_wr(pid)%lon_rho_a(icount,jcount)+dx/pgr*mm
            ngrd_wr(ma)%lat_rho_a(i,j+mm)=                              &
     &                 ngrd_wr(pid)%lat_rho_a(icount,jcount)+dy/pgr*mm
          enddo

!  right
          dx=ngrd_wr(pid)%lon_rho_a(icount+1,jcount+1)-                 &
     &       ngrd_wr(pid)%lon_rho_a(icount+1,jcount)
          dy=ngrd_wr(pid)%lat_rho_a(icount+1,jcount+1)-                 &
     &       ngrd_wr(pid)%lat_rho_a(icount+1,jcount)
          do mm=1,pgr-1
            ngrd_wr(ma)%lon_rho_a(i+pgr,j+mm)=                          &
     &               ngrd_wr(pid)%lon_rho_a(icount+1,jcount)+dx/pgr*mm
            ngrd_wr(ma)%lat_rho_a(i+pgr,j+mm)=                          &
     &               ngrd_wr(pid)%lat_rho_a(icount+1,jcount)+dy/pgr*mm
          enddo

!  insides
          do nn=1,pgr-1                ! y-dir
            dx=ngrd_wr(ma)%lon_rho_a(i+pgr,j+nn)-                       &
     &         ngrd_wr(ma)%lon_rho_a(i    ,j+nn)
            do mm=1,pgr-1              ! x-dir
              ngrd_wr(ma)%lon_rho_a(i+mm,j+nn)=                         &
     &                         ngrd_wr(ma)%lon_rho_a(i,j+nn)+dx/pgr*mm
            enddo
          enddo
          do mm=1,pgr-1                ! x-dir
            dy=ngrd_wr(ma)%lat_rho_a(i+mm,j+pgr)-                       &
     &                         ngrd_wr(ma)%lat_rho_a(i+mm,j)
            do nn=1,pgr-1              ! y-dir
              ngrd_wr(ma)%lat_rho_a(i+mm,j+nn)=                         &
     &                         ngrd_wr(ma)%lat_rho_a(i+mm,j)+dy/pgr*nn
            enddo
          enddo
        end do
      end do
!  place data in temp array
      do i=1,nx-2
        do j=1,ny-2
          lon_rho_t(i,j)=ngrd_wr(ma)%lon_rho_a(i+1,j+1)
          lat_rho_t(i,j)=ngrd_wr(ma)%lat_rho_a(i+1,j+1)
        end do
      end do

! Set perimeter
      call create_extra_rho_grid(nx-2, ny-2,                            &
     &                           lon_rho_t,                             & 
     &                           lat_rho_t,                             &
     &                           ngrd_wr(ma)%lon_rho_a,                 &
     &                           ngrd_wr(ma)%lat_rho_a)

      deallocate( lon_rho_t, lat_rho_t )

      end subroutine create_wrf_moving_grid

!***********************************************************************
!

      end module read_wrf
 

