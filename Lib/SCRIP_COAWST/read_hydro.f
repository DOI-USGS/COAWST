!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines reading WRF Hydro grids  
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

      module read_hydro

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use scripwrap_mod ! defines input file variable names 
      use netcdf_mod    !netCDF include file and a netcdf error handling routine
      use create_fullgrid ! adds additional psi points on all four sides of mesh

      contains

      subroutine load_hydro_grid( MyComm )

      implicit none

      integer (kind=int_kind), intent(in) :: MyComm

      ! local variables 
      integer(int_kind) :: i, j, iunit, icount, jcount
      integer(int_kind) :: t, nx, ny, ma, time_size, nx2, ny2
      integer(int_kind) :: ncstat, nc_file_id, nc_grdsize_id,           &
     &         nc_grdlat_id, nc_grdlon_id, nc_grdmsk_id,                &
     &         nc_grdcrnrlat_id, nc_grdcrnrlon_id

      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max
      real(dbl_kind), allocatable :: lon_2drho_h(:,:)
      real(dbl_kind), allocatable :: lat_2drho_h(:,:)
      real(dbl_kind), allocatable :: lon_psi_h(:,:)
      real(dbl_kind), allocatable :: lat_psi_h(:,:)
      real(dbl_kind), allocatable :: lat_temp(:,:)

      allocate(ngrd_hy(Ngrids_hyd))

      do ma=1,Ngrids_hyd
!     Open the file.
          ncstat=nf_open(hydro_grids(ma),nf_nowrite,nc_file_id)
          call netcdf_error_handler(ncstat)

!     Read dimension id, these are the 'rho' point sizes
          ncstat=nf_inq_dimid(nc_file_id,'x',nc_grdsize_id)
          call netcdf_error_handler(ncstat)
!     Get the grid size in each direction
          ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,                &
     &                                    ngrd_hy(ma)%we_size)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_dimid(nc_file_id,'y',nc_grdsize_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,                &
     &                                    ngrd_hy(ma)%sn_size)
          call netcdf_error_handler(ncstat)

!      Read variable for rho points and masking id 
          ncstat=nf_inq_varid(nc_file_id,'LONGITUDE',nc_grdlon_id)
          call netcdf_error_handler(ncstat)
          ncstat=nf_inq_varid(nc_file_id,'LATITUDE',nc_grdlat_id)
          call netcdf_error_handler(ncstat)

          nx2=ngrd_hy(ma)%we_size
          ny2=ngrd_hy(ma)%sn_size

!      Scale the full domain by 10
          ngrd_hy(ma)%we_size=ngrd_hy(ma)%we_size/10
          ngrd_hy(ma)%sn_size=ngrd_hy(ma)%sn_size/10
          nx=ngrd_hy(ma)%we_size
          ny=ngrd_hy(ma)%sn_size

!      Allocate arrays
          allocate( lon_2drho_h(nx2,ny2) )
          allocate( lat_2drho_h(nx2,ny2) )

!     Get rho points and masking values
          ncstat=nf_get_var_double(nc_file_id,nc_grdlon_id,lon_2drho_h)
          call netcdf_error_handler(ncstat)
          ncstat=nf_get_var_double(nc_file_id,nc_grdlat_id,lat_2drho_h)
          call netcdf_error_handler(ncstat)

!      Allocate 2d arrays

          allocate( ngrd_hy(ma)%lon_rho_h(nx,ny) )
          allocate( ngrd_hy(ma)%lat_rho_h(nx,ny) )
          allocate( ngrd_hy(ma)%mask_rho_h(nx,ny) )

!     Here we use every 10th point.
          icount=-5
          jcount=-5
          do j=1,ny
            icount=-5
            jcount=jcount+10
            do i=1,nx
              icount=icount+10
              ngrd_hy(ma)%lon_rho_h(i,j)=lon_2drho_h(icount,jcount)
              ngrd_hy(ma)%lat_rho_h(i,j)=lat_2drho_h(icount,jcount)
!      Set hydro mask to = 1
              ngrd_hy(ma)%mask_rho_h(i,j)=1
            end do
          end do
          deallocate(lon_2drho_h,lat_2drho_h)

! flip up-dpwn the latitude
          allocate(lat_temp(nx,ny))
          jcount=ny
          do j=1,ny
            do i=1,nx
              lat_temp(i,j)=ngrd_hy(ma)%lat_rho_h(i,jcount)
            end do
            jcount=jcount-1
          end do
          do j=1,ny
            do i=1,nx
              ngrd_hy(ma)%lat_rho_h(i,j)=lat_temp(i,j)
            end do
          end do
          deallocate(lat_temp)
!
      end do
!
      do ma=1,Ngrids_hyd
        nx=ngrd_hy(ma)%we_size
        ny=ngrd_hy(ma)%sn_size
!
!       Need to create lat_psi and lon_psi points
!       because WRF input does not carry these values.
!       nx and ny are the number of 'rho' points.
!       First make lon_psi and lat_psi as local arrays
!       with num points =num_rho_pts-1.  This is 2 less than needed
!       but first step is to just interpolate these interior psi points
!       from the rho points.
        allocate(lon_psi_h(nx-1,ny-1), lat_psi_h(nx-1,ny-1))
        allocate(ngrd_hy(ma)%x_full_grid(nx+1,ny+1))
        allocate(ngrd_hy(ma)%y_full_grid(nx+1,ny+1))

        call create_psimesh(nx, ny, ngrd_hy(ma)%lon_rho_h,              &
     &                              ngrd_hy(ma)%lat_rho_h,              &
     &                               lon_psi_h, lat_psi_h)

        call create_extra_rho_grid(nx-1, ny-1,                          &
     &                             lon_psi_h, lat_psi_h,                &
     &                             ngrd_hy(ma)%x_full_grid,             &
     &                             ngrd_hy(ma)%y_full_grid)

        deallocate(lon_psi_h, lat_psi_h)
      end do
!
!     Find the parent grid indices w.r.t child grid
!
      do ma=1,Ngrids_hyd-1
!       allocate(ngrd_wr(ma)%istr_a,ngrd_wr(ma)%jstr_a,                 &
!    &           ngrd_wr(ma)%iend_a,ngrd_wr(ma)%jend_a)

        nx=ngrd_hy(ma)%we_size
        ny=ngrd_hy(ma)%sn_size
        dist_max=10e6
        xx2=ngrd_hy(ma+1)%lon_rho_h(2,2)
        yy2=ngrd_hy(ma+1)%lat_rho_h(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_hy(ma)%lon_rho_h(i,j)-xx2)**2+             &
     &                 (ngrd_hy(ma)%lat_rho_h(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_hy(ma)%istr_h=i
              ngrd_hy(ma)%jstr_h=j
            endif
          enddo
        enddo

        dist_max=10e6
        xxend_1=ngrd_hy(ma+1)%lon_rho_h(ngrd_hy(ma+1)%we_size-1,        &
     &                                  ngrd_hy(ma+1)%sn_size-1)
        yyend_1=ngrd_hy(ma+1)%lat_rho_h(ngrd_hy(ma+1)%we_size-1,        &
     &                                    ngrd_hy(ma+1)%sn_size-1)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_hy(ma)%lon_rho_h(i,j)-xxend_1)**2+         &
     &                 (ngrd_hy(ma)%lat_rho_h(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_hy(ma)%iend_h=i
              ngrd_hy(ma)%jend_h=j
            endif
          enddo
        enddo
         print*, "WRF starting i & j index of parent w.r.t child grid--"
         print*,ngrd_hy(ma)%istr_h, ngrd_hy(ma)%jstr_h
         print*, "ending i & j index of parent w.r.t child grid--"
         print*,ngrd_hy(ma)%iend_h, ngrd_hy(ma)%jend_h
      enddo

      end subroutine load_hydro_grid
!
!======================================================================
!
      subroutine create_psimesh(nx, ny, lon_rho,lat_rho,                &
     &                      lon_psi, lat_psi) 

      implicit none 
!     This subroutine refines the incoming mesh by two times
!     by using bilinear interpolation formula. 
!     Then take every other point of the resulting mesh so that
!     the resulting grid points are at corners of the mesh.  

      integer (int_kind), intent(in) :: nx, ny
      real (dbl_kind), intent(in) :: lon_rho(nx,ny), lat_rho(nx,ny)
      real (dbl_kind), intent(out) :: lon_psi(nx-1,ny-1),               &
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


      end module read_hydro
 

