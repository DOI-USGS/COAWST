!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines reading WW3 grids  
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

      module read_ww3

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use constants
      use create_fullgrid ! adds additional psi points on all four sides of mesh
      use scripwrap_mod ! defines input file variable names for SCRIP-COAWST wrapper

      implicit none

      contains

      subroutine load_ww3_grid()

!     local variables 
      integer(int_kind) :: i, j, n, m
      integer(int_kind) :: mw, nx, ny, grdsize

      real(dbl_kind) :: bf, scale
      real(dbl_kind), allocatable :: zz(:)
      real(dbl_kind), allocatable :: grd(:)

      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max
      real(dbl_kind), allocatable :: lon_psi_w(:,:), lat_psi_w(:,:)

!     Allocate arrays
      allocate(ngrd_w3(Ngrids_ww3))

      do mw=1,Ngrids_ww3
         ngrd_w3(mw)%Numx_ww3=ww3_numx(mw)
         ngrd_w3(mw)%Numy_ww3=ww3_numy(mw)
      end do

      do mw=1,Ngrids_ww3
         nx=ngrd_w3(mw)%Numx_ww3
         ny=ngrd_w3(mw)%Numy_ww3
         allocate(ngrd_w3(mw)%xx(nx,ny))
         allocate(ngrd_w3(mw)%mask_rho_w(nx,ny))
         open(unit=iunit,file=ww3_bath(mw),status='old')
         do j=1,ny
           read(iunit,*) (ngrd_w3(mw)%xx(i,j),i=1,nx)
           do i=1,nx
             if(ngrd_w3(mw)%xx(i,j)==9999)then
               ngrd_w3(mw)%mask_rho_w(i,j)=0
             else
               ngrd_w3(mw)%mask_rho_w(i,j)=1
             endif
            end do
         end do
         close(iunit)
      end do
!----------------------------------------------------------------------
!     Work on the WW3_COORD input files
!----------------------------------------------------------------------
      do mw=1,Ngrids_ww3
        nx=ngrd_w3(mw)%Numx_ww3
        ny=ngrd_w3(mw)%Numy_ww3

        allocate(ngrd_w3(mw)%lon_rho_w(nx,ny))
        allocate(ngrd_w3(mw)%lat_rho_w(nx,ny))

        open(unit=iunit,file=ww3_xcoord(mw),status='old')
        do j=1,ny
          read(iunit,*) (ngrd_w3(mw)%lon_rho_w(i,j),i=1,nx)
        end do
        open(unit=iunit,file=ww3_ycoord(mw),status='old')
        do j=1,ny
          read(iunit,*) (ngrd_w3(mw)%lat_rho_w(i,j),i=1,nx)
        end do
      end do

!   Create a full grid of lat and lon at psi locations 
      do mw=1,Ngrids_ww3
        nx=ngrd_w3(mw)%Numx_ww3
        ny=ngrd_w3(mw)%Numy_ww3

!       First Need to create lat_psi and lon_psi points
!       because SWAN input does not carry these values
!       Making lon_psi and lat_psi as local arrays 
!       Local arrays need to deallocated, no_psi_pts=no_rho_pts-1
        allocate(lon_psi_w(nx-1,ny-1), lat_psi_w(nx-1,ny-1))

        allocate(ngrd_w3(mw)%x_full_grid(nx+1,ny+1))
        allocate(ngrd_w3(mw)%y_full_grid(nx+1,ny+1))

        call create_psimesh(nx, ny, ngrd_w3(mw)%lon_rho_w,              &
     &                              ngrd_w3(mw)%lat_rho_w,              &
     &                              lon_psi_w, lat_psi_w) 
 
        call create_extra_rho_grid(nx-1, ny-1,                          &
     &                             lon_psi_w, lat_psi_w,                &
     &                             ngrd_w3(mw)%x_full_grid,             &
     &                             ngrd_w3(mw)%y_full_grid)

        deallocate(lon_psi_w, lat_psi_w)
      end do

!     Find the indices of child grid with respect to parent grid
!
      do mw=1,Ngrids_ww3-1
        nx=ngrd_w3(mw)%Numx_ww3
        ny=ngrd_w3(mw)%Numy_ww3
        dist_max=10e6
        xx2=ngrd_w3(mw+1)%lon_rho_w(2,2)
        yy2=ngrd_w3(mw+1)%lat_rho_w(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_w3(mw)%lon_rho_w(i,j)-xx2)**2+             &
     &                 (ngrd_w3(mw)%lat_rho_w(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_w3(mw)%istr_w=i
              ngrd_w3(mw)%jstr_w=j
             endif
          enddo
        enddo

        dist_max=10e6
        xxend_1=ngrd_w3(mw+1)%lon_rho_w(ngrd_w3(mw+1)%Numx_ww3-1,       &
     &                                  ngrd_w3(mw+1)%Numy_ww3-1)
        yyend_1=ngrd_w3(mw+1)%lat_rho_w(ngrd_w3(mw+1)%Numx_ww3-1,       &
     &                                  ngrd_w3(mw+1)%Numy_ww3-1)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_w3(mw)%lon_rho_w(i,j)-xxend_1)**2+         &
     &                 (ngrd_w3(mw)%lat_rho_w(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_w3(mw)%iend_w=i
              ngrd_w3(mw)%jend_w=j
            endif
          enddo
        enddo
         print*,"WW3 starting i & j index of parent w.r.t child grid--"
         print*,ngrd_w3(mw)%istr_w, ngrd_w3(mw)%jstr_w
         print*,"ending i & j index of parent w.r.t child grid--"
         print*,ngrd_w3(mw)%iend_w, ngrd_w3(mw)%jend_w
      enddo


      end subroutine load_ww3_grid
!======================================================================

      subroutine create_psimesh(nx, ny, lon_rho,lat_rho,                &
     &                          lon_psi, lat_psi) 

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
      
      ! Save the longitudnal direction first
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
           
!         now middle o
          lon_psi(i,j)=(x1+x2)*half
          lat_psi(i,j)=(y1+y2)*half
        end do 
      end do 
      end subroutine create_psimesh
    
      end module read_ww3

!***********************************************************************

