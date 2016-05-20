!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines reading SWAN grids  
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

      module read_swan

!-----------------------------------------------------------------------
      use kinds_mod     ! defines common data types 
      use constants
      use create_fullgrid ! adds additional psi points on all four sides of mesh
      use scripwrap_mod ! defines input file variable names for SCRIP-COAWST wrapper

      implicit none

      contains

      subroutine load_swan_grid()

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
      allocate(ngrd_sw(Ngrids_swan))

      do mw=1,Ngrids_swan
         ngrd_sw(mw)%Numx_swan=swan_numx(mw)
         ngrd_sw(mw)%Numy_swan=swan_numy(mw)
      end do

      do mw=1,Ngrids_swan
! Notice that the indices are reversed in order for these files to be read 
         nx=ngrd_sw(mw)%Numy_swan
         ny=ngrd_sw(mw)%Numx_swan
         allocate(ngrd_sw(mw)%xx(nx,ny))
         allocate(ngrd_sw(mw)%mask_rho_w(ny,nx))
         open(unit=iunit,file=swan_bath(mw),status='old')
         do i=1,nx
           read(iunit,*) (ngrd_sw(mw)%xx(i,j),j=1,ny)
           do j=1,ny
!            ngrd_sw(mw)%mask_rho_w(j,i)=ngrd_sw(mw)%xx(i,j)
             if(ngrd_sw(mw)%xx(i,j)==9999)then
               ngrd_sw(mw)%mask_rho_w(j,i)=0
             else
               ngrd_sw(mw)%mask_rho_w(j,i)=1
             endif
            end do
         end do
         close(iunit)
      end do
!----------------------------------------------------------------------
!     Work on the SWAN_COORD input file
!----------------------------------------------------------------------
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan

        allocate(ngrd_sw(mw)%lon_rho_w(nx,ny))
        allocate(ngrd_sw(mw)%lat_rho_w(nx,ny))

        open(unit=iunit,file=swan_coord(mw),status='old')
        n=0
        do
          read(iunit,*,end=1)
          n=n+1
        end do
1       rewind(iunit)

!  why do we need to recompute grid size.  it should be Nxswan_nyswan*2
        grdsize=n/2

        allocate(grd(n),zz(grdsize))

        do i=1,n
          read(iunit,*) bf
          grd(i)=bf
        end do
        do i=1,grdsize
          zz(i)=grd(i)
        end do

!  every point will have lon_rho,lat_rho
!  set a scale = 1 for spehrical, 1.rad earth for cartesian
        if (cartesian(mw).eq.1) then
          scale=1./6371000.
        else
          scale=1.
        endif

        do j=1,ny
          do i=1,nx
            ngrd_sw(mw)%lon_rho_w(i,j)=zz(nx*(j-1)+i)*scale
          end do
        end do

        do i=grdsize+1,n
          zz(i-grdsize)=grd(i)
        end do

        do j=1,ny
          do i=1,nx
            ngrd_sw(mw)%lat_rho_w(i,j)=zz(nx*(j-1)+i)*scale
          end do
        end do

        deallocate(grd, zz)

      end do

!   Create a full grid of lat and lon at psi locations 
      do mw=1,Ngrids_swan
        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan

!       First Need to create lat_psi and lon_psi points
!       because SWAN input does not carry these values
!       Making lon_psi and lat_psi as local arrays 
!       Local arrays need to deallocated, no_psi_pts=no_rho_pts-1
        allocate(lon_psi_w(nx-1,ny-1), lat_psi_w(nx-1,ny-1))

        allocate(ngrd_sw(mw)%x_full_grid(nx+1,ny+1))
        allocate(ngrd_sw(mw)%y_full_grid(nx+1,ny+1))

        call create_psimesh(nx, ny, ngrd_sw(mw)%lon_rho_w,
     &                              ngrd_sw(mw)%lat_rho_w,
     &                               lon_psi_w, lat_psi_w) 
 
        call create_extra_rho_grid(nx-1, ny-1, 
     &                             lon_psi_w, lat_psi_w,
     &                             ngrd_sw(mw)%x_full_grid,  
     &                             ngrd_sw(mw)%y_full_grid)

        deallocate(lon_psi_w, lat_psi_w)
      end do

!     Find the indices of child grid with respect to parent grid
!
      do mw=1,Ngrids_swan-1
        allocate(ngrd_sw(mw)%istr_w,ngrd_sw(mw)%jstr_w,
     &           ngrd_sw(mw)%iend_w,ngrd_sw(mw)%jend_w)

        nx=ngrd_sw(mw)%Numx_swan
        ny=ngrd_sw(mw)%Numy_swan
        dist_max=10e6
        xx2=ngrd_sw(mw+1)%lon_rho_w(2,2)
        yy2=ngrd_sw(mw+1)%lat_rho_w(2,2)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_sw(mw)%lon_rho_w(i,j)-xx2)**2+
     &                 (ngrd_sw(mw)%lat_rho_w(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_sw(mw)%istr_w=i
              ngrd_sw(mw)%jstr_w=j
            endif
          enddo
        enddo

        dist_max=10e6
        xxend_1=ngrd_sw(mw+1)%lon_rho_w(ngrd_sw(mw+1)%Numx_swan-1,
     &                                  ngrd_sw(mw+1)%Numy_swan-1)
        yyend_1=ngrd_sw(mw+1)%lat_rho_w(ngrd_sw(mw+1)%Numx_swan-1,
     &                                  ngrd_sw(mw+1)%Numy_swan-1)

        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_sw(mw)%lon_rho_w(i,j)-xxend_1)**2+
     &                 (ngrd_sw(mw)%lat_rho_w(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_sw(mw)%iend_w=i
              ngrd_sw(mw)%jend_w=j
            endif
          enddo
        enddo
         print*,"SWAN starting i & j index of parent w.r.t child grid--"
         print*,ngrd_sw(mw)%istr_w, ngrd_sw(mw)%jstr_w
         print*,"ending i & j index of parent w.r.t child grid--"
         print*,ngrd_sw(mw)%iend_w, ngrd_sw(mw)%jend_w
      enddo


      end subroutine load_swan_grid
!======================================================================

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
    
      end module read_swan

!***********************************************************************

