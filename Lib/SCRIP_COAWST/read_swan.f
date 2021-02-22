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

      subroutine load_swan_grid( MyComm )

#ifdef MPI
      include 'mpif.h'
#endif
      integer (kind=int_kind), intent(in) :: MyComm

!     local variables 
      integer(int_kind) :: i, j, n, m, ng
      integer(int_kind) :: nxs, nys, nxr, nyr
      integer(int_kind) :: test1, test2, test3
      integer(int_kind) :: mw, mo, nx, ny, grdsize
#ifdef MPI
      integer (kind=int_kind) :: MyError, MyRank
#endif

      real(dbl_kind) :: bf, scale
      real(dbl_kind) :: lons, lats, lond, latd, tol, cff
      real(dbl_kind), allocatable :: zz(:)
      real(dbl_kind), allocatable :: grd(:)

      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max
      real(dbl_kind), allocatable :: lon_psi_w(:,:), lat_psi_w(:,:)

      real (kind=dbl_kind), parameter ::                                &
     &      scaleback = 1.0E-10

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

        call create_psimesh(nx, ny, ngrd_sw(mw)%lon_rho_w,              &
     &                              ngrd_sw(mw)%lat_rho_w,              &
     &                               lon_psi_w, lat_psi_w) 
 
        call create_extra_rho_grid(nx-1, ny-1,                          &
     &                             lon_psi_w, lat_psi_w,                &
     &                             ngrd_sw(mw)%x_full_grid,             &
     &                             ngrd_sw(mw)%y_full_grid)

        deallocate(lon_psi_w, lat_psi_w)

      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Here we determine if the ROMS and SWAN grids are coincident.
!
      roms_swan_samegrid=1
      if (Ngrids_roms.eq.Ngrids_swan) then
        do mw=1,Ngrids_swan
          test1=0
          test2=0
          test3=0
!
          nxs=ngrd_sw(mw)%Numx_swan
          nys=ngrd_sw(mw)%Numy_swan
          nxr=ngrd_rm(mw)%xi_size
          nyr=ngrd_rm(mw)%eta_size
          if ((nxs.eq.nxr).and.(nys.eq.nyr)) test1=1
!
!  Compare swan and roms lon points and see if points are coincident. 
!  The test looks for less than 2.0e-7 rads =  1.1459e-05 degs =~ 1.3m
!  If so, then offset the dst by 2.0e-7 rads
          tol=1.1459e-05*2.
          lons=ngrd_rm(mw)%lon_rho_o(1,1)
          lats=ngrd_rm(mw)%lat_rho_o(1,1)
          lond=ngrd_sw(mw)%lon_rho_w(1,1)
          latd=ngrd_sw(mw)%lat_rho_w(1,1)
          cff=min(abs(lons-lond),abs(lats-latd))
          if (abs(cff).le.tol) then
            test2=1
          end if
! now the other corner
          lons=ngrd_rm(mw)%lon_rho_o(nxr,nyr)
          lats=ngrd_rm(mw)%lat_rho_o(nxr,nyr)
          lond=ngrd_sw(mw)%lon_rho_w(nxs,nys)
          latd=ngrd_sw(mw)%lat_rho_w(nxs,nys)
          cff=min(abs(lons-lond),abs(lats-latd))
          if (abs(cff).le.tol) then
            test3=1
          end if
          roms_swan_samegrid=roms_swan_samegrid*test1*test2*test3
        end do
      else
        roms_swan_samegrid=0
      end if
#ifdef MPI
        CALL mpi_comm_rank (MyComm, MyRank, MyError)
        IF (MyRank.eq.0) THEN
#endif
        write(*,*) 'ROMS and SWAN same grids (0=no or 1=yes) ',         &
     &              roms_swan_samegrid
#ifdef MPI
        END IF
        CALL mpi_barrier (MyComm, MyError)
#endif

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Find the indices of child grid with respect to parent grid
!
      do mw=1,Ngrids_swan-1
        if (roms_swan_samegrid.eq.1) then
          mo=mw
          ngrd_sw(mw+1)%istr_w=ngrd_rm(mo+1)%istr_o
          ngrd_sw(mw+1)%jstr_w=ngrd_rm(mo+1)%jstr_o
          ngrd_sw(mw+1)%iend_w=ngrd_rm(mo+1)%iend_o
          ngrd_sw(mw+1)%jend_w=ngrd_rm(mo+1)%jend_o
          ngrd_sw(mw+1)%ref_fac=ngrd_rm(mo+1)%ref_fac
        else
          nx=ngrd_sw(mw)%Numx_swan
          ny=ngrd_sw(mw)%Numy_swan
          dist_max=10e6
          xx2=ngrd_sw(mw+1)%lon_rho_w(2,2)
          yy2=ngrd_sw(mw+1)%lat_rho_w(2,2)

          do j=1,ny
            do i=1,nx
              dist1=sqrt((ngrd_sw(mw)%lon_rho_w(i,j)-xx2)**2+           &
     &                   (ngrd_sw(mw)%lat_rho_w(i,j)-yy2)**2)
              if(dist1<=dist_max)then
                dist_max=dist1
                ngrd_sw(mw+1)%istr_w=i
                ngrd_sw(mw+1)%jstr_w=j
               endif
            enddo
          enddo

          dist_max=10e6
          xxend_1=ngrd_sw(mw+1)%lon_rho_w(ngrd_sw(mw+1)%Numx_swan-1,    &
     &                                    ngrd_sw(mw+1)%Numy_swan-1)
          yyend_1=ngrd_sw(mw+1)%lat_rho_w(ngrd_sw(mw+1)%Numx_swan-1,    &
     &                                    ngrd_sw(mw+1)%Numy_swan-1)

          do j=1,ny
            do i=1,nx
              dist1=sqrt((ngrd_sw(mw)%lon_rho_w(i,j)-xxend_1)**2+       &
     &                   (ngrd_sw(mw)%lat_rho_w(i,j)-yyend_1)**2)
              if(dist1<=dist_max) then
                dist_max=dist1
                ngrd_sw(mw+1)%iend_w=i
                ngrd_sw(mw+1)%jend_w=j
              endif
            enddo
          enddo
        endif

#ifdef MPI
        CALL mpi_comm_rank (MyComm, MyRank, MyError)
        IF (MyRank.eq.0) THEN
#endif
         print*,"SWAN starting i & j index of parent w.r.t child grid--"
         print*,ngrd_sw(mw+1)%istr_w, ngrd_sw(mw+1)%jstr_w
         print*,"ending i & j index of parent w.r.t child grid--"
         print*,ngrd_sw(mw+1)%iend_w, ngrd_sw(mw+1)%jend_w
#ifdef MPI
        END IF
        CALL mpi_barrier (MyComm, MyError)
#endif

      enddo

      end subroutine load_swan_grid
!======================================================================

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

