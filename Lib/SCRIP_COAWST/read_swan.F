#include "cppdefs.h"
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
      integer(int_kind) :: i, j, n, m, ng, mp, ii, nxp, nxw, nyp, nyw
      integer(int_kind) :: nxs, nys, nxr, nyr
      integer(int_kind) :: test1, test2, test3
      integer(int_kind) :: mw, mo, nx, ny, grdsize
      integer(int_kind) :: k, count, INOUT, ruparent
#ifdef MPI
      integer (kind=int_kind) :: MyError, MyRank
#endif

      real(dbl_kind) :: bf, scale, PX, PY
      real(dbl_kind) :: lons, lats, lond, latd, tol, cff
      real(dbl_kind), allocatable :: zz(:)
      real(dbl_kind), allocatable :: grd(:)
      real(dbl_kind), allocatable :: XX(:), YY(:)

      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max
      real(dbl_kind), allocatable :: lon_psi_w(:,:), lat_psi_w(:,:)

      real (kind=dbl_kind), parameter ::                                &
     &      scaleback = 1.0E-10

!     Allocate arrays
      allocate(ngrd_sw(Ngrids_swan))

#ifdef MPI
        CALL mpi_comm_rank (MyComm, MyRank, MyError)
#endif

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
!
!       n is the number of points in the file = numx * numy *2 for lat and lon
!       grdsize is just numx * numy
!
        n=nx*ny*2
        grdsize=n/2

        allocate(grd(n),zz(grdsize))
        do i=1,n
          read(iunit,*) bf
          grd(i)=bf
        end do
!
!       set a scale = 1 for spehrical, 1.rad earth for cartesian
!
        if (swan_cartesian(mw).eq.1) then
          scale=1./6371000.
        else
          scale=1.
        endif
!
!       here we just load lon
!
        do i=1,grdsize
          zz(i)=grd(i)
        end do
!
        do j=1,ny
          do i=1,nx
            ngrd_sw(mw)%lon_rho_w(i,j)=zz(nx*(j-1)+i)*scale
          end do
        end do
!
!       here we just load lat
!
        do i=grdsize+1,n
          zz(i-grdsize)=grd(i)
        end do
!
        do j=1,ny
          do i=1,nx
            ngrd_sw(mw)%lat_rho_w(i,j)=zz(nx*(j-1)+i)*scale
          end do
        end do

        deallocate(grd, zz)
        CLOSE (iunit)

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
!     roms_swan_samegrid=1
      roms_swan_samegrid=0
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
      do mw=2,Ngrids_swan
        do k=1,mw-1
!         Use PNPOLY to check if all points of the grid are inside the potential parent.
          N=(ngrd_sw(k)%Numx_swan+ngrd_sw(k)%Numy_swan+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_sw(k)%Numy_swan+1
            count=count+1
            XX(count)=ngrd_sw(k)%x_full_grid(1,i)
            YY(count)=ngrd_sw(k)%y_full_grid(1,i)
          end do
          do i=2,ngrd_sw(k)%Numx_swan+1
            count=count+1
            XX(count)=ngrd_sw(k)%x_full_grid(i,ngrd_sw(k)%Numy_swan+1)
            YY(count)=ngrd_sw(k)%y_full_grid(i,ngrd_sw(k)%Numy_swan+1)
          end do
          do i=ngrd_sw(k)%Numy_swan,1,-1
            count=count+1
            XX(count)=ngrd_sw(k)%x_full_grid(ngrd_sw(k)%Numx_swan+1,i)
            YY(count)=ngrd_sw(k)%y_full_grid(ngrd_sw(k)%Numx_swan+1,i)
          end do
          do i=ngrd_sw(k)%Numx_swan,2,-1
            count=count+1
            XX(count)=ngrd_sw(k)%x_full_grid(i,1)
            YY(count)=ngrd_sw(k)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          ruparent=1
          do i=1,ngrd_sw(k)%Numx_swan
            do j=1,ngrd_sw(k)%Numy_swan
              PX=ngrd_sw(k)%lon_rho_w(i,j)
              PY=ngrd_sw(k)%lat_rho_w(i,j)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ruparent=ruparent*INOUT
            enddo
          enddo
          deallocate(XX, YY)
!
          if (ruparent.eq.1) then
            swan_parent(mw)=k
          end if
        end do
!           
!       now find the parent indices on the parent grid for this child
        mp=swan_parent(mw)
#ifdef MPI
        IF (MyRank.eq.0) THEN
#endif
          print*,"SWAN grid ", mw, " has a parent grid ", mp
#ifdef MPI
        END IF
#endif
        nx=ngrd_sw(mp)%Numx_swan
        ny=ngrd_sw(mp)%Numy_swan

        dist_max=10e6
        xx2=ngrd_sw(mw)%lon_rho_w(2,2)
        yy2=ngrd_sw(mw)%lat_rho_w(2,2)
        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_sw(mp)%lon_rho_w(i,j)-xx2)**2+             &
     &                 (ngrd_sw(mp)%lat_rho_w(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_sw(mw)%istr_w=i
              ngrd_sw(mw)%jstr_w=j
             endif
          enddo
        enddo

        dist_max=10e6
        xxend_1=ngrd_sw(mw)%lon_rho_w(ngrd_sw(mw)%Numx_swan-1,          &
     &                                ngrd_sw(mw)%Numy_swan-1)
        yyend_1=ngrd_sw(mw)%lat_rho_w(ngrd_sw(mw)%Numx_swan-1,          &
     &                                ngrd_sw(mw)%Numy_swan-1)
        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_sw(mp)%lon_rho_w(i,j)-xxend_1)**2+         &
     &                 (ngrd_sw(mp)%lat_rho_w(i,j)-yyend_1)**2)
            if(dist1<=dist_max) then
              dist_max=dist1
              ngrd_sw(mw)%iend_w=i
              ngrd_sw(mw)%jend_w=j
            endif
          enddo
        enddo
#ifdef MPI
        CALL mpi_comm_rank (MyComm, MyRank, MyError)
        IF (MyRank.eq.0) THEN
#endif
        print*,"SWAN starting i & j index of parent w.r.t child grid--"
        print*,ngrd_sw(mw)%istr_w, ngrd_sw(mw)%jstr_w
        print*,"ending i & j index of parent w.r.t child grid--"
        print*,ngrd_sw(mw)%iend_w, ngrd_sw(mw)%jend_w
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

!                                                                       
         SUBROUTINE PNPOLY( PX, PY, XX, YY, N, INOUT )                     
!
!        https://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
!        Copyright (c) 1970-2003, Wm. Randolph Franklin
!        Permission is hereby granted, free of charge, to any person 
!        obtaining a copy of this software and associated documentation
!        files (the "Software"), to deal in the Software without 
!        restriction, including without limitation the rights to use, 
!        copy, modify, merge, publish, distribute, sublicense, and/or 
!        sell copies of the Software, and to permit persons to whom the
!        Software is furnished to do so, subject to the following 
!        conditions: 
!        1.Redistributions of source code must retain the above 
!          copyright notice, this list of conditions and the following
!          disclaimers. 
!        2.Redistributions in binary form must reproduce the above
!          copyright notice in the documentation and/or other materials
!          provided with the distribution. 
!        3.The name of W. Randolph Franklin may not be used to endorse
!          or promote products derived from this Software without 
!          specific prior written permission. 
!
!        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!        EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!        OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!        NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!        HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!        WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
!        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
!        OTHER DEALINGS IN THE SOFTWARE. 
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!        Modifications:
!          - jcw 12Dec2015 I/O var declarations, 
!                          allow more than 200 maxdim
!     ..................................................................
!                                                                       
!     SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)
!     REAL X(200),Y(200),XX(N),YY(N)
      implicit none
      real(dbl_kind), intent(in) :: PX, PY, XX(1:N), YY(1:N)
      integer, intent(in) :: N
      integer, intent(inout) :: INOUT
!
      integer :: I,J
      real(dbl_kind) :: X(N),Y(N)
      LOGICAL MX,MY,NX,NY
      INTEGER O, MAXDIM
!     OUTPUT UNIT FOR PRINTED MESSAGES
      DATA O/6/
      MAXDIM=N+1
      IF(N.LE.MAXDIM)GO TO 6
      WRITE(O,7)
7     FORMAT('0WARNING:',I5,' TOO GREAT. RESULTS INVALID')
      RETURN
6     DO 1 I=1,N
        X(I)=XX(I)-PX
1       Y(I)=YY(I)-PY
      INOUT=-1
      DO 2 I=1,N
        J=1+MOD(I,N)
        MX=X(I).GE.0.0
        NX=X(J).GE.0.0
        MY=Y(I).GE.0.0
        NY=Y(J).GE.0.0
        IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2
        IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3
        INOUT=-INOUT
        GO TO 2
3       IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4       INOUT=0
        RETURN
5       INOUT=-INOUT
2     CONTINUE
      RETURN
      END SUBROUTINE PNPOLY

      
      end module read_swan

!***********************************************************************

