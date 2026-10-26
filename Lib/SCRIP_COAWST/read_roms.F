#include "cppdefs.h"
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

      subroutine load_roms_grid( MyComm )

      implicit none

      integer(int_kind), intent(in) :: MyComm
      integer(int_kind) :: i, j, iunit
      integer(int_kind) :: nx, ny, mo, k, N, count, INOUT, ruparent
      integer(int_kind)  :: ncstat, nc_file_id, nc_grdsize_id,          &
     &                      nc_grdlat_id, nc_grdlon_id,                 &
     &                      nc_grdcrnrlat_id, nc_grdcrnrlon_id,         &
     &                      nc_grdmsk_id
#ifdef MPI
      integer(int_kind) :: MyRank, MyError
#endif
      real(dbl_kind) :: xx2, yy2, xxend_1, yyend_1
      real(dbl_kind) :: dist1, dist_max, scale, PX, PY
      real(dbl_kind), allocatable :: lon_psi_o(:,:), lat_psi_o(:,:)
      real(dbl_kind), allocatable :: XX(:), YY(:)

      integer(int_kind) :: mp, nxp, nyp, nxo, nyo
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

#ifdef MPI
      CALL mpi_comm_rank (MyComm, MyRank, MyError)
#endif

      do mo=1,Ngrids_roms

!     Open the file. 
        ncstat=nf_open(roms_grids(mo),nf_nowrite,nc_file_id)
        call netcdf_error_handler(ncstat)

!     Read dimension id 
        ncstat=nf_inq_dimid(nc_file_id,'xi_rho',nc_grdsize_id) 
        call netcdf_error_handler(ncstat)

!     Get the grid size in each direction 
        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,                  &                 
     &                            ngrd_rm(mo)%xi_size)
        call netcdf_error_handler(ncstat)

        ncstat=nf_inq_dimid(nc_file_id,'eta_rho',nc_grdsize_id) 
        call netcdf_error_handler(ncstat)
        ncstat=nf_inq_dimlen(nc_file_id,nc_grdsize_id,                  &
     &                           ngrd_rm(mo)%eta_size)
        call netcdf_error_handler(ncstat)

        if (mo.gt.1) then
!         Get the size of any child grid
          ncstat=nf_get_att_int(nc_file_id,NF_GLOBAL,                   &
     &                           'parent_Imin',ngrd_rm(mo)%istr_o)
          call netcdf_error_handler(ncstat)
          ngrd_rm(mo)%istr_o=ngrd_rm(mo)%istr_o+1   !convert psi to rho
!
          ncstat=nf_get_att_int(nc_file_id,NF_GLOBAL,                   &
     &                           'parent_Imax',ngrd_rm(mo)%iend_o)
          call netcdf_error_handler(ncstat)
!
          ncstat=nf_get_att_int(nc_file_id,NF_GLOBAL,                   &
     &                           'parent_Jmin',ngrd_rm(mo)%jstr_o)
          call netcdf_error_handler(ncstat)
          ngrd_rm(mo)%jstr_o=ngrd_rm(mo)%jstr_o+1   !convert psi to rho
!
          ncstat=nf_get_att_int(nc_file_id,NF_GLOBAL,                   &
     &                           'parent_Jmax',ngrd_rm(mo)%jend_o)
          call netcdf_error_handler(ncstat)
!
          ncstat=nf_get_att_int(nc_file_id,NF_GLOBAL,                   &
     &                           'refine_factor',ngrd_rm(mo)%ref_fac)
          call netcdf_error_handler(ncstat)
#ifdef MPI
          IF (MyRank.eq.0) THEN
#endif
            print*,"ROMS starting i j index of parent w.r.t child grid"
            print*,ngrd_rm(mo)%istr_o, ngrd_rm(mo)%jstr_o
            print*,"ending i & j index of parent w.r.t child grid--"
            print*,ngrd_rm(mo)%iend_o, ngrd_rm(mo)%jend_o
#ifdef MPI
          END IF
          CALL mpi_barrier (MyComm, MyError)
#endif
        endif

!     Read spherical variable (it can be True/False or 1/0)
        ncstat=nf_inq_varid(nc_file_id, 'spherical', checksphere_id)
        call netcdf_error_handler(ncstat)

        ncstat=nf_inq_var(nc_file_id, checksphere_id,  RHNAME, RHTYPE,  &
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

!     Use scale to convert m to degrees, if spherical=T/t/1
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
        allocate(ngrd_rm(mo)%                                           &
     &           lon_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%                                           &
     &           lat_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%                                           &
     &           mask_rho_o(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%                                           &
     &           mask_rho_or(ngrd_rm(mo)%xi_size,ngrd_rm(mo)%eta_size))
        allocate(ngrd_rm(mo)%x_full_grid(ngrd_rm(mo)%xi_size+1,         &
     &                                   ngrd_rm(mo)%eta_size+1))
        allocate(ngrd_rm(mo)%y_full_grid(ngrd_rm(mo)%xi_size+1          &
     &                                  ,ngrd_rm(mo)%eta_size+1))
!   
!       Making lon_psi and lat_psi as local arrays 
!       Local arrays need to deallocated, no_psi_pts=no_rho_pts-1
!
        allocate(lon_psi_o(ngrd_rm(mo)%xi_size-1,                       &
     &                     ngrd_rm(mo)%eta_size-1))
        allocate(lat_psi_o(ngrd_rm(mo)%xi_size-1,                       &
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
        ncstat=nf_get_var_double(nc_file_id, nc_grdlon_id,              &
     &                                     ngrd_rm(mo)%lon_rho_o)
        call netcdf_error_handler(ncstat)
        ngrd_rm(mo)%lon_rho_o=ngrd_rm(mo)%lon_rho_o*scale

        ncstat=nf_get_var_double(nc_file_id, nc_grdlat_id,              &
     &                                     ngrd_rm(mo)%lat_rho_o)
        call netcdf_error_handler(ncstat)
        ngrd_rm(mo)%lat_rho_o=ngrd_rm(mo)%lat_rho_o*scale

!       ncstat=nf_get_var_int(nc_file_id, nc_grdmsk_id,                 &
        ncstat=nf_get_var_double(nc_file_id, nc_grdmsk_id,              &
     &                                     ngrd_rm(mo)%mask_rho_or)
        call netcdf_error_handler(ncstat)

        do i=1,ngrd_rm(mo)%xi_size
          do j=1,ngrd_rm(mo)%eta_size
          ngrd_rm(mo)%mask_rho_o(i,j)=int(ngrd_rm(mo)%mask_rho_or(i,j))
          end do
        end do

        ncstat=nf_get_var_double(nc_file_id, nc_grdcrnrlon_id,          &
     &                                       lon_psi_o)
        call netcdf_error_handler(ncstat)
        lon_psi_o=lon_psi_o*scale

        ncstat=nf_get_var_double(nc_file_id, nc_grdcrnrlat_id,          &
     &                                       lat_psi_o)
        call netcdf_error_handler(ncstat)
        lat_psi_o=lat_psi_o*scale
 
        call create_extra_rho_grid(ngrd_rm(mo)%xi_size-1,               &
     &                             ngrd_rm(mo)%eta_size-1,              &
     &                             lon_psi_o, lat_psi_o,                &
     &                             ngrd_rm(mo)%x_full_grid,             &
     &                             ngrd_rm(mo)%y_full_grid)

        deallocate(lon_psi_o, lat_psi_o)
      end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Find the indices of child grid with respect to parent grid
!     For each child grid, go lower in grid number to find its parent.
!     For the grid i am on (mo), the ngrd_rm(mo) istr iend jstr jend _o are the
!     start and end indices of the parent grid on the parent grid.
      do mo=2,Ngrids_roms
!
!     check lower left and upper right corners
!
        do k=1,mo-1
!
!    Use PNPOLY to check if all points of the grid are inside the potential parent.
!
          N=(ngrd_rm(k)%xi_size+ngrd_rm(k)%eta_size+2)*2-4
          allocate(XX(N), YY(N))
          count=0
          do i=1,ngrd_rm(k)%eta_size+1
            count=count+1
            XX(count)=ngrd_rm(k)%x_full_grid(1,i)
            YY(count)=ngrd_rm(k)%y_full_grid(1,i)
          end do
          do i=2,ngrd_rm(k)%xi_size+1
            count=count+1
            XX(count)=ngrd_rm(k)%x_full_grid(i,ngrd_rm(k)%eta_size+1)
            YY(count)=ngrd_rm(k)%y_full_grid(i,ngrd_rm(k)%eta_size+1)
          end do
          do i=ngrd_rm(k)%eta_size,1,-1
            count=count+1
            XX(count)=ngrd_rm(k)%x_full_grid(ngrd_rm(k)%xi_size+1,i)
            YY(count)=ngrd_rm(k)%y_full_grid(ngrd_rm(k)%xi_size+1,i)
          end do
          do i=ngrd_rm(k)%xi_size,2,-1
            count=count+1
            XX(count)=ngrd_rm(k)%x_full_grid(i,1)
            YY(count)=ngrd_rm(k)%y_full_grid(i,1)
          end do
!
!  For each dst pt call the pnpoly to see if it is inside the src grd.
!
          ruparent=1
          do i=1,ngrd_rm(k)%xi_size
            do j=1,ngrd_rm(k)%eta_size
              PX=ngrd_rm(k)%lon_rho_o(i,j)
              PY=ngrd_rm(k)%lat_rho_o(i,j)
              INOUT=0
              call PNPOLY( PX, PY, XX, YY, N, INOUT )
              ruparent=ruparent*INOUT
            enddo
          enddo
          deallocate(XX, YY)
!
          if (ruparent.eq.1) then
            roms_parent(mo)=k
          end if
        end do
!       now find the parent indices on the parent grid for this child
        mp=roms_parent(mo)
#ifdef MPI
        IF (MyRank.eq.0) THEN
#endif
          print*,"ROMS grid ", mo, " has a parent grid ", mp
#ifdef MPI
        END IF
#endif
        nx=ngrd_rm(mp)%xi_size
        ny=ngrd_rm(mp)%eta_size
        dist_max=10e6
        xx2=ngrd_rm(mo)%lon_rho_o(2,2)
        yy2=ngrd_rm(mo)%lat_rho_o(2,2)
!       lower left
        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_rm(mp)%lon_rho_o(i,j)-xx2)**2+             &
     &                 (ngrd_rm(mp)%lat_rho_o(i,j)-yy2)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_rm(mo)%istr_o=i
              ngrd_rm(mo)%jstr_o=j
             endif
          enddo
        enddo
!       upper right
        dist_max=10e6
        xxend_1=ngrd_rm(mo)%lon_rho_o(ngrd_rm(mo)%xi_size-1,            &
     &                                ngrd_rm(mo)%eta_size-1)
        yyend_1=ngrd_rm(mo)%lat_rho_o(ngrd_rm(mo)%xi_size-1,            &
     &                                ngrd_rm(mo)%eta_size-1)
        do j=1,ny
          do i=1,nx
            dist1=sqrt((ngrd_rm(mp)%lon_rho_o(i,j)-xxend_1)**2+         &
     &                 (ngrd_rm(mp)%lat_rho_o(i,j)-yyend_1)**2)
            if(dist1<=dist_max)then
              dist_max=dist1
              ngrd_rm(mo)%iend_o=i
              ngrd_rm(mo)%jend_o=j
            endif
          enddo
        enddo
#ifdef MPI
        IF (MyRank.eq.0) THEN
#endif
          print*,"ROMS starting i & j index of parent w.r.t child grid-"
          print*,ngrd_rm(mo)%istr_o, ngrd_rm(mo)%jstr_o
          print*,"ending i & j index of parent w.r.t child grid--"
          print*,ngrd_rm(mo)%iend_o, ngrd_rm(mo)%jend_o
#ifdef MPI
        END IF
        CALL mpi_barrier (MyComm, MyError)
#endif
      enddo
      end subroutine load_roms_grid
!
!======================================================================
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
!======================================================================
     
      end module read_roms 



