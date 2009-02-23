#include "cppdefs.h"
      SUBROUTINE get_gridcoords (ng, model, ncid, ncvarid,              &
     &                          Nx, Ny, Xmin, Xmax, X, Ymin, Ymax, Y,   &
     &                          rectangular, Xname, Yname)
!
!svn $Id: get_varcoords.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads the spatial positions of any NetCDF variable     !
!  having the attribute  "coordinates",  as specified by CF rules.     !
!  For example, in CDL syntax:                                         !
!                                                                      !
!       float my_var(time, lat, lon) ;                                 !
!               my_var:long_name = "my variable long name" ;           !
!               my_var:units = "my variable units" ;                   !
!               my_var:coordinates = "lon lat" ;                       !
!               my_var:time = "my_var_time" ;                          !
!                                                                      !
!  Notice that the associated coordinate names "lon" and "lat" are     !
!  separated by a single blank space.  Both "lon" and "lat" can be     !
!  1D or 2D arrays. If 1D array, the positions are rectangular and     !
!  and full 2D arrays are filled with the same values.                 !
!                                                                      !
!  It also determines the rectangular switch  which indicates that     !
!  the spatial positions have a plaid distribution.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
      USE mod_ncparam
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid
      integer, intent(in) :: Nx, Ny

      character (len=20), intent(in) :: Xname, Yname

      logical, intent(out) :: rectangular

      real(r8), intent(out) :: Xmin, Xmax, Ymin, Ymax

      real(r8), intent(out) :: X(Nx,Ny)
      real(r8), intent(out) :: Y(Nx,Ny)

!
!  Local variable declarations
!
      integer :: i, ic, j, jc
      integer :: alen, atype, blank, ndims, status, varid

      integer, dimension(2) :: start, total

      real(r8), dimension(Nx) :: Xwrk
      real(r8), dimension(Ny) :: Ywrk

      character (len=20) :: name   !, Xname, Yname
      character (len=40) :: text
      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Get coarse variable coordinates.
!-----------------------------------------------------------------------
!
!  Get names of variable coordinates.
!
!      IF (InpThread) THEN
!        status=nf90_inquire_attribute(ncid, ncvarid, 'coordinates',     &
!     &                                atype, alen)
!        IF (status.eq.nf90_noerr) THEN
!          status=nf90_get_att(ncid, ncvarid, 'coordinates',             &
!     &                        text(1:alen))
!          IF (status.eq.nf90_noerr) THEN
!            blank=INDEX(text(1:alen),' ')
!            Xname=text(1:blank-1)
!            Yname=text(blank+1:alen)
!          ELSE
!            status=nf90_inquire_variable(ncid, ncvarid, name)
!            WRITE (stdout,10) 'coordinates', TRIM(name)
!            exit_flag=4
!            ioerror=status
!            RETURN
!          END IF       
!        ELSE
!          WRITE (stdout,20) 'coordinates', TRIM(name)
!          exit_flag=4
!          ioerror=status
!          RETURN
!        END IF
!        rectangular=.FALSE.
!
!  Read in X-coordinates.
!

!
!  Open grid NetCDF file for reading.
!
        ncname=GRDname(ng)
        IF (ncGRDid(ng).eq.-1) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, ncGRDid(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
!
        status=nf90_inq_varid(ncid, TRIM(Xname), varid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Xname)
          exit_flag=4
          ioerror=status
        END IF
        status=nf90_inquire_variable(ncid, varid, ndims = ndims)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,40) TRIM(Xname)
          exit_flag=4
          ioerror=status
        END IF
        IF (ndims.eq.1) THEN
          start(1)=1
          total(1)=Nx
          status=nf90_get_var(ncid, varid, Xwrk, start, total)
        ELSE
          start(1)=1
          total(1)=Nx
          start(2)=1
          total(2)=Ny
          status=nf90_get_var(ncid, varid, X, start, total)
        END IF        
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(Xname)
          exit_flag=4
          ioerror=status
        END IF
        IF (ndims.eq.1) THEN
          rectangular=.TRUE.
          jc=0
          DO j=1,Ny
            DO i=1,Nx
              X(i,j)=Xwrk(i)
            END DO
          END DO
        ELSE
          jc=1
          DO j=2,Ny
            IF (X(1,j).eq.X(1,1)) THEN
              jc=jc+1
            END IF
          END DO
        END IF
!
!  Read in Y-coordinates.
!
        status=nf90_inq_varid(ncid, TRIM(Yname), varid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Yname)
          exit_flag=4
          ioerror=status
        END IF
        status=nf90_inquire_variable(ncid, varid, ndims = ndims)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,40) TRIM(Yname)
          exit_flag=4
          ioerror=status
        END IF
        IF (ndims.eq.1) THEN
          start(1)=1
          total(1)=Ny
          status=nf90_get_var(ncid, varid, Ywrk, start, total)
        ELSE
          start(1)=1
          total(1)=Nx
          start(2)=1
          total(2)=Ny
          status=nf90_get_var(ncid, varid, Y, start, total)
        END IF        
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(Yname)
          exit_flag=4
          ioerror=status
        END IF
        IF (ndims.eq.1) THEN
          rectangular=.TRUE.
          ic=0
          DO j=1,Ny
            DO i=1,Nx
              Y(i,j)=Ywrk(j)
            END DO
          END DO
        ELSE
          ic=1
          DO i=2,Nx
            IF (Y(i,1).eq.Y(1,1)) THEN
              ic=ic+1
            END IF
          END DO
        END IF
!
!  Determine "rectangular" switch.
!
        IF (((ic.ne.0).and.(ic.eq.Nx)).and.                             &
     &      ((jc.ne.0).and.(jc.eq.Ny))) THEN
          rectangular=.TRUE.
        END IF
!
!  Determine minimum and maximum positions.
!
        Xmin=1.0E+35_r8
        Xmax=-1.0E+35_r8
        Ymin=1.0E+35_r8
        Ymax=-1.0E+35_r8
        DO j=1,Ny
          DO i=1,Nx
            Xmin=MIN(Xmin,X(i,j))
            Xmax=MAX(Xmax,X(i,j))
            Ymin=MIN(Ymin,Y(i,j))
            Ymax=MAX(Ymax,Y(i,j))
          END DO
        END DO
!      END IF

!
! Close GRID NetCDF file.
!
      status=nf90_close(ncGRDid(ng))
      ncGRDid(ng)=-1
!

 10   FORMAT (/,' GET_GRIDCOORDS - error while reading attribute: ', a, &
     &          ' for variable: ', a)
 20   FORMAT (/,' GET_GRIDCOORDS - error while inquiring attribute: ',  &
     &        a,' for variable: ', a)
 30   FORMAT (/,' GET_GRIDCOORDS - cannot inquire ID for variable: ', a)
 40   FORMAT (/,' GET_GRIDCOORDS - error while inquiring dimensions',   &
     &          ' for variable: ', a)
 50   FORMAT (/,' GET_GRIDCOORDS - error while reading variable: ', a)

      RETURN
      END SUBROUTINE get_gridcoords
