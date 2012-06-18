#include "cppdefs.h"
      SUBROUTINE get_gridcoords (ngp, ncid, Ilen, Jlen, variable_name, A)
!
!svn $Id: get_varcoords.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads the spatial positions of the NetCDF variables.   !
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
      integer, intent(in) :: ngp, ncid, Ilen, Jlen

      character (len=20), intent(in) :: variable_name

      real(r8), intent(out) :: A(1:Ilen*Jlen)

!
!  Local variable declarations
!
      integer :: i,j, ic
      integer :: ndims, status, varid

      integer, dimension(2) :: start, total

      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Get coarse variable coordinates.
!-----------------------------------------------------------------------
!
!  Open grid NetCDF file for reading.
!
      ncname=GRDname(ngp)
      IF (InpThread.and.(ncGRDid(ngp).eq.-1)) THEN
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncGRDid(ngp))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Read in variable.
!
        status=nf90_inq_varid(ncid, TRIM(variable_name), varid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(variable_name)
          exit_flag=4
          ioerror=status
        END IF
        start(1)=1
        total(1)=Ilen
        start(2)=1
        total(2)=Jlen
        status=nf90_get_var(ncid, varid, A, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(variable_name)
          exit_flag=4
          ioerror=status
        END IF
!
! Close GRID NetCDF file.
!
        status=nf90_close(ncGRDid(ngp))
        ncGRDid(ngp)=-1
      END IF
!

 10   FORMAT (/,' GET_GRIDCOORDS - error while reading attribute: ', a, &
     &          ' for variable: ', a)
 20   FORMAT (/,' GET_GRIDCOORDS - error while inquiring attribute: ',  &
     &        a,' for variable: ', a)
 30   FORMAT (/,' GET_GRIDCOORDS - cannot inquire ID for variable: ', a)
 40   FORMAT (/,' GET_GRIDCOORDS - error while inquiring dimensions',   &
     &          ' for variable: ', a)
 50   FORMAT (/,' GET_GRIDCOORDS - error while reading variable: ', a,  &
     &          ' Check sizes of Lm and Mm in the input file.')

      RETURN
      END SUBROUTINE get_gridcoords
