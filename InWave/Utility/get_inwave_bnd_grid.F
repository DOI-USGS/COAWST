#include "cppdefs.h"
#ifdef INWAVE_MODEL
!
!************************************************************************
      SUBROUTINE get_inwave_bnd_grid(ng, nc_name)
!************************************************************************
!
!svn $Id: get_inwave_grid.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!======================================================================!
!                                                                      !
!  This routine reads the inwave boundary grid and returns bin ........!
!                         directions at the boundary.                  !
!                                                                      !
!======================================================================!
!
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
      USE mod_ncparam
      USE mod_inwave_params
      USE mod_inwave_bound
      USE inwave_iounits
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      character (len=80), intent(in) :: nc_name

!
!  Local variable declarations
!
      integer :: i, ndims, status, varid, ncid, numelem

      integer, dimension(1) :: start, total
      integer, dimension(nf90_max_var_dims) :: dimIDs
      real(r8) :: cnvrad

!
!-----------------------------------------------------------------------
!  Open Boundary NetCDF file for reading.
!-----------------------------------------------------------------------
!

!  Hard code this for now, since swan does not init it.

      ncSPSid=-1
      IF (ncSPSid.eq.-1) THEN
        status=nf90_open(TRIM(nc_name), nf90_nowrite, ncSPSid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,5) TRIM(nc_name)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
      ncid=ncSPSid

!
!-----------------------------------------------------------------------
! Determine number of directions.
!-----------------------------------------------------------------------
!
      status=nf90_inq_varid(ncSPSid,'energy_angle_c', varid)
      status=nf90_inquire_variable(ncSPSid,varid,dimids = dimIDs)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) TRIM('num energy bins at the boundary')
        exit_flag=4
        ioerror=status
      END IF
      status=nf90_inquire_dimension(ncid, dimIDs(1), len=ndims)

      WAVEB(ng) % ND_bnd = ndims

      allocate(WAVEB(ng) % WD_bnd(ndims) )

!
!-----------------------------------------------------------------------
! Get the angles.
!-----------------------------------------------------------------------
!
      start(1)=1
      total(1)=ndims
      status=nf90_get_var(ncid, varid, WAVEB(ng) % WD_BND, start, total)
      
      cnvrad=pi/180.0_r8
      
      DO i=1,ndims
        WAVEB(ng) % WD_BND(i)=cnvrad*WAVEB(ng) % WD_BND(i)
      ENDDO

!
!-----------------------------------------------------------------------
! Close BND NetCDF file.
!-----------------------------------------------------------------------
!
      status=nf90_close(ncSPSid)
      ncSPSid=-1
!

  5   FORMAT (/,' GET_INWAVE_BND - error while opening file: ', a)
 10   FORMAT (/,' GET_INWAVE_BND - error while reading attribute: ', a,&
     &          ' for variable: ', a)
 20   FORMAT (/,' GET_INWAVE_BND - error while inquiring attribute: ', &
     &        a,' for variable: ', a)
 30   FORMAT (/,' GET_INWAVE_BND - cannot inquire ID for variable: ',a)
 40   FORMAT (/,' GET_INWAVE_BND - error while inquiring dimensions',  &
     &          ' for variable: ', a)
 50   FORMAT (/,' GET_INWAVE_BND - error while reading variable: ', a)
      RETURN
      END SUBROUTINE get_inwave_bnd_grid
#else
      SUBROUTINE get_inwave_bnd_grid (ng, tile)
      RETURN
      END SUBROUTINE get_inwave_bnd_grid
#endif
