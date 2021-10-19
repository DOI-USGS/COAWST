#include "cppdefs.h"
#ifdef INWAVE_MODEL
!
!************************************************************************
      SUBROUTINE get_inwave_grid(ng, nc_name, WD)
!************************************************************************
!
!svn $Id: get_inwave_grid.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!=======================================================================
!                                                                      !
!  This routine reads the inwave grid and returns bin directions.      !
!                                                                      !
!=======================================================================
!
      USE inwave_iounits
      USE mod_iounits
      USE mod_netcdf
      USE mod_parallel
      USE mod_scalars
      USE mod_ncparam
      USE mod_inwave_params
      USE strings_mod,        ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      character (len=80), intent(in) :: nc_name
      real(r8), intent(inout) :: WD(ND)
!
!  Local variable declarations
!
      integer ::  ndims, status, varid, ncid, numelem
      integer ::  model, my_Nangles
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Open grid NetCDF file for reading.
!-----------------------------------------------------------------------
!
      model=iNLM
      IF (ncSPSid.eq.-1) THEN
        CALL netcdf_open (ng, model, TRIM(nc_name), 0, ncSPSid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,5) TRIM(nc_name)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
! Determine number of directions.
!-----------------------------------------------------------------------
!
      CALL netcdf_inq_var (ng, model, TRIM(nc_name),                    &
     &                     ncid = ncSPSid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Check if the number of nested grids is correct.
!
      CALL netcdf_get_dim (ng, model, TRIM(nc_name),                    &
     &                     ncid = ncSPSid,                              &
     &                     DimName = 'ND',                              &
     &                     DimSize = my_Nangles)
! catch to get old way
      IF (exit_flag.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,60) 'cant find ND, search energy_angle_centers'
 60               FORMAT (/,' get_inwave_grid - ', a)
        END IF
        CALL netcdf_get_dim (ng, model, TRIM(nc_name),                  &
     &                       ncid = ncSPSid,                            &
     &                       DimName = 'energy_angle_centers',          &
     &                       DimSize = my_Nangles)
        IF (exit_flag.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,61) 'found energy_angle_centers!!  ALL GOOD.'
 61         FORMAT (/,' get_inwave_grid - ', a)
          END IF
        END IF
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (my_Nangles.ne.ND) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'inconsistent num energy bins ND = ',       &
     &                      ND, my_Nangles
 10             FORMAT (/,' get_inwave_grid - ', a, i4, 2x, i4,         &
     &            /,16x,'in input file:'2x,a)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!-----------------------------------------------------------------------
! Get the angles.
!-----------------------------------------------------------------------
!
      CALL netcdf_get_fvar (ng, model, NGCname, 'energy_angle',         &
     &                      WD,                                         &
     &                      ncid = ncSPSid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
! Close GRID NetCDF file.
!-----------------------------------------------------------------------
!
      CALL netcdf_close (ng, model, ncSPSid, TRIM(nc_name), .FALSE.)
!

  5   FORMAT (/,' GET_INWAVE_GRID - error while opening file: ', a)
 20   FORMAT (/,' GET_INWAVE_GRID - error while inquiring attribute: ', &
     &        a,' for variable: ', a)
 30   FORMAT (/,' GET_INWAVE_GRID - cannot inquire ID for variable: ',a)
 35   FORMAT (/,' GET_INWAVE_GRID - wrong number of bins: ',a)
 40   FORMAT (/,' GET_INWAVE_GRID - error while inquiring dimensions',  &
     &          ' for variable: ', a)
 50   FORMAT (/,' GET_INWAVE_GRID - error while reading variable: ', a)
      RETURN
      END SUBROUTINE get_inwave_grid
#else
      SUBROUTINE get_inwave_grid (ng, tile)
      RETURN
      END SUBROUTINE get_inwave_grid
#endif
