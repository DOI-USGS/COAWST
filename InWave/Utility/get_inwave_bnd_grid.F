#include "cppdefs.h"
#ifdef INWAVE_MODEL
!
!************************************************************************
      SUBROUTINE get_inwave_bnd_grid(ng, nfiles, S)
!************************************************************************
!
!svn $Id: get_inwave_bnd_grid.F 1336 2008-01-24 02:45:56Z jcwarner $
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
      integer, intent(in) :: ng, nfiles
      TYPE(T_IO), intent(in) :: S(nfiles)
!
!  Local variable declarations
!
      logical :: foundit
      integer :: i, ndims, status, varid, ncid, ifile, NDmax, NDmin
      integer :: nvatt, nvdim, Vid
      integer, dimension(1) :: start, total
      integer, dimension(nf90_max_var_dims) :: dimIDs
      real(r8) :: cnvrad
      character (len=256) :: Fname
!
!-----------------------------------------------------------------------
! Determine number of directions.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      QUERY: DO ifile=1,nfiles
        Fname=S(ifile)%name
!
!  Inquire about requested variable.
!
        CALL netcdf_inq_var (ng, iNLM, Fname,                           &
     &                       MyVarName = 'energy_angle_c',              &
     &                       SearchVar = foundit,                       &
     &                       VarID = Vid,                               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (exit_flag.ne.NoError) RETURN
!
!  Determine if gridded or point data.  Set variable dimensions.
!
        IF (foundit) THEN
          status=nf90_open(TRIM(Fname), nf90_nowrite, ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,5) TRIM(Fname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
          status=nf90_inq_varid(ncid,'energy_angle_c', Vid)
          status=nf90_inquire_variable(ncid, Vid, dimids=dimIDs)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM('num energy bins at the boundary')
            exit_flag=4
            ioerror=status
          END IF

          status=nf90_inquire_dimension(ncid, dimIDs(1), len=ndims)
          WAVEB(ng)%ND_bnd=ndims
          allocate(WAVEB(ng)%WD_bnd(ndims))
!
!-----------------------------------------------------------------------
! Get the angles.
!-----------------------------------------------------------------------
!
          CALL netcdf_get_fvar (ng, iNLM, Fname, 'energy_angle_c',      &
     &                          WAVEB(ng)%WD_bnd,                       &
     &                          ncid)
          CALL netcdf_close (ng, iNLM, ncid, Fname, .FALSE.)
!
          cnvrad=pi/180.0_r8
          DO i=1,ndims
            WAVEB(ng) % WD_bnd(i)=cnvrad*WAVEB(ng) % WD_bnd(i)
          ENDDO
          EXIT QUERY
        END IF
      END DO QUERY
!
  5   FORMAT (/,'GET_INWAVE_BND_GRID - error while opening file: ', a)
 10   FORMAT (/,'GET_INWAVE_BND_GRID - error while reading attribute: ' &
     &         , a, ' for variable: ', a)
 20   FORMAT (/,'GET_INWAVE_BND_GRID - error while inquiring attribute:'&
     &         , a,' for variable: ', a)
 30   FORMAT (/,'GET_INWAVE_BND_GRID - cannot inquire ID for variable: '&
     &         ,a)
 40   FORMAT (/,'GET_INWAVE_BND_GRID - error while inquiring dims',     &
     &          ' for variable: ', a)
 50   FORMAT (/,'GET_INWAVE_BND_GRID - error while reading var: ', a)
      RETURN
      END SUBROUTINE get_inwave_bnd_grid
#else
      SUBROUTINE get_inwave_bnd_grid (ng, tile)
      RETURN
      END SUBROUTINE get_inwave_bnd_grid
#endif
