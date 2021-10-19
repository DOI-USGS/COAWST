#include "cppdefs.h"
#ifdef INWAVE_MODEL

!
!************************************************************************
      SUBROUTINE get_inwave_ini(ng, LBi, UBi, LBj, UBj, IniRec,         &
                                ncname, ac_ini, ta_ini)
!************************************************************************
!
! svn $Id: get_inwave_ini.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!=======================================================================!
!                                                                       !
!          This routine reads the inwave initial condition              !
!                                                                       !
!=======================================================================!
!
      USE inwave_iounits
      USE mod_param
      USE mod_grid
      USE mod_iounits
      USE mod_netcdf
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_parallel
      USE strings_mod,        ONLY : FoundError
      USE nf_fread3d_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, IniRec
      integer, intent(in) :: LBi, UBi, LBj, UBj
      character (len=80), intent(in) :: ncname
      real(r8), intent(inout) :: ac_ini(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: ta_ini(LBi:UBi,LBj:UBj,ND)
!
!  Local variable declarations.
!
      real(r8) :: Ta
      real(r8) :: Amin, Amax
      real(r8) , parameter :: Ascl=1.0_r8

      integer :: ndims, status
      integer :: Vsize(4)
      integer, dimension(nf90_max_var_dims) :: dimIDs

      integer :: i, j, k, varid, gtype
      integer :: model

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
      character (len=6 ) :: string
!
      SourceFile=MyFile
      string=' NLM: '                    ! nonlinear model, restart
      model=iNLM
!
!-----------------------------------------------------------------------
!  Read in initial condition data from INITIAL INWAVE NetCDF file.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.
!
      IF (ncIWINIid(ng).eq.-1) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, ncIWINIid(ng))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,5) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read the absolute period
!-----------------------------------------------------------------------
!
      CALL netcdf_inq_var (ng, model, TRIM(ncname),                     &
     &                     ncid = ncIWINIid(ng))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL netcdf_inq_varid(ng, model, TRIM(ncname),                    &
     &                      TRIM(Vname(1,idACtp)), ncIWINIid(ng), varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      gtype=r3dvar
      status=nf_fread3d(ng, model, TRIM(ncname), ncIWINIid(ng),         &
     &                  TRIM(Vname(1,idACtp)), varid,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj, 1, ND,                      &
     &                  Ascl, Amin, Amax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
# ifdef CHECKSUM
     &                  ta_ini,                                         &
     &                  checksum = Fhash)
# else
     &                  ta_ini)
# endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,60) string, TRIM(Vname(1,idACtp)),              &
     &                      IniRec, TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
# ifdef CHECKSUM
          WRITE (stdout,70) TRIM(Vname(2,idACtp)), Amin, Amax,          &
     &                      Fhash
# else
          WRITE (stdout,70) TRIM(Vname(2,idACtp)), Amin, Amax
# endif
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read initial wave energy
!-----------------------------------------------------------------------
!
      CALL netcdf_inq_varid(ng, model, TRIM(ncname),                    &
     &                      TRIM(Vname(1,idACen)), ncIWINIid(ng), varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      gtype=r3dvar
      status=nf_fread3d(ng, model, TRIM(ncname), ncIWINIid(ng),         &
     &                  TRIM(Vname(1,idACen)), varid,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj, 1, ND,                      &
     &                  Ascl, Amin, Amax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
# ifdef CHECKSUM
     &                  ta_ini,                                         &
     &                  checksum = Fhash)
# else
     &                  AC_ini)
# endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,60) string, TRIM(Vname(1,idACen)),              &
     &                      IniRec, TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
# ifdef CHECKSUM
          WRITE (stdout,70) TRIM(Vname(2,idACen)), Amin, Amax,          &
     &                      Fhash
# else
          WRITE (stdout,70) TRIM(Vname(2,idACen)), Amin, Amax
# endif
        END IF
      END IF
!
!-----------------------------------------------------------------------
! Close GRID NetCDF file.
!-----------------------------------------------------------------------
!
      CALL netcdf_close (ng, model, ncIWINIid(ng), TRIM(ncname),        &
     &                   .FALSE.)

 5    FORMAT (/,' GET_INWAVE_INI - error while opening file: ', a)
 60   FORMAT (/,a,'GET_STATE - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,18x,'in input NetCDF file: ',a)
#ifdef CHECKSUM
 70   FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,' CheckSum = ',i0,')')
#else
 70   FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
#endif

      RETURN
      END SUBROUTINE get_inwave_ini
#else
      SUBROUTINE get_inwave_ini (ng, LBi, UBi, LBj, UBj, nc_name,       &
     &                           Ac_ini, Cx_ini, Cy_ini, Ct_ini, Ta)
      RETURN
      END SUBROUTINE get_inwave_ini
#endif

