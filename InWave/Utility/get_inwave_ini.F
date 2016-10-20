#include "cppdefs.h"
#ifdef INWAVE_MODEL

!
!************************************************************************
      SUBROUTINE get_inwave_ini(ng, LBi, UBi, LBj, UBj, ncname, AC_ini, &
     &                          Ta)
!************************************************************************
!
!svn $Id: get_inwave_ini.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!=======================================================================!
!                                                                       !
!          This routine reads the inwave initial condition              !
!                                                                       !
!=======================================================================!
!
      USE mod_param
      USE mod_grid
      USE mod_iounits
      USE mod_netcdf
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE nf_fread3d_mod
      USE inwave_iounits

!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      character (len=80), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj

      real(r8) :: ac_ini(LBi:UBi,LBj:UBj,ND)
      real(r8) :: Ta
      real(r8) :: Amin,Amax
      
      real(r8) , parameter :: Ascl=1.0_r8

      integer :: status
      integer, dimension(4) :: start, total
      integer :: Vsize(4)

      integer :: i, varid, tindex, gtype
      integer :: nxnc,nync,ndnc,nxunc,nyvnc



!# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
!      CALL wclock_on (ng, iNLM, 3)
!
!# endif

!
!-----------------------------------------------------------------------
!  Read in initial condition data from INITIAL INWAVE NetCDF file.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.

      status=nf90_open(TRIM(ncname), nf90_nowrite, ncIWINIid(ng))
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,10)TRIM(ncname)
        exit_flag=2
         ioerror=status
      ENDIF
!
!-----------------------------------------------------------------------
!  Read the absolute period
!-----------------------------------------------------------------------
!
      status=nf90_inq_varid(ncIWINIid(ng),TRIM(Vname(1,idACtp)),varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,20)TRIM(Vname(1,idACtp))
        exit_flag=2
         ioerror=status
      ENDIF

      status=nf90_get_var(ncIWINIid(ng),varid,Ta)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30)TRIM(Vname(1,idACtp))
        exit_flag=2
         ioerror=status
      ENDIF
!
!-----------------------------------------------------------------------
!  Read dimensions
!-----------------------------------------------------------------------
!
      status=nf90_inq_dimid(ncIWINIid(ng),'xi_rho', varid)
      status=nf90_inquire_dimension(ncIWINIid(ng), varid, len=nxnc)
      status=nf90_inq_dimid(ncIWINIid(ng),'xi_u', varid)
      status=nf90_inquire_dimension(ncIWINIid(ng), varid, len=nxunc)
      status=nf90_inq_dimid(ncIWINIid(ng),'eta_rho', varid)
      status=nf90_inquire_dimension(ncIWINIid(ng), varid, len=nync)
      status=nf90_inq_dimid(ncIWINIid(ng),'eta_v', varid)
      status=nf90_inquire_dimension(ncIWINIid(ng), varid, len=nyvnc)
      status=nf90_inq_dimid(ncIWINIid(ng),'energy_angle_c', varid)
      status=nf90_inquire_dimension(ncIWINIid(ng), varid, len=ndnc)

!
!-----------------------------------------------------------------------
!  Read initial wave energy
!-----------------------------------------------------------------------
!
      status=nf90_inq_varid(ncIWINIid(ng),TRIM(Vname(1,idACen)),varid)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,20)TRIM(Vname(1,idACen))
        exit_flag=2
         ioerror=status
      ENDIF
!     gtype=var_flag(varid)*r3dvar
      gtype=r3dvar

      status= nf_fread3d (ng, iNLM, TRIM(ncname), ncIWINIid(ng),        &
     &            TRIM(Vname(1,idACen)),varid, tindex, gtype,           &
     &            Vsize, LBi, UBi, LBj, UBj, 1, ND, Ascl,               &
     &            Amin, Amax,                                           &
#ifdef MASKING
     &            GRID (ng)% rmask,                                     &
#endif
     &            AC_ini)

      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30)TRIM(Vname(1,idACen))
        exit_flag=2
         ioerror=status
      ENDIF
!
!
!-----------------------------------------------------------------------
! Close GRID NetCDF file.
!-----------------------------------------------------------------------
!
      status=nf90_close(ncIWINIid(ng))
      ncIWINIid(ng)=-1

!# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
!      CALL wclock_off (ng, iNLM, 3)
!# endif
  10  FORMAT (/'GET_INWAVE_INI - unable to open input NetCDF file: ',a)
  20  FORMAT (/'GET_INWAVE_INI - unable to find the array: ',a)
  30  FORMAT (/'GET_INWAVE_INI - unable to read the array: ',a)


      RETURN
      END SUBROUTINE get_inwave_ini
#else
      SUBROUTINE get_inwave_ini (ng, LBi, UBi, LBj, UBj, nc_name, Ac_ini,&
     &   Cx_ini, Cy_ini, Ct_ini, Ta)
      RETURN
      END SUBROUTINE get_inwave_ini
#endif

