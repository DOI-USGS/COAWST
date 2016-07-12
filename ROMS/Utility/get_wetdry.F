#include "cppdefs.h"
#ifdef WET_DRY
      SUBROUTINE get_wetdry (ng, model, IniRec)
!
!svn $Id: get_wetdry.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine reads wetting and drying masks from restart NetCDF  !
!  file.
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, IniRec
!
!  Local variable declarations.
!
      integer :: tile, LBi, UBi, LBj, UBj
      integer :: gtype, i, status
      integer :: varid_pmask, varid_rmask, varid_umask, varid_vmask
      integer :: Vsize(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin

      character (len=256) :: ncname
!
      SourceFile='get_wetdry.F'
!
!-----------------------------------------------------------------------
!  Inquire about the contents of restart NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=INI(ng)%name
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idPwet)),        &
     &                     varid_pmask)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idRwet)),        &
     &                     varid_rmask)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idUwet)),        &
     &                     varid_umask)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idVwet)),        &
     &                     varid_vmask)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Open grid NetCDF file for reading.
!
      IF (INI(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, INI(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,20) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read in wet/dry mask from rstart file.
!-----------------------------------------------------------------------
!
!  Set 2D arrays bounds.
!
# ifdef DISTRIBUTE
      tile=MyRank
# else
      tile=-1
# endif
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in wet/dry mask at PSI-points.
!
      gtype=p2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idPwet), varid_pmask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % pmask,                               &
     &                  GRID(ng) % pmask_wet)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in wet/dry mask at RHO-points.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idRwet), varid_rmask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % rmask,                               &
     &                  GRID(ng) % rmask_wet)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in wet/dry mask at U-points.
!
      gtype=u2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idUwet), varid_umask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % umask_wet)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in wet/dry mask at V-points.
!
      gtype=v2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idVwet), varid_vmask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % vmask,                               &
     &                  GRID(ng) % vmask_wet)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % pmask_wet)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % rmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % umask_wet)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % vmask_wet)
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % pmask_wet,                         &
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet)
# endif
!
  10  FORMAT (/,' GET_WETDRY - unable to find grid variable: ',a,       &
     &        /,14x,'in grid NetCDF file: ',a)
  20  FORMAT (/,' GET_WETDRY - unable to open grid NetCDF file: ',a)
  30  FORMAT (/,' GET_WETDRY - error while reading variable: ',a,       &
     &        /,12x,'in grid NetCDF file: ',a)
#else
      SUBROUTINE get_wetdry
#endif
      RETURN
      END SUBROUTINE get_wetdry
