#include "cppdefs.h"
      MODULE get_wetdry_mod
#ifdef WET_DRY
!
!git $Id$
!svn $Id: get_wetdry.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  If restart, this module reads wetting and drying masks from restart !
!  file using either the standard NetCDF library or the Parallel-IO    !
!  (PIO) library.                                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE nf_fread2d_mod,  ONLY : nf_fread2d
      USE strings_mod,     ONLY : find_string
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: get_wetdry
      PRIVATE :: get_wetdry_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: get_wetdry_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_wetdry (ng, tile, model, IniRec)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, IniRec
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Read in wetting and drying mask according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (INI(ng)%IOtype)
        CASE (io_nf90)
          CALL get_wetdry_nf90 (ng, tile, model, IniRec,                &
     &                          LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_wetdry_pio (ng, tile, model, IniRec,                 &
     &                         LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) INI(ng)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_WETDRY - Illegal input file type, io_type = ',i0,   &
     &        /,14x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_wetdry
!
!***********************************************************************
      SUBROUTINE get_wetdry_nf90 (ng, tile, model, IniRec,              &
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, IniRec
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: gtype, i, status, vindex
      integer :: Vsize(4)
!
      real(dp), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_wetdry_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of restart NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=INI(ng)%name
!
!  Open restart (initial) NetCDF for read/write.
!
      IF (INI(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 1, INI(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idPwet)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      INI(ng)%Vid(idPwet)=Vindex
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idRwet)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      INI(ng)%Vid(idRwet)=Vindex
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idUwet)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      INI(ng)%Vid(idUwet)=Vindex
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idVwet)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      INI(ng)%Vid(idVwet)=Vindex
!
!-----------------------------------------------------------------------
!  Read in wet/dry mask from rstart file.
!-----------------------------------------------------------------------
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
     &                  Vname(1,idPwet), INI(ng)%Vid(idPwet),           &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % pmask,                               &
     &                  GRID(ng) % pmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idPwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at RHO-points.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idRwet), INI(ng)%Vid(idRwet),           &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % rmask,                               &
     &                  GRID(ng) % rmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idRwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at U-points.
!
      gtype=u2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idUwet), INI(ng)%Vid(idUwet),           &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % umask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at V-points.
!
      gtype=v2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idVwet), INI(ng)%Vid(idVwet),           &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % vmask,                               &
     &                  GRID(ng) % vmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
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
  10  FORMAT (/,' GET_WETDRY_NF90 - unable to open grid NetCDF',        &
     &        ' file: ',a)
  20  FORMAT (/,' GET_WETDRY_NF90 - unable to find grid variable: ',a,  &
     &        /,19x,'in NetCDF file: ',a)
  30  FORMAT (/,' GET_WETDRY_NF90 - error while reading variable: ',a,  &
     &        /,19x,'in NetCDF file: ',a)
  40  FORMAT (2x,'GET_WETDRY_NF90    - ',a,/,23x,                       &
     &        '(Grid = ',i2.2,', File: ',a,')',/,23x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
!
      RETURN
      END SUBROUTINE get_wetdry_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_wetdry_pio (ng, tile, model, IniRec,               &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, IniRec
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, status
      integer :: index_pmask, index_rmask, index_umask, index_vmask
      integer :: Vsize(4)
!
      real(dp), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_wetdry_pio"
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
      TYPE (My_VarDesc), pointer :: pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of restart NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=INI(ng)%name
!
!  Open restart (initial) NetCDF for read/write.
!
      IF (INI(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 1, INI(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL pio_netcdf_check_dim (ng, model, ncname,                     &
     &                           pioFile = INI(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     pioFile = INI(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idPwet)),        &
     &                     index_pmask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idRwet)),        &
     &                     index_rmask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idUwet)),        &
     &                     index_umask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idVwet)),        &
     &                     index_vmask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Read in wet/dry mask from rstart file.
!-----------------------------------------------------------------------
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
      pioVar%vd => var_desc(index_pmask)
      IF (KIND(GRID(ng)%pmask_wet).eq.8) THEN
        pioVar%dkind=PIO_double
        ioDesc => ioDesc_dp_p2dvar(ng)
      ELSE
        pioVar%dkind=PIO_real
        ioDesc => ioDesc_sp_p2dvar(ng)
      END IF
      pioVar%gtype=p2dvar
!
      status=nf_fread2d(ng, model, ncname, INI(ng)%pioFile,             &
     &                  Vname(1,idPwet), pioVar,                        &
     &                  IniRec, ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % pmask,                               &
     &                  GRID(ng) % pmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idPwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at RHO-points.
!
      pioVar%vd => var_desc(index_rmask)
      IF (KIND(GRID(ng)%rmask_wet).eq.8) THEN
        pioVar%dkind=PIO_double
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        pioVar%dkind=PIO_real
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
      pioVar%gtype=r2dvar
!
      status=nf_fread2d(ng, model, ncname, INI(ng)%pioFile,             &
     &                  Vname(1,idRwet), pioVar,                        &
     &                  IniRec, ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % rmask,                               &
     &                  GRID(ng) % rmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idRwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at U-points.
!
      pioVar%vd => var_desc(index_umask)
      IF (KIND(GRID(ng)%umask_wet).eq.8) THEN
        pioVar%dkind=PIO_double
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        pioVar%dkind=PIO_real
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
      pioVar%gtype=u2dvar
!
      status=nf_fread2d(ng, model, ncname, INI(ng)%pioFile,             &
     &                  Vname(1,idUwet), pioVar,                        &
     &                  IniRec, ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % umask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at V-points.
!
      pioVar%vd => var_desc(index_vmask)
      IF (KIND(GRID(ng)%vmask_wet).eq.8) THEN
        pioVar%dkind=PIO_double
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        pioVar%dkind=PIO_real
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
      pioVar%gtype=v2dvar
!
      status=nf_fread2d(ng, model, ncname, INI(ng)%pioFile,             &
     &                  Vname(1,idVwet), pioVar,                        &
     &                  IniRec, ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % vmask,                               &
     &                  GRID(ng) % vmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
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
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % pmask_wet,                         &
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet)
!
  10  FORMAT (/,' GET_WETDRY_PIO - unable to open grid NetCDF',         &
     &        ' file: ',a)
  20  FORMAT (/,' GET_WETDRY_PIO - unable to find grid variable: ',a,   &
     &        /,18x,'in NetCDF file: ',a)
  30  FORMAT (/,' GET_WETDRY_PIO - error while reading variable: ',a,   &
     &        /,18x,'in NetCDF file: ',a)
  40  FORMAT (2x,'GET_WETDRY_PIO     - ',a,/,23x,                       &
     &        '(Grid = ',i2.2,', File: ',a,')',/,23x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
!
      RETURN
      END SUBROUTINE get_wetdry_pio
# endif
#endif
      END MODULE get_wetdry_mod

