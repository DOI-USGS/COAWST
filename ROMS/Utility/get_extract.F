#include "cppdefs.h"
      MODULE get_extract_mod

#ifdef GRID_EXTRACT
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module reads extraction grid information from input file using !
!  either standard NetCDF library or the Parallel-IO (PIO) library.    !
!                                                                      !
!  This grid geometry will be used in field extraction by decimation   !
!  or interpolation, which are written to output NetCDF file.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_extract
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
!
      USE exchange_2d_xtr_mod
      USE extract_field_mod,   ONLY : interp_coords
# ifdef DISTRIBUTE
      USE mp_exchange_mod,     ONLY : mp_exchange2d_xtr
# endif
      USE nf_fread2d_xtr_mod,  ONLY : nf_fread2d_xtr
      USE strings_mod,         ONLY : FoundError, find_string
!
      implicit none
!
      PUBLIC  :: get_extract
      PRIVATE :: get_extract_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: get_extract_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_extract (ng, tile, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Read in GRID NetCDF file according to IO type.
!-----------------------------------------------------------------------
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (GRX(ng)%IOtype)
        CASE (io_nf90)
          CALL get_extract_nf90 (ng, tile, model,                      &
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_extract_pio (ng, tile, model,                       &
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) GRX(ng)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Set packed, global geometry arrays needed for field extraction by
!  decimation or interpolation.
!-----------------------------------------------------------------------
!
      IF (ExtractFlag(ng).ge.1) THEN
        CALL interp_coords (ng, tile, model, p2dvar,                    &
     &                      GRID(ng)%Gangle_psi,                        &
# ifdef MASKING
     &                      GRID(ng)%Gmask_psi,                         &
# endif
     &                      GRID(ng)%Gx_psi,                            &
     &                      GRID(ng)%Gy_psi,                            &
# ifdef MASKING
     &                      EXTRACT(ng)%Gmask_psi,                      &
# endif
     &                      EXTRACT(ng)%Gx_psi,                         &
     &                      EXTRACT(ng)%Gy_psi,                         &
     &                      EXTRACT(ng)%Iout_psi,                       &
     &                      EXTRACT(ng)%Jout_psi)
!
        CALL interp_coords (ng, tile, model, r2dvar,                    &
     &                      GRID(ng)%Gangle_rho,                        &
# ifdef MASKING
     &                      GRID(ng)%Gmask_rho,                         &
# endif
     &                      GRID(ng)%Gx_rho,                            &
     &                      GRID(ng)%Gy_rho,                            &
# ifdef MASKING
     &                      EXTRACT(ng)%Gmask_rho,                      &
# endif
     &                      EXTRACT(ng)%Gx_rho,                         &
     &                      EXTRACT(ng)%Gy_rho,                         &
     &                      EXTRACT(ng)%Iout_rho,                       &
     &                      EXTRACT(ng)%Jout_rho)
!
        CALL interp_coords (ng, tile, model, u2dvar,                    &
     &                      GRID(ng)%Gangle_u,                          &
# ifdef MASKING
     &                      GRID(ng)%Gmask_u,                           &
# endif
     &                      GRID(ng)%Gx_u,                              &
     &                      GRID(ng)%Gy_u,                              &
# ifdef MASKING
     &                      EXTRACT(ng)%Gmask_u,                        &
# endif
     &                      EXTRACT(ng)%Gx_u,                           &
     &                      EXTRACT(ng)%Gy_u,                           &
     &                      EXTRACT(ng)%Iout_u,                         &
     &                      EXTRACT(ng)%Jout_u)
!
        CALL interp_coords (ng, tile, model, v2dvar,                    &
     &                      GRID(ng)%Gangle_v,                          &
# ifdef MASKING
     &                      GRID(ng)%Gmask_v,                           &
# endif
     &                      GRID(ng)%Gx_v,                              &
     &                      GRID(ng)%Gy_v,                              &
# ifdef MASKING
     &                      EXTRACT(ng)%Gmask_v,                        &
# endif
     &                      EXTRACT(ng)%Gx_v,                           &
     &                      EXTRACT(ng)%Gy_v,                           &

     &                      EXTRACT(ng)%Iout_v,                         &
     &                      EXTRACT(ng)%Jout_v)
      END IF
!
  10  FORMAT (' GET_EXTRACT - Illegal input file type, io_type = ',i0,  &
     &        /,12x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_extract
!
!***********************************************************************
      SUBROUTINE get_extract_nf90 (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: cr, gtype, i, status, vindex
      integer :: Vsize(4)
# ifdef CHECKSUM
      integer(i8b) :: Fhash
# endif
!
      real(dp), parameter :: Fscl = 1.0_dp

      real(r8) :: Fmax, Fmin
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_extract_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of extraction grid NetCDF file: Inquire
!  about the dimensions and variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=GRX(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (GRX(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, GRX(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = GRX(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'spherical',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'h',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pm',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pn',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# if (defined CURVGRID && defined UV_ADV)
      IF (.not.find_string(var_name,n_var,'dndx',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dndx', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dmde',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dmde', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef CURVGRID
      IF (.not.find_string(var_name,n_var,'angle',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'angle', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef MASKING
      IF (.not.find_string(var_name,n_var,'mask_rho',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_u',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_v',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_psi',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL netcdf_get_lvar (ng, model, ncname, 'spherical',           &
     &                        spherical,                                &
     &                        ncid = GRX(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid variables.
!-----------------------------------------------------------------------
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Scan the variable list and read in needed variables.
!
      IF (Master) WRITE (stdout,'(1x)')
!
      DO i=1,n_var

        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in bathymetry.
!
          CASE ('h')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % h,                      &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % h)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'bathymetry at RHO-points: h',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % h)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % h)
# endif
# ifdef MASKING
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % rmask,                  &
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % rmask,                  &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % rmask)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on RHO-points: mask_rho',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % rmask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % rmask)
#  endif
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            gtype=u2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % umask,                  &
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % umask,                  &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % umask)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on U-points: mask_u',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % umask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % umask)
#  endif
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            gtype=v2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % vmask,                  &
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % vmask,                  &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % vmask)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on V-points: mask_v',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % vmask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % vmask)
#  endif
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            gtype=p2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % pmask,                  &
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % pmask,                  &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % pmask)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on PSI-points: mask_psi',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_p2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pmask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pmask)
#  endif
# endif
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % f,                      &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % f)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Coriolis parameter at RHO-points: f',&
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % f)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % f)
# endif
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % pm,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % pm)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal XI-grid spacing: pm',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pm)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pm)
# endif
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % pn,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % pn)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal ETA-grid spacing: pn',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pn)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pn)
# endif
# if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % dmde,                   &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % dmde)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'ETA-derivative of inverse metric '// &
     &                            'factor pm: dmde',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % dmde)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % dmde)
#  endif
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % dndx,                   &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % dndx)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'XI-derivative of inverse metric '//  &
     &                            'factor pn: dndx',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % dndx)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % dndx)
#  endif
# endif
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            gtype=p2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % pmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % xp,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % xp)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of PSI-points: x_psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xp)
# endif
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            gtype=p2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % pmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % yp,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % yp)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of PSI-points: y-psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yp)
# endif
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % xr,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % xr)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of RHO-points: x-rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (.not.spherical) THEN
              EXTRACT(ng)%LonMin(ng)=Fmin
              EXTRACT(ng)%LonMax(ng)=Fmax
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xr)
# endif
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % yr,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % yr)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of RHO-points: y_rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (.not.spherical) THEN
              EXTRACT(ng)%LatMin(ng)=Fmin
              EXTRACT(ng)%LatMax(ng)=Fmax
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yr)
# endif
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            gtype=u2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % umask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % xu,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % xu)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of U-points: x_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xu)
# endif
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            gtype=u2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % umask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % yu,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % yu)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of U-points: y_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yu)
# endif
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            gtype=v2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % vmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % xv,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % xv)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of V-points: x_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xv)
# endif
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            gtype=v2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % vmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % yv,                     &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % yv)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of V-points: y_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yv)
# endif
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % pmask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % lonp,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % lonp)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of PSI-points: lon_psi', &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonp)
# endif
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % pmask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % latp,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % latp)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of PSI-points lat_psi',   &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latp)
# endif
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % rmask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % lonr,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % lonr)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of RHO-points: lon_rho', &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
              EXTRACT(ng)%LonMin(ng)=Fmin
              EXTRACT(ng)%LonMax(ng)=Fmax

# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonr)
# endif
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % rmask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % latr,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % latr)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of RHO-points lat_rho',   &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
              EXTRACT(ng)%LatMin(ng)=Fmin
              EXTRACT(ng)%LatMax(ng)=Fmax
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latr)
# endif
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % umask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % lonu,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % lonu)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of U-points: lon_u',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonu)
# endif
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % umask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % latu,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % latu)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of U-points: lat_u',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latu)
# endif
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % vmask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % lonv,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % lonv)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of V-points: lon_v',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonv)
# endif
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,    &
     &                              var_name(i), var_id(i),             &
     &                              0, gtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
# ifdef MASKING
     &                              EXTRACT(ng) % vmask,                &
# endif
# ifdef CHECKSUM
     &                              EXTRACT(ng) % latv,                 &
     &                              checksum = Fhash)
# else
     &                              EXTRACT(ng) % latv)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of V-points: lat_v',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latv)
# endif
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            gtype=r2dvar
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%ncid,      &
     &                            var_name(i), var_id(i),               &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
# ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
# endif
# ifdef CHECKSUM
     &                            EXTRACT(ng) % angler,                 &
     &                            checksum = Fhash)
# else
     &                            EXTRACT(ng) % angler)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'angle between XI-axis and EAST: '//  &
     &                            'angler',                             &
     &                            ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % angler)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % angler)
# endif
        END SELECT
      END DO
!
! Close GRID NetCDF file.
!
      CALL netcdf_close (ng, model, GRX(ng)%ncid, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' GET_EXTRACT_NF90 - unable to open grid NetCDF file:', &
     &        1x,a)
  20  FORMAT (/,' GET_EXTRACT_NF90 - unable to find grid variable: ',   &
     &        a,/,20x,'in grid NetCDF file: ',a)
  30  FORMAT (2x,'GET_EXTRACT_NF90 - ',a,/,22x,                         &
     &        '(Grid = ',i2.2,', File: ',a,')',/,22x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  40  FORMAT (/,' GET_EXTRACT_NF90 - error while reading variable: ',   &
     &        a,/,20x,'in grid NetCDF file: ',a)
  50  FORMAT (/,2x,'GET_EXTRACT_NF90 - Reading adjoint sensitivity',    &
     &        ' scope arrays from file:',/22x,a,/)
# ifdef CHECKSUM
  60  FORMAT (22x,'(CheckSum = ',i0,')')
# endif
!
      RETURN
      END SUBROUTINE get_extract_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_extract_pio (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: cr, i, status, vindex
      integer :: Vsize(4)
#  ifdef CHECKSUM
      integer(i8b) :: Fhash
#  endif
!
      real(dp), parameter :: Fscl = 1.0_dp

      real(r8) :: Fmax, Fmin
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_extract_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
      TYPE (My_VarDesc) :: pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=GRX(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (GRX(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, GRX(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, model, ncname,                       &
     &                         pioFile = GRX(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'spherical',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'h',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pm',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pn',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined CURVGRID && defined UV_ADV)
      IF (.not.find_string(var_name,n_var,'dndx',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dndx', TRIM(ncname)
!       exit_flag=2
!       RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dmde',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dmde', TRIM(ncname)
!       exit_flag=2
!       RETURN
      END IF
#  endif
#  ifdef CURVGRID
      IF (.not.find_string(var_name,n_var,'angle',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'angle', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef MASKING
      IF (.not.find_string(var_name,n_var,'mask_rho',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_u',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_v',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_psi',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL pio_netcdf_get_lvar (ng, model, ncname,                    &
     &                            'spherical', spherical,               &
     &                            pioFile = GRX(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid variables.
!-----------------------------------------------------------------------
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Scan the variable list and read in needed variables.
!
      IF (Master) WRITE (stdout,'(1x)')
!
      DO i=1,n_var

        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in bathymetry.
!
          CASE ('h')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%h).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % h,                      &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % h)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
#  ifdef SINGLE_PRECISION
              EXTRACT(ng)%Hmin(ng)=REAL(Fmin,dp)
              EXTRACT(ng)%Hmax(ng)=REAL(Fmax,dp)
#  else
              EXTRACT(ng)%Hmin(ng)=Fmin
              EXTRACT(ng)%Hmax(ng)=Fmax
#  endif
              IF (Master) THEN
                WRITE (stdout,30) 'bathymetry at RHO-points: h',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % h)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % h)
#  endif
#  ifdef MASKING
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%rmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % rmask,                  &
#   ifdef CHECKSUM
     &                            EXTRACT(ng) % rmask,                  &
     &                            checksum = Fhash)
#   else
     &                            EXTRACT(ng) % rmask)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on RHO-points: mask_rho',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % rmask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % rmask)
#   endif
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(EXTRACT(ng)%umask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % umask,                  &
#   ifdef CHECKSUM
     &                            EXTRACT(ng) % umask,                  &
     &                            checksum = Fhash)
#   else
     &                            EXTRACT(ng) % umask)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on U-points: mask_u',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % umask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % umask)
#   endif
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(EXTRACT(ng)%vmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % vmask,                  &
#   ifdef CHECKSUM
     &                            EXTRACT(ng) % vmask,                  &
     &                            checksum = Fhash)
#   else
     &                            EXTRACT(ng) % vmask)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on V-points: mask_v',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % vmask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % vmask)
#   endif
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(EXTRACT(ng)%pmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
     &                            EXTRACT(ng) % pmask,                  &
#   ifdef CHECKSUM
     &                            EXTRACT(ng) % pmask,                  &
     &                            checksum = Fhash)
#   else
     &                            EXTRACT(ng) % pmask)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on PSI-points: mask_psi',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_p2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pmask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pmask)
#   endif
#  endif
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%pn).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % f,                      &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % f)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Coriolis parameter at RHO-points: f',&
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pn)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pn)
#  endif
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%pn).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % pm,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % pm)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal XI-grid spacing: pm',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pn)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pn)
#  endif
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%pn).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % pn,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % pn)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal ETA-grid spacing: pn',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % pn)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % pn)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%dmde).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#   endif
#   ifdef CHECKSUM
     &                            EXTRACT(ng) % dmde,                   &
     &                            checksum = Fhash)
#   else
     &                            EXTRACT(ng) % dmde)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'ETA-derivative of inverse metric '// &
     &                            'factor pm: dmde',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % dmde)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % dmde)
#   endif
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%dndx).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#   endif
#   ifdef CHECKSUM
     &                            EXTRACT(ng) % dndx,                   &
     &                            checksum = Fhash)
#   else
     &                            EXTRACT(ng) % dndx)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'XI-derivative of inverse metric '//  &
     &                            'factor pn: dndx',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % dndx)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % dndx)
#   endif
#  endif
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(EXTRACT(ng)%xp).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % pmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % xp,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % xp)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of PSI-points: x_psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xp)
#  endif
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(EXTRACT(ng)%yp).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % pmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % yp,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % yp)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of PSI-points: y-psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yp)
#  endif
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%xr).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % xr,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % xr)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of RHO-points: x-rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (.not.spherical) THEN
              EXTRACT(ng)%LonMin(ng)=Fmin
              EXTRACT(ng)%LonMax(ng)=Fmax
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xr)
#  endif
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%yr).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % yr,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % yr)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of RHO-points: y_rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (.not.spherical) THEN
              EXTRACT(ng)%LatMin(ng)=Fmin
              EXTRACT(ng)%LatMax(ng)=Fmax
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yr)
#  endif
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(EXTRACT(ng)%xu).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % umask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % xu,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % xu)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of U-points: x_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xu)
#  endif
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(EXTRACT(ng)%yu).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % umask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % yu,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % yu)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of U-points: y_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yu)
#  endif
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(EXTRACT(ng)%xv).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % vmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % xv,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % xv)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of V-points: x_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % xv)
#  endif
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(EXTRACT(ng)%yv).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % vmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % yv,                     &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % yv)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of V-points: y_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              .FALSE., .FALSE.,                   &
     &                              EXTRACT(ng) % yv)
#  endif
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=p2dvar
              IF (KIND(EXTRACT(ng)%lonp).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % pmask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % lonp,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % lonp)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of PSI-points: lon_psi', &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonp)
#  endif
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=p2dvar
              IF (KIND(EXTRACT(ng)%latp).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % pmask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % latp,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % latp)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of PSI-points lat_psi',   &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latp)
#  endif
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(EXTRACT(ng)%lonr).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % rmask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % lonr,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % lonr)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of RHO-points: lon_rho', &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
              EXTRACT(ng)%LonMin(ng)=Fmin
              EXTRACT(ng)%LonMax(ng)=Fmax
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonr)
#  endif
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(EXTRACT(ng)%latr).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % rmask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % latr,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % latr)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of RHO-points lat_rho',   &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
              EXTRACT(ng)%LatMin(ng)=Fmin
              EXTRACT(ng)%LatMax(ng)=Fmax
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latr)
#  endif
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=u2dvar
              IF (KIND(EXTRACT(ng)%lonu).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % umask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % lonu,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % lonu)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of U-points: lon_u',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonu)
#  endif
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=u2dvar
              IF (KIND(EXTRACT(ng)%latu).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % umask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % latu,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % latu)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of U-points: lat_u',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latu)
#  endif
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=v2dvar
              IF (KIND(EXTRACT(ng)%lonv).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % vmask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % lonv,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % lonv)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of V-points: lon_v',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % lonv)
#  endif
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=v2dvar
              IF (KIND(EXTRACT(ng)%latv).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile, &
     &                              var_name(i), pioVar,                &
     &                              0, ioDesc, Vsize,                   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Fscl, Fmin, Fmax,                   &
#  ifdef MASKING
     &                              EXTRACT(ng) % vmask,                &
#  endif
#  ifdef CHECKSUM
     &                              EXTRACT(ng) % latv,                 &
     &                              checksum = Fhash)
#  else
     &                              EXTRACT(ng) % latv)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of V-points: lat_v',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d_xtr (ng, tile, model, 1,               &
     &                                LBi, UBi, LBj, UBj,               &
     &                                NghostPoints,                     &
     &                                .FALSE., .FALSE.,                 &
     &                                EXTRACT(ng) % latv)
#  endif
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(EXTRACT(ng)%angler).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d_xtr(ng, model, ncname, GRX(ng)%pioFile,   &
     &                            var_name(i), pioVar,                  &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            EXTRACT(ng) % rmask,                  &
#  endif
#  ifdef CHECKSUM
     &                            EXTRACT(ng) % angler,                 &
     &                            checksum = Fhash)
#  else
     &                            EXTRACT(ng) % angler)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'angle between XI-axis and EAST: '//  &
     &                            'angler',                             &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_xtr_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    EXTRACT(ng) % angler)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_xtr (ng, tile, model, 1,                 &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              EXTRACT(ng) % angler)
#  endif
        END SELECT
      END DO
!
! Close GRID NetCDF file.
!
      CALL pio_netcdf_close (ng, model, GRX(ng)%pioFile, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' GET_EXTRACT_PIO - unable to open grid NetCDF file:',  &
     &        1x,a)
  20  FORMAT (/,' GET_EXTRACT_PIO - unable to find grid variable: ',a,  &
     &        /,19x,'in grid NetCDF file: ',a)
  30  FORMAT (2x,'GET_EXTRACT_PIO - ',a,/,22x,                          &
     &        '(Grid = ',i2.2,', File: ',a,')',/,22x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  40  FORMAT (/,' GET_EXTRACT_PIO - error while reading variable: ',a,  &
     &        /,12x,'in grid NetCDF file: ',a)
  50  FORMAT (/,2x,'GET_EXTRACT_PIO - Reading adjoint sensitivity',     &
     &        ' scope arrays from file:',/22x,a,/)
#  ifdef CHECKSUM
  60  FORMAT (22x,'(CheckSum = ',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE get_extract_pio
# endif
#endif
      END MODULE get_extract_mod
