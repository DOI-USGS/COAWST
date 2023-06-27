#include "cppdefs.h"
      MODULE get_grid_mod

#ifndef ANA_GRID
!
!git $Id$
!svn $Id: get_grid.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module reads grid information from input file using either     !
!  the standard NetCDF library or the Parallel-IO (PIO) library.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# ifdef NESTING
      USE nesting_mod,     ONLY : fill_contact
# endif
      USE nf_fread2d_mod,  ONLY : nf_fread2d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      PUBLIC  :: get_grid
      PRIVATE :: get_grid_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: get_grid_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_grid (ng, tile, model)
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
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (GRD(ng)%IOtype)
        CASE (io_nf90)
          CALL get_grid_nf90 (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_grid_pio (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) GRD(ng)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_GRID - Illegal input file type, io_type = ',i0,     &
     &        /,12x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_grid
!
!***********************************************************************
      SUBROUTINE get_grid_nf90 (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
     defined SO_SEMI
#  ifndef OBS_SPACE
      logical :: GotScope(6)
!
#  endif
# endif
      integer :: cr, gtype, i, status, vindex
# if defined UV_DRAG_GRID && !defined ANA_DRAG
      integer :: varid_dragL, varid_dragQ, varid_ZoBL
# endif
#if defined UV_WAVEDRAG
      integer :: varid_dragW
#endif
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
     &  __FILE__//", get_grid_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=GRD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (GRD(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, GRD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency.
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = GRD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = GRD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef NESTING
!
!  Determine contact region index "cr" for which nested grid "ng" is
!  the receiver grid.
!
      DO i=1,Ncontact
        IF (Rcontact(i)%receiver_grid.eq.ng) THEN
          cr=i
          EXIT
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'xl',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'el',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
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
# ifdef ICESHELF
      IF (.not.find_string(var_name,n_var,'zice',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'zice', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
      IF (.not.find_string(var_name,n_var,'f',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
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
# if defined WTYPE_GRID && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE) && \
    !defined ANA_WTYPE
      IF (.not.find_string(var_name,n_var,'wtype_grid',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'wtype_grid', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifndef ANA_SPONGE
      IF (LuvSponge(ng)) THEN
        IF (.not.find_string(var_name,n_var,'visc_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'visc_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#  ifdef SOLVE3D
      IF (ANY(LtracerSponge(:,ng))) THEN
        IF (.not.find_string(var_name,n_var,'diff_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'diff_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#  endif
# endif
# if defined UV_DRAG_GRID && !defined ANA_DRAG
#  if defined UV_LOGDRAG || defined BBL_MODEL
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idZoBL)),        &
     &                     varid_ZoBL)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idZoBL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef UV_LDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idragL)),        &
     &                     varid_dragL)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef UV_QDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idragQ)),        &
     &                     varid_dragQ)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragQ)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
# endif
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL netcdf_get_lvar (ng, model, ncname, 'spherical',           &
     &                        spherical,                                &
     &                        ncid = GRD(ng)%ncid)
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
!  Read in basin X-length.
!
          CASE ('xl')
            CALL netcdf_get_fvar (ng, model, ncname,                    &
     &                            'xl', xl(ng),                         &
     &                            ncid = GRD(ng)%ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) EXIT
!
!  Read in basin Y-length.
!
          CASE ('el')
            CALL netcdf_get_fvar (ng, model, ncname,                    &
     &                            'el', el(ng),                         &
     &                            ncid = GRD(ng)%ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) EXIT
!
!  Read in bathymetry.
!
          CASE ('h')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % h,                             &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % h)
# endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
# ifdef SINGLE_PRECISION
              Hmin(ng)=REAL(Fmin,dp)
              Hmax(ng)=REAL(Fmax,dp)
# else
              Hmin(ng)=Fmin
              Hmax(ng)=Fmax
# endif
              IF (Master) THEN
                WRITE (stdout,30) 'bathymetry at RHO-points: h',        &
     &                            ng, TRIM(ncname), hmin(ng), hmax(ng)
# ifdef CHECKSUM
                WRITE (stdout,60) Fhash
# endif
              END IF
            END IF
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % h,                   &
     &                        GRID(ng) % h)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % h)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % h)
# endif
# ifdef MASKING
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
#  ifdef CHECKSUM
     &                        GRID(ng) % rmask,                         &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % rmask)
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
#  ifdef NESTING
#   if !defined NOFILL_NESTING_MASK
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, 'rmask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % rmask,               &
     &                        GRID(ng) % rmask)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % rmask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % rmask)
#  endif
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
#  ifdef CHECKSUM
     &                        GRID(ng) % umask,                         &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % umask)
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
#  ifdef NESTING
#   if !defined NOFILL_NESTING_MASK
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        u2dvar, 'umask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % umask,               &
     &                        GRID(ng) % umask)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % umask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % umask)
#  endif
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
#  ifdef CHECKSUM
     &                        GRID(ng) % vmask,                         &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % vmask)
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
#  ifdef NESTING
#   if !defined NOFILL_NESTING_MASK
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        v2dvar, 'vmask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % vmask,               &
     &                        GRID(ng) % vmask)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % vmask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % vmask)
#  endif
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
#  ifdef CHECKSUM
     &                        GRID(ng) % pmask,                         &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % pmask)
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
              CALL exchange_p2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pmask)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pmask)
#  endif
# endif
# ifdef ICESHELF
!
!  Read in ice shelf thickness.
!
          CASE ('zice')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % zice,                          &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % zice)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'ice shelf thickness: zice',          &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % zice)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % zice)
#  endif
# endif
# if defined WTYPE_GRID && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE) && \
    !defined ANA_WTYPE
!
!  Read in Jerlov water type.
!
          CASE ('wtype_grid')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        MIXING(ng) % Jwtype,                      &
     &                        checksum = Fhash)
#  else
     &                        MIXING(ng) % Jwtype)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Jerlov water type: wtype_grid',      &
     &                            ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                MIXING(ng) % Jwtype)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          MIXING(ng) % Jwtype)
#  endif
# endif
# ifndef ANA_SPONGE
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  viscosity (nondimensional) in specific areas of the domain.
!
          CASE ('visc_factor')
            IF (LuvSponge(ng)) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          MIXING(ng) % visc_factor,               &
     &                          checksum = Fhash)
#  else
     &                          MIXING(ng) % visc_factor)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal viscosity sponge '//    &
     &                              'factor: visc_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % visc_factor)
              END IF
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % visc_factor)
#  endif
            END IF

#  ifdef SOLVE3D
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  diffusivity (nondimensional) in specific areas of the domain.
!
          CASE ('diff_factor')
            IF (ANY(LtracerSponge(:,ng))) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          MIXING(ng) % diff_factor,               &
     &                          checksum = Fhash)
#   else
     &                          MIXING(ng) % diff_factor)
#   endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal diffusivity sponge '//  &
     &                              'factor: diff_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % diff_factor)
              END IF
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % diff_factor)
#   endif
            END IF
#  endif
# endif
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % f,                             &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % f)
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
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % f,                   &
     &                        GRID(ng) % f)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % f)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % f)
# endif
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % pm,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % pm)
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
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % pm,                  &
     &                        GRID(ng) % pm)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pm)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pm)
# endif
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % pn,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % pn)
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
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % pn,                  &
     &                        GRID(ng) % pn)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pn)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pn)
# endif
# if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % dmde,                          &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % dmde)
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
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % dmde,                &
     &                        GRID(ng) % dmde)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dmde)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dmde)
#  endif
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % dndx,                          &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % dndx)
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
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % dndx,                &
     &                        GRID(ng) % dndx)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dndx)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dndx)
#  endif
# endif
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % pmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % xp,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % xp)
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
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xp)
# endif
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % pmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % yp,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % yp)
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
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yp)
# endif
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % xr,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % xr)
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
              LonMin(ng)=Fmin
              LonMax(ng)=Fmax
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xr,                &
     &                          GRID(ng) % xr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xr)
# endif
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % yr,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % yr)
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
              LatMin(ng)=Fmin
              LatMax(ng)=Fmax
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yr,                &
     &                          GRID(ng) % yr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yr)
# endif
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % xu,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % xu)
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
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xu,                &
     &                          GRID(ng) % xu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xu)
# endif
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % yu,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % yu)
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
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yu,                &
     &                          GRID(ng) % yu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yu)
# endif
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % xv,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % xv)
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
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xv,                &
     &                          GRID(ng) % xv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xv)
# endif
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % yv,                            &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % yv)
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
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yv,                &
     &                          GRID(ng) % yv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yv)
# endif
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % pmask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % lonp,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % lonp)
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
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonp)
# endif
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % pmask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % latp,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % latp)
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
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latp)
# endif
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LonMin(ng), LonMax(ng),           &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % lonr,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % lonr)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of RHO-points: lon_rho', &
     &                              ng, TRIM(ncname),                   &
     &                              LonMin(ng), LonMax(ng)
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xr,                &
     &                          GRID(ng) % lonr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonr)
# endif
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LatMin(ng), LatMax(ng),           &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % latr,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % latr)
# endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of RHO-points lat_rho',   &
     &                              ng, TRIM(ncname),                   &
     &                              LatMin(ng), LatMax(ng)
# ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
# endif
                END IF
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yr,                &
     &                          GRID(ng) % latr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latr)
# endif
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % lonu,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % lonu)
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
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xu,                &
     &                          GRID(ng) % lonu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonu)
# endif
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % latu,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % latu)
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
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yu,                &
     &                          GRID(ng) % latu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latu)
# endif
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % lonv,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % lonv)
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
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xv,                &
     &                          GRID(ng) % lonv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonv)
# endif
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef CHECKSUM
     &                          GRID(ng) % latv,                        &
     &                          checksum = Fhash)
# else
     &                          GRID(ng) % latv)
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
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yv,                &
     &                          GRID(ng) % latv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latv)
# endif
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef CHECKSUM
     &                        GRID(ng) % angler,                        &
     &                        checksum = Fhash)
# else
     &                        GRID(ng) % angler)
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
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, 'angler', spval_check,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % angler,              &
     &                        GRID(ng) % angler)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % angler)
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % angler)
# endif
# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
     defined SO_SEMI
#  ifndef OBS_SPACE
!
!  Read in adjoint sensitivity spatial scope masking at RHO-points.
!
          CASE ('scope_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % Rscope,                        &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % Rscope)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on RHO-points: scope_rho',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            GotScope(1)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Rscope)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Rscope)
#   endif
!
!  Read in adjoint sensitivity spatial scope masking at U-points.
!
          CASE ('scope_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % umask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % Uscope,                        &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % Uscope)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on U-points: scope_u',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            GotScope(2)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Uscope)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Uscope)
#   endif
!
!  Read in adjoint sensitivity spatial scope masking at V-points.
!
          CASE ('scope_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % Vscope,                        &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % Vscope)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on V-points: scope_v',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            GotScope(3)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Vscope)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Vscope)
#   endif
#  endif
# endif
        END SELECT
      END DO
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF

# if defined UV_DRAG_GRID && !defined ANA_DRAG
#  if defined UV_LOGDRAG || defined BBL_MODEL
!
!  Read in spacially varying bottom roughness length (m).
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,                &
     &                  Vname(1,idZoBL), varid_ZoBL,                    &
     &                  0, gtype, Vsize,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
#   ifdef CHECKSUM
     &                  GRID(ng) % ZoBot,                               &
     &                  checksum = Fhash)
#   else
     &                  GRID(ng) % ZoBot)
#   endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZoBL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) 'time invariant, bottom roughness '//       &
     &                      'length scale: ZoBot',                      &
     &                      ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
          WRITE (stdout,60) Fhash
#   endif
        END IF
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % ZoBot)
      END IF
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % ZoBot)
#   endif
#  endif
#  ifdef UV_LDRAG
!
!  Read in spacially varying linear drag coefficients (m/s).
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,                &
     &                  Vname(1,idragL), varid_dragL,                   &
     &                  0, gtype, Vsize,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
#   ifdef CHECKSUM
     &                  GRID(ng) % rdrag,                               &
     &                  checksum = Fhash)
#   else
     &                  GRID(ng) % rdrag)
#   endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idragL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) 'linear bottom drag coefficient: rdrag',    &
     &                      ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
        END IF
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % rdrag)
      END IF
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % rdrag)
#   endif
#  endif
#  ifdef UV_QDRAG
!
!  Read in spacially varying quadratic drag coefficients.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,                &
     &                  Vname(1,idragQ), varid_dragQ,                   &
     &                  0, gtype, Vsize,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
#   ifdef CHECKSUM
     &                  GRID(ng) % rdrag2,                              &
     &                  checksum = Fhash)
#   else
     &                  GRID(ng) % rdrag2)
#   endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idragQ)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) 'quadratic bottom drag coefficient: rdrag2',&
     &                      ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
        END IF
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % rdrag2)
      END IF
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % rdrag2)
#   endif
#  endif
# endif
# if defined UV_WAVEDRAG
!
!  Read in spacially varying linear drag coefficients (m/s).
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,                &
     &                  Vname(1,idragW), varid_dragW,                   &
     &                  0, gtype, Vsize,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  GRID(ng) % wavedrag)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idragW)),        &
     &                                TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % wavedrag)
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, model, 1,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        GRID(ng) % wavedrag)
#  endif
# endif
!
! Close GRID NetCDF file.
!
      CALL netcdf_close (ng, model, GRD(ng)%ncid, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
     defined SO_SEMI
#  ifndef OBS_SPACE
!
!-----------------------------------------------------------------------
!  Inquire adjoint sensitivity forcing file. Read scope arrays again.
!  These fields take precedence
!-----------------------------------------------------------------------
!
      ncname=ADS(ng)%name
!
!  Open adjoint sensitivity NetCDF file for reading.
!
      IF (ADS(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, ADS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = ADS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ADS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Check if the adjoint sensitivity scope arrays are available.
!
      GotScope(4)=find_string(var_name,n_var,'scope_rho',vindex)
      GotScope(5)=find_string(var_name,n_var,'scope_u',vindex)
      GotScope(6)=find_string(var_name,n_var,'scope_v',vindex)
!
      IF ((.not.GotScope(1)).and.(.not.GotScope(4))) THEN
        IF (Master) WRITE (stdout,20) 'scope_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((.not.GotScope(2)).and.(.not.GotScope(5))) THEN
        IF (Master) WRITE (stdout,20) 'scope_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((.not.GotScope(3)).and.(.not.GotScope(6))) THEN
        IF (Master) WRITE (stdout,20) 'scope_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (Master) THEN
        IF (GotScope(4)) THEN
          WRITE (stdout,50) TRIM(ADS(ng)%name)
        ELSE
          WRITE (stdout,50) TRIM(GRD(ng)%name)
        END IF
      END IF
!
!  Scan adjoint sensitivity variables.
!
      DO i=1,n_var

        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in adjoint sensitivity spatial scope masking at RHO-points.
!
          CASE ('scope_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ADS(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % Rscope,                        &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % Rscope)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on RHO-points: scope_rho',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Rscope)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Rscope)
#   endif
!
!  Read in adjoint sensitivity spatial scope masking at U-points.
!
          CASE ('scope_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, ADS(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % umask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % Uscope,                        &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % Uscope)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on U-points: scope_u',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Uscope)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Uscope)
#   endif
!
!  Read in adjoint sensitivity spatial scope masking at V-points.
!
          CASE ('scope_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, ADS(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % Vscope,                        &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % Vscope)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on V-points: scope_v',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Vscope)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Vscope)
#   endif
        END SELECT
      END DO
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF
#  endif
# endif
!
  10  FORMAT (/,' GET_GRID_NF90 - unable to open grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID_NF90 - unable to find grid variable: ',a,    &
     &        /,12x,'in grid NetCDF file: ',a)
  30  FORMAT (2x,'GET_GRID_NF90    - ',a,/,22x,                         &
     &        '(Grid = ',i2.2,', File: ',a,')',/,22x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  40  FORMAT (/,' GET_GRID_NF90 - error while reading variable: ',a,    &
     &        /,12x,'in grid NetCDF file: ',a)
  50  FORMAT (/,2x,'GET_GRID_NF90    - Reading adjoint sensitivity',    &
     &        ' scope arrays from file:',/22x,a,/)

# ifdef CHECKSUM
  60  FORMAT (22x,'(CheckSum = ',i0,')')
# endif
!
      RETURN
      END SUBROUTINE get_grid_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_grid_pio (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
     defined SO_SEMI
#   ifndef OBS_SPACE
      logical :: GotScope(6)
!
#   endif
#  endif
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
     &  __FILE__//", get_grid_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
      TYPE (My_VarDesc) :: pioVar
#  if defined UV_DRAG_GRID && !defined ANA_DRAG
      TYPE (My_VarDesc) :: pioVar_dragL, pioVar_dragQ, pioVar_ZoBL
#  endif

!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=GRD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (GRD(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, GRD(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency.
!
      CALL pio_netcdf_check_dim (ng, model, ncname,                     &
     &                           pioFile = GRD(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, model, ncname,                       &
     &                         pioFile = GRD(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef NESTING
!
!  Determine contact region index "cr" for which nested grid "ng" is
!  the receiver grid.
!
      DO i=1,Ncontact
        IF (Rcontact(i)%receiver_grid.eq.ng) THEN
          cr=i
          EXIT
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'xl',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'el',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
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
#  ifdef ICESHELF
      IF (.not.find_string(var_name,n_var,'zice',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'zice', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
      IF (.not.find_string(var_name,n_var,'f',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
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
#  if defined WTYPE_GRID && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE) && \
    !defined ANA_WTYPE
      IF (.not.find_string(var_name,n_var,'wtype_grid',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'wtype_grid', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifndef ANA_SPONGE
      IF (LuvSponge(ng)) THEN
        IF (.not.find_string(var_name,n_var,'visc_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'visc_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#   ifdef SOLVE3D
      IF (ANY(LtracerSponge(:,ng))) THEN
        IF (.not.find_string(var_name,n_var,'diff_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'diff_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#   endif
#  endif
#  if defined UV_DRAG_GRID && !defined ANA_DRAG
#   ifdef UV_LOGDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idZoBL)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idZoBL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      ELSE
        pioVar_ZoBL%vd=var_desc(vindex)
        pioVar_ZoBL%gtype=r2dvar
      END IF
#   endif
#   ifdef UV_LDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idragL)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      ELSE
        pioVar_dragL%vd=var_desc(vindex)
        pioVar_dragL%gtype=r2dvar
      END IF
#   endif
#   ifdef UV_QDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idragQ)),        &
     &                     vindex)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragQ)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      ELSE
        pioVar_dragQ%vd=var_desc(vindex)
        pioVar_dragQ%gtype=r2dvar
      END IF
#   endif
#  endif
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL pio_netcdf_get_lvar (ng, model, ncname,                    &
     &                            'spherical', spherical,               &
     &                            pioFile = GRD(ng)%pioFile)
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
!  Read in basin X-length.
!
          CASE ('xl')
            CALL pio_netcdf_get_fvar (ng, model, ncname,                &
     &                                'xl', xl(ng),                     &
     &                                pioFile = GRD(ng)%pioFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) EXIT
!
!  Read in basin Y-length.
!
          CASE ('el')
            CALL pio_netcdf_get_fvar (ng, model, ncname,                &
     &                                'el', el(ng),                     &
     &                                pioFile = GRD(ng)%pioFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) EXIT
!
!  Read in bathymetry.
!
          CASE ('h')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%h).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % h,                             &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % h)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
#  ifdef SINGLE_PRECISION
              Hmin(ng)=REAL(Fmin,dp)
              Hmax(ng)=REAL(Fmax,dp)
#  else
              Hmin(ng)=Fmin
              Hmax(ng)=Fmax
#  endif
              IF (Master) THEN
                WRITE (stdout,30) 'bathymetry at RHO-points: h',        &
     &                            ng, TRIM(ncname), hmin(ng), hmax(ng)
#  ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#  endif
              END IF
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % h,                   &
     &                        GRID(ng) % h)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % h)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % h)
#  endif
#  ifdef MASKING
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%rmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
#   ifdef CHECKSUM
     &                        GRID(ng) % rmask,                         &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % rmask)
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
#   ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, 'rmask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % rmask,               &
     &                        GRID(ng) % rmask)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % rmask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % rmask)
#   endif
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%umask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
#   ifdef CHECKSUM
     &                        GRID(ng) % umask,                         &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % umask)
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
#   ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        u2dvar, 'umask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % umask,               &
     &                        GRID(ng) % umask)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % umask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % umask)
#   endif
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%vmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
#   ifdef CHECKSUM
     &                        GRID(ng) % vmask,                         &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % vmask)
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
#   ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        v2dvar, 'vmask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % vmask,               &
     &                        GRID(ng) % vmask)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % vmask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % vmask)
#   endif
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(GRID(ng)%pmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
#   ifdef CHECKSUM
     &                        GRID(ng) % pmask,                         &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % pmask)
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
              CALL exchange_p2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pmask)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pmask)
#   endif
#  endif
#  ifdef ICESHELF
!
!  Read in ice shelf thickness.
!
          CASE ('zice')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%zice).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % zice,                          &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % zice)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'ice shelf thickness: zice',          &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % zice)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % zice)
#   endif
#  endif
#  if defined WTYPE_GRID && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE) && \
    !defined ANA_WTYPE
!
!  Read in Jerlov water type.
!
          CASE ('wtype_grid')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%Jwtype).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        MIXING(ng) % Jwtype,                      &
     &                        checksum = Fhash)
#   else
     &                        MIXING(ng) % Jwtype)
#   endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Jerlov water type: wtype_grid',      &
     &                            ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                MIXING(ng) % Jwtype)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          MIXING(ng) % Jwtype)
#   endif
#  endif
#  ifndef ANA_SPONGE
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  viscosity (nondimensional) in specific areas of the domain.
!
          CASE ('visc_factor')
            IF (LuvSponge(ng)) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(MIXING(ng)%visc_factor).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          MIXING(ng) % visc_factor,               &
     &                          checksum = Fhash)
#   else
     &                          MIXING(ng) % visc_factor)
#   endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal viscosity sponge '//    &
     &                              'factor: visc_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
#   ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#   endif
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % visc_factor)
              END IF
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % visc_factor)
#   endif
            END IF

#   ifdef SOLVE3D
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  diffusivity (nondimensional) in specific areas of the domain.
!
          CASE ('diff_factor')
            IF (ANY(LtracerSponge(:,ng))) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(MIXING(ng)%diff_factor).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#    ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#    endif
#    ifdef CHECKSUM
     &                          MIXING(ng) % diff_factor,               &
     &                          checksum = Fhash)
#    else
     &                          MIXING(ng) % diff_factor)
#    endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal diffusivity sponge '//  &
     &                              'factor: diff_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % diff_factor)
              END IF
#    ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % diff_factor)
#    endif
            END IF
#   endif
#  endif
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%f).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % f,                             &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % f)
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
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % f,                   &
     &                        GRID(ng) % f)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % f)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % f)
#  endif
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%pm).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % pm,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % pm)
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
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % pm,                  &
     &                        GRID(ng) % pm)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pm)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pm)
#  endif
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%pn).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % pn,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % pn)
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
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % pn,                  &
     &                        GRID(ng) % pn)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pn)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pn)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%dmde).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % dmde,                          &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % dmde)
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
#   ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % dmde,                &
     &                        GRID(ng) % dmde)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dmde)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dmde)
#   endif
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%dndx).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        GRID(ng) % dndx,                          &
     &                        checksum = Fhash)
#   else
     &                        GRID(ng) % dndx)
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
#   ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % dndx,                &
     &                        GRID(ng) % dndx)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dndx)
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dndx)
#   endif
#  endif
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(GRID(ng)%xp).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % pmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % xp,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % xp)
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
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xp)
#  endif
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(GRID(ng)%yp).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % pmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % yp,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % yp)
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
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yp)
#  endif
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%xr).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % xr,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % xr)
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
#  ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xr,                &
     &                          GRID(ng) % xr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xr)
#  endif
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%yr).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % yr,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % yr)
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
#  ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yr,                &
     &                          GRID(ng) % yr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yr)
#  endif
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%xu).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % xu,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % xu)
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
#  ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xu,                &
     &                          GRID(ng) % xu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xu)
#  endif
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%yu).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % yu,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % yu)
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
#  ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yu,                &
     &                          GRID(ng) % yu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yu)
#  endif
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%xv).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % xv,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % xv)
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
#  ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xv,                &
     &                          GRID(ng) % xv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xv)
#  endif
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%yv).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % yv,                            &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % yv)
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
#  ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yv,                &
     &                          GRID(ng) % yv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yv)
#  endif
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=p2dvar
              IF (KIND(GRID(ng)%lonp).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % pmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % lonp,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % lonp)
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
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonp)
#  endif
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=p2dvar
              IF (KIND(GRID(ng)%latp).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % pmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % latp,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % latp)
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
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latp)
#  endif
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(GRID(ng)%lonr).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LonMin(ng), LonMax(ng),           &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % lonr,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % lonr)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of RHO-points: lon_rho', &
     &                              ng, TRIM(ncname),                   &
     &                              LonMin(ng), LonMax(ng)
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xr,                &
     &                          GRID(ng) % lonr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonr)
#  endif
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(GRID(ng)%latr).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LatMin(ng), LatMax(ng),           &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % latr,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % latr)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of RHO-points lat_rho',   &
     &                              ng, TRIM(ncname),                   &
     &                              LatMin(ng), LatMax(ng)
#  ifdef CHECKSUM
                  WRITE (stdout,60) Fhash
#  endif
                END IF
              END IF
#  ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yr,                &
     &                          GRID(ng) % latr)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latr)
#  endif
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=u2dvar
              IF (KIND(GRID(ng)%lonu).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % umask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % lonu,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % lonu)
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
#  ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xu,                &
     &                          GRID(ng) % lonu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonu)
#  endif
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=u2dvar
              IF (KIND(GRID(ng)%latu).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % umask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % latu,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % latu)
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
#  ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yu,                &
     &                          GRID(ng) % latu)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latu)
#  endif
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=v2dvar
              IF (KIND(GRID(ng)%lonv).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % lonv,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % lonv)
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
#  ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xv,                &
     &                          GRID(ng) % lonv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonv)
#  endif
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=v2dvar
              IF (KIND(GRID(ng)%latv).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#  endif
#  ifdef CHECKSUM
     &                          GRID(ng) % latv,                        &
     &                          checksum = Fhash)
#  else
     &                          GRID(ng) % latv)
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
#  ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yv,                &
     &                          GRID(ng) % latv)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latv)
#  endif
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%angler).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        GRID(ng) % angler,                        &
     &                        checksum = Fhash)
#  else
     &                        GRID(ng) % angler)
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
#   ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, 'angler', spval_check,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % angler,              &
     &                        GRID(ng) % angler)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % angler)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % angler)
#  endif
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
     defined SO_SEMI
#   ifndef OBS_SPACE
!
!  Read in adjoint sensitivity spatial scope masking at RHO-points.
!
          CASE ('scope_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%Rscope).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#    endif
#    ifdef CHECKSUM
     &                        GRID(ng) % Rscope,                        &
     &                        checksum = Fhash)
#    else
     &                        GRID(ng) % Rscope)
#    endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on RHO-points: scope_rho',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
              END IF
            END IF
            GotScope(1)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Rscope)
            END IF
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Rscope)
#    endif
!
!  Read in adjoint sensitivity spatial scope masking at U-points.
!
          CASE ('scope_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%Uscope).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask,                         &
#    endif
#    ifdef CHECKSUM
     &                        GRID(ng) % Uscope,                        &
     &                        checksum = Fhash)
#    else
     &                        GRID(ng) % Uscope)
#    endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on U-points: scope_u',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
              END IF
            END IF
            GotScope(2)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Uscope)
            END IF
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Uscope)
#    endif
!
!  Read in adjoint sensitivity spatial scope masking at V-points.
!
          CASE ('scope_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%Vscope).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#    endif
#    ifdef CHECKSUM
     &                        GRID(ng) % Vscope,                        &
     &                        checksum = Fhash)
#    else
     &                        GRID(ng) % Vscope)
#    endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on V-points: scope_v',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
              END IF
            END IF
            GotScope(3)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Vscope)
            END IF
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Vscope)
#    endif
#   endif
#  endif
        END SELECT
      END DO
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF

#  if defined UV_DRAG_GRID && !defined ANA_DRAG
#   ifdef UV_LOGDRAG
!
!  Read in spacially varying bottom roughness length (m).
!
      IF (KIND(GRID(ng)%ZoBot).eq.8) THEN
        pioVar_ZoBL%dkind=PIO_double
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        pioVar_ZoBL%dkind=PIO_real
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,             &
     &                  Vname(1,idZoBL), pioVar_ZoBL,                   &
     &                  0, ioDesc, Vsize,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
#    ifdef CHECKSUM
     &                  GRID(ng) % ZoBot,                               &
     &                  checksum = Fhash)
#    else
     &                  GRID(ng) % ZoBot)
#    endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZoBL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) 'time invariant, bottom roughness '//       &
     &                      'length scale: ZoBot',                      &
     &                      ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
          WRITE (stdout,60) Fhash
#    endif
        END IF
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % ZoBot)
      END IF
#    ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % ZoBot)
#    endif
#   endif
#   ifdef UV_LDRAG
!
!  Read in spacially varying linear drag coefficients (m/s).
!
      IF (KIND(GRID(ng)%rdrag).eq.8) THEN
        pioVar_dragL%dkind=PIO_double
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        pioVar_dragL%dkind=PIO_real
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,             &
     &                  Vname(1,idragL), pioVar_dragL,                  &
     &                  0, ioDesc, Vsize,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
#    ifdef CHECKSUM
     &                  GRID(ng) % rdrag,                               &
     &                  checksum = Fhash)
#    else
     &                  GRID(ng) % rdrag)
#    endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idragL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) 'linear bottom drag coefficient: rdrag',    &
     &                      ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
        END IF
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % rdrag)
      END IF
#    ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % rdrag)
#    endif
#   endif
#   ifdef UV_QDRAG
!
!  Read in spacially varying quadratic drag coefficients.
!
      IF (KIND(GRID(ng)%rdrag2).eq.8) THEN
        pioVar_dragQ%dkind=PIO_double
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        pioVar_dragQ%dkind=PIO_real
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,             &
     &                  Vname(1,idragQ), pioVar_dragQ,                  &
     &                  0, ioDesc, Vsize,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
#    ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#    endif
#    ifdef CHECKSUM
     &                  GRID(ng) % rdrag2,                              &
     &                  checksum = Fhash)
#    else
     &                  GRID(ng) % rdrag2)
#    endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idragQ)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) 'quadratic bottom drag coefficient: rdrag2',&
     &                      ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
        END IF
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % rdrag2)
      END IF
#    ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % rdrag2)
#    endif
#   endif
#  endif
!
! Close GRID NetCDF file.
!
      CALL pio_netcdf_close (ng, model, GRD(ng)%pioFile, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
     defined SO_SEMI
#   ifndef OBS_SPACE
!
!-----------------------------------------------------------------------
!  Inquire adjoint sensitivity forcing file. Read scope arrays again.
!  These fields take precedence
!-----------------------------------------------------------------------
!
      ncname=ADS(ng)%name
!
!  Open adjoint sensitivity NetCDF file for reading.
!
      IF (ADS(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, ADS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL pio_netcdf_check_dim (ng, model, ncname,                     &
     &                           pioFile = ADS(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, model, ncname,                       &
     &                         pioFile = ADS(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Check if the adjoint sensitivity scope arrays are available.
!
      GotScope(4)=find_string(var_name,n_var,'scope_rho',vindex)
      GotScope(5)=find_string(var_name,n_var,'scope_u',vindex)
      GotScope(6)=find_string(var_name,n_var,'scope_v',vindex)
!
      IF ((.not.GotScope(1)).and.(.not.GotScope(4))) THEN
        IF (Master) WRITE (stdout,20) 'scope_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((.not.GotScope(2)).and.(.not.GotScope(5))) THEN
        IF (Master) WRITE (stdout,20) 'scope_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((.not.GotScope(3)).and.(.not.GotScope(6))) THEN
        IF (Master) WRITE (stdout,20) 'scope_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (Master) THEN
        IF (GotScope(4)) THEN
          WRITE (stdout,50) TRIM(ADS(ng)%name)
        ELSE
          WRITE (stdout,50) TRIM(GRD(ng)%name)
        END IF
      END IF
!
!  Scan adjoint sensitivity variables.
!
      DO i=1,n_var

        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in adjoint sensitivity spatial scope masking at RHO-points.
!
          CASE ('scope_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%Rscope).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, ADS(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#    endif
#    ifdef CHECKSUM
     &                        GRID(ng) % Rscope,                        &
     &                        checksum = Fhash)
#    else
     &                        GRID(ng) % Rscope)
#    endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on RHO-points: scope_rho',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Rscope)
            END IF
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Rscope)
#    endif
!
!  Read in adjoint sensitivity spatial scope masking at U-points.
!
          CASE ('scope_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%Uscope).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, ADS(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask,                         &
#    endif
#    ifdef CHECKSUM
     &                        GRID(ng) % Uscope,                        &
     &                        checksum = Fhash)
#    else
     &                        GRID(ng) % Uscope)
#    endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on U-points: scope_u',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Uscope)
            END IF
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Uscope)
#    endif
!
!  Read in adjoint sensitivity spatial scope masking at V-points.
!
          CASE ('scope_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%Vscope).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, ADS(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#    endif
#    ifdef CHECKSUM
     &                        GRID(ng) % Vscope,                        &
     &                        checksum = Fhash)
#    else
     &                        GRID(ng) % Vscope)
#    endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'adjoint sensitivity spatial '//      &
     &                            'scope on V-points: scope_v',         &
     &                            ng, TRIM(ncname), Fmin, Fmax
#    ifdef CHECKSUM
                WRITE (stdout,60) Fhash
#    endif
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Vscope)
            END IF
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Vscope)
#    endif
        END SELECT
      END DO
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF
#   endif
#  endif
!
  10  FORMAT (/,' GET_GRID_PIO - unable to open grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID_PIO - unable to find grid variable: ',a,     &
     &        /,16x,'in grid NetCDF file: ',a)
  30  FORMAT (2x,'GET_GRID_PIO     - ',a,/,22x,                         &
     &        '(Grid = ',i2.2,', File: ',a,')',/,22x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  40  FORMAT (/,' GET_GRID_PIO - error while reading variable: ',a,     &
     &        /,12x,'in grid NetCDF file: ',a)
  50  FORMAT (/,2x,'GET_GRID_PIO     - Reading adjoint sensitivity',    &
     &        ' scope arrays from file:',/22x,a,/)
#  ifdef CHECKSUM
  60  FORMAT (22x,'(CheckSum = ',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE get_grid_pio
# endif
#endif
      END MODULE get_grid_mod
