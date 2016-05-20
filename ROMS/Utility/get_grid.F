#include "cppdefs.h"
#ifndef ANA_GRID
      SUBROUTINE get_grid (ng, model)
!
!svn $Id: get_grid.F 787 2008-10-14 00:37:50Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine reads grid information from GRID NetCDF file.       !
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
      USE mod_scalars
# ifdef WTYPE_GRID
      USE mod_mixing
# endif
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# ifdef NESTING
      USE nesting_mod, ONLY : fill_contact
# endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
      logical :: GotScope(6)
# endif
      integer :: tile, LBi, UBi, LBj, UBj
      integer :: cr, gtype, i, status, vindex
# if defined UV_DRAG_GRID && !defined ANA_DRAG
      integer :: varid_dragL, varid_dragQ, varid_ZoBL
# endif
#if defined UV_WAVEDRAG
      integer :: varid_dragW
#endif
      integer :: Vsize(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin

      character (len=256) :: ncname
!
      SourceFile='get_grid.F'
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=GRD(ng)%name
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
        IF (Master) WRITE (stdout,10) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'el',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'spherical',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'h',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# ifdef ICESHELF
      IF (.not.find_string(var_name,n_var,'zice',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'zice', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef WTYPE_GRID
      IF (.not.find_string(var_name,n_var,'wtype_grid',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'wtype_grid', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
      IF (.not.find_string(var_name,n_var,'f',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'f', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pm',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pn',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# if (defined CURVGRID && defined UV_ADV)
      IF (.not.find_string(var_name,n_var,'dndx',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'dndx', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dmde',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'dmde', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef CURVGRID
      IF (.not.find_string(var_name,n_var,'angle',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'angle', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef MASKING
      IF (.not.find_string(var_name,n_var,'mask_rho',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_u',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_v',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_psi',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# if defined WTYPE_GRID && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE) && \
    !defined ANA_WTYPE
      IF (.not.find_string(var_name,n_var,'wtype_grid',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'wtype_grid', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifndef ANA_SPONGE
      IF (LuvSponge(ng)) THEN
        IF (.not.find_string(var_name,n_var,'visc_factor',vindex)) THEN
          IF (Master) WRITE (stdout,10) 'visc_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#  ifdef SOLVE3D
      IF (ANY(LtracerSponge(:,ng))) THEN
        IF (.not.find_string(var_name,n_var,'diff_factor',vindex)) THEN
          IF (Master) WRITE (stdout,10) 'diff_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#  endif
# endif
# if defined UV_DRAG_GRID && !defined ANA_DRAG
#  ifdef UV_LOGDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idZoBL)),        &
     &                     varid_ZoBL)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idZoBL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef UV_LDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idragL)),        &
     &                     varid_dragL)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef UV_QDRAG
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idragQ)),        &
     &                     varid_dragQ)) THEN
        IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragQ)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
# endif
!
!  Open grid NetCDF file for reading.
!
      IF (GRD(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, GRD(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,20) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL netcdf_get_lvar (ng, model, ncname, 'spherical',           &
     &                        spherical,                                &
     &                        ncid = GRD(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid variables.
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
!  Scan the variable list and read in needed variables.
!
      DO i=1,n_var

        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in basin X-length.
!
          CASE ('xl')
            CALL netcdf_get_fvar (ng, model, ncname, 'xl',              &
     &                            xl(ng),                               &
     &                            ncid = GRD(ng)%ncid)
            IF (exit_flag.ne.NoError) EXIT
!
!  Read in basin Y-length.
!
          CASE ('el')
            CALL netcdf_get_fvar (ng, model, ncname, 'el',              &
     &                            el(ng),                               &
     &                            ncid = GRD(ng)%ncid)
            IF (exit_flag.ne.NoError) EXIT
!
!  Read in bathymetry.
!
          CASE ('h')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, hmin(ng), hmax(ng),                 &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        GRID(ng) % h)
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) 'h', TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % h,                   &
     &                        GRID(ng) % h)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % rmask)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, 'rmask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % rmask,               &
     &                        GRID(ng) % rmask)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % umask)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        u2dvar, 'umask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % umask,               &
     &                        GRID(ng) % umask)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % vmask)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        v2dvar, 'vmask', spval_check,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % vmask,               &
     &                        GRID(ng) % vmask)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % pmask)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
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
!  Read in ice shelf thicknesses.
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
     &                        GRID(ng) % zice)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
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
     &                        var_name(i), var_id(i), 0,                &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        MIXING(ng) % Jwtype)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
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
     &                          var_name(i), var_id(i), 0,              &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          MIXING(ng) % visc_factor)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
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
     &                          var_name(i), var_id(i), 0,              &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
     &                          MIXING(ng) % diff_factor)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
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
     &                        GRID(ng) % f)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % f,                   &
     &                        GRID(ng) % f)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % pm)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % pm,                  &
     &                        GRID(ng) % pm)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % pn)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % pn,                  &
     &                        GRID(ng) % pn)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % dmde)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % dmde,                &
     &                        GRID(ng) % dmde)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % dndx)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, var_name(i), spval_check,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % dndx,                &
     &                        GRID(ng) % dndx)
            IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % xp)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
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
     &                        GRID(ng) % yp)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
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
     &                        GRID(ng) % xr)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xr,                &
     &                          GRID(ng) % xr)
              IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % yr)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yr,                &
     &                          GRID(ng) % yr)
              IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % xu)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xu,                &
     &                          GRID(ng) % xu)
              IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % yu)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yu,                &
     &                          GRID(ng) % yu)
              IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % xv)
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xv,                &
     &                          GRID(ng) % xv)
              IF (exit_flag.ne.NoError) RETURN
            END IF
# endif
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
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
     &                        GRID(ng) % yv)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef NESTING
            IF (.not.spherical) THEN
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yv,                &
     &                          GRID(ng) % yv)
              IF (exit_flag.ne.NoError) RETURN
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
     &                          GRID(ng) % lonp)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
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
     &                          GRID(ng) % latp)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
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
     &                          GRID(ng) % lonr)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xr,                &
     &                          GRID(ng) % lonr)
              IF (exit_flag.ne.NoError) RETURN
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
     &                          GRID(ng) % latr)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          r2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yr,                &
     &                          GRID(ng) % latr)
              IF (exit_flag.ne.NoError) RETURN
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
     &                          GRID(ng) % lonu)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xu,                &
     &                          GRID(ng) % lonu)
              IF (exit_flag.ne.NoError) RETURN
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
     &                          GRID(ng) % latu)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          u2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yu,                &
     &                          GRID(ng) % latu)
              IF (exit_flag.ne.NoError) RETURN
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
     &                          GRID(ng) % lonv)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Xv,                &
     &                          GRID(ng) % lonv)
              IF (exit_flag.ne.NoError) RETURN
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
     &                          GRID(ng) % latv)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef NESTING
              CALL fill_contact(ng, model, tile,                        &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          v2dvar, var_name(i), spval_check,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          CONTACT_METRIC(cr) % Yv,                &
     &                          GRID(ng) % latv)
              IF (exit_flag.ne.NoError) RETURN
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
     &                        GRID(ng) % angler)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  ifdef NESTING
            CALL fill_contact(ng, model, tile,                          &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        r2dvar, 'angler', spval_check,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CONTACT_METRIC(cr) % angler,              &
     &                        GRID(ng) % angler)
            IF (exit_flag.ne.NoError) RETURN
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
# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
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
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        GRID(ng) % Rscope)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            GotScope(1)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Rscope)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Rscope)
#  endif
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
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        GRID(ng) % Uscope)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            GotScope(2)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Uscope)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Uscope)
#  endif
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
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        GRID(ng) % Vscope)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            GotScope(3)=.TRUE.
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Vscope)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Vscope)
#  endif
# endif
        END SELECT
      END DO
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,30) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF

# if defined UV_DRAG_GRID && !defined ANA_DRAG
#  ifdef UV_LOGDRAG
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
     &                  GRID(ng) % ZoBot)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idZoBL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
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
     &                  GRID(ng) % rdrag)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idragL)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
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
     &                  GRID(ng) % rdrag2)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idragQ)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
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
          IF (status.ne.nf90_noerr) THEN
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
      IF (exit_flag.ne.NoError) RETURN

# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Inquire adjoint sensitivity forcing file. Read scope arrays again.
!  These fields take precedence
!-----------------------------------------------------------------------
!
      ncname=ADS(ng)%name
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
!  Check if the adjoint sensitivity scope arrays are available.
!
      GotScope(4)=find_string(var_name,n_var,'scope_rho',vindex)
      GotScope(5)=find_string(var_name,n_var,'scope_u',vindex)
      GotScope(6)=find_string(var_name,n_var,'scope_v',vindex)
!
      IF ((.not.GotScope(1)).and.(.not.GotScope(4))) THEN
        IF (Master) WRITE (stdout,10) 'scope_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((.not.GotScope(2)).and.(.not.GotScope(5))) THEN
        IF (Master) WRITE (stdout,10) 'scope_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((.not.GotScope(3)).and.(.not.GotScope(6))) THEN
        IF (Master) WRITE (stdout,10) 'scope_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (Master) THEN
        IF (GotScope(4)) THEN
          WRITE (stdout,40) TRIM(ADS(ng)%name)
        ELSE
          WRITE (stdout,40) TRIM(GRD(ng)%name)
        END IF
      END IF
!
!  Open adjoint sensitivity NetCDF file for reading.
!
      IF (GotScope(4).and.(ADS(ng)%ncid.eq.-1)) THEN
        CALL netcdf_open (ng, model, ncname, 0, ADS(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,20) TRIM(ncname)
          RETURN
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
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        GRID(ng) % Rscope)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Rscope)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Rscope)
#  endif
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
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        GRID(ng) % Uscope)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Uscope)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Uscope)
#  endif
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
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        GRID(ng) % Vscope)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % Vscope)
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % Vscope)
#  endif
        END SELECT
      END DO
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,30) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF
# endif
!
  10  FORMAT (/,' GET_GRID - unable to find grid variable: ',a,         &
     &        /,12x,'in grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID - unable to open grid NetCDF file: ',a)
  30  FORMAT (/,' GET_GRID - error while reading variable: ',a,         &
     &        /,12x,'in grid NetCDF file: ',a)
  40  FORMAT (/,' GET_GRID - Reading adjoint sensitivity scope arrays', &
     &        ' from file:',/12x,a)
#else
      SUBROUTINE get_grid
#endif
      RETURN
      END SUBROUTINE get_grid
