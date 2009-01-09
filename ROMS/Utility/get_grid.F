#include "cppdefs.h"
#ifndef ANA_GRID
      SUBROUTINE get_grid (ng, model)
!
!svn $Id: get_grid.F 787 2008-10-14 00:37:50Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
      logical :: GotScope(6)
# endif
      logical :: find_string

      integer :: tile, LBi, UBi, LBj, UBj
      integer :: gtype, i, status, vindex

      integer :: Vsize(4)

      integer :: nf_fread2d

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin

      character (len=1 ) :: char1
      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=GRDname(ng)
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
        IF (Master) WRITE (stdout,10) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
!
!  Open grid NetCDF file for reading.
!
      IF (InpThread.and.(ncGRDid(ng).eq.-1)) THEN
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncGRDid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=4
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag, 1)
# endif
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      DO i=1,n_var
        IF (TRIM(var_name(i)).eq.'spherical') THEN
          IF (Inpthread) THEN
            status=nf90_get_var(ncGRDid(ng), var_id(i), char1)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
            ELSE
              IF ((char1.eq.'t').or.(char1.eq.'T')) THEN
                spherical=.TRUE.
              END IF
            END IF
          END IF
          EXIT
        END IF
      END DO
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag, 1)
# endif
      IF (exit_flag.ne.NoError) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'spherical', TRIM(ncname)
        END IF
        RETURN
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastl (ng, model, spherical, 1)
# endif
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
            IF (InpThread) THEN
              status=nf90_get_var(ncGRDid(ng), var_id(i), xl(ng))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasti (ng, model, exit_flag, 1)
# endif
            IF (exit_flag.ne.NoError) EXIT
# ifdef DISTRIBUTE
            CALL mp_bcastf (ng, model, xl(ng), 1)
# endif
!
!  Read in basin Y-length.
!
          CASE ('el')
            IF (InpThread) THEN
              status=nf90_get_var(ncGRDid(ng), var_id(i), el(ng))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
              END IF
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasti (ng, model, exit_flag, 1)
# endif
            IF (exit_flag.ne.NoError) EXIT
# ifdef DISTRIBUTE
            CALL mp_bcastf (ng, model, el(ng), 1)
# endif
!
!  Read in bathymetry.
!
          CASE ('h')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, hmin(ng), hmax(ng),                 &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % h(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % h)
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % h)
# endif
# ifdef MASKING
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask(LBi,LBj),                &
     &                        GRID(ng) % rmask(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % rmask)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % rmask)
#  endif
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask(LBi,LBj),                &
     &                        GRID(ng) % umask(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % umask)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % umask)
#  endif
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask(LBi,LBj),                &
     &                        GRID(ng) % vmask(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % vmask)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % vmask)
#  endif
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask(LBi,LBj),                &
     &                        GRID(ng) % pmask(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                             GRID(ng) % pmask)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % pmask)
#  endif
# endif
# ifdef ICESHELF
!
!  Read in ice shelf thicknesses.
!
          CASE ('zice')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % zice(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % zice)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % zice)
#  endif
# endif
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % f(LBi,LBj))
#undef SOLITON_REVERSE
#ifdef SOLITON_REVERSE
             GRID(ng) % f(:,:)=-GRID(ng) % f(:,:)
#endif
#undef SOLITON_REVERSE
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % f)
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % f)
# endif
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % pm(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % pm)
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % pm)
# endif
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % pn(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % pn)
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % pn)
# endif
# if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % dmde(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % dmde)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % dmde)
#  endif
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % dndx(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % dndx)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % dndx)
#  endif
# endif
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % pmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % xp(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % xp)
# endif
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,    &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % pmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % yp(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % yp)
# endif
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % xr(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % xr)
# endif
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,    &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % yr(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % yr)
# endif
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
# endif
     &                        GRID(ng) % xu(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % xu)
# endif
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,    &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
# endif
     &                        GRID(ng) % yu(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % yu)
# endif
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % xv(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % xv)
# endif
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,    &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % yv(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, .FALSE., .FALSE.,         &
     &                          GRID(ng) % yv)
# endif
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,   &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % pmask(LBi,LBj),              &
# endif
     &                          GRID(ng) % lonp(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % lonp)
# endif
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,  &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % pmask(LBi,LBj),              &
# endif
     &                          GRID(ng) % latp(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % latp)
# endif
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,   &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, LonMin(ng), LonMax(ng),           &
# ifdef MASKING
     &                          GRID(ng) % rmask(LBi,LBj),              &
# endif
     &                          GRID(ng) % lonr(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % lonr)
# endif
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,  &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, LatMin(ng), LatMax(ng),           &
# ifdef MASKING
     &                          GRID(ng) % rmask(LBi,LBj),              &
# endif
     &                          GRID(ng) % latr(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % latr)
# endif
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,   &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask(LBi,LBj),              &
# endif
     &                          GRID(ng) % lonu(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % lonu)
# endif
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,  &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask(LBi,LBj),              &
# endif
     &                          GRID(ng) % latu(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % latu)
# endif
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,   &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask(LBi,LBj),              &
# endif
     &                          GRID(ng) % lonv(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % lonv)
# endif
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model,  ncGRDid(ng), var_id(i), 0,  &
     &                          gtype, Vsize, LBi, UBi, LBj, UBj,       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask(LBi,LBj),              &
# endif
     &                          GRID(ng) % latv(LBi,LBj))
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
# ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, .FALSE., .FALSE.,       &
     &                            GRID(ng) % latv)
# endif
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        GRID(ng) % angler(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % angler)
# endif
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % angler)
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
!
!  Read in adjoint sensitivity spatial scope masking at RHO-points.
!
          CASE ('scope_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % Rscope(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            GotScope(1)=.TRUE.
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % Rscope)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % Rscope)
#  endif
!
!  Read in adjoint sensitivity spatial scope masking at U-points.
!
          CASE ('scope_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % Uscope(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            GotScope(2)=.TRUE.
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % Uscope)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % Uscope)
#  endif
!
!  Read in adjoint sensitivity spatial scope masking at V-points.
!
          CASE ('scope_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncGRDid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % Vscope(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            GotScope(3)=.TRUE.
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % Vscope)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % Vscope)
#  endif
# endif
        END SELECT
      END DO
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,30) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF
!
! Close GRID NetCDF file.
!
      IF (InpThread) THEN
        status=nf90_close(ncGRDid(ng))
        ncGRDid(ng)=-1
      END IF

# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Inquire adjoint sensitivity forcing file. Read scope arrays again.
!  These fields take precedence
!-----------------------------------------------------------------------
!
      ncname=ADSname(ng)
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
          WRITE (stdout,40) TRIM(ADSname(ng))
        ELSE
          WRITE (stdout,40) TRIM(GRDname(ng))
        END IF
      END IF
!
!  Open adjoint sensitivity NetCDF file for reading.
!
      IF (InpThread.and.GotScope(4).and.(ncADSid(ng).eq.-1)) THEN
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncADSid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=4
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag, 1)
# endif
      IF (exit_flag.ne.NoError) RETURN
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
            status=nf_fread2d(ng, model, ncADSid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % Rscope(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % Rscope)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % Rscope)
#  endif
!
!  Read in adjoint sensitivity spatial scope masking at U-points.
!
          CASE ('scope_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncADSid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % Uscope(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % Uscope)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          GRID(ng) % Uscope)
#  endif
!
!  Read in adjoint sensitivity spatial scope masking at V-points.
!
          CASE ('scope_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncADSid(ng), var_id(i), 0,     &
     &                        gtype, Vsize, LBi, UBi, LBj, UBj,         &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#  endif
     &                        GRID(ng) % Vscope(LBi,LBj))
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % Vscope)
#  endif
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
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
 40   FORMAT (/,' GET_GRID - Reading adjoint sensitivity scope arrays', &
     &        ' from file:',/12x,a)
#else
      SUBROUTINE get_grid
#endif
      RETURN
      END SUBROUTINE get_grid
