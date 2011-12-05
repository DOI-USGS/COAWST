#include "cppdefs.h"
      SUBROUTINE wrt_rst (ng)
!
!svn $Id: wrt_rst.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes fields into restart NetCDF file.                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
#ifdef ICE_MODEL
      USE mod_ice
#endif
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
#ifdef NCEP_FLUXES
      USE mod_forces
#endif
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# if defined PERFECT_RESTART || defined SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
# if defined PERFECT_RESTART && defined SOLVE3D
      USE nf_fwrite4d_mod, ONLY : nf_fwrite4d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status, varid
# if defined PERFECT_RESTART || defined SOLVE3D
      integer :: ntmp(1)
# endif
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile='wrt_rst.F'
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
      gfactor=-1
#else
      gfactor=1
#endif
!
!  Set time record index.
!
      tRSTindx(ng)=tRSTindx(ng)+1
      NrecRST(ng)=NrecRST(ng)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (LcycleRST(ng)) THEN
        tRSTindx(ng)=MOD(tRSTindx(ng)-1,2)+1
      END IF

#ifdef PERFECT_RESTART
!
!  Write out time-stepping indices.
!
# ifdef SOLVE3D
      ntmp(1)=1+MOD((iic(ng)-1)-ntstart(ng),2)
      CALL netcdf_put_ivar (ng, iNLM, RSTname(ng), 'nstp',              &
     &                      ntmp, (/tRSTindx(ng)/), (/1/),              &
     &                      ncid = ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, iNLM, RSTname(ng), 'nrhs',              &
     &                      ntmp, (/tRSTindx(ng)/), (/1/),              &
     &                      ncid = ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN

      ntmp(1)=3-ntmp(1)
      CALL netcdf_put_ivar (ng, iNLM, RSTname(ng), 'nnew',              &
     &                      ntmp, (/tRSTindx(ng)/), (/1/),              &
     &                      ncid = ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
# endif
      CALL netcdf_put_ivar (ng, iNLM, RSTname(ng), 'kstp',              &
     &                      kstp(ng:), (/tRSTindx(ng)/), (/1/),         &
     &                      ncid = ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, iNLM, RSTname(ng), 'krhs',              &
     &                      krhs(ng:), (/tRSTindx(ng)/), (/1/),         &
     &                      ncid = ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, iNLM, RSTname(ng), 'knew',              &
     &                      knew(ng:), (/tRSTindx(ng)/), (/1/),         &
     &                      ncid = ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, RSTname(ng),                      &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/tRSTindx(ng)/), (/1/),                    &
     &                      ncid = ncRSTid(ng),                         &
     &                      varid = rstVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN

#if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-dependent bathymetry (m)
!
      IF (Hout(idbath,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idbath,ng),    &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng) % h,                                &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idRwet,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRwet)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idUwet,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUwet)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idVwet,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVwet)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
!
!  Write out free-surface (m).
!
      scale=1.0_r8
#ifdef PERFECT_RESTART
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idFsur,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 3, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef WET_DRY
     &                   OCEAN(ng) % zeta,                              &
     &                   SetFillVal = .FALSE.)
# else
     &                   OCEAN(ng) % zeta)
# endif
#else
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idFsur,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef WET_DRY
     &                   OCEAN(ng) % zeta(:,:,KOUT),                    &
     &                   SetFillVal = .FALSE.)
# else
     &                   OCEAN(ng) % zeta(:,:,KOUT))
# endif
#endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#ifdef PERFECT_RESTART
!
!  Write out RHS of free-surface equation.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idRzet,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 2, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng) % rzeta)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRzet)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
#ifdef PERFECT_RESTART
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idUbar,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 3, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   OCEAN(ng) % ubar)
#else
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idUbar,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   OCEAN(ng) % ubar(:,:,KOUT))
#endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#ifdef PERFECT_RESTART
!
!  Write out RHS of 2D momentum equation in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idRu2d,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 2, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   OCEAN(ng) % rubar)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRu2d)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
#ifdef PERFECT_RESTART
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVbar,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 3, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   OCEAN(ng) % vbar)
#else
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idVbar,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   OCEAN(ng) % vbar(:,:,KOUT))
#endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#ifdef PERFECT_RESTART
!
!  Write out RHS of 2D momentum equation in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idRv2d,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 2, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   OCEAN(ng) % rvbar)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRv2d)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
#ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
# ifdef PERFECT_RESTART
      status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstVid(idUvel,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 2, scale,     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   OCEAN(ng) % u)
# else
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idUvel,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   OCEAN(ng) % u(:,:,:,NOUT))
# endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef PERFECT_RESTART
!
!  Write out RHS of 3D momentum equation in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstVid(idRu3d,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   OCEAN(ng) % ru)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRu3d)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
# ifdef PERFECT_RESTART
      status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstVid(idVvel,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 2, scale,     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   OCEAN(ng) % v)
# else
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVvel,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   OCEAN(ng) % v(:,:,:,NOUT))
# endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef PERFECT_RESTART
!
!  Write out RHS of 3D momentum equation in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstVid(idRv3d,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   OCEAN(ng) % rv)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRv3d)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
# ifdef PERFECT_RESTART
        status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstTid(itrc,ng),      &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 1, 2, scale,   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % t(:,:,:,:,itrc))
# else
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstTid(itrc,ng),      &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % t(:,:,:,NOUT,itrc))
# endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idDano,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng) % rho)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idHsbl,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % hsbl)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idHbbl,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % hbbl)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined BVF_MIXING  || defined GLS_MIXING || \
     defined MY25_MIXING || defined LMD_MIXING
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVvis,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Akv)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idTdif,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Akt(:,:,:,itemp))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idSdif,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akt(:,:,:,isalt))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
# endif
# ifdef ICE_MODEL
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idUice,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % umask,                                &
#endif
     &                 ICE (ng) % ui(:,:,IUOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUice)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idvice,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#endif
     &                 ICE (ng) % vi(:,:,IUOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVice)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice concentration
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idAice,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % ai(:,:,IOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idAice)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice average thickness
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idHice,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % hi(:,:,IOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHice)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out snow average thickness
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idHsno,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % hsn(:,:,IOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHsno)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface water thickness (on ice)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idSfwat,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % sfwat(:,:,IOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSfwat)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice age
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idAgeice,ng),    &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % ageice(:,:,IOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idAgeice)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice/snow surface temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idTice,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % tis)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTice)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice interior temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idTimid,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % ti(:,:,IOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTimid)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out internal ice stress component 11
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idSig11,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % sig11(:,:,IEOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig11)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out internal ice stress component 22
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idSig22,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % sig22(:,:,IEOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig22)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out internal ice stress component 12
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idSig12,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % sig12(:,:,IEOUT))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig12)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-ocean friction velocity
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idTauiw,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % utau_iw)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTauiw)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-ocean momentum transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idChuiw,ng),     &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % chu_iw)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idChuiw)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out salinity in molecular sub-layer under ice
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idS0mk,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % s0mk)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idS0mk)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out temperature in molecular sub-layer under ice
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idT0mk,ng),      &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#endif
     &                 ICE (ng) % t0mk)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idT0mk)), trstindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
!
!  Write out turbulent kinetic energy.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstVid(idMtke,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % tke)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtke)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define turbulent kinetic energy time length scale.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite4d(ng, iNLM, ncRSTid(ng), rstVid(idMtls,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % gls)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtls)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define vertical mixing turbulent length scale.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVmLS,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Lscale)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVmLS)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define turbulent kinetic energy vertical diffusion coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVmKK,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Akk)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVmKK)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  ifdef GLS_MIXING
!
!  Define turbulent length scale vertical diffusion coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVmKP,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Akp)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVmKP)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
# endif
# ifdef NCEP_FLUXES
!
!  Write out NCEP gustiness squared
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idWg2d,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % wg2_d)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idWg2d)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out NCEP momentum transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idCdd,ng),       &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % cd_d)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCdd)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out NCEP sensible heat transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idChd,ng),       &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % ch_d)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idChd)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out NCEP latent heat transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idCed,ng),       &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % ce_d)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCed)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out model gustiness squared
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idWg2m,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % wg2_m)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idWg2m)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out model momentum transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idCdm,ng),       &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % cd_m)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCdm)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out model sensible heat transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idChm,ng),       &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % ch_m)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idChm)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out model latent heat transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idCem,ng),       &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % ce_m)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCem)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out NCEP near-surface air density
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idRhoa,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FORCES(ng) % rhoa_n)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRhoa)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idUbld(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     SEDBED(ng) % bedldu(:,:,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out bed load transport in V-direction.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idVbld(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     SEDBED(ng) % bedldv(:,:,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idfrac(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % bed_frac(:,:,:,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idBmas(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % bed_mass(:,:,:,NOUT,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (i.eq.itauc) THEN
          scale=rho0
        ELSE
          scale=1.0_r8
        END IF
        gtype=gfactor*b3dvar
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idSbed(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % bed(:,:,:,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties. Notice that only the
!  first four properties (mean grain diameter, mean grain density,
!  mean settling velocity, mean critical erosion stress,
!  ripple length and ripple height) are written.
!
      DO i=1,6
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idBott(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBED(ng) % bottom(:,:,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tRSTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
#endif
#ifdef WEC
!
!  Write out 2D U-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idU2Sd,ng),    &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   OCEAN(ng) % ubar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),          &
     &                                  tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 2D V-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idV2Sd,ng),    &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   OCEAN(ng) % vbar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),          &
     &                                  tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# ifdef SOLVE3D
!
!  Write out 3D U-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idU3Sd,ng),    &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   OCEAN(ng) % u_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF(Master) WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),           &
     &                                 tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 3D V-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idV3Sd,ng),    &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   OCEAN(ng) % v_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),          &
     &                                  tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, RSTname(ng), ncRSTid(ng))
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      IF (Master) WRITE (stdout,20) KOUT, NOUT, tRSTindx(ng)
#else
      IF (Master) WRITE (stdout,20) KOUT, tRSTindx(ng)
#endif
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
#ifdef SOLVE3D
  20  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
#else
  20  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
#endif

      RETURN
      END SUBROUTINE wrt_rst
