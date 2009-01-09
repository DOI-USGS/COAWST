#include "cppdefs.h"
      SUBROUTINE wrt_rst (ng)
!
!svn $Id: wrt_rst.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
      USE mod_stepping
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
      integer :: ntmp(1), start(4), total(4)

      integer :: nf_fwrite2d
# if defined PERFECT_RESTART || defined SOLVE3D
      integer :: nf_fwrite3d
# endif
# if defined PERFECT_RESTART && defined SOLVE3D
      integer :: nf_fwrite4d
# endif
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
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
#if defined WRITE_WATER && defined MASKING
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
      start(1)=tRSTindx(ng)
      total(1)=1

#ifdef PERFECT_RESTART
!
!  Write out time-stepping indices.
!
# ifdef SOLVE3D
      ntmp(1)=1+MOD((iic(ng)-1)-ntstart(ng),2)
      status=nf90_inq_varid(ncRSTid(ng), 'nstp', varid)
      status=nf90_put_var(ncRSTid(ng), varid, ntmp, start, total)

      status=nf90_inq_varid(ncRSTid(ng), 'nrhs', varid)
      status=nf90_put_var(ncRSTid(ng), varid, ntmp, start, total)

      ntmp(1)=3-ntmp(1)
      status=nf90_inq_varid(ncRSTid(ng), 'nnew', varid)
      status=nf90_put_var(ncRSTid(ng), varid, ntmp, start, total)
# endif
      status=nf90_inq_varid(ncRSTid(ng), 'kstp', varid)
      status=nf90_put_var(ncRSTid(ng), varid, kstp(ng:), start, total)

      status=nf90_inq_varid(ncRSTid(ng), 'krhs', varid)
      status=nf90_put_var(ncRSTid(ng), varid, krhs(ng:), start, total)

      status=nf90_inq_varid(ncRSTid(ng), 'knew', varid)
      status=nf90_put_var(ncRSTid(ng), varid, knew(ng:), start, total)
#endif
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        status=nf90_put_var(ncRSTid(ng), rstVid(idtime,ng),             &
     &                      time(ng:), start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
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
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     GRID(ng) % h(LBi,LBj))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % zeta(LBi,LBj,1))
#else
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idFsur,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % zeta(LBi,LBj,KOUT))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rzeta(LBi,LBj,1))
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
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % ubar(LBi,LBj,1))
#else
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idUbar,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % ubar(LBi,LBj,KOUT))
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
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rubar(LBi,LBj,1))
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
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % vbar(LBi,LBj,1))
#else
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idVbar,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % vbar(LBi,LBj,KOUT))
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
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rvbar(LBi,LBj,1))
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
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % u(LBi,LBj,1,1))
# else
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idUvel,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % u(LBi,LBj,1,NOUT))
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
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % ru(LBi,LBj,0,1))
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
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % v(LBi,LBj,1,1))
# else
      status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstVid(idVvel,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % v(LBi,LBj,1,NOUT))
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
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % rv(LBi,LBj,0,1))
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
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % t(LBi,LBj,1,1,itrc))
# else
        status=nf_fwrite3d(ng, iNLM, ncRSTid(ng), rstTid(itrc,ng),      &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % t(LBi,LBj,1,NOUT,itrc))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rho(LBi,LBj,1))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % hsbl(LBi,LBj))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % hbbl(LBi,LBj))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Akv(LBi,LBj,0))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Akt(LBi,LBj,0,itemp))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   MIXING(ng) % Akt(LBi,LBj,0,isalt))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % tke(LBi,LBj,0,1))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % gls(LBi,LBj,0,1))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Lscale(LBi,LBj,0))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Akk(LBi,LBj,0))
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
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Akp(LBi,LBj,0))
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
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % bedldu(LBi,LBj,i))
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
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % bedldv(LBi,LBj,i))
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
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bed_frac(LBi,LBj,1,i))
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
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bed_mass(LBi,LBj,1,NOUT,i))
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
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bed(LBi,LBj,1,i))
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
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bottom(LBi,LBj,i))
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
#ifdef WET_DRY
!
!  Write out wet mask at RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idWmsk,ng),      &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   GRID(ng) % rmask_wet(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idWmsk)), tRSTindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
#ifdef NEARSHORE_MELLOR
!
!  Write out 2D U-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncRSTid(ng), rstVid(idU2Sd,ng),    &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % ubar_stokes(LBi,LBj))
        IF (OutThread.and.(status.ne.nf90_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),tRSTindx(ng)
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
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % vbar_stokes(LBi,LBj))
        IF (OutThread.and.(status.ne.nf90_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),tRSTindx(ng)
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
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % u_stokes(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf90_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),tRSTindx(ng)
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
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % v_stokes(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf90_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),tRSTindx(ng)
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
      IF (OutThread) THEN
        status=nf90_sync(ncRSTid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef SOLVE3D
        WRITE (stdout,30) KOUT, NOUT, tRSTindx
#else
        WRITE (stdout,30) KOUT, tRSTindx
#endif
      END IF
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_RST - unable to synchronize restart NetCDF to ',  &
     &        'disk.')
#ifdef SOLVE3D
  30  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
#else
  30  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
#endif
      RETURN
      END SUBROUTINE wrt_rst
