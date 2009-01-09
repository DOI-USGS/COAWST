#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE wrt_avg (ng)
!
!svn $Id: wrt_avg.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged fields into averages     !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
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
      integer :: gfactor, gtype, i, itrc, status
      integer :: start(4), total(4)

      integer :: nf_fwrite2d
# ifdef SOLVE3D
      integer :: nf_fwrite3d
# endif
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
        gfactor=-1
# else
        gfactor=1
# endif
!
!  Set time record index.
!
      tAVGindx(ng)=tAVGindx(ng)+1
      NrecAVG(ng)=NrecAVG(ng)+1
!
!  Write out averaged time.
!
      IF (OutThread) THEN
        start(1)=tAVGindx(ng)
        total(1)=1
        status=nf90_put_var(ncAVGid(ng), avgVid(idtime,ng),             &
     &                      AVGtime(ng:), start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idFsur,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     AVERAGE(ng) % avgzeta(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided free-surface (m).
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idFsuD,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   TIDES(ng) % zeta_detided(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsuD)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUbar,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
# endif
     &                     AVERAGE(ng) % avgu2d(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idu2dD,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   TIDES(ng) % ubar_detided(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu2dD)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVbar,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
# endif
     &                     AVERAGE(ng) % avgv2d(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idv2dD,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   TIDES(ng) % vbar_detided(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idv2dD)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Write out time-averaged mass fluxes for 3D momentum coupling.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUfx1,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgDU_avg1(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx1)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUfx2,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgDU_avg2(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx2)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVfx1,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgDV_avg1(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx1)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVfx2,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgDV_avg2(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx2)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef AVERAGES_NEARSHORE
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Hout(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idW2xx,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgSxx2d(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Hout(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idW2xy,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgSxy2d(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Hout(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idW2yy,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgSyy2d(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the XI-direction.
!
      IF (Hout(idU2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idU2rs,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgu2RS(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the ETA-direction
!
      IF (Hout(idV2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idV2rs,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgv2RS(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the XI-direction.
!
      IF (Hout(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idU2Sd,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgu2Sd(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idV2Sd,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgv2Sd(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idUvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgu3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idu3dD,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#   endif
     &                   TIDES(ng) % u_detided(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu3dD)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idVvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgv3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idv3dD,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#   endif
     &                   TIDES(ng) % v_detided(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idv3dD)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_NEARSHORE
!
!  Write out 3D radiation stress, Sxx-component.
!
      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idW3xx,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgSxx3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xx)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Sxy-component.
!
      IF (Hout(idW3xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idW3xy,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgSxy3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xy)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Syy-component.
!
      IF (Hout(idW3yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idW3yy,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgSyy3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3yy)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szx-component.
!
      IF (Hout(idW3zx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idW3zx,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgSzx3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zx)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szy-component.
!
      IF (Hout(idW3zy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idW3zy,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgSzy3d(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zy)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress in the XI-direction, u_Rstress.
!
      IF (Hout(idU3rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idU3rs,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgu3RS(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3rs)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress in the ETA-direction, v_Rstress.
!
      IF (Hout(idV3rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idV3rs,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgv3RS(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3rs)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D stokes momentum component (m/s) in the XI-direction.
!
      IF (Hout(idU3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idU3Sd,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgu3Sd(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3Sd)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idV3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idV3Sd,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgv3Sd(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3Sd)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idOvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgw3d(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idWvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgwvel(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgTid(itrc,ng),    &
     &                       tAVGindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#  endif
     &                       AVERAGE(ng) % avgt(LBi,LBj,1,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                        tAVGindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idDano,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgrho(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idHsbl,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avghsbl(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idHbbl,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avghbbl(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
!
!  Write out quadratic <zeta*zeta> term.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idZZav,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgZZ(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idZZav)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idU2av,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgU2(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idU2av)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idV2av,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgV2(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idV2av)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
!
!  Write u-volume flux.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idHUav,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgHuon(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHUav)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write v-volume flux.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idHVav,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgHvom(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHVav)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <u*u> term.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idUUav,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgUU(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUUav)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <u*v> term.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idUVav,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgUV(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUVav)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <v*v> term.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idVVav,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgVV(LBi,LBj,1))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVVav)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO i=1,NAT
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idTTav(i),ng), &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgTT(LBi,LBj,1,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTTav(i))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out active tracer volume fluxes.
!
      DO i=1,NAT
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(iHUTav(i),ng), &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgHuonT(LBi,LBj,1,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,iHUTav(i))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(iHVTav(i),ng), &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgHvomT(LBi,LBj,1,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,iHVTav(i))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO i=1,NAT
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idUTav(i),ng), &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgUT(LBi,LBj,1,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUTav(i))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idVTav(i),ng), &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgVT(LBi,LBj,1,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVTav(i))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
#  ifdef AVERAGES_AKV
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idVvis,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgAKv(LBi,LBj,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKT
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idTdif,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgAKt(LBi,LBj,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKS
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idSdif,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgAKs(LBi,LBj,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Write out surface net heat flux.
!
      IF (Hout(idTsur(itemp),ng)) THEN
#   ifdef ADJOINT
#    if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
        defined OPT_OBSERVATIONS
        scale=1.0_r8/(rho0*Cp)
#    else
        scale=1.0_r8
#    endif
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng),                       &
     &                     avgVid(idTsur(itemp),ng),                    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgstf(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Hout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng),                       &
     &                     avgVid(idTsur(isalt),ng),                    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     AVERAGE(ng) % avgswf(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idLhea,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     AVERAGE(ng) % avglhf(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idShea,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     AVERAGE(ng) % avgshf(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idLrad,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     AVERAGE(ng) % avglrf(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idevap,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#     endif
     &                     AVERAGE(ng) % avgevap(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idrain,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#     endif
     &                     AVERAGE(ng) % avgrain(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#   ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idSrad,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     AVERAGE(ng) % avgsrf(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tAVGindx(ng)
          END IF          
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Write out surface u-momentun stress.
!
      IF (Hout(idUsms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUsms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgsus(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tAVGindx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentun stress.
!
      IF (Hout(idVsms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVsms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgsvs(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentun stress.
!
      IF (Hout(idUbms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUbms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgbus(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tAVGindx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentun stress.
!
      IF (Hout(idVbms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVbms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgbvs(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Hout(idUbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVGid(ng),                     &
     &                       avgVid(idUbld(itrc),ng),                   &
     &                       tAVGindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask(LBi,LBj),                 &
#  endif
     &                       AVERAGE(ng) % avgbedldu(LBi,LBj,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))),            &
     &                          tAVGindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Hout(idVbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVGid(ng),                     &
     &                       avgVid(idVbld(itrc),ng),                   &
     &                       tAVGindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask(LBi,LBj),                 &
#  endif
     &                       AVERAGE(ng) % avgbedldv(LBi,LBj,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))),            &
     &                          tAVGindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncAVGid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,30) tAVGindx(ng)
      END IF
!
  10  FORMAT (/,' WRT_AVG - error while writing variable: ',a,/,11x,    &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_AVG - unable to synchronize averages NetCDF to ', &
     &        'disk.')
  30  FORMAT (6x,'WRT_AVG   - wrote averaged fields into time ',        &
     &        'record =',t72,i7.7)
#else
      SUBROUTINE wrt_avg
#endif
      RETURN
      END SUBROUTINE wrt_avg
