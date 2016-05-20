#include "cppdefs.h"
#ifdef AVERAGES2
#undef AVERAGES_DETIDE
      SUBROUTINE wrt_avg2 (ng)
!
!svn $Id: wrt_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      USE mod_average2
      USE mod_forces
      USE mod_grid
#ifdef ICE_MODEL
      USE mod_ice
#endif
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
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
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
      integer :: Fcount, gfactor, gtype, i, itrc, status

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
      AVG2(ng)%Rindex=AVG2(ng)%Rindex+1
      Fcount=AVG2(ng)%Fcount
      AVG2(ng)%Nrec(Fcount)=AVG2(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, AVG2(ng)%name,                    &
     &                      TRIM(Vname(idtime,ng)), AVG2time(ng:),      &
     &                      (/AVG2(ng)%Rindex/), (/1/),                 &
     &                      ncid = AVG2(ng)%ncid,                       &
     &                      varid = AVG2(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Aout2(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idFsur),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgzeta)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), AVG2(ng)%Rindex
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
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                       &
     &                   AVG2(ng)%Vid(idFsuD),                          &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask_avg2,                         &
#  endif
     &                   TIDES(ng) % zeta_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsuD)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout2(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUbar),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgu2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), AVG2(ng)%Rindex
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
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                       &
     &                   AVG2(ng)%Vid(idu2dD),                          &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_avg2,                         &
#  endif
     &                   TIDES(ng) % ubar_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu2dD)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout2(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVbar),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgv2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), AVG2(ng)%Rindex
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
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idv2dD), &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_avg2,                         &
#  endif
     &                   TIDES(ng) % vbar_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idv2dD)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout2(idu2dE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idu2dE),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgu2dE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout2(idv2dN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idv2dN),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgv2dN)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Write out time-averaged mass fluxes for 3D momentum coupling.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idUfx1), &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_avg2,                         &
#  endif
     &                   AVERAGE2(ng) % avgDU_avg1)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx1)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idUfx2), &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_avg2,                         &
#  endif
     &                   AVERAGE2(ng) % avgDU_avg2)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx2)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idVfx1), &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_avg2,                         &
#  endif
     &                   AVERAGE2(ng) % avgDV_avg1)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx1)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idVfx2), &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_avg2,                         &
#  endif
     &                   AVERAGE2(ng) % avgDV_avg2)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx2)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef NEARSHORE_MELLOR
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Aout2(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idW2xx),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgSxx2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Aout2(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idW2xy),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgSxy2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Aout2(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idW2yy),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgSyy2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the XI-direction.
!
      IF (Aout2(idU2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idU2rs),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgu2RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the ETA-direction
!
      IF (Aout2(idV2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idV2rs),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgv2RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the XI-direction.
!
      IF (Aout2(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idU2Sd),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgu2Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Aout2(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idV2Sd),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgv2Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), AVG2(ng)%Rindex
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
      IF (Aout2(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUvel),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgu3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), AVG2(ng)%Rindex
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
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idu3dD), &
     &                   AVG2(ng)%Rindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % umask_avg2,                         &
#   endif
     &                   TIDES(ng) % u_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu3dD)), AVG2(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout2(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVvel),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgv3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout2(idu3dE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idu3dE),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgu3dE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout2(idv3dN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idv3dN),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgv3dN)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), AVG2(ng)%Rindex
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
        IF (Aout2(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                   &
     &                       AVG2(ng)%Tid(itrc),                        &
     &                       AVG2(ng)%Rindex, gtype,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask_avg2,                     &
#  endif
     &                       AVERAGE2(ng) % avgt(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                        AVG2(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# ifdef BIO_GOANPZ
!---------------------------------------------
!  Write out stationary tracer type variables.
!---------------------------------------------
      DO itrc=1,NTS(ng)
        !print *, 'XXX',AVG2(ng)%ncid,avg2TSid(itrc,ng),AVG2(ng)%Rindex
        !print *, 'xxx',Aout2(idTSvar(itrc),ng)
        IF (Aout2(idTSvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                   &
     &                       AVG2(ng)%TSid(itrc),                       &
     &                       AVG2(ng)%Rindex, gtype,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask_avg2,                     &
#  endif
     &                       AVERAGE2(ng) % avgst(:,:,1,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTSvar(itrc))),           &
     &                        AVG2(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Write out density anomaly.
!
      IF (Aout2(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idDano),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgrho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), AVG2(ng)%Rindex
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
      IF (Aout2(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idHsbl),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avghsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), AVG2(ng)%Rindex
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
      IF (Aout2(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idHbbl),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avghbbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef ICE_MODEL
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout2(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUice),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avguice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout2(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVice),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgvice)
       IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward ice momentum component (m/s) at RHO-points.
!
      IF (Aout2(idUiceE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUiceE),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avguiceE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUiceE)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward ice momentum component (m/s) at RHO-points.
!
      IF (Aout2(idViceN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idViceN),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgviceN)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idViceN)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

!
! Write out ice concentration.
!
      IF (Aout2(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idAice),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgaice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice thickness.
!
      IF (Aout2(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idHice),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avghice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice/snow surface temperature.
!
      IF (Aout2(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idTice),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgtice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice internal temperature.
!
      IF (Aout2(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idTimid),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgtimid)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out snow thickness.
!
      IF (Aout2(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idHsno),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avghsno)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef MELT_PONDS
!
! Write out surface water fraction (on ice).
!
      IF (Aout2(idApond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idApond),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgapond)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idApond)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out surface water thickness (on ice).
!
      IF (Aout2(idHpond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idHpond),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avghpond)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHpond)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#endif
!
! Write out ice age.
!
      IF (Aout2(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idAgeice),                      &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgageice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-ocean mass flux
!
      IF (Aout2(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idIomflx),                      &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgiomflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 11
!
      IF (Aout2(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idSig11),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgsig11)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 12
!
      IF (Aout2(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idSig12),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgsig12)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 22
!
      IF (Aout2(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idSig22),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                   AVERAGE2(ng) % avgsig22)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out temperature of molecular sublayer under ice.
!
      IF (Aout2(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idT0mk),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgT0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out salinity of molecular sublayer under ice.
!
      IF (Aout2(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idS0mk),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgS0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water friction velocity
!
      IF (Aout2(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idTauiw),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgutau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water momentum transfer coefficient
!
      IF (Aout2(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idChuiw),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgchu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef SOLVE3D
!
!  Write out surface net heat flux.
!
      IF (Aout2(idTsur(itemp),ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
        defined OPT_OBSERVATIONS
        scale=1.0_r8/(rho0*Cp)
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0*Cp
#  endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idTsur(itemp)),                 &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgstf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),             &
     &                        AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Aout2(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idTsur(isalt)),                 &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgswf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))),             &
     &                        AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Aout2(idLhea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idLhea),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avglhf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout2(idShea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idShea),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgshf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Aout2(idLrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idLrad),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avglrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-wind.
!
      IF (Aout2(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUair),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avguwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-wind.
!
      IF (Aout2(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVair),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#  endif
     &                     AVERAGE2(ng) % avgvwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward surface wind (m/s) at RHO-points.
!
      IF (Aout2(idUairE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUairE),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avguwindE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUairE)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward surface wind (m/s) at RHO-points.
!
      IF (Aout2(idVairN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVairN),                       &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgvwindN)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVairN)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Aout2(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idevap),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#    endif
     &                     AVERAGE2(ng) % avgevap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Aout2(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idrain),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#    endif
     &                     AVERAGE2(ng) % avgrain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Aout2(idSrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idSrad),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg2,                       &
#   endif
     &                     AVERAGE2(ng) % avgsrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out surface u-momentum stress.
!
      IF (Aout2(idUsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUsms),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgsus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout2(idVsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVsms),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgsvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout2(idUbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idUbms),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgbus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout2(idVbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                     &
     &                     AVG2(ng)%Vid(idVbms),                        &
     &                     AVG2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_avg2,                       &
# endif
     &                     AVERAGE2(ng) % avgbvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), AVG2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Aout2(idUbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                   &
     &                       AVG2(ng)%Vid(idUbld(itrc)),                &
     &                       AVG2(ng)%Rindex, gtype,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask_avg2,                     &
#  endif
     &                       AVERAGE2(ng) % avgbedldu(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))),            &
     &                          AVG2(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Aout2(idVbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, AVG2(ng)%ncid,                   &
     &                       AVG2(ng)%Vid(idVbld(itrc)),                &
     &                       AVG2(ng)%Rindex, gtype,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask_avg2,                     &
#  endif
     &                       AVERAGE2(ng) % avgbedldv(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))),            &
     &                          AVG2(ng)%Rindex
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
      CALL netcdf_sync (ng, iNLM, AVG2(ng)%name, AVG2(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) AVG2(ng)%Rindex
!
  10  FORMAT (/,' WRT_AVG2 - error while writing variable: ',a,/,11x,   &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_AVG2   - wrote averaged fields into time ',       &
     &        'record =',t72,i7.7)
#else
      SUBROUTINE wrt_avg2
#endif
      RETURN
      END SUBROUTINE wrt_avg2
