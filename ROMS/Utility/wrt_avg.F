#include "cppdefs.h"
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      SUBROUTINE wrt_avg (ng)
!
!svn $Id: wrt_avg.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged fields into averages     !
!  NetCDF file.                                                        !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average
      USE mod_forces
      USE mod_grid
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
!
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
      SourceFile='wrt_avg.F'
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
      AVG(ng)%Rindex=AVG(ng)%Rindex+1
      Fcount=AVG(ng)%Fcount
      AVG(ng)%Nrec(Fcount)=AVG(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, AVG(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), AVGtime(ng:),       &
     &                      (/AVG(ng)%Rindex/), (/1/),                  &
     &                      ncid = AVG(ng)%ncid,                        &
     &                      varid = AVG(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idFsur), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgzeta)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), AVG(ng)%Rindex
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
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idFsuD), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     TIDES(ng) % zeta_detided)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbar), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgu2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), AVG(ng)%Rindex
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
      IF (Aout(idu2dD,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idu2dD), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     TIDES(ng) % ubar_detided)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbar), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgv2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), AVG(ng)%Rindex
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
      IF (Aout(idv2dD,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idv2dD), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     TIDES(ng) % vbar_detided)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu2dE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idu2dE), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgu2dE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv2dN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idv2dN), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgv2dN)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), AVG(ng)%Rindex
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
      IF (Aout(idUfx1,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUfx1), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDU_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idUfx2,ng)) THEN
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUfx2), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDU_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVfx1,ng)) THEN
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVfx1), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDV_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVfx2,ng)) THEN
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVfx2), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDV_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WEC_MELLOR
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Aout(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW2xx), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgSxx2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Aout(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW2xy), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgSxy2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Aout(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW2yy), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgSyy2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WEC
!
!  Write out 2D wec stress in the XI-direction.
!
      IF (Aout(idU2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idU2rs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgu2RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D wec stress in the ETA-direction
!
      IF (Aout(idV2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idV2rs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgv2RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the XI-direction.
!
      IF (Aout(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idU2Sd), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgu2Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idV2Sd), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgv2Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WEC_VF
!
!  Write out quasi-static sea level adjustment.
!
      IF (Aout(idWztw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWztw), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWztw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWztw)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quasi-static pressure.
!
      IF (Aout(idWqsp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWqsp), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWqsp)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWqsp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Bernoulli head.
!
      IF (Aout(idWbeh,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWbeh), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWbeh)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbeh)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_HEIGHT
      IF (Aout(idWamp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWamp), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWamp)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_LENGTH
      IF (Aout(idWlen,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWlen), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWlen)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_LENGTHP
      IF (Aout(idWlep,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWlep), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWlep)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlep)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_DIR
      IF (Aout(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWdir), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdir)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_TOP_PERIOD
      IF (Aout(idWptp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWptp), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWptp)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_BOT_PERIOD
      IF (Aout(idWpbt,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWpbt), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWpbt)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef BBL_MODEL
      IF (Aout(idWorb,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWorb), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWorb)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      IF (Aout(idWdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWdif), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdif)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR
      IF (Aout(idWdib,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWdib), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdib)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdib)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Aout(idWdiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWdiw), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdiw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdiw)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ROLLER_SVENDSEN
      IF (Aout(idWbrk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWbrk), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWbrk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbrk)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WEC_ROLLER
      IF (Aout(idWdis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWdis), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdis)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdis)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ROLLER_RENIERS
      IF (Aout(idWrol,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWrol), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWrol)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWrol)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef UV_KIRBY
      IF (Aout(idUwav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUwav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUwav)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Aout(idVwav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVwav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVwav)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), AVG(ng)%Rindex
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
      IF (Aout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUvel), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgu3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), AVG(ng)%Rindex
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
      IF (Aout(idu3dD,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idu3dD), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     TIDES(ng) % u_detided)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVvel), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgv3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), AVG(ng)%Rindex
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
      IF (Aout(idv3dD,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idv3dD), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#   endif
     &                     TIDES(ng) % v_detided)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu3dE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idu3dE), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgu3dE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv3dN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idv3dN), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgv3dN)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef AVERAGES_UV_KIRBY
!
!  Write out coupling current velocity (x component).
!
      IF (Aout(idUwav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUwav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avguWave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out coupling current velocity (y component).
!
      IF (Aout(idVwav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVwav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgvWave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined WEC_MELLOR
!
!  Write out 3D radiation stress, Sxx-component.
!
      IF (Aout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3xx), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                           &
#   endif
     &                     AVERAGE(ng) % avgSxx3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xx)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Sxy-component.
!
      IF (Aout(idW3xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3xy), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgSxy3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xy)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Syy-component.
!
      IF (Aout(idW3yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3yy), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgSyy3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3yy)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szx-component.
!
      IF (Aout(idW3zx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3zx), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgSzx3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zx)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szy-component.
!
      IF (Aout(idW3zy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3zy), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgSzy3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zy)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WEC
!
!  Write out 3D wec stress in the XI-direction, u_Rstress.
!
      IF (Aout(idU3rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idU3rs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     AVERAGE(ng) % avgu3RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3rs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D wec stress in the ETA-direction, v_Rstress.
!
      IF (Aout(idV3rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idV3rs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     AVERAGE(ng) % avgv3RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3rs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D stokes momentum component (m/s) in the XI-direction.
!
      IF (Aout(idU3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idU3Sd), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     AVERAGE(ng) % avgu3Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3Sd)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idV3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idV3Sd), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     AVERAGE(ng) % avgv3Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3Sd)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega-stokes vertical velocity (m/s).
!
      IF (Aout(idW3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3Sd), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgW3Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3Sd)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out z-coordinate w-stokes vertical velocity (m/s).
!
      IF (Aout(idW3St,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idW3St), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgW3St)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3St)), AVG(ng)%Rindex
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
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idOvel), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgw3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Aout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWvel), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgwvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
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
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Tid(itrc), &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgt(:,:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  ifdef BEST_NPZ
#   ifdef BIOFLUX
!---------------------------------------------
!  Write out bioflux arrays.
!---------------------------------------------
      DO itrc=1,NT(ng)
        DO itrc2=1,NT(ng)
!g        !print *, 'XXX',AVG(ng)%ncid,AVG(ng)%TSid(itrc),AVG(ng)%Rindex
!g       !print *, 'xxx',Aout(idTSvar(itrc),ng)
        IF (Aout(idTBFvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar

          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%TBFid(itrc),                       &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       1,NT(ng),1, NT(ng), scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgbflx(:,:))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTBFvar(itrc))),          &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status

            RETURN
          END IF
        END IF
      END DO
        END DO
#   endif

#   ifdef STATIONARY
!---------------------------------------------
!  Write out 3D stationary variable array.
!---------------------------------------------
      DO itrc=1,NTS(ng)
        !print *, 'XXX',AVG(ng)%ncid,AVG(ng)%TSid(itrc),AVG(ng)%Rindex
        !print *, 'xxx',Aout(idTSvar(itrc),ng)
        IF (Aout(idTSvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%TSid(itrc),&
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgst(:,:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTSvar(itrc))),           &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#   ifdef STATIONARY2
!---------------------------------------------
!  Write out 2D stationary variable array.
!---------------------------------------------
      DO itrc=1,NTS2(ng)
        IF (Aout(idTS2var(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%TS2id(itrc),                       &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgst2(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTS2var(itrc))),          &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif

#   ifdef PROD3
!---------------------------------------------
!  Write 3D production variables.
!---------------------------------------------
      DO itrc=1,NPT3(ng)
        !print *,
	'XXX',AVG(ng)%ncid,AVG(ng)%PT3id(itrc),AVG(ng)%Rindex
        !print *, 'xxx',Aout(idPT3var(itrc),ng)
        IF (Aout(idPT3var(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%PT3id(itrc),                       &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgpt3(:,:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
             WRITE (stdout,10) TRIM(Vname(1,idPT3var(itrc))),           &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif

#   ifdef PROD2
!---------------------------------------------
!  Write 2D production variables.
!---------------------------------------------
      DO itrc=1,NPT2(ng)
        IF (Aout(idPT2var(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,AVG(ng)%PT2id(itrc),&
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgpt2(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
             WRITE (stdout,10) TRIM(Vname(1,idPT2var(itrc))),           &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
!-----------------------------
!  Write out benthic variables.
!-----------------------------
#   ifdef BENTHIC
      DO itrc=1,NBeT(ng)
        IF (Aout(idBeTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Bid(itrc), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgbt(:,:,1,itrc))
! switch if have more than one benthic level

!          gtype=gfactor*r3dvar
!          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Bid(itrc),&
!     &                       AVG(ng)%Rindex, gtype,                    &
!     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,      &
!#  ifdef MASKING
!     &                       GRID(ng) % rmask,                          &
!#  endif
!     &                       AVERAGE(ng) % avgbt(LBi,LBj,1,itrc))

          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBeTvar(itrc))),          &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif

!-----------------------------
!  Write out ice bio variables.
!-----------------------------
#   ifdef ICE_BIO
#    ifdef CLIM_ICE_1D
!  Write out ice variables.

      DO itrc=1,NIceT(ng)
        IF (Aout(idIceBvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                    &
     &                     AVG(ng)%IceBid(itrc),                        &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgit(:,:,itrc))
! switch if have more than one ice level

!          gtype=gfactor*r3dvar
!          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
!     &                       AVG(ng)%IceBid(itrc),                      &
!     &                       AVG(ng)%Rindex, gtype,                     &
!     &                       LBi, UBi, LBj, UBj, 1, NIL(ng), scale,     &
!#  ifdef MASKING
!     &                       GRID(ng) % rmask,                          &
!#  endif
!     &                       AVERAGE(ng) % avgit(LBi,LBj,itrc))

          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idIceBvar(itrc))),            &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#    else
      IF (Aout(idIcePhL,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                   AVG(ng)%Vid(idIcePhL),                         &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgIPhL)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIcePhL)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!

      IF (Aout(idIceNO3,ng)) THEN
!       IF (Hout(idIceBvar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                   AVG(ng)%Vid(idIceNO3),                         &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgINO3)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIceNO3)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!

      IF (Aout(idIceNH4,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                   AVG(ng)%Vid(idIceNH4),                         &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgINH4)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIceNH4)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
        IF (Hout(idIceLog,ng)) THEN
!       IF (Hout(idIceBvar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                   AVG(ng)%Vid(idIceLog),                         &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgILog)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIceLog)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#    endif
#   endif
#  endif

#  ifdef BIO_GOANPZ
!---------------------------------------------
!  Write out stationary tracer type variables.
!---------------------------------------------
      DO itrc=1,NTS(ng)
        !print *, 'XXX',AVG(ng)%ncid,AVG(ng)%TSid(itrc),AVG(ng)%Rindex
        !print *, 'xxx',Aout(idTSvar(itrc),ng)
        IF (Aout(idTSvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%TSid(itrc),&
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgst(:,:,1,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTSvar(itrc))),            &
     &                        AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      DO itrc=1,NAT
        IF (Aout(idTrcD(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(idTrcD(itrc)),                 &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       TIDES(ng) % t_detided(:,:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTrcD(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idDano), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgrho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef NEMURO_SED1
!
!  Write out PON in sediment.
!
      IF (Aout(idPONsed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idPONsed),                       &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgPONsed)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPONsed)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out OPAL in sediment.
!
      IF (Aout(idOPALsed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idOPALsed),                      &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgOPALsed)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOPALsed)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out DENIT in sediment.
!
      IF (Aout(idDENITsed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idDENITsed),                     &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgDENITsed)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDENITsed)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out buried PON in sediment.
!
      IF (Aout(idPONbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idPONbur),                       &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgPONbur)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPONbur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out buried OPAL in sediment.
!
      IF (Aout(idOPALbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idOPALbur),                      &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgOPALbur)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOPALbur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef PRIMARY_PROD
!
!  Write out Net primary production.
!
      IF (Aout(idNPP,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idNPP),                          &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgNPP)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idNPP)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Aout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHsbl), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avghsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), AVG(ng)%Rindex
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
      IF (Aout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHbbl), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avghbbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef SSSFLX
!
!  Write out sea surface salinity correction flux
!
      IF (Aout(idSSSf,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idSSSf), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgsssflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSSSf)), AVG(ng)%Rindex
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
      IF (Aout(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUice), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % umask_avg,                          &
#   endif
     &                   AVERAGE(ng) % avguice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVice), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % vmask_avg,                          &
#   endif
     &                   AVERAGE(ng) % avgvice)
       IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward ice momentum component (m/s) at RHO-points.
!
      IF (Aout(idUiceE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idUiceE),                        &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avguiceE)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUiceE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward ice momentum component (m/s) at RHO-points.
!
      IF (Aout(idViceN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idViceN),                        &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgviceN)
       IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idViceN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice concentration.
!
      IF (Aout(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idAice), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgaice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice thickness.
!
      IF (Aout(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHice), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avghice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice/snow surface temperature.
!
      IF (Aout(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idTice), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgtice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice internal temperature.
!
      IF (Aout(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idTimid),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgtimid)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out snow thickness.
!
      IF (Aout(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHsno), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avghsno)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef MELT_PONDS
!
! Write out surface water fraction (on ice).
!
      IF (Aout(idApond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idApond),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgapond)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idApond)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out surface water thickness (on ice).
!
      IF (Aout(idHpond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHpond),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avghpond)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHpond)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#endif
!
! Write out ice age.
!
      IF (Aout(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                   AVG(ng)%Vid(idAgeice),                         &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgageice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-ocean mass flux
!
      IF (Aout(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                   AVG(ng)%Vid(idIomflx),                         &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgiomflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 11
!
      IF (Aout(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idSig11),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgsig11)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 12
!
      IF (Aout(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idSig12),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgsig12)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 22
!
      IF (Aout(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idSig22),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgsig22)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out temperature of molecular sublayer under ice.
!
      IF (Aout(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idT0mk), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgT0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out salinity of molecular sublayer under ice.
!
      IF (Aout(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idS0mk), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgS0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out frazil ice growth rate.
!
      IF (Aout(idWfr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWfr),  &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgWfr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWfr)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice growth/melt rate.
!
      IF (Aout(idWai,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWai),  &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgWai)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWai)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice growth/melt rate.
!
      IF (Aout(idWao,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWao),  &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgWao)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWao)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice growth/melt rate.
!
      IF (Aout(idWio,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWio),  &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgWio)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWio)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice melt runoff rate.
!
      IF (Aout(idWro,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWro),  &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgWro)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWro)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice divergence rate.
!
      IF (Aout(idWdiv,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idWdiv), &
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgWdiv)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdiv)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water friction velocity
!
      IF (Aout(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idTauiw),&
     &                   AVG(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgutau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water momentum transfer coefficient
!
      IF (Aout(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idChuiw),&
     &                   AVG(ng)%Rindex, gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   AVERAGE(ng) % avgchu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out 2D potential vorticity.
!
      IF (Aout(id2dPV,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*p2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id2dPV), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
     &                     AVERAGE(ng) % avgpvor2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id2dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D relative vorticity.
!
      IF (Aout(id2dRV,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*p2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id2dRV), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
     &                     AVERAGE(ng) % avgrvor2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id2dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D potential vorticity.
!
      IF (Aout(id3dPV,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id3dPV), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     AVERAGE(ng) % avgpvor3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id3dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D relative vorticity.
!
      IF (Aout(id3dRV,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id3dRV), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     AVERAGE(ng) % avgrvor3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id3dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out quadratic <zeta*zeta> term.
!
      IF (Aout(idZZav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idZZav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgZZ)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idZZav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      IF (Aout(idU2av,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idU2av), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     AVERAGE(ng) % avgU2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      IF (Aout(idV2av,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idV2av), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     AVERAGE(ng) % avgV2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write u-volume flux.
!
      IF (Aout(idHUav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHUav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgHuon)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write v-volume flux.
!
      IF (Aout(idHVav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idHVav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgHvom)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*u> term.
!
      IF (Aout(idUUav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUUav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgUU)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*v> term.
!
      IF (Aout(idUVav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUVav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUV)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <v*v> term.
!
      IF (Aout(idVVav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVVav), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVV)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO i=1,NT(ng)
        IF (Aout(idTTav(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(idTTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgTT(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out active tracer volume fluxes.
!
      DO i=1,NT(ng)
        IF (Aout(iHUTav(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(iHUTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       AVERAGE(ng) % avgHuonT(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,iHUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(iHVTav(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(iHVTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       AVERAGE(ng) % avgHvomT(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,iHVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO i=1,NT(ng)
        IF (Aout(idUTav(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(idUTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       AVERAGE(ng) % avgUT(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(idVTav(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(idVTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       AVERAGE(ng) % avgVT(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVvis), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgAKv,                        &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idTdif), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgAKt,                        &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idSdif), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgAKs,                        &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined RUNOFF
!
!  Write out surface air pressure.
!
      IF (Aout(idRunoff,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idRunoff),                       &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgRunoff)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRunoff)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Aout(idPair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idPair), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgPair)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Aout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUair), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgUwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

      IF (Aout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVair), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgVwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
       defined OPT_OBSERVATIONS
       scale=1.0_r8/(rho0*Cp)
#   else
        scale=1.0_r8
#   endif
#  else
#   ifdef SOGLOBEC
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
#  endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idTsur(itemp)),                  &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgstf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                      &
     &                     AVG(ng)%Vid(idTsur(isalt)),                  &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgswf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))),             &
     &                        AVG(ng)%Rindex
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
      IF (Aout(idLhea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idLhea), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avglhf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idShea), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgshf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idLrad), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avglrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), AVG(ng)%Rindex
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
      IF (Aout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idevap), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgevap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Aout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idrain), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgrain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), AVG(ng)%Rindex
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
      IF (Aout(idSrad,ng)) THEN
#    if defined ADJOINT || defined SOGLOBEC
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idSrad), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgsrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), AVG(ng)%Rindex
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
      IF (Aout(idUsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUsms), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     AVERAGE(ng) % avgsus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVsms), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     AVERAGE(ng) % avgsvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout(idUbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbms), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     AVERAGE(ng) % avgbus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbms), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     AVERAGE(ng) % avgbvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef BBL
!
!  Write out bottom u-momentum stress at rho points.
!
      IF (Aout(idUbrs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbrs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbrs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress at rho points.
!
      IF (Aout(idVbrs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbrs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbrs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-wave stress at rho points.
!
      IF (Aout(idUbws,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbws), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                             &
# endif
     &                     AVERAGE(ng) % avgUbws)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-wave stress at rho points.
!
      IF (Aout(idVbws,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbws), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbws)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom max u-wave curr stress at rho points.
!
      IF (Aout(idUbcs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbcs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbcs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out max bottom v-wave curr stress at rho points.
!
      IF (Aout(idVbcs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbcs), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbcs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out max bottom wave curr stress magnitude at rho points.
!
      IF (Aout(idUVwc,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUVwc), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUVwc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-wave orbital vel at rho points.
!
      IF (Aout(idUbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbot), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-wave orbital vel at rho points.
!
      IF (Aout(idVbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbot), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-cur above bed at rho points.
!
      IF (Aout(idUbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUbur), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbur)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-cur above the bed at rho points.
!
      IF (Aout(idVbvr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVbvr), &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbvr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif

# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Aout(idUbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(idUbld(itrc)),                 &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       SEDBED(ng) % avgbedldu(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Aout(idVbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, AVG(ng)%ncid,                    &
     &                       AVG(ng)%Vid(idVbld(itrc)),                 &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       SEDBED(ng) % avgbedldv(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))),            &
     &                          AVG(ng)%Rindex
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
      CALL netcdf_sync (ng, iNLM, AVG(ng)%name, AVG(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef NESTING
      IF (Master) WRITE (stdout,20) AVG(ng)%Rindex, ng
#else
      IF (Master) WRITE (stdout,20) AVG(ng)%Rindex
#endif
!
  10  FORMAT (/,' WRT_AVG - error while writing variable: ',a,/,11x,    &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_AVG   - wrote averaged fields into time ',        &
# ifdef NESTING
     &        'record =',t72,i7.7,2x,i2.2)
# else
     &        'record =',t72,i7.7)
# endif
#else
      SUBROUTINE wrt_avg
#endif
      RETURN
      END SUBROUTINE wrt_avg
