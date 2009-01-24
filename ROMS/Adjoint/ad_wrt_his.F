#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_wrt_his (ng)
!
!svn $Id: ad_wrt_his.F 677 2008-08-05 20:17:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested adjoint model fields into adjoint     !
!  history NetCDF file.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT_NOT_YET || defined BBL_MODEL_NOT_YET
      USE mod_sediment
# endif
      USE mod_stepping

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
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
      integer :: i, j, gfactor, gtype, status
      integer :: kout
# ifdef SOLVE3D
      integer :: itrc, k, nout
# endif

      integer :: nf_fwrite2d
# ifdef SOLVE3D
      integer :: nf_fwrite3d
# endif

      real(r8) :: scale
      real(r8) :: Tval(1)
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out adjoint fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
!  Determine time index to write.
!
# ifdef SOLVE3D
      kout=kstp(ng)
# else
      kout=kstp(ng)
# endif
# ifdef SOLVE3D
      IF (iic(ng).ne.ntend(ng)) THEN
        nout=nnew(ng)
      ELSE
        nout=nstp(ng)
      END IF
# endif
!
!  Set time record index.
!
      tADJindx(ng)=tADJindx(ng)+1
      NrecADJ(ng)=NrecADJ(ng)+1
!
!  If requested, set time index to recycle time records in the adjoint
!  file.
!
      IF (LcycleADJ(ng)) THEN
        tADJindx(ng)=MOD(tADJindx(ng)-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtTime(ng)) THEN
        IF (LwrtPER(ng)) THEN
          Tval(1)=REAL(Nrun,r8)*day2sec
        ELSE
          Tval(1)=time(ng)
        END IF
        CALL netcdf_put_fvar (ng, iADM, ADJname(ng),                    &
     &                        TRIM(Vname(1,idtime)), tval,              &
     &                        (/tADJindx(ng)/), (/1/),                  &
     &                        ncid = ncADJid(ng),                       &
     &                        varid = adjVid(idtime,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=rho0                            ! m2/s2 to Pa
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idUsms,ng),      &
     &                   tADJindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % ad_ustr(LBi,LBj,1,Lfout(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=rho0                            ! m2/s2 to Pa
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idVsms,ng),      &
     &                   tADJindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % ad_vstr(LBi,LBj,1,Lfout(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net heat flux. Notice that different tracer fluxes
!  are written at their own fixed time-dimension (of size Nfrec) to
!  allow 4DVAR adjustments at other times in addition to initial time.
!
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          scale=rho0*Cp                     ! Celsius m/s to W/m2
        ELSE
          scale=1.0_r8
        END IF
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng),                       &
     &                     adjVid(idTsur(itrc),ng),                     &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     FORCES(ng) % ad_tflux(LBi,LBj,1,             &
     &                                           Lfout(ng),itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))), Lfout(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
!
!  Write out bathymetry.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idbath,ng),      &
     &                   tADJindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   GRID(ng)% ad_h(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idbath)), tADJindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        IF (LwrtState2d(ng)) THEN
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idFsur,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       OCEAN(ng)% ad_zeta(LBi,LBj,kout))
        ELSE
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idFsur,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       OCEAN(ng)% ad_zeta_sol(LBi,LBj))
        ENDIF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        IF (LwrtState2d(ng)) THEN
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUbar,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % umask(LBi,LBj),                 &
# endif
     &                       OCEAN(ng) % ad_ubar(LBi,LBj,kout))
        ELSE
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUbar,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % umask(LBi,LBj),                 &
# endif
     &                       OCEAN(ng) % ad_ubar_sol(LBi,LBj))
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        IF (LwrtState2d(ng)) THEN
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVbar,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
      &                      GRID(ng) % vmask(LBi,LBj),                 &
# endif
     &                       OCEAN(ng) % ad_vbar(LBi,LBj,kout))
        ELSE
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVbar,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
      &                      GRID(ng) % vmask(LBi,LBj),                 &
# endif
     &                       OCEAN(ng) % ad_vbar_sol(LBi,LBj))
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idUvel,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % ad_u(LBi,LBj,1,nout))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idVvel,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % ad_v(LBi,LBj,1,nout))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tADJindx(ng)
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
          status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjTid(itrc,ng),    &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#  endif
     &                       OCEAN(ng) % ad_t(LBi,LBj,1,nout,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          tADJindx(ng)
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
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idDano,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % ad_rho(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idVvis,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     MIXING(ng) % ad_Akv(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idTdif,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     MIXING(ng) % ad_Akt(LBi,LBj,0,itemp))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idSdif,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     MIXING(ng) % ad_Akt(LBi,LBj,0,isalt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifndef ADJUST_STFLUX
!
!  Write out net surface active tracer fluxes.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTsur(itrc),ng)) THEN
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
          IF (itrc.eq.itemp) THEN
!!          scale=rho0*Cp
            scale=1.0_r8/(rho0*Cp)
          ELSE
            scale=1.0_r8
          END IF
#   else
          scale=1.0_r8
#   endif
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iADM, ncADJid(ng),                     &
     &                       adjVid(idTsur(itrc),ng),                   &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#   endif
     &                       FORCES(ng) % ad_stflx(LBi,LBj,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
# ifndef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS
!!      scale=rho0
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUsms,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     FORCES(ng) % ad_sustr(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
!!      scale=rho0
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVsms,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     FORCES(ng) % ad_svstr(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
!!      scale=-rho0
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUbms,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
# endif
     &                     FORCES(ng) % ad_bustr_sol(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
!!      scale=-rho0
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVbms,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
# endif
     &                     FORCES(ng) % ad_bvstr_sol(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize adjoint history NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncADJid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ADJname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag, 1)
# endif
# ifdef SOLVE3D
      IF (Master) WRITE (stdout,30) kout, nout, tADJindx(ng)
# else
      IF (Master) WRITE (stdout,30) kout, tADJindx(ng)
# endif
!
  10  FORMAT (/,' AD_WRT_HIS - error while writing variable: ',a,/,14x, &
     &        'into adjoint NetCDF file for time record: ',i4)
  20  FORMAT (/,' AD_WRT_HIS - unable to synchronize to disk adjoint ', &
     &        'file:',/,14x,a)
# ifdef SOLVE3D
  30  FORMAT (3x,'AD_WRT_HIS   - wrote adjoint  fields (Index=', i1,    &
     &        ',',i1,') into time record = ',i7.7)
# else
  30  FORMAT (3x,'AD_WRT_HIS   - wrote adjoint  fields (Index=', i1,    &
     &        ') into time record = ',i7.7)
# endif
      RETURN
      END SUBROUTINE ad_wrt_his
#else
      SUBROUTINE ad_wrt_his
      RETURN
      END SUBROUTINE ad_wrt_his
#endif
