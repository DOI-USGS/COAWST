#include "cppdefs.h"
#if defined TANGENT || defined TL_IOMS
      SUBROUTINE tl_wrt_his (ng)
!
!svn $Id: tl_wrt_his.F 677 2008-08-05 20:17:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine writes requested tangent linear model fields into      !
!  into tangent history NetCDF file.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_iounits
# ifdef SOLVE3D
      USE mod_mixing
# endif
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
      integer :: gfactor, gtype, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k, tile
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
!  Write out tangent linear fields.
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
!  Set time record index.
!
      tTLMindx(ng)=tTLMindx(ng)+1
      NrecTLM(ng)=NrecTLM(ng)+1
!
!  If requested, set time index to recycle time records in the tangent
!  linear file.
!
      IF (LcycleTLM(ng)) THEN
        tTLMindx(ng)=MOD(tTLMindx(ng)-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtPER(ng)) THEN
        Tval(1)=REAL(Nrun,r8)*day2sec
      ELSE
        Tval(1)=time(ng)
      END IF
      CALL netcdf_put_fvar (ng, iTLM, TLMname(ng),                      &
     &                      TRIM(Vname(1,idtime)), tval,                &
     &                      (/tTLMindx(ng)/), (/1/),                    &
     &                      ncid = ncTLMid(ng),                         &
     &                      varid = tlmVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=rho0                            ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idUsms,ng),      &
     &                   tTLMindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % tl_ustr(LBi,LBj,1,Lfout(ng)))
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
      scale=rho0                            ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idVsms,ng),      &
     &                   tTLMindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % tl_vstr(LBi,LBj,1,Lfout(ng)))
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
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          scale=rho0*Cp                     ! Celsius m/s to W/m2
        ELSE
          scale=1.0_r8
        END IF
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iTLM, ncTLMid(ng),                       &
     &                     tlmVid(idTsur(itrc),ng),                     &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     FORCES(ng) % tl_tflux(LBi,LBj,1,             &
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
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idFsur,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % tl_zeta(LBi,LBj,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idRzet,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % tl_rzeta(LBi,LBj,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRzet)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idUbar,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % tl_ubar(LBi,LBj,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idRu2d,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % tl_rubar(LBi,LBj,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu2d)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif 
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idRuct,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#    endif
     &                     COUPLING(ng) % tl_rufrc(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRuct)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idUfx1,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     COUPLING(ng) % tl_DU_avg1(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idUfx2,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     COUPLING(ng) % tl_DU_avg2(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idVbar,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % tl_vbar(LBi,LBj,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idRv2d,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % tl_rvbar(LBi,LBj,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv2d)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idRvct,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#    endif
     &                     COUPLING(ng) % tl_rvfrc(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRvct)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idVfx1,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     COUPLING(ng) % tl_DV_avg1(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iTLM, ncTLMid(ng), tlmVid(idVfx2,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     COUPLING(ng) % tl_DV_avg2(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idUvel,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % tl_u(LBi,LBj,1,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idRu3d,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % tl_ru(LBi,LBj,1,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu3d)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idVvel,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % tl_v(LBi,LBj,1,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idRv3d,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % tl_rv(LBi,LBj,1,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv3d)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmTid(itrc,ng),    &
     &                       tTLMindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#  endif
     &                       OCEAN(ng) % tl_t(LBi,LBj,1,NOUT,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          tTLMindx(ng)
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
        status=nf_fwrite3d(ng, iTLM, ncTLMid(ng), tlmVid(idDano,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % tl_rho(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined FORWARD_MIXING     && \
     (defined BVF_MIXING     || defined GLS_MIXING  || \
      defined LMD_MIXING     || defined MY25_MIXING)
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idVvis,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     MIXING(ng) % Akv(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tTLMindx(ng)
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
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idTdif,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     MIXING(ng) % Akt(LBi,LBj,0,itemp))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idSdif,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     MIXING(ng) % Akt(LBi,LBj,0,isalt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idMtke,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     MIXING(ng) % tl_tke(LBi,LBj,0,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idVmKK,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     MIXING(ng) % tl_Akk(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKK)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idMtls,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     MIXING(ng) % tl_gls(LBi,LBj,0,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idVmLS,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#    endif
     &                     MIXING(ng) % tl_Lscale(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmLS)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    ifdef GSL_MIXING
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncTLMid(ng), tlmVid(idVmKP,ng),    &
     &                     tTLMindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#     endif
     &                     MIXING(ng) % tl_Akp(LBi,LBj,0))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKP)), tTLMindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncTLMid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(TLMname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag, 1)
# endif
# ifdef SOLVE3D
      IF (Master) WRITE (stdout,30) KOUT, NOUT, tTLMindx(ng)
# else
      IF (Master)  WRITE (stdout,30) KOUT, tTLMindx(ng)
# endif
!
  10  FORMAT (/,' TL_WRT_HIS - error while writing variable: ',a,/,14x, &
     &        'into tangent NetCDF file for time record: ',i4)
  20  FORMAT (/,' TL_WRT_HIS - unable to synchronize to disk tangent ', &
     &        'linear file:',/,14x,a)
# ifdef SOLVE3D
  30  FORMAT (3x,'TL_WRT_HIS   - wrote tangent  fields (Index=', i1,    &
     &        ',',i1,') into time record = ',i7.7)
# else
  30  FORMAT (3x,'TL_WRT_HIS   - wrote tangent  fields (Index=', i1,    &
     &        ') into time record = ',i7.7)
# endif
      RETURN
      END SUBROUTINE tl_wrt_his
#else
      SUBROUTINE tl_wrt_his
      RETURN
      END SUBROUTINE tl_wrt_his
#endif
