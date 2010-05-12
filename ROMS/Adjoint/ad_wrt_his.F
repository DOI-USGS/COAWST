#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_wrt_his (ng)
!
!svn $Id: ad_wrt_his.F 431 2009-12-26 20:36:20Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
      USE mod_forces
# ifdef WEAK_CONSTRAINT
      USE mod_fourdvar
# endif
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
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
# endif
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
#  ifdef ADJUST_BOUNDARY
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
#  endif
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
# ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
# endif
      integer :: i, j, gfactor, gtype, status
      integer :: kout
# ifdef WEAK_CONSTRAINT
      integer :: kfout
# endif
# ifdef SOLVE3D
      integer :: itrc, k, nout
# endif
      real(r8) :: scale
      real(r8) :: Tval(1)
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
# ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
# endif
!
      SourceFile='ad_wrt_his.F'
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
# ifdef WEAK_CONSTRAINT
      kfout=2
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
# ifdef WEAK_CONSTRAINT
          Tval(1)=ForceTime(ng)
# else
          Tval(1)=time(ng)
# endif
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
      scale=1.0_r8                          ! m2/s2
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idUsms,ng),      &
     &                   tADJindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Lfout(ng)))
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
      scale=1.0_r8                          ! m2/s2
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idVsms,ng),      &
     &                   tADJindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Lfout(ng)))
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
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_r8                      ! kinematic flux units
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iADM, ncADJid(ng),                     &
     &                       adjVid(idTsur(itrc),ng),                   &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng)% ad_tflux(:,:,:,Lfout(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))), Lfout(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
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
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng)% ad_h,                                &
     &                   SetFillVal = .FALSE.)
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
# ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idFsur,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng)% f_zetaG(:,:,kfout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
# endif
          scale=1.0_r8
          gtype=gfactor*r2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idFsur,ng),&
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
# ifdef WET_DRY
     &                         OCEAN(ng)% ad_zeta(:,:,kout),            &
     &                         SetFillVal = .FALSE.)
# else
     &                         OCEAN(ng)% ad_zeta(:,:,kout))
# endif
          ELSE
            status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idFsur,ng),&
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
# ifdef WET_DRY
     &                         OCEAN(ng)% ad_zeta_sol,                  &
     &                         SetFillVal = .FALSE.)
# else
     &                         OCEAN(ng)% ad_zeta_sol)
# endif
          ENDIF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
# ifdef WEAK_CONSTRAINT
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iADM, ADJname(ng), ncADJid(ng),     &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          adjVid(idSbry(isFsur),ng),              &
     &                          tADJindx(ng), r2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,   &
     &                                                     Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
# ifdef WEAK_CONSTRAINT
#  ifndef SOLVE3D
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUbar,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       OCEAN(ng) % f_ubarG(:,:,kfout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
# endif
          scale=1.0_r8
          gtype=gfactor*u2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUbar,ng),&
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % umask,                        &
# endif
     &                         OCEAN(ng) % ad_ubar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idUbar,ng),&
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % umask,                        &
# endif
     &                         OCEAN(ng) % ad_ubar_sol)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
# ifdef WEAK_CONSTRAINT
#  ifndef SOLVE3D
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iADM, ADJname(ng), ncADJid(ng),     &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          adjVid(idSbry(isUbar),ng),              &
     &                          tADJindx(ng), u2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,   &
     &                                                     Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
# ifdef WEAK_CONSTRAINT
#  ifndef SOLVE3D
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVbar,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
      &                      GRID(ng) % vmask,                          &
#   endif
     &                       OCEAN(ng) % f_vbarG(:,:,kfout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
# endif
          scale=1.0_r8
          gtype=gfactor*v2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVbar,ng),&
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
      &                        GRID(ng) % vmask,                        &
# endif
     &                         OCEAN(ng) % ad_vbar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ncADJid(ng), adjVid(idVbar,ng),&
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
      &                        GRID(ng) % vmask,                        &
# endif
     &                         OCEAN(ng) % ad_vbar_sol)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
# ifdef WEAK_CONSTRAINT
#  ifndef SOLVE3D
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iADM, ADJname(ng), ncADJid(ng),     &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          adjVid(idSbry(isVbar),ng),              &
     &                          tADJindx(ng), v2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,   &
     &                                                     Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idUvel,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       OCEAN(ng) % f_uG(:,:,:,kfout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          scale=1.0_r8
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idUvel,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       OCEAN(ng) % ad_u(:,:,:,nout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iADM, ADJname(ng), ncADJid(ng),     &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          adjVid(idSbry(isUvel),ng),              &
     &                          tADJindx(ng), u3dvar,                   &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,    &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idVvel,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       OCEAN(ng) % f_vG(:,:,:,kfout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          scale=1.0_r8
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idVvel,ng),  &
     &                       tADJindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       OCEAN(ng) % ad_v(:,:,:,nout))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), tADJindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iADM, ADJname(ng), ncADJid(ng),     &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          adjVid(idSbry(isVvel),ng),              &
     &                          tADJindx(ng), v3dvar,                   &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,    &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        tADJindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
#  ifdef WEAK_CONSTRAINT
          IF (WRTforce(ng)) THEN
            scale=1.0_r8
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjTid(itrc,ng),  &
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         OCEAN(ng) % f_tG(:,:,:,kfout,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            tADJindx(ng)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
#  endif
            scale=1.0_r8
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjTid(itrc,ng),  &
     &                         tADJindx(ng), gtype,                     &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         OCEAN(ng) % ad_t(:,:,:,nout,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            tADJindx(ng)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
#  ifdef WEAK_CONSTRAINT
          END IF
#  endif
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_r8
          status=nf_fwrite3d_bry (ng, iADM, ADJname(ng), ncADJid(ng),   &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            adjVid(idSbry(isTvar(itrc)),ng),      &
     &                            tADJindx(ng), r3dvar,                 &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,  &
     &                                                 Lbout(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          tADJindx(ng)
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
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iADM, ncADJid(ng), adjVid(idDano,ng),    &
     &                     tADJindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % ad_rho)
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
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % ad_Akv)
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
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % ad_Akt(:,:,:,itemp))
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
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % ad_Akt(:,:,:,isalt))
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
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
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
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_stflx(:,:,itrc))
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
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
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
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % ad_sustr)
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
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
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
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % ad_svstr)
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
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % ad_bustr_sol)
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
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % ad_bvstr_sol)
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
      CALL netcdf_sync (ng, iADM, ADJname(ng), ncADJid(ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      IF (Master) WRITE (stdout,20) kout, nout, tADJindx(ng)
# else
      IF (Master) WRITE (stdout,20) kout, tADJindx(ng)
# endif
!
  10  FORMAT (/,' AD_WRT_HIS - error while writing variable: ',a,/,14x, &
     &        'into adjoint NetCDF file for time record: ',i4)
# ifdef SOLVE3D
  20  FORMAT (3x,'AD_WRT_HIS   - wrote adjoint  fields (Index=', i1,    &
     &        ',',i1,') into time record = ',i7.7)
# else
  20  FORMAT (3x,'AD_WRT_HIS   - wrote adjoint  fields (Index=', i1,    &
     &        ') into time record = ',i7.7)
# endif
      RETURN
      END SUBROUTINE ad_wrt_his
#else
      SUBROUTINE ad_wrt_his
      RETURN
      END SUBROUTINE ad_wrt_his
#endif
