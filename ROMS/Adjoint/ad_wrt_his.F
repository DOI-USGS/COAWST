#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_wrt_his (ng)
!
!svn $Id: ad_wrt_his.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested adjoint model fields into adjoint     !
!  history NetCDF file.                                                !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
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
      integer :: Fcount, i, j, gfactor, gtype, status
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
# if defined WEAK_CONSTRAINT && !defined TIME_CONV
      kfout=2
# endif
# ifdef SOLVE3D
      IF (iic(ng).ne.ntend(ng)) THEN
        nout=nnew(ng)
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        LwrtState3d(ng)=.FALSE.
# endif
      ELSE
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        LwrtState3d(ng)=.TRUE.
# endif
        nout=nstp(ng)
      END IF
# endif
!
!  Set time record index.
!
      ADM(ng)%Rindex=ADM(ng)%Rindex+1
      Fcount=ADM(ng)%Fcount
      ADM(ng)%Nrec(Fcount)=ADM(ng)%Nrec(Fcount)+1
# if defined WEAK_CONSTRAINT && defined TIME_CONV
      kfout=NrecTC(ng)-ADM(ng)%Rindex+1
# endif
!
!  If requested, set time index to recycle time records in the adjoint
!  file.
!
      IF (LcycleADJ(ng)) THEN
        ADM(ng)%Rindex=MOD(ADM(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtTime(ng)) THEN
        IF (LwrtPER(ng)) THEN
          Tval(1)=REAL(ADM(ng)%Rindex,r8)*day2sec
        ELSE
# ifdef WEAK_CONSTRAINT
#  ifdef TIME_CONV
          IF (ADM(ng)%Rindex.le.NrecTC(ng)) THEN
            Tval(1)=dstart*day2sec+                                     &
     &              REAL(kfout-1,r8)*nADJ(ng)*dt(ng)
          ELSE
            Tval(1)=dstart*day2sec
          END IF
#  else
          Tval(1)=ForceTime(ng)
#  endif
# else
          Tval(1)=time(ng)
# endif
        END IF
        CALL netcdf_put_fvar (ng, iADM, ADM(ng)%name,                   &
     &                        TRIM(Vname(1,idtime)), tval,              &
     &                        (/ADM(ng)%Rindex/), (/1/),                &
     &                        ncid = ADM(ng)%ncid,                      &
     &                        varid = ADM(ng)%Vid(idtime))
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
      status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idUsms),   &
     &                   ADM(ng)%Rindex, gtype,                         &
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
      status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVsms),   &
     &                   ADM(ng)%Rindex, gtype,                         &
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
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idTsur(itrc)),                 &
     &                       ADM(ng)%Rindex, gtype,                     &
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
      IF (Hout(idbath,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idbath), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng)% ad_h,                              &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
# ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idFsur),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
#  ifdef TIME_CONV
     &                       OCEAN(ng)% f_zetaS(:,:,kfout))
#  else
     &                       OCEAN(ng)% f_zetaG(:,:,kfout))
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), ADM(ng)%Rindex
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
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idFsur),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
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
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idFsur),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
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
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), ADM(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, ADM(ng)%name, ADM(ng)%ncid,   &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          ADM(ng)%Vid(idSbry(isFsur)),            &
     &                          ADM(ng)%Rindex, r2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,   &
     &                                                     Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        ADM(ng)%Rindex
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
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idUbar),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#  endif
#  ifdef TIME_CONV
     &                       OCEAN(ng) % f_ubarS(:,:,kfout))
#  else
     &                       OCEAN(ng) % f_ubarG(:,:,kfout))
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
# endif
          scale=1.0_r8
          gtype=gfactor*u2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idUbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
# endif
     &                         OCEAN(ng) % ad_ubar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idUbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
# endif
     &                         OCEAN(ng) % ad_ubar_sol)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), ADM(ng)%Rindex
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
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iADM, ADM(ng)%name, ADM(ng)%ncid,   &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          ADM(ng)%Vid(idSbry(isUbar)),            &
     &                          ADM(ng)%Rindex, u2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,   &
     &                                                     Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        ADM(ng)%Rindex
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
        IF (WRTforce(ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idVbar),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#  endif
#  ifdef TIME_CONV
     &                       OCEAN(ng) % f_vbarS(:,:,kfout))
#  else
     &                       OCEAN(ng) % f_vbarG(:,:,kfout))
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
# endif
          scale=1.0_r8
          gtype=gfactor*v2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idVbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         OCEAN(ng) % ad_vbar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idVbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         OCEAN(ng) % ad_vbar_sol)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), ADM(ng)%Rindex
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
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iADM, ADM(ng)%name, ADM(ng)%ncid,   &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          ADM(ng)%Vid(idSbry(isVbar)),            &
     &                          ADM(ng)%Rindex, v2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,   &
     &                                                     Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        ADM(ng)%Rindex
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
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idUvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#   endif
#  ifdef TIME_CONV
     &                       OCEAN(ng) % f_uS(:,:,:,kfout))
#  else
     &                       OCEAN(ng) % f_uG(:,:,:,kfout))
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          scale=1.0_r8
          gtype=gfactor*u3dvar
#  if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
          IF (LwrtState3d(ng)) THEN
#  endif
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idUvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#    endif
     &                       OCEAN(ng) % ad_u(:,:,:,nout))
#  if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
          ELSE
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idUvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#    endif
     &                       OCEAN(ng) % ad_u_sol(:,:,:))
          ENDIF
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), ADM(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, iADM, ADM(ng)%name, ADM(ng)%ncid,   &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          ADM(ng)%Vid(idSbry(isUvel)),            &
     &                          ADM(ng)%Rindex, u3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,    &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        ADM(ng)%Rindex
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
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idVvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
#  ifdef TIME_CONV
     &                       OCEAN(ng) % f_vS(:,:,:,kfout))
#  else
     &                       OCEAN(ng) % f_vG(:,:,:,kfout))
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          scale=1.0_r8
          gtype=gfactor*v3dvar
#  if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
          IF (LwrtState3d(ng)) THEN
#  endif
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idVvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       OCEAN(ng) % ad_v(:,:,:,nout))
#  if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
          ELSE
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idVvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       OCEAN(ng) % ad_v_sol(:,:,:))
          END IF
#  endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), ADM(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, iADM, ADM(ng)%name, ADM(ng)%ncid,   &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          ADM(ng)%Vid(idSbry(isVvel)),            &
     &                          ADM(ng)%Rindex, v3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,    &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        ADM(ng)%Rindex
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
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Tid(itrc),                       &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
#  ifdef TIME_CONV
     &                         OCEAN(ng) % f_tS(:,:,:,kfout,itrc))
#  else
     &                         OCEAN(ng) % f_tG(:,:,:,kfout,itrc))
#  endif
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            ADM(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
#  endif
            scale=1.0_r8
            gtype=gfactor*r3dvar
#  if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
            IF (LwrtState3d(ng)) THEN
#  endif
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Tid(itrc),                       &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         OCEAN(ng) % ad_t(:,:,:,nout,itrc))
#  if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
            ELSE
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Tid(itrc),                       &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         OCEAN(ng) % ad_t_sol(:,:,:,itrc))
            END IF
#  endif
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            ADM(ng)%Rindex
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
          status=nf_fwrite3d_bry (ng, iADM, ADM(ng)%name, ADM(ng)%ncid, &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            ADM(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                            ADM(ng)%Rindex, r3dvar,               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,  &
     &                                                 Lbout(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          ADM(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idDano), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % ad_rho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), ADM(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVvis), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % ad_Akv,                         &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), ADM(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idTdif), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % ad_Akt(:,:,:,itemp),            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), ADM(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idSdif), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % ad_Akt(:,:,:,isalt),            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), ADM(ng)%Rindex
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
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idTsur(itrc)),                 &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          ADM(ng)%Rindex
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
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idUsms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % ad_sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), ADM(ng)%Rindex
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
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVsms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % ad_svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), ADM(ng)%Rindex
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
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idUbms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % ad_bustr_sol)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), ADM(ng)%Rindex
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
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVbms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % ad_bvstr_sol)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), ADM(ng)%Rindex
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
      CALL netcdf_sync (ng, iADM, ADM(ng)%name, ADM(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      IF (Master) WRITE (stdout,20) kout, nout, ADM(ng)%Rindex
# else
      IF (Master) WRITE (stdout,20) kout, ADM(ng)%Rindex
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
