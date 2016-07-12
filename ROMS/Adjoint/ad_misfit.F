#include "cppdefs.h"
      MODULE ad_misfit_mod
#if  defined ADJOINT && defined FOUR_DVAR  && defined OBSERVATIONS && \
    !defined WEAK_CONSTRAINT
!
!svn $Id: ad_misfit.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the model minus observations adjoint misfit   !
!  forcing for each state variable.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_misfit (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_misfit_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
# ifdef SOLVE3D
     &                     nstp(ng),                                    &
# endif
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef SOLVE3D
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_v,                              &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
     &                     OCEAN(ng) % ad_t,                            &
# endif
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
     &                     OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE ad_misfit
!
!***********************************************************************
      SUBROUTINE ad_misfit_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
# ifdef SOLVE3D
     &                           nstp,                                  &
# endif
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef SOLVE3D
     &                           z_r, z_v,                              &
     &                           ad_u, ad_v, ad_t,                      &
# endif
     &                           ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
      USE ad_extract_obs_mod, ONLY : ad_extract_obs2d
# ifdef SOLVE3D
      USE ad_extract_obs_mod, ONLY : ad_extract_obs3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp
# ifdef SOLVE3D
      integer, intent(in) :: nstp
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: ObsSum, ObsVoid
      integer :: i, ie, iobs, is, kfrc
# ifdef SOLVE3D
      integer :: itrc, j, k, nfrc
# endif
      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute model minus observations adjoint misfit forcing.
# ifdef DISTRIBUTE
!  Notice that ad_mp_exchange*d calls are not required here
!  (see weak constraint routine ad_htobs.F) because we are
!  using the adjoint state arrays directly which will be
!  exchanged during the initialization step.
# endif
!-----------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
# ifndef BGQC
!
!  Set observation scale (ObsScale). The scale factor is used for
!  screenning  or normalization of the observations. If only screening,
!  the scale is one for good observations and zero for bad observations.
!  If normalizing, the scale is terms of screenning and some specified
!  norm. In incremental 4DVAR only screenning is carried out.
!
        DO iobs=1,Mobs
          ObsScale(iobs)=IniVal
        END DO
# endif
!
!  Compute adjoint forcing terms at observation locations.
!
        DO iobs=1,Nobs(ng)
          ADmodVal(iobs)=ObsErr(iobs)*                                  &
     &                   (NLmodVal(iobs)+TLmodVal(iobs)-ObsVal(iobs))
        END DO

# ifdef BGQC
!
!  Reject observations that fail background quality control check.
!
        DO iobs=1,Nobs(ng)
          ADmodVal(iobs)=ObsScale(iobs)*ADmodVal(iobs)
        END DO
# endif
!
!  Set adjoint time index to update.
!
        kfrc=kstp
# ifdef SOLVE3D
        nfrc=nstp
# endif
!
!  Free-surface.
!
        IF (FOURDVAR(ng)%ObsCount(isFsur).gt.0) THEN
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           isFsur,                                &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           rXmin(ng), rXmax(ng),                  &
     &                           rYmin(ng), rYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs,                      &
     &                           ad_zeta(:,:,kfrc),                     &
# ifdef MASKING
     &                           rmask,                                 &
# endif
     &                           ADmodVal)
        END IF
!
!  2D u-momentum component.
!
        IF (FOURDVAR(ng)%ObsCount(isUbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           isUbar,                                &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs,                      &
     &                           ad_ubar(:,:,kfrc),                     &
# ifdef MASKING
     &                           umask,                                 &
# endif
     &                           ADmodVal)
        END IF
!
!  2D v-momentum component.
!
        IF (FOURDVAR(ng)%ObsCount(isVbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           isVbar,                                &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs,                      &
     &                           ad_vbar(:,:,kfrc),                     &
# ifdef MASKING
     &                           vmask,                                 &
# endif
     &                           ADmodVal)
        END IF
# ifdef SOLVE3D
!
!  3D u-momentum component.
!
        IF (FOURDVAR(ng)%ObsCount(isUvel).gt.0) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                z_v(i,j,k)=0.5_r8*(z_r(i-1,j,k)+                        &
     &                             z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL ad_extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           isUvel,                                &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           ad_u(:,:,:,nfrc),                      &
     &                           z_v,                                   &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ADmodVal)
        END IF
!
!  3D v-momentum component.
!
        IF (FOURDVAR(ng)%ObsCount(isVvel).gt.0) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                z_v(i,j,k)=0.5_r8*(z_r(i,j-1,k)+                        &
     &                             z_r(i,j  ,k))
              END DO
            END DO
          END DO
          CALL ad_extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           isVvel,                                &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           ad_v(:,:,:,nfrc),                      &
     &                           z_v,                                   &
#  ifdef MASKING
     &                           vmask,                                 &
#  endif
     &                           ADmodVal)
        END IF
!
!  Tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0) THEN
            CALL ad_extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,        &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             isTvar(itrc),                        &
     &                             Mobs, 1, Nobs(ng),                   &
     &                             rXmin(ng), rXmax(ng),                &
     &                             rYmin(ng), rYmax(ng),                &
     &                             time(ng), dt(ng),                    &
     &                             ObsType, ObsScale,                   &
     &                             Tobs, Xobs, Yobs, Zobs,              &
     &                             ad_t(:,:,:,nfrc,itrc),               &
     &                             z_r,                                 &
#  ifdef MASKING
     &                             rmask,                               &
#  endif
     &                             ADmodVal)
          END IF
        END DO
# endif
# if defined DISTRIBUTE && !defined BGQC
!
!-----------------------------------------------------------------------
!  Collect all observations scaling factors.
!-----------------------------------------------------------------------
!
        CALL mp_collect (ng, model, Mobs, IniVal, ObsScale)
# endif
!
!-----------------------------------------------------------------------
!  Set counters for number of rejected observations for each state
!  variable.
!-----------------------------------------------------------------------
!
        DO iobs=1,Nobs(ng)
          IF (ObsScale(iobs).lt.1.0) THEN
            IF  (ObsType(iobs).eq.isFsur) THEN
              FOURDVAR(ng)%ObsReject(isFsur)=                           &
     &                              FOURDVAR(ng)%ObsReject(isFsur)+1
            ELSE IF (ObsType(iobs).eq.isUbar) THEN
              FOURDVAR(ng)%ObsReject(isUbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUbar)+1
            ELSE IF (ObsType(iobs).eq.isVbar) THEN
              FOURDVAR(ng)%ObsReject(isVbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVbar)+1
# ifdef SOLVE3D
            ELSE IF (ObsType(iobs).eq.isUvel) THEN
              FOURDVAR(ng)%ObsReject(isUvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUvel)+1
            ELSE IF (ObsType(iobs).eq.isVvel) THEN
              FOURDVAR(ng)%ObsReject(isVvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVvel)+1
            ELSE
              DO itrc=1,NT(ng)
                IF (ObsType(iobs).eq.isTvar(itrc)) THEN
                  i=isTvar(itrc)
                  FOURDVAR(ng)%ObsReject(i)=FOURDVAR(ng)%ObsReject(i)+1
                END IF
              END DO
# endif
            END IF
          END IF
        END DO
!
!  Load total available and rejected observations into structure
!  array.
!
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%ObsCount(0)=FOURDVAR(ng)%ObsCount(0)+            &
     &                             FOURDVAR(ng)%ObsCount(i)
          FOURDVAR(ng)%ObsReject(0)=FOURDVAR(ng)%ObsReject(0)+          &
     &                              FOURDVAR(ng)%ObsReject(i)
        END DO
!
!  Report.
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            ObsSum=0
            ObsVoid=0
            is=NstrObs(ng)
            DO i=1,NstateVar(ng)
              IF (FOURDVAR(ng)%ObsCount(i).gt.0) THEN
                ie=is+FOURDVAR(ng)%ObsCount(i)-1
                WRITE (stdout,10) TRIM(Vname(1,idSvar(i))), is, ie,     &
     &                            ie-is+1, FOURDVAR(ng)%ObsReject(i)
                is=ie+1
                ObsSum=ObsSum+FOURDVAR(ng)%ObsCount(i)
                ObsVoid=ObsVoid+FOURDVAR(ng)%ObsReject(i)
              END IF
            END DO
            WRITE (stdout,20) ObsSum, ObsVoid,                          &
     &                        FOURDVAR(ng)%ObsCount(0),                 &
     &                        FOURDVAR(ng)%ObsReject(0)
            WRITE (stdout,30) time_code(ng), NstrObs(ng), NendObs(ng),  &
     &                        iic(ng)
  10        FORMAT (10x,a,t25,4(1x,i10))
  20        FORMAT (/,10x,'Total',t47,2(1x,i10),                        &
     &              /,10x,'Obs Tally',t47,2(1x,i10),/)
  30        FORMAT (3x,' AD_MISFIT   - Added observations misfit ',     &
     &              'forcing,',t62,'t = ',a,/,19x,'(Observation ',      &
     &              'records = ',i7.7,' - ',i7.7,', iic = ',i7.7,')')
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE ad_misfit_tile
#endif
      END MODULE ad_misfit_mod
