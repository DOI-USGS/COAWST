#include "cppdefs.h"
      MODULE ad_misfit_mod

#if  defined ADJOINT && defined FOUR_DVAR  && defined OBSERVATIONS && \
    !defined WEAK_CONSTRAINT
!
!git $Id$
!svn $Id: ad_misfit.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the model minus observations adjoint misfit   !
!  forcing for each state variable.                                    !
!                                                                      !
!  The observation screening and quality control variable "ObsScale"   !
!  is not modified in this routine. Their values are the ones set in   !
!  obs_write.F during the running of the nonlinear model.              !
!                                                                      !
!=======================================================================
!
      implicit none
!
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
      integer :: i, ic, ie, iobs, is, kfrc
# ifdef SOLVE3D
      integer :: itrc, j, k, nfrc
# endif
      real(r8) :: angle

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8) :: ad_uradial(Mobs), ad_vradial(Mobs)

# include "set_bounds.h"
!
!=======================================================================
!  Compute model minus observations adjoint misfit forcing.
# ifdef DISTRIBUTE
!  Notice that ad_mp_exchange*d calls are not required here
!  (see weak constraint routine ad_htobs.F) because we are
!  using the adjoint state arrays directly which will be
!  exchanged during the initialization step.
# endif
!
!  The processing flag used to reject (ObsVetting=0) or accept
!  (ObsVetting=1) observations is computed here but it is never
!  used. The observation screening and quality control variable
!  (ObsScale) is only computed in routine obs_write.
!=======================================================================
!
      IF (ProcessObs(ng)) THEN
!
!  Set adjoint time index to update.
!
        kfrc=kstp
# ifdef SOLVE3D
        nfrc=nstp
# endif
!
!  Initialize observation reject/accept processing flag.
!
        DO iobs=1,Mobs
          ObsVetting(iobs)=IniVal
        END DO
!
!  Compute adjoint forcing terms at observation locations.
!
        DO iobs=1,Nobs(ng)
          ADmodVal(iobs)=ObsErr(iobs)*                                  &
     &                   (NLmodVal(iobs)+TLmodVal(iobs)-ObsVal(iobs))
        END DO
!
!-----------------------------------------------------------------------
!  Load observation reject and accept flag into screening variable.
!-----------------------------------------------------------------------
!
!  In the primal formulation, the current time window (or survey)
!  observations are loaded to the working arrays using local indexes
!  (array elements 1:Nobs) as opposed to global indexes in the dual
!  formulation. Recall that the screening variable ObsScale is computed
!  only once to facilitate Background Quality Control on the first pass
!  of the NLM. Therefore, we need to load its values from the saved
!  global arrays values, ObsScaleGlobal.
!
        ic=0
        DO iobs=NstrObs(ng),NendObs(ng)
          ic=ic+1
          ObsScale(ic)=ObsScaleGlobal(iobs)
        END DO

# ifdef BGQC
!
!-----------------------------------------------------------------------
!  Reject observations that fail background quality control check.
!-----------------------------------------------------------------------
!
        DO iobs=1,Nobs(ng)
          ADmodVal(iobs)=ObsScale(iobs)*ADmodVal(iobs)
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Free-surface.
!-----------------------------------------------------------------------
!
        IF (FOURDVAR(ng)%ObsCount(isFsur).gt.0) THEN
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isFsur),                 &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           rXmin(ng), rXmax(ng),                  &
     &                           rYmin(ng), rYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs,                      &
     &                           ad_zeta(:,:,kfrc),                     &
# ifdef MASKING
     &                           rmask,                                 &
# endif
     &                           ADmodVal)
        END IF
!
!-----------------------------------------------------------------------
!  2D u-momentum component.
!-----------------------------------------------------------------------
!
        IF (FOURDVAR(ng)%ObsCount(isUbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isUbar),                 &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs,                      &
     &                           ad_ubar(:,:,kfrc),                     &
# ifdef MASKING
     &                           umask,                                 &
# endif
     &                           ADmodVal)
        END IF
!
!-----------------------------------------------------------------------
!  2D v-momentum component.
!-----------------------------------------------------------------------
!
        IF (FOURDVAR(ng)%ObsCount(isVbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isVbar),                 &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs,                      &
     &                           ad_vbar(:,:,kfrc),                     &
# ifdef MASKING
     &                           vmask,                                 &
# endif
     &                           ADmodVal)
        END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D u-momentum component.
!-----------------------------------------------------------------------
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
     &                           ObsState2Type(isUvel),                 &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           ad_u(:,:,:,nfrc),                      &
     &                           z_v,                                   &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ADmodVal)
        END IF
!
!-----------------------------------------------------------------------
!  3D v-momentum component.
!-----------------------------------------------------------------------
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
     &                           ObsState2Type(isVvel),                 &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           ad_v(:,:,:,nfrc),                      &
     &                           z_v,                                   &
#  ifdef MASKING
     &                           vmask,                                 &
#  endif
     &                           ADmodVal)
        END IF
!
!-----------------------------------------------------------------------
!  Radial Velocity.  The observations are in terms of radial speed and
!  angle (stored in obs_meta).  The observation angle converts the
!  velocity components to geographical EAST and North components.
#  ifdef RADIAL_ANGLE_CCW_EAST
!  The radial velocity observations are processed as magnitude and
!  heading angle (obs_meta; radians) in the math convention: an
!  azimuth that is counterclockwise from TRUE East.
!
!  In curvilinear coordinates, the radial forward problem is:
!
!     radial = u * COS(obs_meta - angler) + v * SIN(obs_meta - angler)
!
!  In the adjoint, we get:
!
!        ad_vradial = ad_vradial + ADmodVal * SIN(obs_meta - angler)
!        ad_uradial = ad_uradial + ADmodVal * COS(obs_meta - angler)
#  else
!  By default, the radial velocity observations are processed as
!  magnitude and heading angle (obs_meta; radians) in the navigation
!  convention: an azimuth that is clockwise from TRUE North.
!
!  In curvilinear coordinates, the radial forward problem is:
!
!     radial = u * SIN(obs_meta + angler) + v * COS(obs_meta + angler)
!
!  In the adjoint, we get:
!
!        ad_vradial = ad_vradial + ADmodVal * COS(obs_meta + angler)
!        ad_uradial = ad_uradial + ADmodVal * SIN(obs_meta + angler)
#  endif
!-----------------------------------------------------------------------
!
        IF (FOURDVAR(ng)%ObsCount(isRadial).gt.0) THEN
          DO iobs=1,Nobs(ng)
            ad_uradial(iobs)=IniVal
            ad_vradial(iobs)=IniVal
          END DO
          DO iobs=1,Nobs(ng)
            IF (ObsType(iobs).eq.ObsState2Type(isRadial)) THEN
#  ifdef RADIAL_ANGLE_CCW_EAST
#   ifdef CURVGRID
              angle=ObsMeta(iobs)-ObsAngler(iobs)
              ad_uradial(iobs)=ad_uradial(iobs)+                        &
     &                         ADmodVal(iobs)*COS(angle)
              ad_vradial(iobs)=ad_vradial(iobs)+                        &
     &                         ADmodVal(iobs)*SIN(angle)
#   else
              ad_uradial(iobs)=ad_uradial(iobs)+                        &
     &                         ADmodVal(iobs)*COS(ObsMeta(iobs))
              ad_vradial(iobs)=ad_vradial(iobs)+                        &
     &                         ADmodVal(iobs)*SIN(ObsMeta(iobs))
#   endif
#  else
#   ifdef CURVGRID
              angle=ObsMeta(iobs)+ObsAngler(iobs)
              ad_uradial(iobs)=ad_uradial(iobs)+                        &
     &                         ADmodVal(iobs)*SIN(angle)
              ad_vradial(iobs)=ad_vradial(iobs)+                        &
     &                         ADmodVal(iobs)*COS(angle)
#   else
              ad_uradial(iobs)=ad_uradial(iobs)+                        &
     &                         ADmodVal(iobs)*SIN(ObsMeta(iobs))
              ad_vradial(iobs)=ad_vradial(iobs)+                        &
     &                         ADmodVal(iobs)*COS(ObsMeta(iobs))
#   endif
#  endif
            END IF
          END DO
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
     &                           ObsState2Type(isRadial),               &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           ad_v(:,:,:,nfrc), z_v,                 &
#  ifdef MASKING
     &                           vmask,                                 &
#  endif
     &                           ad_vradial)
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
     &                           ObsState2Type(isRadial),               &
     &                           Mobs, 1, Nobs(ng),                     &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           ad_u(:,:,:,nfrc), z_v,                 &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ad_uradial)
        END IF
!
!-----------------------------------------------------------------------
!  Tracer type variables.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
          IF (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0) THEN
            CALL ad_extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,        &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ObsState2Type(isTvar(itrc)),         &
     &                             Mobs, 1, Nobs(ng),                   &
     &                             rXmin(ng), rXmax(ng),                &
     &                             rYmin(ng), rYmax(ng),                &
     &                             time(ng), dt(ng),                    &
     &                             ObsType, ObsVetting,                 &
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
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  For debugging purposes, collect all observations reject/accept
!  processing flag.
!-----------------------------------------------------------------------
!
        CALL mp_collect (ng, model, Mobs, IniVal, ObsVetting)
# endif
!
!-----------------------------------------------------------------------
!  Set counters for the number of rejected observations for each state
!  variable. Although unnecessary, the counters are recomputed here to
!  check if "ObsScale" changed from its initial values.
!-----------------------------------------------------------------------
!
        DO iobs=1,Nobs(ng)
          IF (ObsScale(iobs).lt.1.0) THEN
            IF (ObsType(iobs).eq.ObsState2Type(isFsur)) THEN
              FOURDVAR(ng)%ObsReject(isFsur)=                           &
     &                              FOURDVAR(ng)%ObsReject(isFsur)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isUbar)) THEN
              FOURDVAR(ng)%ObsReject(isUbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUbar)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isVbar)) THEN
              FOURDVAR(ng)%ObsReject(isVbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVbar)+1
# ifdef SOLVE3D
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isUvel)) THEN
              FOURDVAR(ng)%ObsReject(isUvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUvel)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isVvel)) THEN
              FOURDVAR(ng)%ObsReject(isVvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVvel)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isRadial)) THEN
              FOURDVAR(ng)%ObsReject(isRadial)=                         &
     &                              FOURDVAR(ng)%ObsReject(isRadial)+1
            ELSE
              DO itrc=1,NT(ng)
                IF (ObsType(iobs).eq.ObsState2Type(isTvar(itrc))) THEN
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
        DO i=1,NobsVar(ng)
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
            DO i=1,NobsVar(ng)
              IF (FOURDVAR(ng)%ObsCount(i).gt.0) THEN
                ie=is+FOURDVAR(ng)%ObsCount(i)-1
                WRITE (stdout,10) TRIM(ObsName(i)), is, ie,             &
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
  30        FORMAT (2x,' AD_MISFIT       - Added observations misfit ', &
     &              'forcing,',t75,a,/,22x,'(Observation ',             &
     &              'records = ',i0,' - ',i0,', iic = ',i0,')')
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE ad_misfit_tile
#endif
      END MODULE ad_misfit_mod
