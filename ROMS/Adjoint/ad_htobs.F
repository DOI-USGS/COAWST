#include "cppdefs.h"
      MODULE ad_htobs_mod
#if defined ADJOINT && defined OBSERVATIONS && defined WEAK_CONSTRAINT
!
!git $Id$
!svn $Id: ad_htobs.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the adjoint observation operator,             !
!                                                                      !
!       transpose(H) * X                                               !
!                                                                      !
!  which loads the observation vector into the adjoint forcing arrays. !
!                                                                      !
!  The observation screening and quality control variable "ObsScale"   !
!  is not modified in this routine. Their values are the ones set in   !
!  obs_write.F during the running of the nonlinear model.              !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_htobs (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_htobs_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef SOLVE3D
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_v,                               &
     &                    OCEAN(ng) % f_u,                              &
     &                    OCEAN(ng) % f_v,                              &
     &                    OCEAN(ng) % f_t,                              &
# endif
     &                    OCEAN(ng) % f_ubar,                           &
     &                    OCEAN(ng) % f_vbar,                           &
     &                    OCEAN(ng) % f_zeta)
      RETURN
      END SUBROUTINE ad_htobs
!
!***********************************************************************
      SUBROUTINE ad_htobs_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef SOLVE3D
     &                          z_r, z_v,                               &
     &                          f_u, f_v, f_t,                          &
# endif
     &                          f_ubar, f_vbar, f_zeta)
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
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d
#  endif
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
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Mstr, Mend, ObsSum, ObsVoid
# ifdef DISTRIBUTE
      integer :: Ncollect
# endif
      integer :: i, ie, iobs, is, j

# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: angle

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8) :: ad_uradial(Mobs), ad_vradial(Mobs)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute model minus observations adjoint misfit forcing. The
!  representer coefficients (or its approximation PSI) has been
!  already loaded into vector ADmodVal in the conjugate gradient
!  or read in.
!
# ifdef DISTRIBUTE
!  The adjoint of this operator is tricky in parallel (tile partitions)
!  because we need to avoid adding observation contributions to ghost
!  points. In the case that an observation is located between neighbor
!  tiles, both tiles need to process it and the contribution to f_var
!  (an adjoint variable) is only done in the tile where (i,j) or
!  (i,j,k) is not a ghost point.
!
!  Alternatively, only one tile process such observation and the
!  ad_mp_exchange*d routine is used to add the contribution to the
!  correct (i,j) or (i,j,k) point. This is the strategy used here.
!
# endif
!  The processing flag used to reject (ObsVetting=0) or accept
!  (ObsVetting=1) observations is computed here but it is never
!  used. The observation screening and quality control variable
!  (ObsScale) is only computed in routine obs_write.
!-----------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
!
!  Set starting and ending indices of representer coefficient vector to
!  proccess. The adjoint forcing is only computed for current time
!  survey observations.
!
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
!
!  Initialize observation reject/accept processing flag.
!
        DO iobs=Mstr,Mend
          ObsVetting(iobs)=IniVal
        END DO

# if defined OBS_SPACE && defined RBL4DVAR_FCT_SENSITIVITY
!
!  Define ADmodVal to be the impact forcing in observation space.
!
#  ifndef OBS_IMPACT
        IF (.not.LadjVAR(ng)) THEN
#  endif
          IF (Lobspace(ng)) THEN
            DO iobs=Mstr,Mend
              ADmodVal(iobs)=ObsVal(iobs)
              ObsScale(iobs)=1.0_r8
            END DO
          ELSE
            DO iobs=Mstr,Mend
              ADmodVal(iobs)=0.0_r8
              ObsScale(iobs)=1.0_r8
            END DO
          END IF
#  ifndef OBS_IMPACT
        END IF
#  endif
# endif

# ifdef BGQC
!
!  Reject observation that fail background quality control check.
!
        DO iobs=Mstr,Mend
          ADmodVal(iobs)=ObsScale(iobs)*ADmodVal(iobs)
        END DO
# endif
!
!  Free-surface.
!
        DO i=LBi,UBi
          DO j=LBj,UBj
            f_zeta(i,j)=0.0_r8
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isFsur).gt.0) THEN
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isFsur),                 &
     &                           Mobs, Mstr, Mend,                      &
     &                           rXmin(ng), rXmax(ng),                  &
     &                           rYmin(ng), rYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs,                      &
     &                           f_zeta,                                &
# ifdef MASKING
     &                           rmask,                                 &
# endif
     &                           ADmodVal)
# ifdef DISTRIBUTE
          CALL ad_mp_exchange2d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           f_zeta)
# endif
        END IF
!
!  2D u-momentum component.
!
        DO i=LBi,UBi
          DO j=LBj,UBj
            f_ubar(i,j)=0.0_r8
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isUbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isUbar),                 &
     &                           Mobs, Mstr, Mend,                      &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs,                      &
     &                           f_ubar,                                &
# ifdef MASKING
     &                           umask,                                 &
# endif
     &                           ADmodVal)
# ifdef DISTRIBUTE
          CALL ad_mp_exchange2d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           f_ubar)
# endif
        END IF
!
!  2D v-momentum component.
!
        DO i=LBi,UBi
          DO j=LBj,UBj
            f_vbar(i,j)=0.0_r8
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isVbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ObsState2Type(isVbar),                 &
     &                           Mobs, Mstr, Mend,                      &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs,                      &
     &                           f_vbar,                                &
# ifdef MASKING
     &                           vmask,                                 &
# endif
     &                           ADmodVal)
# ifdef DISTRIBUTE
          CALL ad_mp_exchange2d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           f_vbar)
# endif
        END IF

# ifdef SOLVE3D
!
!  3D u-momentum component.
!
        DO k=1,N(ng)
          DO i=LBi,UBi
            DO j=LBj,UBj
              f_u(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
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
     &                           Mobs, Mstr, Mend,                      &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           f_u, z_v,                              &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ADmodVal)
#  ifdef DISTRIBUTE
          CALL ad_mp_exchange3d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           f_u)
#  endif
        END IF
!
!  3D v-momentum component.
!
        DO k=1,N(ng)
          DO i=LBi,UBi
            DO j=LBj,UBj
              f_v(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
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
     &                           Mobs, Mstr, Mend,                      &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           f_v, z_v,                              &
#  ifdef MASKING
     &                           vmask,                                 &
#  endif
     &                           ADmodVal)
#  ifdef DISTRIBUTE
          CALL ad_mp_exchange3d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           f_v)
#  endif
        END IF
!
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
!  In the adjoint, u and v are given by:
!
!        f_v = f_v + ADmodVal * SIN(obs_meta - angler)
!        f_u = f_u + ADmodVal * COS(obs_meta - angler)
#  else
!  By default, the radial velocity observations are processed as
!  magnitude and heading angle (obs_meta; radians) in the navigation
!  convention: an azimuth that is clockwise from TRUE North.
!
!  In curvilinear coordinates, the radial forward problem is:
!
!     radial = u * SIN(obs_meta + angler) + v * COS(obs_meta + angler)
!
!  In the adjoint, u and v are given by:
!
!        f_v = f_v + ADmodVal * COS(obs_meta + angler)
!        f_u = f_u + ADmodVal * SIN(obs_meta + angler)
#  endif
!
        IF (FOURDVAR(ng)%ObsCount(isRadial).gt.0) THEN
          DO iobs=Mstr,Mend
            ad_uradial(iobs)=IniVal
            ad_vradial(iobs)=IniVal
          END DO
          DO iobs=Mstr,Mend
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
              ADmodVal(iobs)=0.0_r8
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
     &                           Mobs, Mstr, Mend,                      &
     &                           vXmin(ng), vXmax(ng),                  &
     &                           vYmin(ng), vYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           f_v, z_v,                              &
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
     &                           Mobs, Mstr, Mend,                      &
     &                           uXmin(ng), uXmax(ng),                  &
     &                           uYmin(ng), uYmax(ng),                  &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsVetting,                   &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           f_u, z_v,                              &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ad_uradial)
        END IF

#  ifdef DISTRIBUTE
!
!  Exchange adjoint velocites forcing terms, after all the 3D velocity
!  observations are processed.
!
        IF ((FOURDVAR(ng)%ObsCount(isUvel).gt.0).or.                    &
     &      (FOURDVAR(ng)%ObsCount(isVvel).gt.0).or.                    &
     &      (FOURDVAR(ng)%ObsCount(isRadial).gt.0)) THEN
          CALL ad_mp_exchange3d (ng, tile, iADM, 2,                     &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints, EWperiodic(ng),          &
     &                           NSperiodic(ng), f_u, f_v)
        END IF
#  endif
!
!  Tracer type variables.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=LBi,UBi
              DO j=LBj,UBj
                f_t(i,j,k,itrc)=0.0_r8
              END DO
            END DO
          END DO
          IF (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0) THEN
            CALL ad_extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,        &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ObsState2Type(isTvar(itrc)),         &
     &                             Mobs, Mstr, Mend,                    &
     &                             rXmin(ng), rXmax(ng),                &
     &                             rYmin(ng), rYmax(ng),                &
     &                             time(ng), dt(ng),                    &
     &                             ObsType, ObsVetting,                 &
     &                             Tobs, Xobs, Yobs, Zobs,              &
     &                             f_t(:,:,:,itrc), z_r,                &
#  ifdef MASKING
     &                             rmask,                               &
#  endif
     &                             ADmodVal)
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange3d (ng, tile, iADM, 1,                   &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             NghostPoints,                        &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             f_t(:,:,:,itrc))
#  endif

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
        Ncollect=Mend-Mstr+1
        CALL mp_collect (ng, model, Ncollect, IniVal,                   &
     &                   ObsVetting(Mstr:))
# endif
!
!-----------------------------------------------------------------------
!  Set counters for the number of rejected observations for each state
!  variable. Although unnecessary, the counters are recomputed here to
!  check if "ObsScale" changed from its initial values.
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          IF (ObsScale(iobs).lt.1.0) THEN
            IF  (ObsType(iobs).eq.ObsState2Type(isFsur)) THEN
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
  30        FORMAT (3x,' AD_HTOBS    - Computed adjoint observations ', &
     &              'forcing,',t68,a,/,19x,'(Observation ',      &
     &              'records = ',i7.7,' - ',i7.7,', iic = ',i7.7,')')
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE ad_htobs_tile
#endif
      END MODULE ad_htobs_mod
