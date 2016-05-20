#include "cppdefs.h"
#undef NEW_CODE
      MODULE hmixing_mod

#if (defined DIFF_3DCOEF || defined VISC_3DCOEF) && defined SOLVE3D
!
!svn $Id: hmixing.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes time-dependent 3D horizontal mixing           !
!  coefficients.                                                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Smagorinsky, J, 1963: General circulation experiments with        !
!      the primitive equations: I. The basic experiment, Mon.          !
!      Wea. Rev., 91, 99-164.                                          !
!                                                                      !
!    Holland, W.R., J.C. Chow, and F.O. Bryan, 1998: Application       !
!      of a Third-Order Upwind Scheme in the NCAR Ocean Model, J.      !
!      Climate, 11, 1487-1493.                                         !
!                                                                      !
!    Webb, D.J., B.A. De Cuevas, and C.S. Richmond, 1998: Improved     !
!      Advection Schemes for Ocean Models, J. Atmos. Oceanic           !
!      Technol., 15, 1171-1187.                                        !
!                                                                      !
!    Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic Friction       !
!      with a Smagorinsky-like Viscosity for Use in Large-Scale        !
!      Eddy-Permitting Ocean Models, Monthly Weather Review, 128,      !
!      2935-2946.                                                      !
!                                                                      !
!    Marchesiello, P., L. Debreu, and Xavien Couvelard, 2008:          !
!      Spurious diapycnal mixing in terrain-following coordinate       !
!      models" advection problem and solutions, DRAFT.                 !
!                                                                      !
!  This routine was adapted from a routine provided by Patrick         !
!  Marchiesello (April 2008).                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: hmixing

      CONTAINS
!
!***********************************************************************
      SUBROUTINE hmixing (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
# if defined VEGETATION && defined VEG_HMIXING
      USE mod_vegarr
      USE vegetation_hmixing_mod, ONLY:vegetation_hmixing_cal
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
!
# if defined VEGETATION && defined VEG_HMIXING
!
!-----------------------------------------------------------------------
!  Add the effect of vegetation on horizontal viscosity.
!-----------------------------------------------------------------------
!
      CALL vegetation_hmixing_cal (ng, tile)
# endif
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 28)
# endif
      CALL hmixing_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
# endif
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   GRID(ng) % omn,                                &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
# ifdef DIFF_3DCOEF
     &                   MIXING(ng) % Hdiffusion,                       &
#  ifdef TS_U3ADV_SPLIT
     &                   MIXING(ng) % diff3d_u,                         &
     &                   MIXING(ng) % diff3d_v,                         &
#  else
     &                   MIXING(ng) % diff3d_r,                         &
#  endif
# endif
# ifdef VISC_3DCOEF
     &                   MIXING(ng) % Hviscosity,                       &
#  ifdef UV_U3ADV_SPLIT
     &                   MIXING(ng) % Uvis3d_r,                         &
     &                   MIXING(ng) % Vvis3d_r,                         &
#  else
     &                   MIXING(ng) % visc3d_r,                         &
#  endif
# endif
# if defined VEGETATION && defined VEG_HMIXING
     &                   VEG(ng) % visc3d_r_veg,                        &
# endif
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 28)
# endif
      RETURN
      END SUBROUTINE hmixing
!
!***********************************************************************
      SUBROUTINE hmixing_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
     &                         pm, pn, omn, om_u, on_v,                 &
     &                         Hz, z_r,                                 &
# ifdef DIFF_3DCOEF
     &                         Hdiffusion,                              &
#  ifdef TS_U3ADV_SPLIT
     &                         diff3d_u, diff3d_v,                      &
#  else
     &                         diff3d_r,                                &
#  endif
# endif
# ifdef VISC_3DCOEF
     &                         Hviscosity,                              &
#  ifdef UV_U3ADV_SPLIT
     &                         Uvis3d_r, Vvis3d_r,                      &
#  else
     &                         visc3d_r,                                &
#  endif
# endif
# if defined VEGETATION && defined VEG_HMIXING
     &                         visc3d_r_veg,                            &
# endif
     &                         u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef DIFF_3DCOEF
      real(r8), intent(in) :: Hdiffusion(LBi:,LBj:)
#  endif
#  ifdef VISC_3DCOEF
      real(r8), intent(in) :: Hviscosity(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
      real(r8), intent(out) :: diff3d_u(LBi:,LBj:,:)
      real(r8), intent(out) :: diff3d_v(LBi:,LBj:,:)
#   else
      real(r8), intent(out) :: diff3d_r(LBi:,LBj:,:)
#   endif
#  endif
#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
      real(r8), intent(out) :: Uvis3d_r(LBi:,LBj:,:)
      real(r8), intent(out) :: Vvis3d_r(LBi:,LBj:,:)
#   else
      real(r8), intent(out) :: visc3d_r(LBi:,LBj:,:)
#   endif
#  endif
#  if defined VEGETATION && defined VEG_HMIXING 
      real(r8), intent(in) :: visc3d_r_veg(LBi:,LBj:,:)
#  endif
# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIFF_3DCOEF
      real(r8), intent(in) :: Hdiffusion(LBi:UBi,LBj:UBj)
#  endif
#  ifdef VISC_3DCOEF
      real(r8), intent(in) :: Hviscosity(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
      real(r8), intent(out) :: diff3d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: diff3d_v(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(out) :: diff3d_r(LBi:UBi,LBj:UBj,N(ng))
#   endif
#  endif
#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
      real(r8), intent(out) :: Uvis3d_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: Vvis3d_r(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(out) :: visc3d_r(LBi:UBi,LBj:UBj,N(ng))
#   endif
#  endif
#  if defined VEGETATION && defined VEG_HMIXING 
      real(r8), intent(in) :: visc3d_r_veg(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: SmagorCoef = 0.1_r8
      real(r8), parameter :: PecletCoef = 1.0_r8 / 12.0_r8

      real(r8) :: DefRate, cff, clip_diff, clip_scale

# include "set_bounds.h"

# if defined UV_SMAGORINSKY || defined UV_U3ADV_SPLIT || \
     defined TS_SMAGORINSKY
!
!-----------------------------------------------------------------------
!  Compute velocity dependent, Smagorinsky (1963) horizontal
!  viscosity coefficients. This are based on the local deformation
!  rate, DefRate, which has a horizontal tension and horizontal
!  shearing strain terms:
!
!     DefRate = SQRT [ (du/dx)^2 + (dvdy)^2 + 0.5 * (dvdx + dudy)^2 ]
!                            tension                shearing strain
!
!  The harmonic viscosity coefficient is computed as:
!
!     Asmag = SmagorCoef * dx * dy * DefRate
!
!  The biharmonic viscosity coefficient follows Griffies and Hallberg
!  (2000) formulation:
!
!     Bsmag = PecletCoef * (dx * dy)^2 * DefRate
!
#  ifdef TS_SMAGORINSKY
!  It is not clear physically how to compute the Smagorinsky tracer
!  diffusion from momentum viscosity.  The diffusivity is usually
!  computed as a time independent constant which is a function of
!  grid space that it is smaller that the viscosity (Griffies and
!  Hallberg, 2000).
#  endif
!-----------------------------------------------------------------------
!
!  Compute viscosity and clipping scale.
!
      clip_scale=0.01_r8*grdmax(ng)**3

      DO k=1,N(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!
!  Compute local deformation rate at RHO-points.
!
            DefRate=SQRT(((u(i+1,j,k,nrhs)-                             &
                           u(i  ,j,k,nrhs))*pm(i,j))**2+                &
     &                   ((v(i,j+1,k,nrhs)-                             &
     &                     v(i,j  ,k,nrhs))*pn(i,j))**2+                &
     &                   0.5_r8*(0.25_r8*pn(i,j)*                       &
     &                           (u(i  ,j+1,k,nrhs)+                    &
     &                            u(i+1,j+1,k,nrhs)-                    &
     &                            u(i  ,j-1,k,nrhs)-                    &
     &                            u(i+1,j-1,k,nrhs))+                   &
     &                           0.25_r8*pm(i,j)*                       &
     &                           (v(i+1,j  ,k,nrhs)+                    &
     &                            v(i+1,j+1,k,nrhs)-                    &
     &                            v(i-1,j  ,k,nrhs)-                    &
     &                            v(i-1,j+1,k,nrhs)))**2)

#  ifdef UV_SMAGORINSKY
!
!  Smagorinsky viscosity.
!
#   if defined UV_VIS2
            visc3d_r(i,j,k)=Hviscosity(i,j)+                            &
     &                      SmagorCoef*omn(i,j)*DefRate
#   elif defined UV_VIS4
            visc3d_r(i,j,k)=Hviscosity(i,j)+                            &
     &                      PecletCoef*(omn(i,j)**2)*DefRate
#    ifdef MIX_GEO_UV
            visc3d_r(i,j,k)=MIN(clip_scale, visc3d_r(i,j,k))
#    endif
            visc3d_r(i,j,k)=SQRT(visc3d_r(i,j,k))
#   endif
#   ifdef MASKING
            visc3d_r(i,j,k)=visc3d_r(i,j,k)*rmask(i,j)
#   endif

#  elif defined UV_U3ADV_SPLIT
!
!  Compute momentum horizontal viscosity coefficient [SQRT(m4/s)]
!  as the sum of the biharmonic viscosity term and a diapycnal
!  hyperdiffusive term resulting from the splitting into advective
!  and diffusive components the 3rd-order upstream bias horizontal
!  advection.
!
            cff=0.5_r8*pm(i,j)*                                         &
     &          ABS(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
            Uvis3d_r(i,j,k)=Hviscosity(i,j)+                            &
     &                      PecletCoef*(omn(i,j)**2)*                   &
     &                      MAX(DefRate, cff)
#   ifdef MIX_GEO_UV
            Uvis3d_r(i,j,k)=MIN(clip_scale, Uvis3d_r(i,j,k))
#   endif
            Uvis3d_r(i,j,k)=SQRT(Uvis3d_r(i,j,k))
#   ifdef MASKING
            Uvis3d_r(i,j,k)=Uvis3d_r(i,j,k)*rmask(i,j)
#   endif
!
            cff=0.5_r8*pn(i,j)*                                         &
     &          ABS(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            Vvis3d_r(i,j,k)=Hviscosity(i,j)+                            &
     &                      PecletCoef*(omn(i,j)**2)*                   &
     &                      MAX(DefRate, cff)
#   ifdef MIX_GEO_UV
            Vvis3d_r(i,j,k)=MIN(clip_scale, Uvis3d_r(i,j,k))
#   endif
            Vvis3d_r(i,j,k)=SQRT(Vvis3d_r(i,j,k))
#   ifdef MASKING
            Vvis3d_r(i,j,k)=Vvis3d_r(i,j,k)*rmask(i,j)
#   endif
#  endif

#  ifdef TS_SMAGORINSKY
!
!  Smagorinsky diffusion at RHO-points.
!
#   ifdef TS_DIF2
            diff3d_r(i,j,k)=Hdiffusion(i,j)+                            &
     &                      SmagorCoef*omn(i,j)*DefRate
#   elif defined TS_DIF4
            diff3d_r(i,j,k)=Hdiffusion(i,j)+                            &
     &                      PecletCoef*(omn(i,j)**2)*DefRate
#    ifdef MIX_GEO_TS
            diff3d_r(i,j,k)=MIN(clip_scale, diff3d_r(i,j,k))
#    endif
            diff3d_r(i,j,k)=SQRT(diff3d_r(i,j,k))
#   endif
#   ifdef MASKING
            diff3d_r(i,j,k)=diff3d_r(i,j,k)*rmask(i,j)
#   endif
#  endif
          END DO
        END DO
      END DO
# endif

# if defined VEGETATION && defined VEG_HMIXING 

      clip_scale=0.01_r8*grdmax(ng)**3

      DO k=1,N(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!
#  if defined UV_VIS2
            visc3d_r(i,j,k)=visc3d_r(i,j,k)+visc3d_r_veg(i,j,k)  
#  elif defined UV_VIS4
            visc3d_r(i,j,k)=visc3d_r(i,j,k)+visc3d_r_veg(i,j,k)
#   ifdef MIX_GEO_UV
            visc3d_r(i,j,k)=MIN(clip_scale, visc3d_r(i,j,k))
#   endif
            visc3d_r(i,j,k)=SQRT(visc3d_r(i,j,k))
#  endif
#  ifdef MASKING
            visc3d_r(i,j,k)=visc3d_r(i,j,k)*rmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif

# if defined TS_U3ADV_SPLIT
!
!-----------------------------------------------------------------------
!  Compute tracer horizontal diffusion coefficient as the sum of the
!  biharmonic diffusion term and a diapycnal hyperdiffusion term
!  resulting from the splitting into advective and diffusive components
!  the 3rd-order upstream bias horizontal advection.
!-----------------------------------------------------------------------
!
!  Following Holland et al. (1998) and Webb et al. (1998), the 3rd-order
!  upstream bias horizontal advection can be splitted into advective and
!  diffusive terms. The advective term is just the standard 4th-order
!  cenrtered differences operator.  The diffusive term is a 3rd-order
!  diapycnal hyperdiffusive operator which is proportional to the
!  absolute local velocity (Marchiesello et al., 2008):
!
!         B = 1/12 * ABS(u) * (GridSpace ** 3)
!
!  The resulting diapycnak diffusion, B, is scaled and clipped to achive
!  better stability. The squared-root of the total (diapycnal and
!  biharmonic) diffusion coefficient taken since the harmonic operator
!  is applied twice to achieve biharmonic mixing.
!
!  Squared-root biharmonic diffusion coefficients [SQRT(m4/s)] at
!  U-points.
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
            diff3d_u(i,j,k)=0.5_r8*(Hdiffusion(i-1,j)+                  &
     &                              Hdiffusion(i  ,j))+                 &
     &                      PecletCoef*(om_u(i,j)**3)*                  &
     &                      ABS(u(i,j,k,nrhs))
#  ifdef MIX_GEO_TS
#   ifdef NEW_CODE
            clip_scale=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))*om_u(i,j)/        &
     &                 (z_r(i,j,k)-z_r(i-1,j,k))
            clip_diff=0.05_r8*clip_scale**4/dt(ng)
#   else
            clip_scale=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))/                  &
     &                 (z_r(i,j,k)-z_r(i-1,j,k))
            clip_diff=diff3d_u(i,j,k)*                                  &
     &                MIN(1.0_r8, clip_scale*clip_scale)
#   endif
            diff3d_u(i,j,k)=MIN(clip_diff, diff3d_u(i,j,k))
#  endif
            diff3d_u(i,j,k)=SQRT(diff3d_u(i,j,k))
#  ifdef MASKING
            diff3d_u(i,j,k)=diff3d_u(i,j,k)*umask(i,j)
#  endif
          END DO
        END DO
!
!  Squared-root biharmonic diffusion coefficients [SQRT(m4/s)] at
!  V-points.
!
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
            diff3d_v(i,j,k)=0.5_r8*(Hdiffusion(i,j-1)+                  &
     &                              Hdiffusion(i,j  ))+                 &
     &                      PecletCoef*(on_v(i,j)**3)*                  &
     &                      ABS(v(i,j,k,nrhs))
#  ifdef MIX_GEO_TS
#   ifdef NEW_CODE
            clip_scale=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*on_v(i,j)/        &
     &                 (z_r(i,j,k)-z_r(i,j-1,k))
            clip_diff=0.05_r8*clip_scale**4/dt(ng)
#   else
            clip_scale=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))/                  &
     &                 (z_r(i,j,k)-z_r(i,j-1,k))
            clip_diff=diff3d_v(i,j,k)*                                  &
     &                MIN(1.0_r8, clip_scale*clip_scale)
#   endif
            diff3d_v(i,j,k)=MIN(clip_diff, diff3d_v(i,j,k))
#  endif
            diff3d_v(i,j,k)=SQRT(diff3d_v(i,j,k))
#  ifdef MASKING
            diff3d_v(i,j,k)=diff3d_v(i,j,k)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Apply boundary conditions
!-----------------------------------------------------------------------
!
!  East-West gradient boundary conditions.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
# if defined DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
              diff3d_u(Iend+1,j,k)=diff3d_u(Iend,j,k)
              diff3d_v(Iend+1,j,k)=diff3d_v(Iend,j,k)
#  else
              diff3d_r(Iend+1,j,k)=diff3d_r(Iend,j,k)
#  endif
# endif
# if defined VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
              Uvis3d_r(Iend+1,j,k)=Uvis3d_r(Iend,j,k)
              Vvis3d_r(Iend+1,j,k)=Vvis3d_r(Iend,j,k)
#  else
              visc3d_r(Iend+1,j,k)=visc3d_r(Iend,j,k)
#  endif
# endif
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
# if defined DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
              diff3d_u(Istr-1,j,k)=diff3d_u(Istr,j,k)
              diff3d_v(Istr-1,j,k)=diff3d_v(Istr,j,k)
#  else
              diff3d_r(Istr-1,j,k)=diff3d_r(Istr,j,k)
#  endif
# endif
# if defined VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
              Uvis3d_r(Istr-1,j,k)=Uvis3d_r(Istr,j,k)
              Vvis3d_r(Istr-1,j,k)=Vvis3d_r(Istr,j,k)
#  else
              visc3d_r(Istr-1,j,k)=visc3d_r(Istr,j,k)
#  endif
# endif
            END DO
          END DO
        END IF
      END IF
!
!  North-South gradient boundary conditions.
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
              diff3d_u(i,Jend+1,k)=diff3d_u(i,Jend,k)
              diff3d_v(i,Jend+1,k)=diff3d_v(i,Jend,k)
#  else
              diff3d_r(i,Jend+1,k)=diff3d_r(i,Jend,k)
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
              Uvis3d_r(i,Jend+1,k)=Uvis3d_r(i,Jend,k)
              Vvis3d_r(i,Jend+1,k)=Vvis3d_r(i,Jend,k)
#  else
              visc3d_r(i,Jend+1,k)=visc3d_r(i,Jend,k)
#  endif
# endif
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
              diff3d_u(i,Jstr-1,k)=diff3d_u(i,Jstr,k)
              diff3d_v(i,Jstr-1,k)=diff3d_v(i,Jstr,k)
#  else
              diff3d_r(i,Jstr-1,k)=diff3d_r(i,Jstr,k)
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
              Uvis3d_r(i,Jstr-1,k)=Uvis3d_r(i,Jstr,k)
              Vvis3d_r(i,Jstr-1,k)=Vvis3d_r(i,Jstr,k)
#  else
              visc3d_r(i,Jstr-1,k)=visc3d_r(i,Jstr,k)
#  endif
# endif
            END DO
          END DO
        END IF
      END IF
!
!  Boundary corners.
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
            diff3d_u(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_u(Istr  ,Jstr-1,k)+       &
     &                                 diff3d_u(Istr-1,Jstr  ,k))
            diff3d_v(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_v(Istr  ,Jstr-1,k)+       &
     &                                 diff3d_v(Istr-1,Jstr  ,k))
#  else
            diff3d_r(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_r(Istr  ,Jstr-1,k)+       &
     &                                 diff3d_r(Istr-1,Jstr  ,k))
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
            Uvis3d_r(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (Uvis3d_r(Istr  ,Jstr-1,k)+       &
     &                                 Uvis3d_r(Istr-1,Jstr  ,k))
            Vvis3d_r(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (Vvis3d_r(Istr  ,Jstr-1,k)+       &
     &                                 Vvis3d_r(Istr-1,Jstr  ,k))
#  else
            visc3d_r(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (visc3d_r(Istr  ,Jstr-1,k)+       &
     &                                 visc3d_r(Istr-1,Jstr  ,k))
#  endif
# endif
          END DO
        END IF

        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
            diff3d_u(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_u(Iend  ,Jstr-1,k)+       &
     &                                 diff3d_u(Iend+1,Jstr  ,k))
            diff3d_v(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_v(Iend  ,Jstr-1,k)+       &
     &                                 diff3d_v(Iend+1,Jstr  ,k))
#  else
            diff3d_r(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_r(Iend  ,Jstr-1,k)+       &
     &                                 diff3d_r(Iend+1,Jstr  ,k))
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
            Uvis3d_r(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (Uvis3d_r(Iend  ,Jstr-1,k)+       &
     &                                 Uvis3d_r(Iend+1,Jstr  ,k))
            Vvis3d_r(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (Vvis3d_r(Iend  ,Jstr-1,k)+       &
     &                                 Vvis3d_r(Iend+1,Jstr  ,k))
#  else
            visc3d_r(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (visc3d_r(Iend  ,Jstr-1,k)+       &
     &                                 visc3d_r(Iend+1,Jstr  ,k))
#  endif
# endif
          END DO
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile).and.                         &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
            diff3d_u(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_u(Istr  ,Jend+1,k)+       &
     &                                 diff3d_u(Istr-1,Jend  ,k))
            diff3d_v(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_v(Istr  ,Jend+1,k)+       &
     &                                 diff3d_v(Istr-1,Jend  ,k))
#  else
            diff3d_r(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_r(Istr  ,Jend+1,k)+       &
     &                                 diff3d_r(Istr-1,Jend  ,k))
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
            Uvis3d_r(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (Uvis3d_r(Istr  ,Jend+1,k)+       &
     &                                 Uvis3d_r(Istr-1,Jend  ,k))
            Vvis3d_r(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (Vvis3d_r(Istr  ,Jend+1,k)+       &
     &                                 Vvis3d_r(Istr-1,Jend  ,k))
#  else
            visc3d_r(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (visc3d_r(Istr  ,Jend+1,k)+       &
     &                                 visc3d_r(Istr-1,Jend  ,k))
#  endif
# endif
          END DO
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile).and.                         &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
            diff3d_u(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_u(Iend  ,Jend+1,k)+       &
     &                                 diff3d_u(Iend+1,Jend  ,k))
            diff3d_v(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_v(Iend  ,Jend+1,k)+       &
     &                                 diff3d_v(Iend+1,Jend  ,k))
#  else
            diff3d_r(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_r(Istr  ,Jend+1,k)+       &
     &                                 diff3d_r(Istr-1,Jend  ,k))
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
            Uvis3d_r(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (Uvis3d_r(Iend  ,Jend+1,k)+       &
     &                                 Vvis3d_r(Iend+1,Jend  ,k))
            Vvis3d_r(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (Vvis3d_r(Iend  ,Jend+1,k)+       &
     &                                 Vvis3d_r(Iend+1,Jend  ,k))
#  else
            visc3d_r(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (visc3d_r(Iend  ,Jend+1,k)+       &
     &                                 visc3d_r(Iend+1,Jend  ,k))
#  endif
# endif
          END DO
        END IF
      END IF
!
!  Periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          diff3d_u)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          diff3d_v)
#  else
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          diff3d_r)
#  endif
# endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Uvis3d_r)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Vvis3d_r)
#  else
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          visc3d_r)
#  endif
# endif
      END IF

# ifdef DISTRIBUTE
!
!  Exhange boundary data.
!
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff3d_u, diff3d_v)
#   else
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff3d_r)
#   endif
#  endif
#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Uvis3d_r, Vvis3d_r)
#   else
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    visc3d_r)
#   endif
#  endif
# endif

      RETURN
      END SUBROUTINE hmixing_tile
#endif
      END MODULE hmixing_mod
