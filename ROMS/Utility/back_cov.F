#include "cppdefs.h" 
      MODULE back_cov_mod
#ifdef S4DVAR
!
!svn $Id: back_cov.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the background error covariance (B) which     !
!  is used as preconditioning for the descent algorithm.  When the     !
!  background term is added to the cost function,  the appropriate     !
!  scaling for preconditioning is B^(-1/2) for the model variables     !
!  and transpose[B^(1/2)] for the gradient.                            !
!                                                                      !
!  However,  if the cost function contains no background term,  at     !
!  preconditioning the model variables are scaled by  L^(1/2)  and     !
!  the gradient by L^(-1/2).  Here, L is a norm operator, itself a     !
!  diagonal matrix.                                                    !
!                                                                      !
!======================================================================!
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE back_cov (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL back_cov_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nstp(ng), nnew(ng), Lold(ng),                 &
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % omn,                               &
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
# ifdef SOLVE3D
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
#  ifdef ICESHELF
     &                    GRID(ng) % zice,                              &
#  endif
     &                    COUPLING(ng) % Zt_avg1,                       &
#  if defined SEdIMENT && defined SED_MORPH
     &                    GRID(ng) % bed_thick,                         &
#  endif
     &                    OCEAN(ng) % b_u,                              &
     &                    OCEAN(ng) % b_v,                              &
     &                    OCEAN(ng) % b_t,                              &
# endif
     &                    OCEAN(ng) % b_ubar,                           &
     &                    OCEAN(ng) % b_vbar,                           &
     &                    OCEAN(ng) % b_zeta)
      RETURN
      END SUBROUTINE back_cov
!
!***********************************************************************
      SUBROUTINE back_cov_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nstp, nnew, Lold,                       &
     &                          h, omn, om_u, om_v, on_u, on_v,         &
# ifdef SOLVE3D
     &                          Hz, z_r, z_w,                           &
#  ifdef ICESHELF
     &                          zice,                                   &
#  endif
     &                          Zt_avg1,                                &
#  if defined SEdIMENT && defined SED_MORPH
     &                          bed_thick,                              &
#  endif
     &                          b_u, b_v, b_t,                          &
# endif
     &                          b_ubar, b_vbar, b_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars

# ifdef SOLVE3D
!
      USE set_depth_mod, ONLY : set_depth_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  ifdef SOLVE3D
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEdIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed_thick(LBi:,LBj:,:)
#   endif
#  endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_w(LBi:,LBj:,:)
      real(r8), intent(out) :: b_u(LBi:,LBj:,:)
      real(r8), intent(out) :: b_v(LBi:,LBj:,:)
      real(r8), intent(out) :: b_t(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(out) :: b_ubar(LBi:,LBj:)
      real(r8), intent(out) :: b_vbar(LBi:,LBj:)
      real(r8), intent(out) :: b_zeta(LBi:,LBj:)
# else
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEdIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
#  endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_w(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: b_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: b_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: b_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
      real(r8), intent(out) :: b_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: b_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: b_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef INPUT_COVARIANCE
      logical :: k_bound, r_bound, u_bound, v_bound
      integer :: iobs
# endif
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k, order

      real(r8), dimension(NT(ng)) :: Tnorm, TrcCoef
# endif
      real(r8) :: cff, cff1, cff2, fac
!
# include "set_bounds.h"

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set tracer norm coefficients.
!-----------------------------------------------------------------------
!
#  ifdef N2NORM_PROFILE
      cff=1.0
#  else
      cff=1.0/bvf_bak
#  endif
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          TrcCoef(itrc)=Tcoef(ng)*Tcoef(ng)*g*g*cff
          Tnorm(itrc)=TrcCoef(itrc)
        ELSE IF (itrc.eq.isalt) THEN
          TrcCoef(itrc)=Scoef(ng)*Scoef(ng)*g*g*cff
          Tnorm(itrc)=TrcCoef(itrc)
        ELSE
          TrcCoef(itrc)=1.0_r8
          Tnorm(itrc)=TrcCoef(itrc)
        END IF
      END DO
# endif
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses for preconditioning.
!-----------------------------------------------------------------------
! 
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
#  endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
# endif
# ifndef BACKGROUND
!
!-----------------------------------------------------------------------
!  If the cost function contains no background term, set background
!  error covariance arrays with norm diagonal matrix (L). Scale the
!  model variables by L^(1/2) and the gradient by L^(-1/2).
!-----------------------------------------------------------------------
!
#  ifdef INPUT_COVARIANCE
      DO iobs=1,Nobs(ng)
        i=INT(Xobs(iobs))
        j=INT(Yobs(iobs))
#   ifdef SOLVE3D
        k=INT(Zobs(iobs))
        k_bound=((1.le.k).and.(k.le.N(ng)))
#   endif
        r_bound=((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))
        u_bound=((Istr .le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))
        v_bound=((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((Jstr .le.j).and.(j.le.JendR))
        cff1=SQRT(ObsErr(iobs))
        IF (r_bound.and.                                                &
     &      (ObsType(iobs).eq.isFsur)) THEN
          b_zeta(i,j)=cff1
        ELSE IF (u_bound.and.                                           &
     &           (ObsType(iobs).eq.isUbar)) THEN
          b_ubar(i,j)=cff1
        ELSE IF (v_bound.and.                                           &
     &           (ObsType(iobs).eq.isVbar)) THEN
          b_vbar(i,j)=cff1
#   ifdef SOLVE3D
        ELSE IF (u_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isUvel)) THEN
          b_u(i,j,k)=cff1
        ELSE IF (v_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isVvel)) THEN
          b_v(i,j,k)=cff1
        ELSE
          DO itrc=1,NT(ng)
            IF (r_bound.and.k_bound.and.                                &
     &          (ObsType(iobs).eq.isTvar(itrc))) THEN
              b_t(i,j,k,itrc)=cff1
            END IF
          END DO
#   endif
        END IF
      END DO
#  else
      fac=0.5_r8*rho0/TotVolume
!
!  Free-surface preconditioning.
!
      cff1=1.0_r8/SQRT(FOURDVAR(ng)%ObsVar(isFsur))
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#   if defined ENERGY1_NORM
          cff=fac*omn(i,j)*g/h(i,j)
          b_zeta(i,j)=SQRT(cff)
#   elif defined ENERGY2_NORM
          cff=fac*omn(i,j)*g
          b_zeta(i,j)=SQRT(cff)
#   elif defined ENERGY3_NORM
          cff=g/h(i,j)
          b_zeta(i,j)=SQRT(cff)
#   else
          b_zeta(i,j)=cff1
#   endif
        END DO
      END DO
!
!  2D momentum precondtioning.
!
      cff1=1.0_r8/SQRT(FOURDVAR(ng)%ObsVar(isUbar))
      cff2=1.0_r8/SQRT(FOURDVAR(ng)%ObsVar(isVbar))
      DO j=JstrR,JendR
        DO i=Istr,IendR
#   if defined ENERGY1_NORM
          cff=fac*om_u(i,j)*on_u(i,j)
          b_ubar(i,j)=SQRT(cff)
#   elif defined ENERGY2_NORM
          cff=fac*om_u(i,j)*on_u(i,j)*0.5_r8*(h(i-1,j)+h(i,j))
          b_ubar(i,j)=SQRT(cff)
#   elif defined ENERGY3_NORM
          b_ubar(i,j)=1.0_r8
#   else
          b_ubar(i,j)=cff1
#   endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
#   if defined ENERGY1_NORM
          cff=fac*om_v(i,j)*on_v(i,j)
          b_vbar(i,j)=SQRT(cff)
#   elif defined ENERGY2_NORM
          cff=fac*om_v(i,j)*on_v(i,j)*0.5_r8*(h(i,j-1)+h(i,j))
          b_vbar(i,j)=SQRT(cff)
#   elif defined ENERGY3_NORM
          b_vbar(i,j)=1.0_r8
#   else
          b_vbar(i,j)=cff2
#   endif
        END DO
      END DO
#   ifdef SOLVE3D
!
!  3D momentum preconditioning.
!
      cff1=1.0_r8/SQRT(FOURDVAR(ng)%ObsVar(isUvel))
      cff2=1.0_r8/SQRT(FOURDVAR(ng)%ObsVar(isVvel))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
#    if defined ENERGY1_NORM
            cff=fac*om_u(i,j)*on_u(i,j)*                                &
     &          (Hz(i-1,j,k)+Hz(i,j,k))/(h(i-1,j)+h(i,j))
            b_u(i,j,k)=SQRT(cff)
#    elif defined ENERGY2_NORM
            cff=fac*om_u(i,j)*on_u(i,j)*                                &
     &          0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            b_u(i,j,k)=SQRT(cff)
#    elif defined ENERGY3_NORM
            b_u(i,j,k)=1.0_r8
#    else
            b_u(i,j,k)=cff1
#    endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
#    if defined ENERGY1_NORM
            cff=fac*om_v(i,j)*on_v(i,j)*                                &
     &          (Hz(i,j-1,k)+Hz(i,j,k))/(h(i,j-1)+h(i,j))
            b_v(i,j,k)=SQRT(cff)
#    elif defined ENERGY2_NORM
            cff=fac*om_v(i,j)*on_v(i,j)*                                &
     &          0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            b_v(i,j,k)=SQRT(cff)
#    elif defined ENERGY3_NORM
            b_v(i,j,k)=1.0_r8
#    else
            b_v(i,j,k)=cff2
#    endif
          END DO
        END DO
      END DO
!
!  Tracers preconditioning.
!
      DO itrc=1,NT(ng)
        cff1=1.0_r8/SQRT(FOURDVAR(ng)%ObsVar(isTvar(itrc)))
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
#    ifdef N2NORM_PROFILE
              cff2=pcoef_N2(npN2)
              DO order=npN2-1,0,-1
                cff2=z_r(i,j,k)*cff2+pcoef_N2(order)
              END DO
!!            cff2=1.0E-7_r8+1.0E-5_r8*                                 &
!!   &             (1.0_r8+TANH(z_r(i,j,k)*0.002_r8))
              TrcCoef(itrc)=Tnorm(itrc)/ABS(cff2)
#    endif
#    if defined ENERGY1_NORM
              cff=fac*TrcCoef(itrc)*omn(i,j)*Hz(i,j,k)/h(i,j)
              b_t(i,j,k,itrc)=SQRT(cff)
#    elif defined ENERGY2_NORM
              cff=fac*TrcCoef(itrc)*omn(i,j)*Hz(i,j,k)
              b_t(i,j,k,itrc)=SQRT(cff)
#    elif defined ENERGY3_NORM
              cff=TrcCoef(itrc)
              b_t(i,j,k,itrc)=SQRT(cff)
#    else
              b_t(i,j,k,itrc)=cff1
#    endif
            END DO
          END DO
        END DO
      END DO
#   endif
#  endif
# endif
      RETURN
      END SUBROUTINE back_cov_tile
#endif
      END MODULE back_cov_mod
