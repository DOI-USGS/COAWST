#include "cppdefs.h"
#ifdef SOLVE3D
      SUBROUTINE set_scoord (ng)
!
!svn $Id: set_scoord.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine sets and initializes relevant variables associated     !
!  with the vertical terrain-following coordinates transformation.     !
!                                                                      !
!  Definitions:                                                        !
!                                                                      !
!    N(ng) : Number of vertical levels for each nested grid.           !
!                                                                      !
!     zeta : time-varying free-surface, zeta(x,y,t), (m)               !
!                                                                      !
!        h : bathymetry, h(x,y), (m, positive, maybe time-varying)     !
!                                                                      !
!       hc : critical (thermocline, pycnocline) depth (m, positive)    !
!                                                                      !
!        z : vertical depths, z(x,y,s,t), meters, negative             !
!              z_w(x,y,0:N(ng))      at   W-points  (top/bottom cell)  !
!              z_r(z,y,1:N(ng))      at RHO-points  (cell center)      !
!                                                                      !
!              z_w(x,y,0    ) = -h(x,y)                                !
!              z_w(x,y,N(ng)) = zeta(x,y,t)                            !
!                                                                      !
!        s : nondimensional stretched vertical coordinate,             !
!             -1 <= s <= 0                                             !
!                                                                      !
!              s = 0   at the free-surface, z(x,y, 0,t) = zeta(x,y,t)  !
!              s = -1  at the bottom,       z(x,y,-1,t) = - h(x,y,t)   !
!                                                                      !
!              sc_w(k) = (k-N(ng))/N(ng)       k=0:N,    W-points      !
!              sc_r(k) = (k-N(ng)-0.5)/N(ng)   k=1:N,  RHO-points      !
!                                                                      !
!        C : nondimensional vertical stretching function, C(s),        !
!              -1 <= C(s) <= 0                                         !
!                                                                      !
!              C(s) = 0    for s = 0,  at the free-surface             !
!              C(s) = -1   for s = -1, at the bottom                   !
!                                                                      !
!              Cs_w(k) = F(s,theta_s,theta_b)  k=0:N,    W-points      !
!              Cs_r(k) = C(s,theta_s,theta_b)  k=1:N,  RHO-points      !
!                                                                      !
!       Zo : vertical transformation functional, Zo(x,y,s):            !
!                                                                      !
!              Zo(x,y,s) = H(x,y)C(s)      separable functions         !
!                                                                      !
!                                                                      !
!  Two vertical transformations are supported, z => z(x,y,s,t):        !
!                                                                      !
!  (1) Original transformation (Shchepetkin and McWilliams, 2005): In  !
!      ROMS since 1999 (version 1.8):                                  !
!                                                                      !
!        z(x,y,s,t) = Zo(x,y,s) + zeta(x,y,t) * [1 + Zo(x,y,s)/h(x,y)] !
!                                                                      !
!      where                                                           !
!                                                                      !
!        Zo(x,y,s) = hc * s + [h(x,y) - hc] * C(s)                     !
!                                                                      !
!        Zo(x,y,s) = 0         for s = 0,  C(s) = 0,  at the surface   !
!        Zo(x,y,s) = -h(x,y)   for s = -1, C(s) = -1, at the bottom    !
!                                                                      !
!  (2) New transformation: In UCLA-ROMS since 2005:                    !
!                                                                      !
!        z(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t) + h(x,y)] * Zo(x,y,s) !
!                                                                      !
!      where                                                           !
!                                                                      !
!        Zo(x,y,s) = [hc * s(k) + h(x,y) * C(k)] / [hc + h(x,y)]       !
!                                                                      !
!        Zo(x,y,s) = 0         for s = 0,  C(s) = 0,  at the surface   !
!        Zo(x,y,s) = -1        for s = -1, C(s) = -1, at the bottom    !
!                                                                      !
!      At the rest state, corresponding to zero free-surface, this     !
!      transformation yields the following unperturbed depths, zhat:   !
!                                                                      !
!        zhat = z(x,y,s,0) = h(x,y) * Zo(x,y,s)                        !
!                                                                      !
!             = h(x,y) * [hc * s(k) + h(x,y) * C(k)] / [hc + h(x,y)]   !
!                                                                      !
!      and                                                             !
!                                                                      !
!        d(zhat) = ds * h(x,y) * hc / [hc + h(x,y)]                    !
!                                                                      !
!      As a consequence, the uppermost grid box retains very little    !
!      dependency from bathymetry in the areas where hc << h(x,y),     !
!      that is deep areas. For example, if hc=250 m, and  h(x,y)       !
!      changes from 2000 to 6000 meters, the uppermost grid box        !
!      changes only by a factor of 1.08 (less than 10%).               !
!                                                                      !
!      Notice that:                                                    !
!                                                                      !
!      * Regardless of the design of C(s), transformation (2) behaves  !
!        like equally-spaced sigma-coordinates in shallow areas, where !
!        h(x,y) << hc.  This is advantageous because high vertical     !
!        resolution and associated CFL limitation is avoided in these  !
!        areas.                                                        !
!                                                                      !
!      * Near-surface refinement is close to geopotential coordinates  !
!        in deep areas (level thickness do not depend or weakly-depend !
!        on the bathymetry).  Contrarily,  near-bottom refinement is   !
!        like sigma-coordinates with thicknesses roughly proportional  !
!        to depth reducing high r-factors in these areas.              !
!                                                                      !
!                                                                      !
!  This generic transformation design facilitates numerous vertical    !
!  stretching functions, C(s).  These functions are set-up in this     !
!  routine in terms of several stretching parameters specified in      !
!  the standard input file.                                            !
!                                                                      !
!  C(s) vertical stretching function properties:                       !
!                                                                      !
!  * a nonlinear, monotonic function                                   !
!  * a continuous differentiable function, or                          !
!  * a piecewise function with smooth transition and differentiable    !
!  * must be constrained by -1 <= C(s) <= 0, with C(0)=0 at the        !
!    free-surface and C(-1)=-1 at the bottom (bathymetry).             !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic !
!         modeling system (ROMS): a split-explicit, free-surface,      !
!         topography-following-coordinate oceanic model, Ocean         !
!         Modelling, 9, 347-404.                                       !
!                                                                      !
!    Song, Y. and D. Haidvogel, 1994: A semi-implicit ocean            !
!         circulation model using a generalized topography-            !
!         following coordinate system,  J.  Comp.  Physics,            !
!         115, 228-244.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: k

      real(r8) :: Aweight, Bweight, Cweight, Cbot, Csur, Hscale
      real(r8) :: ds, exp_bot, exp_sur, sc_r, sc_w
      real(r8) :: cff, cff1, cff2
      real(r8) :: zhc, z1, z2, z3
!
!-----------------------------------------------------------------------
!  Set thickness controlling vertical coordinate stretching.
!-----------------------------------------------------------------------
!
!  Set hc <= hmin, in the original formulation (Vtransform=1) to avoid
!  [h(x,y)-hc] to be negative which results in dz/ds to be negative.
!  Notice that this restriction is REMOVED in the new transformation
!  (Vtransform=2): hc can be any value. It works for both hc < hmin
!  and hc > hmin.
!
      IF (Vtransform(ng).eq.1) THEN
# if defined WET_DRY
        hc(ng)=MIN(MAX(hmin(ng),0.0_r8),Tcline(ng))
# else
        hc(ng)=MIN(hmin(ng),Tcline(ng))
# endif
      ELSE IF (Vtransform(ng).eq.2) THEN
        hc(ng)=Tcline(ng)
      END IF
!
!-----------------------------------------------------------------------
!  Original vertical strectching function, Song and Haidvogel (1994).
!-----------------------------------------------------------------------
!
      IF (Vstretching(ng).eq.1) THEN
!
!  This vertical stretching function is defined as:
!
!      C(s) = (1 - b) * [SINH(s * a) / SINH(a)] +
!
!             b * [-0.5 + 0.5 * TANH(a * (s + 0.5)) / TANH(0.5 * a)]
!
!  where the stretching parameters (a, b) are specify at input:
!
!         a = theta_s               0 <  theta_s <= 8
!         b = theta_b               0 <= theta_b <= 1
!
!  If theta_b=0, the refinement is surface intensified as theta_s is
!  increased.
!  If theta_b=1, the refinement is both bottom ans surface intensified
!  as theta_s is increased.
!
        IF (theta_s(ng).ne.0.0_r8) THEN
          cff1=1.0_r8/SINH(theta_s(ng))
          cff2=0.5_r8/TANH(0.5_r8*theta_s(ng))
        END IF
        SCALARS(ng)%sc_w(0)=-1.0_r8
        SCALARS(ng)%Cs_w(0)=-1.0_r8
        ds=1.0_r8/REAL(N(ng),r8)
        DO k=1,N(ng)
          SCALARS(ng)%sc_w(k)=ds*REAL(k-N(ng),r8)
          SCALARS(ng)%sc_r(k)=ds*(REAL(k-N(ng),r8)-0.5_r8)
          IF (theta_s(ng).ne.0.0_r8) THEN
            SCALARS(ng)%Cs_w(k)=(1.0_r8-theta_b(ng))*                   &
     &                          cff1*SINH(theta_s(ng)*                  &
     &                                    SCALARS(ng)%sc_w(k))+         &
     &                          theta_b(ng)*                            &
     &                          (cff2*TANH(theta_s(ng)*                 &
     &                                     (SCALARS(ng)%sc_w(k)+        &
     &                                      0.5_r8))-                   &
     &                           0.5_r8)
            SCALARS(ng)%Cs_r(k)=(1.0_r8-theta_b(ng))*                   &
     &                          cff1*SINH(theta_s(ng)*                  &
     &                                    SCALARS(ng)%sc_r(k))+         &
     &                          theta_b(ng)*                            &
     &                          (cff2*TANH(theta_s(ng)*                 &
     &                                     (SCALARS(ng)%sc_r(k)+        &
     &                                      0.5_r8))-                   &
     &                           0.5_r8)
          ELSE
            SCALARS(ng)%Cs_w(k)=SCALARS(ng)%sc_w(k)
            SCALARS(ng)%Cs_r(k)=SCALARS(ng)%sc_r(k)
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  A. Shchepetkin vertical stretching function. This function was
!  improved further to allow bottom refiment (see Vstretching=4).
!-----------------------------------------------------------------------
!
      ELSE IF (Vstretching(ng).eq.2) THEN
!
!  This vertical stretching function is defined, in the simplest form,
!  as:
!
!      C(s) = [1.0 - COSH(theta_s * s)] / [COSH(theta_s) - 1.0]
!
!  it is similar in meaning to the original vertical stretcing function
!  (Song and Haidvogel, 1994), but note that hyperbolic functions are
!  COSH, and not SINH.
!
!  Note that the above definition results in
!
!         -1 <= C(s) <= 0
!
!  as long as
!
!         -1 <= s <= 0
!
!  and, unlike in any previous definition
!
!         d[C(s)]/ds  -->  0      if  s -->  0
!
!  For the purpose of bottom boundary layer C(s) is further modified
!  to allow near-bottom refinement.  This is done by blending it with
!  another function.
!
        Aweight=1.0_r8
        Bweight=1.0_r8
        ds=1.0_r8/REAL(N(ng),r8)
!
        SCALARS(ng)%sc_w(N(ng))=0.0_r8
        SCALARS(ng)%Cs_w(N(ng))=0.0_r8
        DO k=N(ng)-1,1,-1
          sc_w=ds*REAL(k-N(ng),r8)
          SCALARS(ng)%sc_w(k)=sc_w
          IF (theta_s(ng).gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s(ng)*sc_w))/                       &
     &           (COSH(theta_s(ng))-1.0_r8)
            IF (theta_b(ng).gt.0.0_r8) THEN
              Cbot=SINH(theta_b(ng)*(sc_w+1.0_r8))/                     &
     &             SINH(theta_b(ng))-1.0_r8
              Cweight=(sc_w+1.0_r8)**Aweight*                           &
     &                (1.0_r8+(Aweight/Bweight)*                        &
     &                        (1.0_r8-(sc_w+1.0_r8)**Bweight))
              SCALARS(ng)%Cs_w(k)=Cweight*Csur+(1.0_r8-Cweight)*Cbot
            ELSE
              SCALARS(ng)%Cs_w(k)=Csur
            END IF
          ELSE
            SCALARS(ng)%Cs_w(k)=sc_w
          END IF
        END DO
        SCALARS(ng)%sc_w(0)=-1.0_r8
        SCALARS(ng)%Cs_w(0)=-1.0_r8
!
        DO k=1,N(ng)
          sc_r=ds*(REAL(k-N(ng),r8)-0.5_r8)
          SCALARS(ng)%sc_r(k)=sc_r
          IF (theta_s(ng).gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s(ng)*sc_r))/                       &
     &           (COSH(theta_s(ng))-1.0_r8)
            IF (theta_b(ng).gt.0.0_r8) THEN
              Cbot=SINH(theta_b(ng)*(sc_r+1.0_r8))/                     &
     &             SINH(theta_b(ng))-1.0_r8
              Cweight=(sc_r+1.0_r8)**Aweight*                           &
     &                (1.0_r8+(Aweight/Bweight)*                        &
     &                        (1.0_r8-(sc_r+1.0_r8)**Bweight))
              SCALARS(ng)%Cs_r(k)=Cweight*Csur+(1.0_r8-Cweight)*Cbot
            ELSE
              SCALARS(ng)%Cs_r(k)=Csur
            END IF
          ELSE
            SCALARS(ng)%Cs_r(k)=sc_r
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  R. Geyer stretching function for high bottom boundary layer
!  resolution.
!-----------------------------------------------------------------------
!
      ELSE IF (Vstretching(ng).eq.3) THEN
!
!  This stretching function is intended for very shallow coastal
!  applications, like gravity sediment flows.
!
!  At the surface, C(s=0)=0
!
!      C(s) = - LOG(COSH(Hscale * ABS(s) ** alpha)) /
!               LOG(COSH(Hscale))
!
!  At the bottom, C(s=-1)=-1
!
!      C(s) = LOG(COSH(Hscale * (s + 1) ** beta)) /
!             LOG(COSH(Hscale)) - 1
!
!  where
!
!       Hscale : scale value for all hyperbolic functions
!                  Hscale = 3.0    set internally here
!        alpha : surface stretching exponent
!                  alpha = 0.65   minimal increase of surface resolution
!                          1.0    significant amplification
!         beta : bottoom stretching exponent
!                  beta  = 0.58   no amplification
!                          1.0    significant amplification
!                          3.0    super-high bottom resolution
!            s : stretched vertical coordinate, -1 <= s <= 0
!                  s(k) = (k-N)/N       k=0:N,    W-points  (s_w)
!                  s(k) = (k-N-0.5)/N   k=1:N,  RHO-points  (s_rho)
!
!  The stretching exponents (alpha, beta) are specify at input:
!
!         alpha = theta_s
!         beta  = theta_b
!
        exp_sur=theta_s(ng)
        exp_bot=theta_b(ng)
        Hscale=3.0_r8
        ds=1.0_r8/REAL(N(ng),r8)
!
        SCALARS(ng)%sc_w(N(ng))=0.0_r8
        SCALARS(ng)%Cs_w(N(ng))=0.0_r8
        DO k=N(ng)-1,1,-1
          sc_w=ds*REAL(k-N(ng),r8)
          SCALARS(ng)%sc_w(k)=sc_w
          Cbot= LOG(COSH(Hscale*(sc_w+1.0_r8)**exp_bot))/               &
     &          LOG(COSH(Hscale))-1.0_r8
          Csur=-LOG(COSH(Hscale*ABS(sc_w)**exp_sur))/                   &
     &          LOG(COSH(Hscale))
          Cweight=0.5_r8*(1.0_r8-TANH(Hscale*(sc_w+0.5_r8)))
          SCALARS(ng)%Cs_w(k)=Cweight*Cbot+(1.0_r8-Cweight)*Csur
        END DO
        SCALARS(ng)%sc_w(0)=-1.0_r8
        SCALARS(ng)%Cs_w(0)=-1.0_r8
!
        DO k=1,N(ng)
          sc_r=ds*(REAL(k-N(ng),r8)-0.5_r8)
          SCALARS(ng)%sc_r(k)=sc_r
          Cbot= LOG(COSH(Hscale*(sc_r+1.0_r8)**exp_bot))/               &
     &          LOG(COSH(Hscale))-1.0_r8
          Csur=-LOG(COSH(Hscale*ABS(sc_r)**exp_sur))/                   &
     &          LOG(COSH(Hscale))
          Cweight=0.5_r8*(1.0_r8-TANH(Hscale*(sc_r+0.5_r8)))
          SCALARS(ng)%Cs_r(k)=Cweight*Cbot+(1.0_r8-Cweight)*Csur
        END DO
!
!-----------------------------------------------------------------------
!  A. Shchepetkin improved double vertical stretching functions with
!  bottom refiment.
!-----------------------------------------------------------------------
!
      ELSE IF (Vstretching(ng).eq.4) THEN
!
!  The range of meaningful values for the control parameters are:
!
!       0 <  theta_s <= 10.0
!       0 <= theta_b <=  3.0
!
!  Users need to pay attention to extreme r-factor (rx1) values near
!  the bottom.
!
!  This vertical stretching function is defined, in the simplest form,
!  as:
!
!      C(s) = [1.0 - COSH(theta_s * s)] / [COSH(theta_s) - 1.0]
!
!  it is similar in meaning to the original vertical stretcing function
!  (Song and Haidvogel, 1994), but note that hyperbolic functions are
!  COSH, and not SINH.
!
!  Note that the above definition results in
!
!         -1 <= C(s) <= 0
!
!  as long as
!
!         -1 <= s <= 0
!
!  and
!
!         d[C(s)]/ds  -->  0      if  s -->  0
!
!  For the purpose of bottom boundary layer C(s) is further modified
!  to allow near-bottom refinement by using a continuous, second
!  stretching function
!
!         C(s) = [EXP(theta_b * C(s)) - 1.0] / [1.0 - EXP(-theta_b)]
!
!  This double transformation is continuous with respect to "theta_s"
!  and "theta_b", as both values approach to zero.
!
        ds=1.0_r8/REAL(N(ng),r8)
!
        SCALARS(ng)%sc_w(N(ng))=0.0_r8
        SCALARS(ng)%Cs_w(N(ng))=0.0_r8
        DO k=N(ng)-1,1,-1
          sc_w=ds*REAL(k-N(ng),r8)
          SCALARS(ng)%sc_w(k)=sc_w
          IF (theta_s(ng).gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s(ng)*sc_w))/                       &
     &           (COSH(theta_s(ng))-1.0_r8)
          ELSE
            Csur=-sc_w**2
          END IF
          IF (theta_b(ng).gt.0.0_r8) THEN
            Cbot=(EXP(theta_b(ng)*Csur)-1.0_r8)/                        &
     &           (1.0_r8-EXP(-theta_b(ng)))
            SCALARS(ng)%Cs_w(k)=Cbot
          ELSE
            SCALARS(ng)%Cs_w(k)=Csur
          END IF
        END DO
        SCALARS(ng)%sc_w(0)=-1.0_r8
        SCALARS(ng)%Cs_w(0)=-1.0_r8
!
        DO k=1,N(ng)
          sc_r=ds*(REAL(k-N(ng),r8)-0.5_r8)
          SCALARS(ng)%sc_r(k)=sc_r
          IF (theta_s(ng).gt.0.0_r8) THEN
            Csur=(1.0_r8-COSH(theta_s(ng)*sc_r))/                       &
     &           (COSH(theta_s(ng))-1.0_r8)
          ELSE
            Csur=-sc_r**2
          END IF
          IF (theta_b(ng).gt.0.0_r8) THEN
            Cbot=(EXP(theta_b(ng)*Csur)-1.0_r8)/                        &
     &           (1.0_r8-EXP(-theta_b(ng)))
            SCALARS(ng)%Cs_r(k)=Cbot
          ELSE
            SCALARS(ng)%Cs_r(k)=Csur
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Report information about vertical transformation.
!-----------------------------------------------------------------------
!
      IF (Master.and.LwrtInfo(ng)) THEN
        WRITE (stdout,10) ng
        cff=0.5_r8*(hmax(ng)+hmin(ng))
        DO k=N(ng),0,-1
          IF (Vtransform(ng).eq.1) THEN
            zhc=hc(ng)*SCALARS(ng)%sc_w(k)
            z1=zhc+(hmin(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
            z2=zhc+(cff     -hc(ng))*SCALARS(ng)%Cs_w(k)
            z3=zhc+(hmax(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
          ELSE IF (Vtransform(ng).eq.2) THEN
            z1=hmin(ng)*(hc(ng)  *SCALARS(ng)%sc_w(k)+                  &
     &                   hmin(ng)*SCALARS(ng)%Cs_w(k))/(hc(ng)+hmin(ng))
            z2=cff     *(hc(ng)  *SCALARS(ng)%sc_w(k)+                  &
     &                   cff     *SCALARS(ng)%Cs_w(k))/(hc(ng)+cff)
            z3=hmax(ng)*(hc(ng)  *SCALARS(ng)%sc_w(k)+                  &
     &                   hmax(ng)*SCALARS(ng)%Cs_w(k))/(hc(ng)+hmax(ng))
            IF (hc(ng).gt.hmax(ng)) THEN
              zhc=z3      ! same as hmax, other values do not make sense
            ELSE
              zhc=0.5_r8*hc(ng)*(SCALARS(ng)%sc_w(k)+                   &
     &                           SCALARS(ng)%Cs_w(k))
            END IF
          END IF
          WRITE (stdout,20) k, SCALARS(ng)%sc_w(k),                     &
     &                         SCALARS(ng)%Cs_w(k), z1, zhc, z2, z3
        END DO
      END IF

  10  FORMAT (/,' Vertical S-coordinate System, Grid ',i2.2,':'/,/,     &
     &          ' level   S-coord     Cs-curve   Z',3x,                 &
     &          'at hmin       at hc    half way     at hmax',/)
  20  FORMAT (i6,2f12.7,1x,4f12.3)

      RETURN
      END SUBROUTINE set_scoord
#else
      SUBROUTINE set_scoord
      RETURN
      END SUBROUTINE set_scoord
#endif

