#include "cppdefs.h"
#define MPDATA_HOT

      MODULE mpdata_adiff_mod
#if defined NONLINEAR && defined TS_MPDATA && defined SOLVE3D
!
!svn $Id: mpdata_adiff.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes anti-diffusive velocities to correct tracer   !
!  advection using MPDATA Recursive method.                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ua      Andi-diffusive velocity in the XI-direction (m/s).       !
!     Va      Anti-diffusive velocity in the ETA-direction (m/s).      !
!     Wa      Anti-diffusive velocity in the S-direction (m3/s).       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Margolin, L. and P.K. Smolarkiewicz, 1998:  Antidiffusive         !
!      velocities for multipass donor cell advection,  SIAM J.         !
!      Sci. Comput., 907-929.                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC :: mpdata_adiff_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE mpdata_adiff_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef WET_DRY
     &                              rmask_wet, umask_wet, vmask_wet,    &
# endif
     &                              pm, pn, omn, om_u, on_v,            &
     &                              z_r, oHz,                           &
     &                              Huon, Hvom, w,                      &
# ifdef WEC_VF
     &                              W_stokes,                           &
# endif
     &                              t, Ta, Ua, Va, Wa)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: oHz(IminS:,JminS:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:)
      real(r8), intent(in) :: w(LBi:,LBj:,0:)
#  ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
#  endif
      real(r8), intent(inout) :: Ta(IminS:,JminS:,:)

      real(r8), intent(out) :: Ua(IminS:,JminS:,:)
      real(r8), intent(out) :: Va(IminS:,JminS:,:)
      real(r8), intent(out) :: Wa(IminS:,JminS:,0:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: oHz(IminS:ImaxS,JminS:JmaxS,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(inout) :: Ta(IminS:ImaxS,JminS:JmaxS,N(ng))

      real(r8), intent(out) :: Ua(IminS:ImaxS,JminS:JmaxS,N(ng))
      real(r8), intent(out) :: Va(IminS:ImaxS,JminS:JmaxS,N(ng))
      real(r8), intent(out) :: Wa(IminS:ImaxS,JminS:JmaxS,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, is, j, k

      real(r8), parameter :: Large = 1.0E+20_r8
      real(r8), parameter :: eps = 1.0E-18_r8
      real(r8), parameter :: eps2= 1.0E-10_r8

      real(r8) :: A, B, Tmax, Tmin, Um, Vm, X, Y, Z
      real(r8) :: cff, cff1, cff2, sig_alfa, fac
# ifdef MPDATA_HOT
      real(r8) :: AA, BB, CC, AB, AC, BC
      real(r8) :: XX, YY, ZZ, XY, XZ, YZ
      real(r8) :: sig_beta, sig_gama
      real(r8) :: sig_a, sig_b, sig_c, sig_d, sig_e, sig_f
# endif

      real(r8), dimension(IminS:ImaxS,N(ng)) :: C
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Wm

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: mask_dn
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: mask_up

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: beta_dn
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: beta_up
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: odz

# include "set_bounds.h"
# ifdef TS_MPDATA_LIMIT
      fac=0.25_r8
# else
      fac=1.0_r8
# endif
!
!-----------------------------------------------------------------------
!  Compute anti-diffusive MPDATA velocities (Ua,Va,Wa).
!-----------------------------------------------------------------------
!
!  Set boundary conditions to advected tracer, Ta.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrVm2,Jendp2i
              Ta(Istr-1,j,k)=Ta(Istr,j,k)
            END DO
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrVm2,Jendp2i
              Ta(Iend+1,j,k)=Ta(Iend,j,k)
            END DO
          END DO
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrUm2,Iendp2i
              Ta(i,Jstr-1,k)=Ta(i,Jstr,k)
            END DO
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrUm2,Iendp2i
              Ta(i,Jend+1,k)=Ta(i,Jend,k)
            END DO
          END DO
        END IF
      END IF

      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            Ta(Istr-1,Jstr-1,k)=0.5_r8*(Ta(Istr  ,Jstr-1,k)+            &
     &                                  Ta(Istr-1,Jstr  ,k))
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            Ta(Iend+1,Jstr-1,k)=0.5_r8*(Ta(Iend+1,Jstr  ,k)+            &
     &                                  Ta(Iend  ,Jstr-1,k))
          END DO
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            Ta(Istr-1,Jend+1,k)=0.5_r8*(Ta(Istr-1,Jend  ,k)+            &
     &                                  Ta(Istr  ,Jend+1,k))
          END DO
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            Ta(Iend+1,Jend+1,k)=0.5_r8*(Ta(Iend+1,Jend  ,k)+            &
     &                                  Ta(Iend  ,Jend+1,k))
          END DO
        END IF
      END IF
!
!  Compute inverse vertical grid spacing at W-points.
!
      DO k=1,N(ng)-1
        DO j=Jstrm2,Jendp2
          DO i=Istrm2,Iendp2
            odz(i,j,k)=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
          END DO
        END DO
      END DO
      cff=1.0_r8/dt(ng)
!
!  Compute nondimensional U-antidiffusive velocities, Ua. If applicable,
!  retain up to third-order terms of the power series.
!
      DO j=JstrV-1,Jendp1
        k=1
        DO i=IstrUm1,Iendp2
          C(i,k)=0.25_r8*                                               &
     &           ((Ta(i  ,j,k+1)-Ta(i  ,j,k  ))*odz(i  ,j,k  )+         &
     &            (Ta(i-1,j,k+1)-Ta(i-1,j,k  ))*odz(i-1,j,k  ))*        &
     &           (z_r(i  ,j,k+1)-z_r(i  ,j,k)+                          &
     &            z_r(i-1,j,k+1)-z_r(i-1,j,k))/                         &
     &           (Ta(i-1,j,k)+Ta(i,j,k)+eps)
          Wm(i,k)=0.25_r8*dt(ng)*                                       &
     &            (w(i-1,j,k  )*odz(i-1,j,k)*pm(i-1,j)*pn(i-1,j)+       &
     &             w(i  ,j,k  )*odz(i  ,j,k)*pm(i  ,j)*pn(i  ,j))
# ifdef WEC_VF
          Wm(i,k)=Wm(i,k)+                                              &
     &            0.25_r8*dt(ng)*                                       &
     &            (W_stokes(i-1,j,k  )*odz(i-1,j,k)*pm(i-1,j)*pn(i-1,j)+&
     &             W_stokes(i  ,j,k  )*odz(i  ,j,k)*pm(i  ,j)*pn(i  ,j))
# endif
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU-1,Iendp2
            C(i,k)=0.0625_r8*                                           &
     &             ((Ta(i  ,j,k+1)-Ta(i  ,j,k  ))*odz(i  ,j,k  )+       &
     &              (Ta(i  ,j,k  )-Ta(i  ,j,k-1))*odz(i  ,j,k-1)+       &
     &              (Ta(i-1,j,k+1)-Ta(i-1,j,k  ))*odz(i-1,j,k  )+       &
     &              (Ta(i-1,j,k  )-Ta(i-1,j,k-1))*odz(i-1,j,k-1))*      &
     &             (z_r(i  ,j,k+1)-z_r(i  ,j,k-1)+                      &
     &              z_r(i-1,j,k+1)-z_r(i-1,j,k-1))/                     &
     &             (Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25_r8*dt(ng)*                                     &
     &              ((w(i-1,j,k-1)*odz(i-1,j,k-1)+                      &
     &                w(i-1,j,k  )*odz(i-1,j,k  ))*pm(i-1,j)*pn(i-1,j)+ &
     &               (w(i  ,j,k  )*odz(i  ,j,k  )+                      &
     &                w(i  ,j,k-1)*odz(i  ,j,k-1))*pm(i  ,j)*pn(i  ,j))
# ifdef WEC_VF
            Wm(i,k)=Wm(i,k)+                                            &
     &              0.25_r8*dt(ng)*                                     &
     &              ((W_stokes(i-1,j,k-1)*odz(i-1,j,k-1)+               &
     &                W_stokes(i-1,j,k  )*odz(i-1,j,k  ))*              &
     &                                    pm(i-1,j)*pn(i-1,j)+          &
     &               (W_stokes(i  ,j,k  )*odz(i  ,j,k  )+               &
     &                W_stokes(i  ,j,k-1)*odz(i  ,j,k-1))*              &
     &                                    pm(i  ,j)*pn(i  ,j))
# endif
          END DO
        END DO
        k=N(ng)
        DO i=IstrU-1,Iendp2
          C(i,k)=0.25_r8*                                               &
     &           ((Ta(i  ,j,k  )-Ta(i  ,j,k-1))*odz(i  ,j,k-1)+         &
     &            (Ta(i-1,j,k  )-Ta(i-1,j,k-1))*odz(i-1,j,k-1))*        &
     &           (z_r(i  ,j,k  )-z_r(i  ,j,k-1)+                        &
     &            z_r(i-1,j,k  )-z_r(i-1,j,k-1))/                       &
     &           (Ta(i-1,j,k)+Ta(i,j,k)+eps)
          Wm(i,k)=0.25_r8*dt(ng)*                                       &
     &            (w(i-1,j,k-1)*odz(i-1,j,k-1)*pm(i-1,j)*pn(i-1,j)+     &
     &             w(i  ,j,k-1)*odz(i  ,j,k-1)*pm(i  ,j)*pn(i  ,j))
# ifdef WEC_VF
            Wm(i,k)=Wm(i,k)+                                            &
     &              0.25_r8*dt(ng)*                                     &
     &              (W_stokes(i-1,j,k-1)*odz(i-1,j,k-1)*                &
     &                                   pm(i-1,j)*pn(i-1,j)+           &
     &               W_stokes(i  ,j,k-1)*odz(i  ,j,k-1)*                &
     &                                   pm(i  ,j)*pn(i  ,j))
# endif
        END DO
        DO k=1,N(ng)
          DO i=IstrU-1,Iendp2
            IF ((Ta(i-1,j,k).le.0.0_r8).or.                             &
     &          (Ta(i  ,j,k).le.0.0_r8).or.                             &
     &          (ABS(Ta(i-1,j,k)-Ta(i,j,k)).le.eps2)) THEN
              Ua(i,j,k)=0.0_r8
            ELSE
              A=(Ta(i,j,k)-Ta(i-1,j,k))/                                &
     &          (Ta(i,j,k)+Ta(i-1,j,k)+eps)
# ifdef MASKING
              B=0.03125_r8*                                             &
     &          ((Ta(i  ,j+1,k)-Ta(i  ,j  ,k))*                         &
     &           (pn(i  ,j  )+pn(i  ,j+1))*vmask(i  ,j+1)+              &
     &           (Ta(i  ,j  ,k)-Ta(i  ,j-1,k))*                         &
     &           (pn(i  ,j-1)+pn(i  ,j  ))*vmask(i  ,j  )+              &
     &           (Ta(i-1,j+1,k)-Ta(i-1,j  ,k))*                         &
     &           (pn(i-1,j  )+pn(i-1,j+1))*vmask(i-1,j+1)+              &
     &           (Ta(i-1,j  ,k)-Ta(i-1,j-1,k))*                         &
     &           (pn(i-1,j-1)+pn(i-1,j  ))*vmask(i-1,j  ))
# else
              B=0.03125_r8*                                             &
     &          ((Ta(i  ,j+1,k)-Ta(i  ,j  ,k))*                         &
     &           (pn(i  ,j  )+pn(i  ,j+1))+                             &
     &           (Ta(i  ,j  ,k)-Ta(i  ,j-1,k))*                         &
     &           (pn(i  ,j-1)+pn(i  ,j  ))+                             &
     &           (Ta(i-1,j+1,k)-Ta(i-1,j  ,k))*                         &
     &           (pn(i-1,j  )+pn(i-1,j+1))+                             &
     &           (Ta(i-1,j  ,k)-Ta(i-1,j-1,k))*                         &
     &           (pn(i-1,j-1)+pn(i-1,j  )))
# endif
              B=B*(on_v(i  ,j  )+on_v(i  ,j+1)+                         &
     &           on_v(i-1,j  )+on_v(i-1,j+1))/                          &
     &          (Ta(i-1,j,k)+Ta(i,j,k)+eps)
!
              Um=0.125_r8*Huon(i,j,k)*                                  &
     &           dt(ng)*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))*        &
     &           (oHz(i-1,j,k)+oHz(i,j,k))
              Vm=0.03125_r8*dt(ng)*                                     &
     &           (Hvom(i-1,j  ,k)*(pm(i-1,j)+pm(i-1,j-1))*              &
     &                            (pn(i-1,j)+pn(i-1,j-1))*              &
     &                            (oHz(i-1,j,  k)+oHz(i-1,j-1,k))+      &
     &            Hvom(i-1,j+1,k)*(pm(i-1,j+1)+pm(i-1,j))*              &
     &                            (pn(i-1,j+1)+pn(i-1,j))*              &
     &                            (oHz(i-1,j+1,k)+oHz(i-1,j  ,k))+      &
     &            Hvom(i  ,j  ,k)*(pm(i  ,j)+pm(i  ,j-1))*              &
     &                            (pn(i  ,j)+pn(i  ,j-1))*              &
     &                            (oHz(i  ,j  ,k)+oHz(i  ,j-1,k))+      &
     &            Hvom(i  ,j+1,k)*(pm(i  ,j+1)+pm(i  ,j))*              &
     &                            (pn(i  ,j+1)+pn(i  ,j))*              &
     &                            (oHz(i  ,j+1,k)+oHz(i  ,j  ,k)))
!
              X=(ABS(Um)-Um*Um)*A-B*Um*Vm-C(i,k)*Um*Wm(i,k)
              Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-C(i,k)*Vm*Wm(i,k)
              Z=(ABS(Wm(i,k))-Wm(i,k)*Wm(i,k))*C(i,k)-                  &
     &          A*Um*Wm(i,k)-B*Vm*Wm(i,k)

# ifdef MPDATA_HOT
!
              AA=A*A
              BB=B*B
              CC=C(i,k)*C(i,k)
              AB=A*B
              AC=A*C(i,k)
              BC=B*C(i,k)

              XX=X*X
              YY=Y*Y
              ZZ=Z*Z
              XY=X*Y
              XZ=X*Z
              YZ=Y*Z
# endif
!
              sig_alfa=1.0_r8/(1.0_r8-ABS(A)+eps)
# ifdef MPDATA_HOT
              sig_beta=-A/((1.0_r8-ABS(A))*                             &
     &                     (1.0_r8-AA)+eps)
              sig_gama=2.0_r8*ABS(AA*A)/((1.0_r8-ABS(A))*               &
     &                                   (1.0_r8-AA)*                   &
     &                                   (1.0_r8-ABS(AA*A))+eps)
              sig_a=-B/((1.0_r8-ABS(A))*                                &
     &                  (1.0_r8-ABS(AB))+eps)
              sig_b=AB/((1.0_r8-ABS(A))*                                &
     &                  (1.0_r8-AA*ABS(B))+eps)*                        &
     &              (ABS(B)/(1.0_r8-ABS(AB)+eps)+                       &
     &               2.0_r8*A/(1.0_r8-AA+eps))
              sig_c=ABS(A)*BB/((1.0_r8-ABS(A))*                         &
     &                         (1.0_r8-BB*ABS(A))*                      &
     &                         (1.0_r8-ABS(AB))+eps)
              sig_d=-C(i,k)/((1.0_r8-ABS(A))*                           &
     &                       (1.0_r8-ABS(AC))+eps)
              sig_e=AC/((1.0_r8-ABS(A))*                                &
     &                  (1.0_r8-AA*ABS(C(i,k)))+eps)*                   &
     &              (ABS(C(i,k))/(1.0_r8-ABS(AC)+eps)+                  &
     &               2.0_r8*A/(1.0_r8-AA+eps))
              sig_f=ABS(A)*CC/((1.0_r8-ABS(A))*                         &
     &                         (1.0_r8-CC*ABS(A))*                      &
     &                         (1.0_r8-ABS(AC))+eps)

              Ua(i,j,k)=sig_alfa*X+                                     &
     &                  sig_beta*XX+                                    &
     &                  sig_gama*XX*X+                                  &
     &                  sig_a*XY+                                       &
     &                  sig_b*XX*Y+                                     &
     &                  sig_c*X*YY+                                     &
     &                  sig_d*XZ+                                       &
     &                  sig_e*XX*Z+                                     &
     &                  sig_f*X*ZZ
# else
              Ua(i,j,k)=sig_alfa*X
# endif
!
!  Limit by physical velocity.
!
              Ua(i,j,k)=MIN(ABS(Ua(i,j,k)),fac*ABS(Um))*                &
     &                  SIGN(1.0_r8,Ua(i,j,k))
# ifdef MASKING
              Ua(i,j,k)=Ua(i,j,k)*umask(i,j)
# endif
# ifdef WET_DRY
              Ua(i,j,k)=Ua(i,j,k)*umask_wet(i,j)
# endif
            END IF
          END DO
        END DO
      END DO
!
!  Compute nondimensional V-antidiffusive velocities, Va. If applicable,
!  retain up to third-order terms of the power series.
!
      DO j=JstrVm1,Jendp2
        k=1
        DO i=IstrU-1,Iendp1
          C(i,k)=0.25_r8*                                               &
     &           ((Ta(i,j  ,k+1)-Ta(i,j  ,k  ))*odz(i,j  ,k  )+         &
     &            (Ta(i,j-1,k+1)-Ta(i,j-1,k  ))*odz(i,j-1,k  ))*        &
     &           (z_r(i,j  ,k+1)-z_r(i,j  ,k)+                          &
     &            z_r(i,j-1,k+1)-z_r(i,j-1,k))/                         &
     &           (Ta(i,j-1,k)+Ta(i,j,k)+eps)
          Wm(i,k)=0.25_r8*dt(ng)*                                       &
     &            (w(i,j-1,k  )*odz(i,j-1,k  )*pm(i,j-1)*pn(i,j-1)+     &
     &             w(i,j  ,k  )*odz(i,j  ,k  )*pm(i,j  )*pn(i,j  ))
# ifdef WEC_VF
          Wm(i,k)=Wm(i,k)+                                              &
     &            0.25_r8*dt(ng)*                                       &
     &            (W_stokes(i,j-1,k )*odz(i,j-1,k )*pm(i,j-1)*pn(i,j-1)+&
     &             W_stokes(i,j  ,k )*odz(i,j  ,k )*pm(i,j  )*pn(i,j  ))
# endif
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU-1,Iendp1
            C(i,k)=0.0625_r8*                                           &
     &             ((Ta(i,j  ,k+1)-Ta(i,j  ,k  ))*odz(i,j  ,k  )+       &
     &              (Ta(i,j  ,k  )-Ta(i,j  ,k-1))*odz(i,j  ,k-1)+       &
     &              (Ta(i,j-1,k+1)-Ta(i,j-1,k  ))*odz(i,j-1,k  )+       &
     &              (Ta(i,j-1,k  )-Ta(i,j-1,k-1))*odz(i,j-1,k-1))*      &
     &             (z_r(i,j  ,k+1)-z_r(i,j  ,k-1)+                      &
     &              z_r(i,j-1,k+1)-z_r(i,j-1,k-1))/                     &
     &             (Ta(i,j-1,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25_r8*dt(ng)*                                     &
     &              ((w(i,j-1,k-1)*odz(i,j-1,k-1)+                      &
     &                w(i,j-1,k  )*odz(i,j-1,k  ))*pm(i,j-1)*pn(i,j-1)+ &
     &               (w(i,j  ,k  )*odz(i,j  ,k  )+                      &
     &                w(i,j  ,k-1)*odz(i,j  ,k-1))*pm(i,j  )*pn(i,j  ))
# ifdef WEC_VF
          Wm(i,k)=Wm(i,k)+                                              &
     &            0.25_r8*dt(ng)*                                       &
     &              ((W_stokes(i,j-1,k-1)*odz(i,j-1,k-1)+               &
     &                W_stokes(i,j-1,k  )*odz(i,j-1,k  ))*              &
     &                                    pm(i,j-1)*pn(i,j-1)+          &
     &               (W_stokes(i,j  ,k  )*odz(i,j  ,k  )+               &
     &                W_stokes(i,j  ,k-1)*odz(i,j  ,k-1))*              &
     &                                    pm(i,j  )*pn(i,j  ))
# endif
          END DO
        END DO
        k=N(ng)
        DO i=IstrU-1,Iendp1
          C(i,k)=0.25_r8*                                               &
     &           ((Ta(i,j  ,k  )-Ta(i,j  ,k-1))*odz(i,j  ,k-1)+         &
     &            (Ta(i,j-1,k  )-Ta(i,j-1,k-1))*odz(i,j-1,k-1))*        &
     &           (z_r(i,j  ,k  )-z_r(i,j  ,k-1)+                        &
     &            z_r(i,j-1,k  )-z_r(i,j-1,k-1))/                       &
     &           (Ta(i,j-1,k)+Ta(i,j,k)+eps)
          Wm(i,k)=0.25_r8*dt(ng)*                                       &
     &            (w(i,j-1,k-1)*odz(i,j-1,k-1)*pm(i,j-1)*pn(i,j-1)+     &
     &             w(i,j  ,k-1)*odz(i,j  ,k-1)*pm(i,j  )*pn(i,j  ))
# ifdef WEC_VF
          Wm(i,k)=Wm(i,k)+                                              &
     &            0.25_r8*dt(ng)*                                       &
     &            (W_stokes(i,j-1,k-1)*odz(i,j-1,k-1)*                  &
     &                                 pm(i,j-1)*pn(i,j-1)+             &
     &             W_stokes(i,j  ,k-1)*odz(i,j  ,k-1)*                  &
     &                                 pm(i,j  )*pn(i,j  ))
# endif
        END DO
        DO k=1,N(ng)
          DO i=IstrU-1,Iendp1
            IF ((Ta(i,j-1,k).le.0.0_r8).or.                             &
     &          (Ta(i,j  ,k).le.0.0_r8).or.                             &
     &          (ABS(Ta(i,j-1,k)-Ta(i,j,k)).le.eps2)) THEN
              Va(i,j,k)=0.0_r8
            ELSE
# ifdef MASKING
              A=0.03125_r8*                                             &
     &          ((Ta(i+1,j  ,k)-Ta(i  ,j  ,k))*                         &
     &           (pm(i+1,j  )+pm(i  ,j  ))*umask(i+1,j  )+              &
     &           (Ta(i  ,j  ,k)-Ta(i-1,j  ,k))*                         &
     &           (pm(i-1,j  )+pm(i  ,j  ))*umask(i  ,j  )+              &
     &           (Ta(i+1,j-1,k)-Ta(i  ,j-1,k))*                         &
     &           (pm(i+1,j-1)+pm(i  ,j-1))*umask(i+1,j-1)+              &
     &           (Ta(i  ,j-1,k)-Ta(i-1,j-1,k))*                         &
     &           (pm(i-1,j-1)+pm(i  ,j-1))*umask(i  ,j-1))
# else
              A=0.03125_r8*                                             &
     &          ((Ta(i+1,j  ,k)-Ta(i  ,j  ,k))*                         &
     &           (pm(i+1,j  )+pm(i  ,j  ))+                             &
     &           (Ta(i  ,j  ,k)-Ta(i-1,j  ,k))*                         &
     &           (pm(i-1,j  )+pm(i  ,j  ))+                             &
     &           (Ta(i+1,j-1,k)-Ta(i  ,j-1,k))*                         &
     &           (pm(i+1,j-1)+pm(i  ,j-1))+                             &
     &           (Ta(i  ,j-1,k)-Ta(i-1,j-1,k))*                         &
     &           (pm(i-1,j-1)+pm(i  ,j-1)))
# endif
              A=A*(om_u(i  ,j  )+om_u(i+1,j  )+                         &
     &            om_u(i  ,j-1)+om_u(i+1,j-1))/                         &
     &            (Ta(i  ,j-1,k)+Ta(i,j,k)+eps)

              B=(Ta(i,j,k)-Ta(i,j-1,k))/                                &
     &          (Ta(i,j,k)+Ta(i,j-1,k)+eps)
!
              Um=0.03125_r8*dt(ng)*                                     &
     &           (Huon(i+1,j  ,k)*(pm(i+1,j)+pm(i,j))*                  &
     &                            (pn(i+1,j)+pn(i,j))*                  &
     &                            (oHz(i+1,j  ,k)+oHz(i,j  ,k))+        &
     &            Huon(i+1,j-1,k)*(pm(i+1,j-1)+pm(i,j-1))*              &
     &                            (pn(i+1,j-1)+pn(i,j-1))*              &
     &                            (oHz(i+1,j-1,k)+oHz(i,j-1,k))+        &
     &            Huon(i  ,j  ,k)*(pm(i-1,j)+pm(i,j))*                  &
     &                            (pn(i-1,j)+pn(i,j))*                  &
     &                            (oHz(i-1,j  ,k)+oHz(i,j  ,k))+        &
     &            Huon(i  ,j-1,k)*(pm(i-1,j-1)+pm(i,j-1))*              &
     &                            (pn(i-1,j-1)+pn(i,j-1))*              &
     &                            (oHz(i-1,j-1,k)+oHz(i,j-1,k)))
              Vm=0.125_r8*Hvom(i,j,k)*                                  &
     &           dt(ng)*(pn(i,j-1)+pn(i,j))*(pm(i,j-1)+pm(i,j))*        &
     &           (oHz(i,j-1,k)+oHz(i,j,k))
!
              X=(ABS(Um)-Um*Um)*A-B*Um*Vm-C(i,k)*Um*Wm(i,k)
              Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-C(i,k)*Vm*Wm(i,k)
              Z=(ABS(Wm(i,k))-Wm(i,k)*Wm(i,k))*C(i,k)-                  &
     &          A*Um*Wm(i,k)-B*Vm*Wm(i,k)

# ifdef MPDATA_HOT
!
              AA=A*A
              BB=B*B
              CC=C(i,k)*C(i,k)
              AB=A*B
              AC=A*C(i,k)
              BC=B*C(i,k)

              XX=X*X
              YY=Y*Y
              ZZ=Z*Z
              XY=X*Y
              XZ=X*Z
              YZ=Y*Z
# endif
!
              sig_alfa=1.0_r8/(1.0_r8-ABS(B)+eps)
# ifdef MPDATA_HOT
              sig_beta=-B/((1.0_r8-ABS(B))*                             &
     &                     (1.0_r8-BB)+eps)
              sig_gama=2.0_r8*ABS(BB*B)/((1.0_r8-ABS(B))*               &
     &                                   (1.0_r8-BB)*                   &
     &                                   (1.0_r8-ABS(BB*B))+eps)
              sig_a=-A/((1.0_r8-ABS(B))*                                &
     &                  (1.0_r8-ABS(AB))+eps)
              sig_b=AB/((1.0_r8-ABS(B))*                                &
     &                  (1.0_r8-BB*ABS(A))+eps)*                        &
     &              (ABS(A)/(1.0_r8-ABS(AB)+eps)+                       &
     &               2.0_r8*B/(1.0_r8-BB+eps))
              sig_c=ABS(B)*AA/((1.0_r8-ABS(B))*                         &
     &                         (1.0_r8-AA*ABS(B))*                      &
     &                         (1.0_r8-ABS(AB))+eps)
              sig_d=-C(i,k)/((1.0_r8-ABS(B))*                           &
     &                       (1.0_r8-ABS(BC))+eps)
              sig_e=BC/((1.0_r8-ABS(B))*                                &
     &                  (1.0_r8-BB*ABS(C(i,k)))+eps)*                   &
     &              (ABS(C(i,k))/(1.0_r8-ABS(BC)+eps)+                  &
     &               2.0_r8*B/(1.0_r8-BB+eps))
              sig_f=ABS(B)*CC/((1.0_r8-ABS(B))*                         &
     &                         (1.0_r8-CC*ABS(B))*                      &
     &                         (1.0_r8-ABS(BC))+eps)

              Va(i,j,k)=sig_alfa*Y+                                     &
     &                  sig_beta*YY+                                    &
     &                  sig_gama*YY*Y+                                  &
     &                  sig_a*XY+                                       &
     &                  sig_b*Y*XX+                                     &
     &                  sig_c*YY*X+                                     &
     &                  sig_d*YZ+                                       &
     &                  sig_e*YY*Z+                                     &
     &                  sig_f*Y*ZZ
# else
              Va(i,j,k)=sig_alfa*Y
# endif
!
!  Limit by physical velocity.
!
              Va(i,j,k)=MIN(ABS(Va(i,j,k)),fac*ABS(Vm))*                &
     &                  SIGN(1.0_r8,Va(i,j,k))
# ifdef MASKING
              Va(i,j,k)=Va(i,j,k)*vmask(i,j)
# endif
# ifdef WET_DRY
              Va(i,j,k)=Va(i,j,k)*vmask_wet(i,j)
# endif
            END IF
          END DO
        END DO
      END DO
!
!  Apply boundary conditions to anti-diffusive velocities.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBu3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                Ua(Istr,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                Ua(Istr,j,k)=Ua(Istr+1,j,k)
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBu3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                Ua(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstrm1,Jendp1
                Ua(Iend+1,j,k)=Ua(Iend,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBv3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                Va(i,Jstr,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                Va(i,Jstr,k)=Va(i,Jstr+1,k)
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,isBv3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                Va(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istrm1,Iendp1
                Va(i,Jend+1,k)=Va(i,Jend,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!  Compute nondimensional W-antidiffusive velocities, Wa. If applicable,
!  retain up to third-order terms of the power series.
!
      DO j=JstrV-1,Jendp1
        DO k=1,N(ng)-1
          DO i=IstrU-1,Iendp1
            IF ((Ta(i,j,k  ).le.0.0_r8).or.                             &
     &          (Ta(i,j,k+1).le.0.0_r8).or.                             &
     &          (ABS(Ta(i,j,k)-Ta(i,j,k+1)).le.eps2)) THEN
              Wa(i,j,k)=0.0_r8
            ELSE
              C(i,k)=(Ta(i,j,k+1)-Ta(i,j,k))/                           &
     &               (Ta(i,j,k+1)+Ta(i,j,k)+eps)
# ifdef MASKING
              A=0.0625_r8*                                              &
     &          ((Ta(i+1,j,k+1)-Ta(i  ,j,k+1))*                         &
     &           (pm(i+1,j  )+pm(i  ,j  ))*umask(i+1,j)+                &
     &           (Ta(i  ,j,k+1)-Ta(i-1,j,k+1))*                         &
     &           (pm(i  ,j  )+pm(i-1,j  ))*umask(i  ,j)+                &
     &           (Ta(i+1,j,k  )-Ta(i  ,j,k  ))*                         &
     &           (pm(i+1,j  )+pm(i  ,j  ))*umask(i+1,j)+                &
     &           (Ta(i  ,j,k  )-Ta(i-1,j,k  ))*                         &
     &           (pm(i  ,j  )+pm(i-1,j  ))*umask(i  ,j))
              B=0.0625_r8*                                              &
     &          ((Ta(i,j+1,k+1)-Ta(i,j  ,k+1))*                         &
     &           (pn(i  ,j+1)+pn(i  ,j  ))*vmask(i,j+1)+                &
     &           (Ta(i,j  ,k+1)-Ta(i,j-1,k+1))*                         &
     &           (pn(i  ,j  )+pn(i  ,j-1))*vmask(i,j  )+                &
     &           (Ta(i,j+1,k  )-Ta(i,j  ,k  ))*                         &
     &           (pn(i  ,j+1)+pn(i  ,j  ))*vmask(i,j+1)+                &
     &           (Ta(i,j  ,k  )-Ta(i,j-1,k  ))*                         &
     &           (pn(i  ,j  )+pn(i  ,j-1))*vmask(i,j  ))
# else
              A=0.0625_r8*                                              &
     &          ((Ta(i+1,j,k+1)-Ta(i  ,j,k+1))*                         &
     &           (pm(i+1,j  )+pm(i  ,j  ))+                             &
     &           (Ta(i  ,j,k+1)-Ta(i-1,j,k+1))*                         &
     &           (pm(i  ,j  )+pm(i-1,j  ))+                             &
     &           (Ta(i+1,j,k  )-Ta(i  ,j,k  ))*                         &
     &           (pm(i+1,j  )+pm(i  ,j  ))+                             &
     &           (Ta(i  ,j,k  )-Ta(i-1,j,k  ))*                         &
     &           (pm(i  ,j  )+pm(i-1,j  )))
              B=0.0625_r8*                                              &
     &          ((Ta(i,j+1,k+1)-Ta(i,j  ,k+1))*                         &
     &           (pn(i  ,j+1)+pn(i  ,j  ))+                             &
     &           (Ta(i,j  ,k+1)-Ta(i,j-1,k+1))*                         &
     &           (pn(i  ,j  )+pn(i  ,j-1))+                             &
     &           (Ta(i,j+1,k  )-Ta(i,j  ,k  ))*                         &
     &           (pn(i  ,j+1)+pn(i  ,j  ))+                             &
     &           (Ta(i,j  ,k  )-Ta(i,j-1,k  ))*                         &
     &           (pn(i  ,j  )+pn(i  ,j-1)))
# endif
              A=A*(om_u(i+1,j)+om_u(i  ,j))/                            &
     &            (Ta(i,j,k+1)+Ta(i,j,k)+eps)
              B=B*(on_v(i,j+1)+on_v(i,j  ))/                            &
     &            (Ta(i,j,k+1)+Ta(i,j,k)+eps)
!
              Um=0.03125_r8*dt(ng)*                                     &
     &            (Huon(i  ,j,k  )*(pm(i,j)+pm(i-1,j))*                 &
     &                             (pn(i,j)+pn(i-1,j))*                 &
     &                             (oHz(i,j,k  )+oHz(i-1,j,k  ))+       &
     &             Huon(i  ,j,k+1)*(pm(i,j)+pm(i-1,j))*                 &
     &                             (pn(i,j)+pn(i-1,j))*                 &
     &                             (oHz(i,j,k+1)+oHz(i-1,j,k+1))+       &
     &             Huon(i+1,j,k  )*(pm(i,j)+pm(i+1,j))*                 &
     &                             (pn(i,j)+pn(i+1,j))*                 &
     &                             (oHz(i,j,k  )+oHz(i+1,j,k  ))+       &
     &             Huon(i+1,j,k+1)*(pm(i,j)+pm(i+1,j))*                 &
     &                             (pn(i,j)+pn(i+1,j))*                 &
     &                             (oHz(i,j,k+1)+oHz(i+1,j,k+1)))
              Vm=0.03125_r8*dt(ng)*                                     &
     &            (Hvom(i,j  ,k  )*(pm(i,j)+pm(i,j-1))*                 &
     &                             (pn(i,j)+pn(i,j-1))*                 &
     &                             (oHz(i,j,k  )+oHz(i,j-1,k  ))+       &
     &             Hvom(i,j  ,k+1)*(pm(i,j)+pm(i,j-1))*                 &
     &                             (pn(i,j)+pn(i,j-1))*                 &
     &                             (oHz(i,j,k+1)+oHz(i,j-1,k+1))+       &
     &             Hvom(i,j+1,k  )*(pm(i,j)+pm(i,j+1))*                 &
     &                             (pn(i,j)+pn(i,j+1))*                 &
     &                             (oHz(i,j,k  )+oHz(i,j+1,k  ))+       &
     &             Hvom(i,j+1,k+1)*(pm(i,j)+pm(i,j+1))*                 &
     &                             (pn(i,j)+pn(i,j+1))*                 &
     &                             (oHz(i,j,k+1)+oHz(i,j+1,k+1)))

              Wm(i,k)=w(i,j,k)*odz(i,j,k)*pm(i,j)*pn(i,j)*dt(ng)
# ifdef WEC_VF
              Wm(i,k)=Wm(i,k)+W_stokes(i,j,k)*odz(i,j,k)*               &
     &                pm(i,j)*pn(i,j)*dt(ng)
# endif
!
              X=(ABS(Um)-Um*Um)*A-B*Um*Vm-C(i,k)*Um*Wm(i,k)
              Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-C(i,k)*Vm*Wm(i,k)
              Z=(ABS(Wm(i,k))-Wm(i,k)*Wm(i,k))*C(i,k)-                  &
     &          A*Um*Wm(i,k)-B*Vm*Wm(i,k)

# ifdef MPDATA_HOT
!
              AA=A*A
              BB=B*B
              CC=C(i,k)*C(i,k)
              AB=A*B
              AC=A*C(i,k)
              BC=B*C(i,k)

              XX=X*X
              YY=Y*Y
              ZZ=Z*Z
              XY=X*Y
              XZ=X*Z
              YZ=Y*Z
# endif
!
              sig_alfa=1.0_r8/(1.0_r8-ABS(C(i,k))+eps)
# ifdef MPDATA_HOT
              sig_beta=-C(i,k)/((1.0_r8-ABS(C(i,k)))*                   &
     &                          (1.0_r8-CC)+eps)
              sig_gama=2.0_r8*ABS(CC*C(i,k))/                           &
     &                 ((1.0_r8-ABS(C(i,k)))*                           &
     &                  (1.0_r8-CC)*                                    &
     &                  (1.0_r8-ABS(CC*C(i,k)))+eps)
              sig_a=-B/((1.0_r8-ABS(C(i,k)))*                           &
     &                  (1.0_r8-ABS(BC))+eps)
              sig_b=BC/((1.0_r8-ABS(C(i,k)))*                           &
     &                  (1.0_r8-CC*ABS(B))+eps)*                        &
     &                  (ABS(B)/(1.0_r8-ABS(BC)+eps)+                   &
     &                   2.0_r8*C(i,k)/(1.0_r8-CC+eps))
              sig_c=ABS(C(i,k))*BB/((1.0_r8-ABS(C(i,k)))*               &
     &                              (1.0_r8-B*B*ABS(C(i,k)))*           &
     &                              (1.0_r8-ABS(BC))+eps)
              sig_d=-A/((1.0_r8-ABS(C(i,k)))*                           &
     &                  (1.0_r8-ABS(AC))+eps)
              sig_e=AC/((1.0_r8-ABS(C(i,k)))*                           &
     &                  (1.0_r8-CC*ABS(A))+eps)*                        &
     &              (ABS(A)/(1.0_r8-ABS(AC)+eps)+                       &
     &               2.0_r8*C(i,k)/(1.0_r8-CC+eps))
              sig_f=ABS(C(i,k))*AA/((1.0_r8-ABS(C(i,k)))*               &
     &                              (1.0_r8-AA*ABS(C(i,k)))*            &
     &                              (1.0_r8-ABS(AC))+eps)

              Wa(i,j,k)=sig_alfa*Z+                                     &
     &                  sig_beta*ZZ+                                    &
     &                  sig_gama*ZZ*Z+                                  &
     &                  sig_a*YZ+                                       &
     &                  sig_b*ZZ*Y+                                     &
     &                  sig_c*Z*YY+                                     &
     &                  sig_d*XZ+                                       &
     &                  sig_e*ZZ*X+                                     &
     &                  sig_f*Z*XX
# else
              Wa(i,j,k)=sig_alfa*Z
# endif
!
!  Limit by physical velocity.
!
              Wa(i,j,k)=MIN(ABS(Wa(i,j,k)),fac*ABS(Wm(i,k)))*           &
     &                  SIGN(1.0_r8,Wa(i,j,k))
# ifdef MASKING
              Wa(i,j,k)=Wa(i,j,k)*rmask(i,j)
# endif
# ifdef WET_DRY
              Wa(i,j,k)=Wa(i,j,k)*rmask_wet(i,j)
# endif
            END IF
          END DO
        END DO
        DO i=IstrU-1,Iendp1
          Wa(i,j,0)=0.0_r8
          Wa(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Supress false oscillations in the solution by imposing appropriate
!  limits on the transport fluxes. Compute the UP and DOWN beta-ratios
!  described in Smolarkiewicz and Grabowski (1990).
!-----------------------------------------------------------------------
!
!  Build special mask array used to limit the UP and DOWN beta-ratios
!  to avoid including land/sea masking when computing limiting Tmin
!  and Tmax values. Notice that a zero Tmin due to land mask is not
!  considered.
!
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
# ifdef MASKING
          mask_up(i,j)=rmask(i,j)
          mask_dn(i,j)=MAX(1.0_r8,MIN(Large,(1.0_r8-rmask(i,j))*Large))
# else
          mask_up(i,j)=1.0_r8
          mask_dn(i,j)=1.0_r8
# endif
        END DO
      END DO
!
!  Compute UP and DOWN beta-ratios.
!
      DO j=JstrV-1,Jendp1
        k=1
        DO i=IstrU-1,Iendp1
          Tmax=MAX(Ta(i-1,j  ,k  )*mask_up(i-1,j  ),                    &
     &             t (i-1,j  ,k  )*mask_up(i-1,j  ),                    &
     &             Ta(i  ,j  ,k  )*mask_up(i  ,j  ),                    &
     &             t (i  ,j  ,k  )*mask_up(i  ,j  ),                    &
     &             Ta(i+1,j  ,k  )*mask_up(i+1,j  ),                    &
     &             t (i+1,j  ,k  )*mask_up(i+1,j  ),                    &
     &             Ta(i  ,j-1,k  )*mask_up(i  ,j-1),                    &
     &             t (i  ,j-1,k  )*mask_up(i  ,j-1),                    &
     &             Ta(i  ,j+1,k  )*mask_up(i  ,j+1),                    &
     &             t (i  ,j+1,k  )*mask_up(i  ,j+1),                    &
     &             Ta(i  ,j  ,k+1)*mask_up(i  ,j  ),                    &
     &             t (i  ,j  ,k+1)*mask_up(i  ,j  ))
          cff1=Ta(i-1,j  ,k  )*MAX(0.0_r8,Ua(i  ,j  ,k  ))-             &
     &         Ta(i+1,j  ,k  )*MIN(0.0_r8,Ua(i+1,j  ,k  ))+             &
     &         Ta(i  ,j-1,k  )*MAX(0.0_r8,Va(i  ,j  ,k  ))-             &
     &         Ta(i  ,j+1,k  )*MIN(0.0_r8,Va(i  ,j+1,k  ))-             &
     &         Ta(i  ,j  ,k+1)*MIN(0.0_r8,Wa(i  ,j  ,k  ))
          beta_up(i,j,k)=(Tmax-Ta(i,j,k))/(cff1+eps)
!
          Tmin=MIN(Ta(i-1,j  ,k  )*mask_dn(i-1,j  ),                    &
     &             t (i-1,j  ,k  )*mask_dn(i-1,j  ),                    &
     &             Ta(i  ,j  ,k  )*mask_dn(i  ,j  ),                    &
     &             t (i  ,j  ,k  )*mask_dn(i  ,j  ),                    &
     &             Ta(i+1,j  ,k  )*mask_dn(i+1,j  ),                    &
     &             t (i+1,j  ,k  )*mask_dn(i+1,j  ),                    &
     &             Ta(i  ,j-1,k  )*mask_dn(i  ,j-1),                    &
     &             t (i  ,j-1,k  )*mask_dn(i  ,j-1),                    &
     &             Ta(i  ,j+1,k  )*mask_dn(i  ,j+1),                    &
     &             t (i  ,j+1,k  )*mask_dn(i  ,j+1),                    &
     &             Ta(i  ,j  ,k+1)*mask_dn(i  ,j  ),                    &
     &             t (i  ,j  ,k+1)*mask_dn(i  ,j  ))
          cff2=Ta(i  ,j  ,k  )*MAX(0.0_r8,Ua(i+1,j  ,k  ))-             &
     &         Ta(i  ,j  ,k  )*MIN(0.0_r8,Ua(i  ,j  ,k  ))+             &
     &         Ta(i  ,j  ,k  )*MAX(0.0_r8,Va(i  ,j+1,k  ))-             &
     &         Ta(i  ,j  ,k  )*MIN(0.0_r8,Va(i  ,j  ,k  ))+             &
     &         Ta(i  ,j  ,k  )*MAX(0.0_r8,Wa(i  ,j  ,k  ))
          beta_dn(i,j,k)=(Ta(i,j,k)-Tmin)/(cff2+eps)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU-1,Iendp1
            Tmax=MAX(Ta(i-1,j  ,k  )*mask_up(i-1,j  ),                  &
     &               t (i-1,j  ,k  )*mask_up(i-1,j  ),                  &
     &               Ta(i  ,j  ,k  )*mask_up(i  ,j  ),                  &
     &               t (i  ,j  ,k  )*mask_up(i  ,j  ),                  &
     &               Ta(i+1,j  ,k  )*mask_up(i+1,j  ),                  &
     &               t (i+1,j  ,k  )*mask_up(i+1,j  ),                  &
     &               Ta(i  ,j-1,k  )*mask_up(i  ,j-1),                  &
     &               t (i  ,j-1,k  )*mask_up(i  ,j-1),                  &
     &               Ta(i  ,j+1,k  )*mask_up(i  ,j+1),                  &
     &               t (i  ,j+1,k  )*mask_up(i  ,j+1),                  &
     &               Ta(i  ,j  ,k-1)*mask_up(i  ,j  ),                  &
     &               t (i  ,j  ,k-1)*mask_up(i  ,j  ),                  &
     &               Ta(i  ,j  ,k+1)*mask_up(i  ,j  ),                  &
     &               t (i  ,j  ,k+1)*mask_up(i  ,j  ))
            cff1=Ta(i-1,j  ,k  )*MAX(0.0_r8,Ua(i  ,j  ,k  ))-           &
     &           Ta(i+1,j  ,k  )*MIN(0.0_r8,Ua(i+1,j  ,k  ))+           &
     &           Ta(i  ,j-1,k  )*MAX(0.0_r8,Va(i  ,j  ,k  ))-           &
     &           Ta(i  ,j+1,k  )*MIN(0.0_r8,Va(i  ,j+1,k  ))+           &
     &           Ta(i  ,j  ,k-1)*MAX(0.0_r8,Wa(i  ,j  ,k-1))-           &
     &           Ta(i  ,j  ,k+1)*MIN(0.0_r8,Wa(i  ,j  ,k  ))
            beta_up(i,j,k)=(Tmax-Ta(i,j,k))/(cff1+eps)
!
            Tmin=MIN(Ta(i-1,j  ,k  )*mask_dn(i-1,j  ),                  &
     &               t (i-1,j  ,k  )*mask_dn(i-1,j  ),                  &
     &               Ta(i  ,j  ,k  )*mask_dn(i  ,j  ),                  &
     &               t (i  ,j  ,k  )*mask_dn(i  ,j  ),                  &
     &               Ta(i+1,j  ,k  )*mask_dn(i+1,j  ),                  &
     &               t (i+1,j  ,k  )*mask_dn(i+1,j  ),                  &
     &               Ta(i  ,j-1,k  )*mask_dn(i  ,j-1),                  &
     &               t (i  ,j-1,k  )*mask_dn(i  ,j-1),                  &
     &               Ta(i  ,j+1,k  )*mask_dn(i  ,j+1),                  &
     &               t (i  ,j+1,k  )*mask_dn(i  ,j+1),                  &
     &               Ta(i  ,j  ,k-1)*mask_dn(i  ,j  ),                  &
     &               t (i  ,j  ,k-1)*mask_dn(i  ,j  ),                  &
     &               Ta(i  ,j  ,k+1)*mask_dn(i  ,j  ),                  &
     &               t (i  ,j  ,k+1)*mask_dn(i  ,j  ))
            cff2=Ta(i  ,j  ,k  )*MAX(0.0_r8,Ua(i+1,j  ,k  ))-           &
     &           Ta(i  ,j  ,k  )*MIN(0.0_r8,Ua(i  ,j  ,k  ))+           &
     &           Ta(i  ,j  ,k  )*MAX(0.0_r8,Va(i  ,j+1,k  ))-           &
     &           Ta(i  ,j  ,k  )*MIN(0.0_r8,Va(i  ,j  ,k  ))+           &
     &           Ta(i  ,j  ,k  )*MAX(0.0_r8,Wa(i  ,j  ,k  ))-           &
     &           Ta(i  ,j  ,k  )*MIN(0.0_r8,Wa(i  ,j  ,k-1))
            beta_dn(i,j,k)=(Ta(i,j,k)-Tmin)/(cff2+eps)
          END DO
        END DO
        k=N(ng)
        DO i=IstrU-1,Iendp1
          Tmax=MAX(Ta(i-1,j  ,k  )*mask_up(i-1,j  ),                    &
     &             t (i-1,j  ,k  )*mask_up(i-1,j  ),                    &
     &             Ta(i  ,j  ,k  )*mask_up(i  ,j  ),                    &
     &             t (i  ,j  ,k  )*mask_up(i  ,j  ),                    &
     &             Ta(i+1,j  ,k  )*mask_up(i+1,j  ),                    &
     &             t (i+1,j  ,k  )*mask_up(i+1,j  ),                    &
     &             Ta(i  ,j-1,k  )*mask_up(i  ,j-1),                    &
     &             t (i  ,j-1,k  )*mask_up(i  ,j-1),                    &
     &             Ta(i  ,j+1,k  )*mask_up(i  ,j+1),                    &
     &             t (i  ,j+1,k  )*mask_up(i  ,j+1),                    &
     &             Ta(i  ,j  ,k-1)*mask_up(i  ,j  ),                    &
     &             t (i  ,j  ,k-1)*mask_up(i  ,j  ))
          cff1=Ta(i-1,j  ,k  )*MAX(0.0_r8,Ua(i  ,j  ,k  ))-             &
     &         Ta(i+1,j  ,k  )*MIN(0.0_r8,Ua(i+1,j  ,k  ))+             &
     &         Ta(i  ,j-1,k  )*MAX(0.0_r8,Va(i  ,j  ,k  ))-             &
     &         Ta(i  ,j+1,k  )*MIN(0.0_r8,Va(i  ,j+1,k  ))+             &
     &         Ta(i  ,j  ,k-1)*MAX(0.0_r8,Wa(i  ,j  ,k-1))
          beta_up(i,j,k)=(Tmax-Ta(i,j,k))/(cff1+eps)
!
          Tmin=MIN(Ta(i-1,j  ,k  )*mask_dn(i-1,j  ),                    &
     &             t (i-1,j  ,k  )*mask_dn(i-1,j  ),                    &
     &             Ta(i  ,j  ,k  )*mask_dn(i  ,j  ),                    &
     &             t (i  ,j  ,k  )*mask_dn(i  ,j  ),                    &
     &             Ta(i+1,j  ,k  )*mask_dn(i+1,j  ),                    &
     &             t (i+1,j  ,k  )*mask_dn(i+1,j  ),                    &
     &             Ta(i  ,j-1,k  )*mask_dn(i  ,j-1),                    &
     &             t (i  ,j-1,k  )*mask_dn(i  ,j-1),                    &
     &             Ta(i  ,j+1,k  )*mask_dn(i  ,j+1),                    &
     &             t (i  ,j+1,k  )*mask_dn(i  ,j+1),                    &
     &             Ta(i  ,j  ,k-1)*mask_dn(i  ,j  ),                    &
     &             t (i  ,j  ,k-1)*mask_dn(i  ,j  ))
          cff2=Ta(i  ,j  ,k  )*MAX(0.0_r8,Ua(i+1,j  ,k  ))-             &
     &         Ta(i  ,j  ,k  )*MIN(0.0_r8,Ua(i  ,j  ,k  ))+             &
     &         Ta(i  ,j  ,k  )*MAX(0.0_r8,Va(i  ,j+1,k  ))-             &
     &         Ta(i  ,j  ,k  )*MIN(0.0_r8,Va(i  ,j  ,k  ))-             &
     &         Ta(i  ,j  ,k  )*MIN(0.0_r8,Wa(i  ,j  ,k-1))
          beta_dn(i,j,k)=(Ta(i,j,k)-Tmin)/(cff2+eps)
        END DO
      END DO
      DO k=1,N(ng)
        DO j=JstrV-1,Jendp1
          DO i=IstrU-1,Iendp1
            IF (mask_up(i,j).eq.0.0_r8) THEN
              beta_up(i,j,k)=2.0_r8
              beta_dn(i,j,k)=2.0_r8
            END IF
          END DO
        END DO
      END DO
!
!  Calculate monotonic velocities. Scale back to dimensional units.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iendp1
            cff1=MIN(beta_dn(i-1,j,k),beta_up(i,j,k),1.0_r8)
            cff2=MIN(beta_up(i-1,j,k),beta_dn(i,j,k),1.0_r8)
            Ua(i,j,k)=(cff1*MAX(0.0_r8,Ua(i,j,k))+                      &
     &                 cff2*MIN(0.0_r8,Ua(i,j,k)))*                     &
     &                 cff*om_u(i,j)
# ifdef MASKING
            Ua(i,j,k)=Ua(i,j,k)*umask(i,j)
# endif
# ifdef WET_DRY
            Ua(i,j,k)=Ua(i,j,k)*umask_wet(i,j)
# endif
          END DO
        END DO
        DO j=JstrV,Jendp1
          DO i=Istr,Iend
            cff1=MIN(beta_dn(i,j-1,k),beta_up(i,j,k),1.0_r8)
            cff2=MIN(beta_up(i,j-1,k),beta_dn(i,j,k),1.0_r8)
            Va(i,j,k)=(cff1*MAX(0.0_r8,Va(i,j,k))+                      &
     &                 cff2*MIN(0.0_r8,Va(i,j,k)))*                     &
     &                 cff*on_v(i,j)
# ifdef MASKING
            Va(i,j,k)=Va(i,j,k)*vmask(i,j)
# endif
# ifdef WET_DRY
            Va(i,j,k)=Va(i,j,k)*vmask_wet(i,j)
# endif
          END DO
        END DO
        IF (k.lt.N(ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff1=MIN(beta_dn(i,j,k),beta_up(i,j,k+1),1.0_r8)
              cff2=MIN(beta_up(i,j,k),beta_dn(i,j,k+1),1.0_r8)
              Wa(i,j,k)=(cff1*MAX(0.0_r8,Wa(i,j,k))+                    &
     &                   cff2*MIN(0.0_r8,Wa(i,j,k)))*                   &
     &                   cff*omn(i,j)*(z_r(i,j,k+1)-z_r(i,j,k))
# ifdef MASKING
              Wa(i,j,k)=Wa(i,j,k)*rmask(i,j)
# endif
# ifdef WET_DRY
              Wa(i,j,k)=Wa(i,j,k)*rmask_wet(i,j)
# endif
            END DO
          END DO
        END IF
      END DO
!
!  Apply boundary conditions to anti-diffusive velocities.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBu3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Ua(Istr,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Ua(Istr,j,k)=Ua(Istr+1,j,k)
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBu3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Ua(Iend+1,j,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Ua(Iend+1,j,k)=Ua(Iend,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBv3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Va(i,Jstr,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istr,Iend
                Va(i,Jstr,k)=Va(i,Jstr+1,k)
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,isBv3d,ng)%closed) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Va(i,Jend+1,k)=0.0_r8
              END DO
            END DO
          ELSE
            DO k=1,N(ng)
              DO i=Istr,Iend
                Va(i,Jend+1,k)=Va(i,Jend,k)
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE mpdata_adiff_tile
#endif
      END MODULE mpdata_adiff_mod
