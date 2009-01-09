#include "cppdefs.h"
#ifdef SOLVE3D
      SUBROUTINE set_scoord (ng)
!
!svn $Id: set_scoord.F 1601 2008-07-25 15:02:43Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine sets and initializes relevant variables associated     !
!  with the vertical transformation to nondimensional S-coordinate.    !
!                                                                      !
!  Reference:                                                          !
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
# if defined FOUR_DVAR && defined N2NORM_PROFILE
      USE mod_fourdvar
# endif
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
      integer :: k, order
      real(r8) :: cff, cff1, cff2, cff3
!
!-----------------------------------------------------------------------
!  Define S-Curves in domain [-1 < sc < 0] at vertical W- and
!  RHO-points.
!-----------------------------------------------------------------------
!
# if defined WET_DRY
      hc(ng)=MIN(MAX(hmin(ng),0.0_r8),Tcline(ng))
# else
      hc(ng)=MIN(hmin(ng),Tcline(ng))
# endif
      IF (theta_s(ng).ne.0.0_r8) THEN
        cff1=1.0_r8/SINH(theta_s(ng))
        cff2=0.5_r8/TANH(0.5_r8*theta_s(ng))
      END IF
      SCALARS(ng)%sc_w(0)=-1.0_r8
      SCALARS(ng)%Cs_w(0)=-1.0_r8
      cff=1.0_r8/REAL(N(ng),r8)
      DO k=1,N(ng)
        SCALARS(ng)%sc_w(k)=cff*REAL(k-N(ng),r8)
        SCALARS(ng)%sc_r(k)=cff*(REAL(k-N(ng),r8)-0.5_r8)
        IF (theta_s(ng).ne.0.0_r8) THEN
          SCALARS(ng)%Cs_w(k)=(1.0_r8-theta_b(ng))*                     &
     &                        cff1*SINH(theta_s(ng)*                    &
     &                                  SCALARS(ng)%sc_w(k))+           &
     &                        theta_b(ng)*                              &
     &                        (cff2*TANH(theta_s(ng)*                   &
     &                                   (SCALARS(ng)%sc_w(k)+0.5_r8))- &
     &                         0.5_r8)
          SCALARS(ng)%Cs_r(k)=(1.0_r8-theta_b(ng))*                     &
     &                        cff1*SINH(theta_s(ng)*                    &
     &                                  SCALARS(ng)%sc_r(k))+           &
     &                        theta_b(ng)*                              &
     &                        (cff2*TANH(theta_s(ng)*                   &
     &                                   (SCALARS(ng)%sc_r(k)+0.5_r8))- &
     &                         0.5_r8)
        ELSE
          SCALARS(ng)%Cs_w(k)=SCALARS(ng)%sc_w(k)
          SCALARS(ng)%Cs_r(k)=SCALARS(ng)%sc_r(k)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report information about vertical S-levels.
!-----------------------------------------------------------------------
!
      IF (Master.and.LwrtInfo(ng)) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' Vertical S-coordinate System: ',/,/,                &
     &            ' level   S-coord     Cs-curve',10x,                  &
     &            'at_hmin  over_slope     at_hmax',/)
        DO k=N(ng),0,-1
          cff1=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (hmin(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
          cff2=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (0.5*(hmin(ng)+hmax(ng))-hc(ng))*SCALARS(ng)%Cs_w(k)
          cff3=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (hmax(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
          WRITE (stdout,20) k, SCALARS(ng)%sc_w(k),                     &
     &                         SCALARS(ng)%Cs_w(k), cff1, cff2, cff3
  20      FORMAT (i6,2f12.7,4x,3f12.3)
        END DO
      END IF
# if defined FOUR_DVAR && defined N2NORM_PROFILE
!
!-----------------------------------------------------------------------
!  Report N2(z) profile for cost function normalization.
!-----------------------------------------------------------------------
!
      IF (Master.and.LwrtInfo(ng)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' N2 (1/s2) vertical profile used for',               &
     &            ' normalization:',/)
        DO k=N(ng),0,-1
          cff1=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (hmax(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
          cff2=pcoef_N2(npN2)
          DO order=npN2-1,0,-1
            cff2=cff1*cff2+pcoef_N2(order)
          END DO
          WRITE (stdout,40) k, cff1, ABS(cff2)
  40      FORMAT (I6,F12.3,3x,1p,e21.14,0p)
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE set_scoord
#else
      SUBROUTINE set_scoord
      RETURN
      END SUBROUTINE set_scoord
#endif

