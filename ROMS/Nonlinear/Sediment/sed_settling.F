#include "cppdefs.h"
#undef NEUMANN
#undef LINEAR_CONTINUATION

      MODULE sed_settling_mod

#if defined NONLINEAR && defined SEDIMENT && defined SUSPLOAD
!
!svn $Id: sed_settling.F 429 2009-12-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes the vertical settling (sinking) of suspended  !
!  sediment via a semi-Lagrangian advective flux algorithm. It uses a  !
!  parabolic,  vertical reconstructuion of the suspended  sediment in  !
!  the water column with PPT/WENO constraints to avoid oscillations.   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Colella, P. and P. Woodward, 1984: The piecewise parabolic method   !
!    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201. !
!                                                                      !
!  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially        !
!    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_settling

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_settling (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
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
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL sed_settling_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % z_w,                           &
     &                        SEDBED(ng) % settling_flux,               &
     &                        OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_settling
!
!***********************************************************************
      SUBROUTINE sed_settling_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              Hz, z_w,                            &
     &                              settling_flux,                      &
     &                              t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: settling_flux(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: settling_flux(LBi:UBi,LBj:UBj,NST(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, indx, ised, j, k, ks

      real(r8) :: cff, cu, cffL, cffR, dltL, dltR

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add sediment vertical sinking (settling) term.
!-----------------------------------------------------------------------
!
!  Compute inverse thicknesses to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Copy concentration of suspended sediment into scratch array "qc"
!  (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for sediment
!  concentration.
!
        SED_LOOP: DO ised=1,NST
          indx=idsed(ised)
          DO k=1,N(ng)
            DO i=Istr,Iend
              qc(i,k)=t(i,j,k,nstp,indx)

            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
            END DO
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              dltR=Hz(i,j,k)*FC(i,k)
              dltL=Hz(i,j,k)*FC(i,k-1)
              cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
              cffR=cff*FC(i,k)
              cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
              IF ((dltR*dltL).le.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
              cff=(dltR-dltL)*Hz_inv3(i,k)
              dltR=dltR-cff*Hz(i,j,k+1)
              dltL=dltL+cff*Hz(i,j,k-1)
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
              WR(i,k)=(2.0_r8*dltR-dltL)**2
              WL(i,k)=(dltR-2.0_r8*dltL)**2
            END DO
          END DO
          cff=1.0E-14_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              dltL=MAX(cff,WL(i,k  ))
              dltR=MAX(cff,WR(i,k+1))
              qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
              qL(i,k+1)=qR(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
# if defined LINEAR_CONTINUATION
            qL(i,N(ng))=qR(i,N(ng)-1)
            qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
# elif defined NEUMANN
            qL(i,N(ng))=qR(i,N(ng)-1)
            qR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
# else
            qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
            qL(i,N(ng))=qc(i,N(ng))         ! conditions
            qR(i,N(ng)-1)=qc(i,N(ng))
# endif
# if defined LINEAR_CONTINUATION
            qR(i,1)=qL(i,2)
            qL(i,1)=2.0_r8*qc(i,1)-qR(i,1)
# elif defined NEUMANN
            qR(i,1)=qL(i,2)
            qL(i,1)=1.5_r8*qc(i,1)-0.5_r8*qR(i,1)
# else
            qL(i,2)=qc(i,1)                 ! bottom grid boxes are
            qR(i,1)=qc(i,1)                 ! re-assumed to be
            qL(i,1)=qc(i,1)                 ! piecewise constant.
# endif
          END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              dltR=qR(i,k)-qc(i,k)
              dltL=qc(i,k)-qL(i,k)
              cffR=2.0_r8*dltR
              cffL=2.0_r8*dltL
              IF ((dltR*dltL).lt.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
            END DO
          END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
          cff=dt(ng)*ABS(Wsed(ised,ng))
          DO k=1,N(ng)
            DO i=Istr,Iend
              FC(i,k-1)=0.0_r8
              WL(i,k)=z_w(i,j,k-1)+cff
              WR(i,k)=Hz(i,j,k)*qc(i,k)
              ksource(i,k)=k
            END DO
          END DO
          DO k=1,N(ng)
            DO ks=k,N(ng)-1
              DO i=Istr,Iend
                IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                  ksource(i,k)=ks+1
                  FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                END IF
              END DO
            END DO
          END DO
!
!  Finalize computation of flux: add fractional part.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              ks=ksource(i,k)
              cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
              FC(i,k-1)=FC(i,k-1)+                                      &
     &                  Hz(i,j,ks)*cu*                                  &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
            END DO
          END DO
          DO i=Istr,Iend
            DO k=1,N(ng)
              t(i,j,k,nnew,indx)=t(i,j,k,nnew,indx)+(FC(i,k)-FC(i,k-1))
            END DO
            settling_flux(i,j,ised)=FC(i,0)
          END DO
        END DO SED_LOOP
      END DO J_LOOP

      RETURN
      END SUBROUTINE sed_settling_tile
#endif
      END MODULE sed_settling_mod
