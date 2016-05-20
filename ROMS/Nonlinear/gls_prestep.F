#include "cppdefs.h"
      MODULE gls_prestep_mod
#if defined NONLINEAR && defined GLS_MIXING && defined SOLVE3D
!
!svn $Id: gls_prestep.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine perfoms the predictor step for turbulent kinetic       !
!  energy prognostic variables, tke and gls. A NON-conservative,       !
!  but constancy preserving, auxiliary advective substep for tke       !
!  gls equations is carried out. The result of this substep will       !
!  be used to compute advective terms in the corrector substep.        !
!  No dissipation terms are included here.                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: gls_prestep

      CONTAINS
!
!***********************************************************************
      SUBROUTINE gls_prestep (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_mixing
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
      CALL wclock_on (ng, iNLM, 19)
# endif
      CALL gls_prestep_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
# ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % Huon,                           &
     &                       GRID(ng) % Hvom,                           &
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       OCEAN(ng) % W,                             &
# ifdef WEC_VF
     &                       OCEAN(ng) % W_stokes,                      &
# endif
     &                       MIXING(ng) % gls,                          &
     &                       MIXING(ng) % tke)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 19)
# endif
      RETURN
      END SUBROUTINE gls_prestep
!
!***********************************************************************
      SUBROUTINE gls_prestep_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             Huon, Hvom, Hz, pm, pn, W,           &
# ifdef WEC_VF
     &                             W_stokes,                            &
# endif
     &                             gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
      USE tkebc_mod, ONLY : tkebc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
# ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
# endif
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
# ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
# endif
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:N(ng),3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:N(ng),3)
# endif
!
!  Local variable declarations.
!
      integer :: i, indx, j, k

      real(r8), parameter :: Gamma = 1.0_r8/6.0_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(IminS:ImaxS,N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: FCL

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Hz_half

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: XF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradL

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Predictor step for advection of turbulent kinetic energy variables.
!-----------------------------------------------------------------------
!
! Start computation of auxiliary time step fields tke(:,:,:,n+1/2) and
! gls(:,:,:,n+1/2) with computation of horizontal advection terms and
! auxiliary grid-box height field Hz_new()=Hz(:,:,k+1/2,n+1/2);
! This is effectivey an LF step with subsequent interpolation of the
! result half step back, using AM3 weights. The LF step and
! interpolation are perfomed as a single operation, which results in
! weights cff1,cff2,cff3 below.
!
! Either centered fourth-order accurate or standard second order
! accurate versions are supported.
!
! At the same time prepare for corrector step for tke,gls: set tke,
! gls(:,:,:,nnew) to  tke,gls(:,:,:,nstp) multiplied by the
! corresponding grid-box height. This needs done at this time because
! array Hz(:,:,:) will overwritten after 2D time stepping with the
! values computed from zeta(:,:,n+1) rather than zeta(:,:,n), so that
! the old-time-step Hz will be no longer awailable.
!

      DO k=1,N(ng)-1
# ifdef K_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            XF(i,j)=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FX (i,j)=XF(i,j)*                                           &
     &               0.5_r8*(tke(i,j,k,nstp)+tke(i-1,j,k,nstp))
            FXL(i,j)=XF(i,j)*                                           &
     &               0.5_r8*(gls(i,j,k,nstp)+gls(i-1,j,k,nstp))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            EF(i,j)=0.5*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FE (i,j)=EF(i,j)*                                           &
     &               0.5*(tke(i,j,k,nstp)+tke(i,j-1,k,nstp))
            FEL(i,j)=EF(i,j)*                                           &
     &               0.5*(gls(i,j,k,nstp)+gls(i,j-1,k,nstp))
          END DO
        END DO
# else
!
!  Fourth-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=Istrm1,Iendp2
            grad (i,j)=(tke(i,j,k,nstp)-tke(i-1,j,k,nstp))
#  ifdef MASKING
            grad (i,j)=grad (i,j)*umask(i,j)
#  endif
            gradL(i,j)=(gls(i,j,k,nstp)-gls(i-1,j,k,nstp))
#  ifdef MASKING
            gradL(i,j)=gradL(i,j)*umask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              grad (Istr-1,j)=grad (Istr,j)
              gradL(Istr-1,j)=gradL(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              grad (Iend+2,j)=grad (Iend+1,j)
              gradL(Iend+2,j)=gradL(Iend+1,j)
            END DO
          END IF
        END IF
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            XF(i,j)=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FX (i,j)=XF(i,j)*                                           &
     &               0.5_r8*(tke(i-1,j,k,nstp)+tke(i,j,k,nstp)-         &
     &                       cff*(grad (i+1,j)-grad (i-1,j)))
            FXL(i,j)=XF(i,j)*                                           &
     &               0.5_r8*(gls(i-1,j,k,nstp)+gls(i,j,k,nstp)-         &
     &                       cff*(gradL(i+1,j)-gradL(i-1,j)))
          END DO
        END DO
!
        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            grad (i,j)=(tke(i,j,k,nstp)-tke(i,j-1,k,nstp))
#  ifdef MASKING
            grad (i,j)=grad (i,j)*vmask(i,j)
#  endif
            gradL(i,j)=(gls(i,j,k,nstp)-gls(i,j-1,k,nstp))
#  ifdef MASKING
            gradL(i,j)=gradL(i,j)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              grad (i,Jstr-1)=grad (i,Jstr)
              gradL(i,Jstr-1)=gradL(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              grad (i,Jend+2)=grad (i,Jend+1)
              gradL(i,Jend+2)=gradL(i,Jend+1)
            END DO
          END IF
        END IF
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            EF(i,j)=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FE (i,j)=EF(i,j)*                                           &
     &               0.5_r8*(tke(i,j-1,k,nstp)+tke(i,j,k,nstp)-         &
     &                       cff*(grad (i,j+1)-grad (i,j-1)))
            FEL(i,j)=EF(i,j)*                                           &
     &               0.5_r8*(gls(i,j-1,k,nstp)+gls(i,j,k,nstp)-         &
     &                       cff*(gradL(i,j+1)-gradL(i,j-1)))
          END DO
        END DO
# endif
!
!  Time-step horizontal advection.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff1=1.0_r8
          cff2=0.0_r8
          cff3=0.5_r8*dt(ng)
          indx=nstp
        ELSE
          cff1=0.5_r8+Gamma
          cff2=0.5_r8-Gamma
          cff3=(1.0_r8-Gamma)*dt(ng)
          indx=3-nstp
        END IF
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
            cff4=cff3*pm(i,j)*pn(i,j)
            Hz_half(i,j,k)=cff-cff4*(XF(i+1,j)-XF(i,j)+                 &
     &                               EF(i,j+1)-EF(i,j))
            tke(i,j,k,3)=cff*(cff1*tke(i,j,k,nstp)+                     &
     &                        cff2*tke(i,j,k,indx))-                    &
     &                   cff4*(FX (i+1,j)-FX (i,j)+                     &
     &                         FE (i,j+1)-FE (i,j))
            gls(i,j,k,3)=cff*(cff1*gls(i,j,k,nstp)+                     &
     &                        cff2*gls(i,j,k,indx))-                    &
     &                   cff4*(FXL(i+1,j)-FXL(i,j)+                     &
     &                         FEL(i,j+1)-FEL(i,j))
            tke(i,j,k,nnew)=cff*tke(i,j,k,nstp)
            gls(i,j,k,nnew)=cff*gls(i,j,k,nstp)
          END DO
        END DO
      END DO
!
! Compute vertical advection term.
!
      DO j=Jstr,Jend
# ifdef K_C2ADVECTION
        DO k=1,N(ng)
          DO i=Istr,Iend
            CF(i,k)=0.5_r8*(W(i,j,k)+W(i,j,k-1))
#  ifdef WEC_VF
            CF(i,k)=CF(i,k)+0.5_r8*(W_stokes(i,j,k)+W_stokes(i,j,k-1))
#  endif
            FC (i,k)=CF(i,k)*                                           &
     &               0.5_r8*(tke(i,j,k-1,nstp)+tke(i,j,k,nstp))
            FCL(i,k)=CF(i,k)*                                           &
     &               0.5_r8*(gls(i,j,k-1,nstp)+gls(i,j,k,nstp))
          END DO
        END DO
# else
        cff1=7.0_r8/12.0_r8
        cff2=1.0_r8/12.0_r8
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            CF(i,k)=0.5_r8*(W(i,j,k)+W(i,j,k-1))
#  ifdef WEC_VF
            CF(i,k)=CF(i,k)+0.5_r8*(W_stokes(i,j,k)+W_stokes(i,j,k-1))
#  endif
            FC (i,k)=CF(i,k)*(cff1*(tke(i,j,k-1,nstp)+                  &
     &                              tke(i,j,k  ,nstp))-                 &
     &                        cff2*(tke(i,j,k-2,nstp)+                  &
     &                              tke(i,j,k+1,nstp)))
            FCL(i,k)=CF(i,k)*(cff1*(gls(i,j,k-1,nstp)+                  &
     &                              gls(i,j,k  ,nstp))-                 &
     &                        cff2*(gls(i,j,k-2,nstp)+                  &
     &                              gls(i,j,k+1,nstp)))
          END DO
        END DO
        cff1=1.0_r8/3.0_r8
        cff2=5.0_r8/6.0_r8
        cff3=1.0_r8/6.0_r8
        DO i=Istr,Iend
          CF(i,1)=0.5*(W(i,j,0)+W(i,j,1))
#  ifdef WEC_VF
          CF(i,1)=CF(i,1)+0.5_r8*(W_stokes(i,j,0)+W_stokes(i,j,1))
#  endif
          FC (i,1)=CF(i,1)*(cff1*tke(i,j,0,nstp)+                       &
     &                      cff2*tke(i,j,1,nstp)-                       &
     &                      cff3*tke(i,j,2,nstp))
          FCL(i,1)=CF(i,1)*(cff1*gls(i,j,0,nstp)+                       &
     &                      cff2*gls(i,j,1,nstp)-                       &
     &                      cff3*gls(i,j,2,nstp))
          CF(i,N(ng))=0.5*(W(i,j,N(ng))+W(i,j,N(ng)-1))
#  ifdef WEC_VF
          CF(i,N(ng))=CF(i,N(ng))+0.5_r8*                               &
     &                (W_stokes(i,j,N(ng))+W_stokes(i,j,N(ng)-1))
#  endif
          FC (i,N(ng))=CF(i,N(ng))*(cff1*tke(i,j,N(ng)  ,nstp)+         &
     &                              cff2*tke(i,j,N(ng)-1,nstp)-         &
     &                              cff3*tke(i,j,N(ng)-2,nstp))
          FCL(i,N(ng))=CF(i,N(ng))*(cff1*gls(i,j,N(ng)  ,nstp)+         &
     &                              cff2*gls(i,j,N(ng)-1,nstp)-         &
     &                              cff3*gls(i,j,N(ng)-2,nstp))
        END DO
# endif
!
!  Time-step vertical advection term.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff3=0.5_r8*dt(ng)
        ELSE
          cff3=(1.0_r8-Gamma)*dt(ng)
        END IF
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff4=cff3*pm(i,j)*pn(i,j)
            Hz_half(i,j,k)=Hz_half(i,j,k)-cff4*(CF(i,k+1)-CF(i,k))
            cff1=1.0_r8/Hz_half(i,j,k)
            tke(i,j,k,3)=cff1*(tke(i,j,k,3)-                            &
     &                         cff4*(FC (i,k+1)-FC (i,k)))
            gls(i,j,k,3)=cff1*(gls(i,j,k,3)-                            &
     &                         cff4*(FCL(i,k+1)-FCL(i,k)))
          END DO
        END DO
      END DO
!
!  Apply lateral boundary conditions.
!
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 3, nstp,                                         &
     &                 gls, tke)

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tke(:,:,:,3))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          gls(:,:,:,3))
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tke(:,:,:,3),                                 &
     &                    gls(:,:,:,3))
# endif

      RETURN
      END SUBROUTINE gls_prestep_tile
#endif
      END MODULE gls_prestep_mod
