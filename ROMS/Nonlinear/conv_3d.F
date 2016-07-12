#include "cppdefs.h"

      MODULE conv_3d_mod

#if defined NONLINEAR && defined FOUR_DVAR && defined SOLVE3D
!
!svn $Id: conv_3d.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines applies the background error covariance to data      !
!  assimilation fields via the space convolution of the diffusion      !
!  equation (filter) for 3D state variables. The diffusion filter      !
!  is solved using an implicit or explicit vertical algorithm.         !
!                                                                      !
!  For Gaussian (bell-shaped) correlations, the space convolution      !
!  of the diffusion operator is an efficient way  to estimate the      !
!  finite domain error covariances.                                    !
!                                                                      !
!  Notice that "z_r" and "Hz" are assumed to be time invariant in      !
!  the spatial convolution. That it, they are not adjointable.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     LBk        K-dimension lower bound.                              !
!     UBk        K-dimension upper bound.                              !
!     Nghost     Number of ghost points.                               !
!     NHsteps    Number of horizontal diffusion integration steps.     !
!     NVsteps    Number of vertical   diffusion integration steps.     !
!     DTsizeH    Horizontal diffusion pseudo time-step size.           !
!     DTsizeV    Vertical   diffusion pseudo time-step size.           !
!     Kh         Horizontal diffusion coefficients.                    !
!     Kv         Vertical   diffusion coefficients.                    !
!     A          3D state variable to diffuse.                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Convolved 3D state variable.                          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    conv_r3d_tile     Nonlinear 3D convolution at RHO-points          !
!    conv_u3d_tile     Nonlinear 3D convolution at U-points            !
!    conv_v3d_tile     Nonlinear 3D convolution at V-points            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PUBLIC
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE conv_r3d_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Nghost, NHsteps, NVsteps,               &
     &                          DTsizeH, DTsizeV,                       &
     &                          Kh, Kv,                                 &
     &                          pm, pn,                                 &
# ifdef GEOPOTENTIAL_HCONV
     &                          on_u, om_v,                             &
# else
     &                          pmon_u, pnom_v,                         &
# endif
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
     &                          Hz, z_r,                                &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY: dabc_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps, NVsteps

      real(r8), intent(in) :: DTsizeH, DTsizeV
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
#  else
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, k, k1, k2, step

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk,2) :: Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac

# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
#  endif
#  if defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
#   endif
#  else
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FS
#  endif
# endif
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 3D state variable
!  at RHO-points.
!-----------------------------------------------------------------------
!
!  Compute metrics factors.  Notice that "z_r" and "Hz" are assumed to
!  be time invariant in the vertical convolution.  Scratch array are
!  used for efficiency.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=DTsizeH*pm(i,j)*pn(i,j)
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
            FC(i,j,k)=-DTsizeV*Kv(i,j,k)/(z_r(i,j,k+1)-z_r(i,j,k))
#   else
            FC(i,j,k)=DTsizeV*Kv(i,j,k)/(z_r(i,j,k+1)-z_r(i,j,k))
#   endif
          END DO
          FC(i,j,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
#  endif
# endif
        END DO
      END DO
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
      CALL dabc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    A)
# endif
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            Awrk(i,j,k,Nold)=A(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps

# ifdef GEOPOTENTIAL_HCONV
!
!  Diffusion along geopotential surfaces: Compute horizontal and
!  vertical gradients.  Notice the recursive blocking sequence.  The
!  vertical placement of the gradients is:
!
!        dAdx,dAde(:,:,k1) k     rho-points
!        dAdx,dAde(:,:,k2) k+1   rho-points
!          FZ,dAdz(:,:,k1) k-1/2   W-points
!          FZ,dAdz(:,:,k2) k+1/2   W-points
!
        k2=1
        K_LOOP : DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#  ifdef MASKING
                cff=cff*umask(i,j)
#  endif
                dZdx(i,j,k2)=cff*(z_r(i  ,j,k+1)-                       &
     &                            z_r(i-1,j,k+1))
#  ifdef MASKING
                dAdx(i,j,k2)=cff*(Awrk(i  ,j,k+1,Nold)*rmask(i  ,j)-    &
     &                            Awrk(i-1,j,k+1,Nold)*rmask(i-1,j))
#  else
                dAdx(i,j,k2)=cff*(Awrk(i  ,j,k+1,Nold)-                 &
     &                            Awrk(i-1,j,k+1,Nold))
#  endif
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#  ifdef MASKING
                cff=cff*vmask(i,j)
#  endif
                dZde(i,j,k2)=cff*(z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j-1,k+1))
#  ifdef MASKING
                dAde(i,j,k2)=cff*(Awrk(i,j  ,k+1,Nold)*rmask(i,j  )-    &
     &                            Awrk(i,j-1,k+1,Nold)*rmask(i,j-1))
#  else
                dAde(i,j,k2)=cff*(Awrk(i,j  ,k+1,Nold)-                 &
     &                            Awrk(i,j-1,k+1,Nold))
#  endif
              END DO
            END DO
          END IF
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                dAdz(i,j,k2)=0.0_r8
                FZ(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
                dAdz(i,j,k2)=cff*(Awrk(i,j,k+1,Nold)-                  &
     &                            Awrk(i,j,k  ,Nold))
#  ifdef MASKING
                dAdz(i,j,k2)=dAdz(i,j,k2)*rmask(i,j)
#  endif
              END DO
            END DO
          END IF
!
!  Compute components of the rotated A flux (A m3/s) along geopotential
!  surfaces.
!
          IF (k.gt.0) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.25_r8*(Kh(i-1,j)+Kh(i-1,j))*on_u(i,j)
                cff1=MIN(dZdx(i,j,k1),0.0_r8)
                cff2=MAX(dZdx(i,j,k1),0.0_r8)
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  (dAdx(i,j,k1)-                                  &
     &                   0.5_r8*(cff1*(dAdz(i-1,j,k1)+                  &
     &                                 dAdz(i  ,j,k2))+                 &
     &                           cff2*(dAdz(i-1,j,k2)+                  &
     &                                 dAdz(i  ,j,k1))))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.25_r8*(Kh(i,j-1)+Kh(i,j))*om_v(i,j)
                cff1=MIN(dZde(i,j,k1),0.0_r8)
                cff2=MAX(dZde(i,j,k1),0.0_r8)
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (dAde(i,j,k1)-                                  &
     &                   0.5_r8*(cff1*(dAdz(i,j-1,k1)+                  &
     &                                 dAdz(i,j  ,k2))+                 &
     &                           cff2*(dAdz(i,j-1,k2)+                  &
     &                                 dAdz(i,j  ,k1))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  cff=0.5_r8*Kh(i,j)
                  cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                  cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                  cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                  cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
                  FZ(i,j,k2)=cff*                                       &
     &                       (cff1*(cff1*dAdz(i,j,k2)-dAdx(i  ,j,k1))+  &
     &                        cff2*(cff2*dAdz(i,j,k2)-dAdx(i+1,j,k2))+  &
     &                        cff3*(cff3*dAdz(i,j,k2)-dAdx(i  ,j,k2))+  &
     &                        cff4*(cff4*dAdz(i,j,k2)-dAdx(i+1,j,k1)))
                  cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                  cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                  cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                  cff4=MAX(dZde(i,j+1,k1),0.0_r8)
                  FZ(i,j,k2)=FZ(i,j,k2)+                                &
     &                       cff*                                       &
     &                       (cff1*(cff1*dAdz(i,j,k2)-dAde(i,j  ,k1))+  &
     &                        cff2*(cff2*dAdz(i,j,k2)-dAde(i,j+1,k2))+  &
     &                        cff3*(cff3*dAdz(i,j,k2)-dAde(i,j  ,k2))+  &
     &                        cff4*(cff4*dAdz(i,j,k2)-dAde(i,j+1,k1)))
                END DO
              END DO
            END IF
!
!  Time-step harmonic, geopotential diffusion term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
                Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                      &
     &                           Hfac(i,j)*                             &
     &                           (FX(i+1,j  )-FX(i,j)+                  &
     &                            FE(i  ,j+1)-FE(i,j))+                 &
     &                           DTsizeH*                               &
     &                           (FZ(i,j,k2)-FZ(i,j,k1))
              END DO
            END DO
          END IF
        END DO K_LOOP

# else

!
!  Diffusion along S-coordinates: compute XI- and ETA-components of
!  diffusive flux.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*           &
     &                (Awrk(i,j,k,Nold)-Awrk(i-1,j,k,Nold))
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*           &
     &                (Awrk(i,j,k,Nold)-Awrk(i,j-1,k,Nold))
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
!
!  Time-step horizontal diffusion equation.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         Hfac(i,j)*                               &
     &                         (FX(i+1,j)-FX(i,j)+                      &
     &                          FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO
# endif
!
!  Apply boundary conditions.
!
        CALL dabc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Awrk(:,:,:,Nnew))
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Nghost,                                     &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Awrk(:,:,:,Nnew))
# endif
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
#   ifdef SPLINES_VCONV
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation implicitly using parabolic
!  splines.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=Jstr,Jend
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Kv(i,j,k)
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         DTsizeV*oHz(i,j,k)*                      &
     &                         (DC(i,k)-DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   else
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation implicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the diffusion equation into DC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,j,k)-FC(i,j,k-1)
              DC(i,k)=Awrk(i,j,k,Nold)*Hz(i,j,k)
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=(DC(i,N(ng))-                                   &
     &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
     &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
#    ifdef MASKING
            Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*rmask(i,j)
#    endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              Awrk(i,j,k,Nnew)=DC(i,k)
#    ifdef MASKING
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*rmask(i,j)
#    endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   endif
#  else
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation explicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Compute vertical diffusive flux.  Notice that "FC" is assumed to
!  be time invariant.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FS(i,k)=FC(i,j,k)*(Awrk(i,j,k+1,Nold)-                    &
     &                           Awrk(i,j,k  ,Nold))
#   ifdef MASKING
              FS(i,k)=FS(i,k)*rmask(i,j)
#   endif
            END DO
          END DO
          DO i=Istr,Iend
            FS(i,0)=0.0_r8
            FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         oHz(i,j,k)*(FS(i,k  )-                   &
     &                                     FS(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Load convolved solution.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            A(i,j,k)=Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO
      CALL dabc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    A)
# endif

      RETURN
      END SUBROUTINE conv_r3d_tile
!
!***********************************************************************
      SUBROUTINE conv_u3d_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Nghost, NHsteps, NVsteps,               &
     &                          DTsizeH, DTsizeV,                       &
     &                          Kh, Kv,                                 &
     &                          pm, pn,                                 &
# ifdef GEOPOTENTIAL_HCONV
     &                          on_r, om_p,                             &
# else
     &                          pmon_r, pnom_p,                         &
# endif
# ifdef MASKING
#  ifdef GEOPOTENTIAL_HCONV
     &                          pmask, rmask, umask, vmask,             &
#  else
     &                          umask, pmask,                           &
#  endif
# endif
     &                          Hz, z_r,                                &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY: dabc_u3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps, NVsteps

      real(r8), intent(in) :: DTsizeH, DTsizeV
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
#  else
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   else
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
#   endif
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   else
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#   endif
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, k, k1, k2, step

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk,2) :: Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
#  endif
#  if defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
#   endif
#  else
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FS
#  endif
# endif
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZde

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_p
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 3D state variable
!  at U-points.
!-----------------------------------------------------------------------
!
!  Compute metrics factors.  Notice that "z_r" and "Hz" are assumed to
!  be time invariant in the vertical convolution.  Scratch array are
!  used for efficiency.
!
      cff=DTsizeH*0.25_r8
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          Hfac(i,j)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
            FC(i,j,k)=-DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                 &
     &                 (z_r(i-1,j,k+1)+z_r(i,j,k+1)-                    &
     &                  z_r(i-1,j,k  )-z_r(i,j,k  ))
#   else
            FC(i,j,k)=DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                  &
     &                (z_r(i-1,j,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i-1,j,k  )-z_r(i,j,k  ))
#   endif
          END DO
          FC(i,j,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=2.0_r8/(Hz(i-1,j,k)+Hz(i,j,k))
          END DO
#  endif
# endif
        END DO
      END DO
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
      CALL dabc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    A)
# endif
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
            Awrk(i,j,k,Nold)=A(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps

# ifdef GEOPOTENTIAL_HCONV
!
!  Diffusion along geopotential surfaces: Compute horizontal and
!  vertical gradients.  Notice the recursive blocking sequence.  The
!  vertical placement of the gradients is:
!
!        dAdx,dAde(:,:,k1) k     rho-points
!        dAdx,dAde(:,:,k2) k+1   rho-points
!          FZ,dAdz(:,:,k1) k-1/2   W-points
!          FZ,dAdz(:,:,k2) k+1/2   W-points
!
        k2=1
        K_LOOP : DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU-1,Iend+1
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#  ifdef MASKING
                cff=cff*umask(i,j)
#  endif
                dZdx(i,j)=cff*(z_r(i  ,j,k+1)-                          &
     &                         z_r(i-1,j,k+1))
              END DO
            END DO
!
            DO j=Jstr,Jend+1
              DO i=IstrU-1,Iend
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#  ifdef MASKING
                cff=cff*vmask(i,j)
#  endif
                dZde(i,j)=cff*(z_r(i,j  ,k+1)-                          &
     &                         z_r(i,j-1,k+1))
              END DO
            END DO
!
            DO j=Jstr,Jend
              DO i=IstrU-1,Iend
#  ifdef MASKING
                dAdx(i,j,k2)=pm(i,j)*                                   &
     &                       (Awrk(i+1,j,k+1,Nold)*umask(i+1,j)-        &
     &                        Awrk(i  ,j,k+1,Nold)*umask(i  ,j))
                dAdx(i,j,k2)=dAdx(i,j,k2)*rmask(i,j)
#  else
                dAdx(i,j,k2)=pm(i,j)*(Awrk(i+1,j,k+1,Nold)-             &
     &                                Awrk(i  ,j,k+1,Nold))
#  endif
                dZdx_r(i,j,k2)=0.5_r8*(dZdx(i  ,j)+                     &
     &                                 dZdx(i+1,j))
              END DO
            END DO
!
            DO j=Jstr,Jend+1
              DO i=IstrU,Iend
                cff=0.25_r8*(pn(i-1,j  )+pn(i,j  )+                     &
     &                       pn(i-1,j-1)+pn(i,j-1))
#  ifdef MASKING
                dAde(i,j,k2)=cff*                                       &
     &                       (Awrk(i,j  ,k+1,Nold)*umask(i,j  )-        &
     &                        Awrk(i,j-1,k+1,Nold)*umask(i,j-1))
                dAde(i,j,k2)=dAde(i,j,k2)*pmask(i,j)
#  else
                dAde(i,j,k2)=cff*(Awrk(i,j  ,k+1,Nold)-                 &
     &                            Awrk(i,j-1,k+1,Nold))
#  endif
                dZde_p(i,j,k2)=0.5_r8*(dZde(i-1,j)+                     &
     &                                 dZde(i  ,j))
              END DO
            END DO
          END IF
!
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                dAdz(i,j,k2)=0.0_r8
                FZ(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
                dAdz(i,j,k2)=cff*(Awrk(i,j,k+1,Nold)-                  &
     &                            Awrk(i,j,k  ,Nold))
#  ifdef MASKING
                dAdz(i,j,k2)=dAdz(i,j,k2)*umask(i,j)
#  endif
              END DO
            END DO
          END IF
!
!  Compute components of the rotated A flux (A m3/s) along geopotential
!  surfaces.
!
          IF (k.gt.0) THEN
            DO j=Jstr,Jend
              DO i=IstrU-1,Iend
                cff=Kh(i,j)*on_r(i,j)
                cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
                cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
                FX(i,j)=cff*                                            &
     &                  Hz(i,j,k)*                                      &
     &                  (dAdx(i,j,k1)-                                  &
     &                   0.5_r8*(cff1*(dAdz(i  ,j,k1)+                  &
     &                                 dAdz(i+1,j,k2))+                 &
     &                           cff2*(dAdz(i  ,j,k2)+                  &
     &                                 dAdz(i+1,j,k1))))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=IstrU,Iend
                cff=0.0625_r8*(Kh(i-1,j-1)+Kh(i-1,j)+                   &
     &                         Kh(i  ,j-1)+Kh(i  ,j))*om_p(i,j)
                cff1=MIN(dZde_p(i,j,k1),0.0_r8)
                cff2=MAX(dZde_p(i,j,k1),0.0_r8)
                FE(i,j)=cff*                                            &
     &                  (Hz(i-1,j-1,k)+Hz(i-1,j,k)+                     &
     &                   Hz(i  ,j-1,k)+Hz(i  ,j,k))*                    &
     &                  (dAde(i,j,k1)-                                  &
     &                   0.5_r8*(cff1*(dAdz(i,j-1,k1)+                  &
     &                                 dAdz(i,j  ,k2))+                 &
     &                           cff2*(dAdz(i,j-1,k2)+                  &
     &                                 dAdz(i,j  ,k1))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  cff=0.25_r8*(Kh(i-1,j)+Kh(i,j))
                  cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                  cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                  cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                  cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                  FZ(i,j,k2)=cff*                                       &
     &                       (cff1*(cff1*dAdz(i,j,k2)-dAdx(i-1,j,k1))+  &
     &                        cff2*(cff2*dAdz(i,j,k2)-dAdx(i  ,j,k2))+  &
     &                        cff3*(cff3*dAdz(i,j,k2)-dAdx(i-1,j,k2))+  &
     &                        cff4*(cff4*dAdz(i,j,k2)-dAdx(i  ,j,k1)))
                  cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                  cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                  cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                  cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                  FZ(i,j,k2)=FZ(i,j,k2)+                                &
     &                       cff*                                       &
     &                       (cff1*(cff1*dAdz(i,j,k2)-dAde(i,j  ,k1))+  &
     &                        cff2*(cff2*dAdz(i,j,k2)-dAde(i,j+1,k2))+  &
     &                        cff3*(cff3*dAdz(i,j,k2)-dAde(i,j  ,k2))+  &
     &                        cff4*(cff4*dAdz(i,j,k2)-dAde(i,j+1,k1)))
                END DO
              END DO
            END IF
!
!  Time-step harmonic, geopotential diffusion term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                      &
     &                           Hfac(i,j)*                             &
     &                           (FX(i,j  )-FX(i-1,j)+                  &
     &                            FE(i,j+1)-FE(i  ,j))+                 &
     &                           DTsizeH*                               &
     &                           (FZ(i,j,k2)-FZ(i,j,k1))
              END DO
            END DO
          END IF
        END DO K_LOOP

# else

!
!  Diffusion along S-coordinates: compute XI- and ETA-components of
!  diffusive flux.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
              FX(i,j)=pmon_r(i,j)*Kh(i,j)*                              &
     &                (Awrk(i+1,j,k,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              FE(i,j)=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+       &
     &                                     Kh(i-1,j-1)+Kh(i,j-1))*      &
     &                (Awrk(i,j,k,Nold)-Awrk(i,j-1,k,Nold))
#  ifdef MASKING
              FE(i,j)=FE(i,j)*pmask(i,j)
#  endif
            END DO
          END DO
!
!  Time-step horizontal diffusion equation.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         Hfac(i,j)*                               &
     &                         (FX(i,j)-FX(i-1,j)+                      &
     &                          FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO
# endif
!
!  Apply boundary conditions.
!
        CALL dabc_u3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Awrk(:,:,:,Nnew))
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Nghost,                                     &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Awrk(:,:,:,Nnew))
# endif
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
#   ifdef SPLINES_VCONV
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation implicitly using parabolic
!  splines.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                             &
     &                         Hz(i  ,j,k))
            END DO
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              FC(i,k)=cff1*Hzk(i,k  )-DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=IstrU,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=IstrU,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
              DC(i,k)=DC(i,k)*Kv(i,j,k)
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         DTsizeV*oHz(i,j,k)*                      &
     &                         (DC(i,k)-DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   else
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation implicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the diffusion equation into DC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
              BC(i,k)=cff-FC(i,j,k)-FC(i,j,k-1)
              DC(i,k)=Awrk(i,j,k,Nold)*cff
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=IstrU,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=IstrU,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=IstrU,Iend
            DC(i,N(ng))=(DC(i,N(ng))-                                   &
     &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
     &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
#    ifdef MASKING
            Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*umask(i,j)
#    endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              Awrk(i,j,k,Nnew)=DC(i,k)
#    ifdef MASKING
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*umask(i,j)
#    endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   endif
#  else
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation explicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Compute vertical diffusive flux.  Notice that "FC" is assumed to
!  be time invariant.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              FS(i,k)=FC(i,j,k)*(Awrk(i,j,k+1,Nold)-                    &
     &                           Awrk(i,j,k  ,Nold))
#   ifdef MASKING
              FS(i,k)=FS(i,k)*umask(i,j)
#   endif
            END DO
          END DO
          DO i=IstrU,Iend
            FS(i,0)=0.0_r8
            FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         oHz(i,j,k)*(FS(i,k  )-                   &
     &                                     FS(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Load convolved solution.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            A(i,j,k)=Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO
      CALL dabc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    A)
# endif

      RETURN
      END SUBROUTINE conv_u3d_tile
!
!***********************************************************************
      SUBROUTINE conv_v3d_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Nghost, NHsteps, NVsteps,               &
     &                          DTsizeH, DTsizeV,                       &
     &                          Kh, Kv,                                 &
     &                          pm, pn,                                 &
# ifdef GEOPOTENTIAL_HCONV
     &                          on_p, om_r,                             &
# else
     &                          pmon_p, pnom_r,                         &
# endif
# ifdef MASKING
#  ifdef GEOPOTENTIAL_HCONV
     &                          pmask, rmask, umask, vmask,             &
#  else
     &                          vmask, pmask,                           &
#  endif
# endif
     &                          Hz, z_r,                                &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY: dabc_v3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps, NVsteps

      real(r8), intent(in) :: DTsizeH, DTsizeV
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
#  else
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   else
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
#   endif
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   else
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#   endif
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, k, k1, k2, step

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk,2) :: Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
#  endif
#  if defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
#   endif
#  else
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FS
#  endif
# endif
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZde

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dAde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_r
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 3D state variable
!  at V-points.
!-----------------------------------------------------------------------
!
!  Compute metrics factors.  Notice that "z_r" and "Hz" are assumed to
!  be time invariant in the vertical convolution.  Scratch array are
!  used for efficiency.
!
      cff=DTsizeH*0.25_r8
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
            FC(i,j,k)=-DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                 &
     &                 (z_r(i,j-1,k+1)+z_r(i,j,k+1)-                    &
     &                  z_r(i,j-1,k  )-z_r(i,j,k  ))
#   else
            FC(i,j,k)=DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                  &
     &                (z_r(i,j-1,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i,j-1,k  )-z_r(i,j,k  ))
#   endif
          END DO
          FC(i,j,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=2.0_r8/(Hz(i,j-1,k)+Hz(i,j,k))
          END DO
#  endif
# endif
        END DO
      END DO
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
      CALL dabc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    A)
# endif
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
            Awrk(i,j,k,Nold)=A(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps

# ifdef GEOPOTENTIAL_HCONV
!
!  Diffusion along geopotential surfaces: Compute horizontal and
!  vertical gradients.  Notice the recursive blocking sequence.  The
!  vertical placement of the gradients is:
!
!        dAdx,dAde(:,:,k1) k     rho-points
!        dAdx,dAde(:,:,k2) k+1   rho-points
!          FZ,dAdz(:,:,k1) k-1/2   W-points
!          FZ,dAdz(:,:,k2) k+1/2   W-points
!
        k2=1
        K_LOOP : DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=JstrV-1,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#  ifdef MASKING
                cff=cff*umask(i,j)
#  endif
                dZdx(i,j)=cff*(z_r(i  ,j,k+1)-                          &
     &                         z_r(i-1,j,k+1))
              END DO
            END DO
!
            DO j=JstrV-1,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#  ifdef MASKING
                cff=cff*vmask(i,j)
#  endif
                dZde(i,j)=cff*(z_r(i,j  ,k+1)-                          &
     &                         z_r(i,j-1,k+1))
              END DO
            END DO
!
            DO j=JstrV,Jend
              DO i=Istr,Iend+1
                cff=0.25_r8*(pm(i-1,j-1)+pm(i-1,j)+                     &
     &                       pm(i  ,j-1)+pm(i  ,j))
#  ifdef MASKING
                dAdx(i,j,k2)=cff*                                       &
     &                       (Awrk(i  ,j,k+1,Nold)*vmask(i  ,j)-        &
     &                        Awrk(i-1,j,k+1,Nold)*vmask(i-1,j))
                dAdx(i,j,k2)=dAdx(i,j,k2)*pmask(i,j)
#  else
                dAdx(i,j,k2)=cff*(Awrk(i  ,j,k+1,Nold)-                 &
     &                            Awrk(i-1,j,k+1,Nold))
#  endif
                dZdx_p(i,j,k2)=0.5_r8*(dZdx(i,j-1)+                     &
     &                                 dZdx(i,j  ))
              END DO
            END DO
!
            DO j=JstrV-1,Jend
              DO i=Istr,Iend
#  ifdef MASKING
                dAde(i,j,k2)=pn(i,j)*                                   &
     &                       (Awrk(i,j+1,k+1,Nold)*vmask(i,j+1)-        &
     &                        Awrk(i,j  ,k+1,Nold)*vmask(i,j  ))
                dAde(i,j,k2)=dAde(i,j,k2)*rmask(i,j)
#  else
                dAde(i,j,k2)=pn(i,j)*(Awrk(i,j+1,k+1,Nold)-             &
     &                                Awrk(i,j  ,k+1,Nold))
#  endif
                dZde_r(i,j,k2)=0.5_r8*(dZde(i,j  )+                     &
     &                                 dZde(i,j+1))
              END DO
            END DO
          END IF
!
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                dAdz(i,j,k2)=0.0_r8
                FZ(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
                dAdz(i,j,k2)=cff*(Awrk(i,j,k+1,Nold)-                  &
     &                            Awrk(i,j,k  ,Nold))
#  ifdef MASKING
                dAdz(i,j,k2)=dAdz(i,j,k2)*vmask(i,j)
#  endif
              END DO
            END DO
          END IF
!
!  Compute components of the rotated A flux (A m3/s) along geopotential
!  surfaces.
!
          IF (k.gt.0) THEN
            DO j=JstrV,Jend
              DO i=Istr,Iend+1
                cff=0.0625_r8*(Kh(i-1,j-1)+Kh(i-1,j)+                   &
     &                         Kh(i  ,j-1)+Kh(i  ,j))*on_p(i,j)
                cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
                cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
                FX(i,j)=cff*                                            &
     &                  (Hz(i-1,j-1,k)+Hz(i-1,j,k)+                     &
     &                   Hz(i  ,j-1,k)+Hz(i  ,j,k))*                    &
     &                  (dAdx(i,j,k1)-                                  &
     &                   0.5_r8*(cff1*(dAdz(i-1,j,k1)+                  &
     &                                 dAdz(i  ,j,k2))+                 &
     &                           cff2*(dAdz(i-1,j,k2)+                  &
     &                                 dAdz(i  ,j,k1))))
              END DO
            END DO
            DO j=JstrV-1,Jend
              DO i=Istr,Iend
                cff=Kh(i,j)*om_r(i,j)
                cff1=MIN(dZde_r(i,j,k1),0.0_r8)
                cff2=MAX(dZde_r(i,j,k1),0.0_r8)
                FE(i,j)=cff*                                            &
     &                  Hz(i,j,k)*                                      &
     &                  (dAde(i,j,k1)-                                  &
     &                   0.5_r8*(cff1*(dAdz(i,j  ,k1)+                  &
     &                                 dAdz(i,j+1,k2))+                 &
     &                           cff2*(dAdz(i,j  ,k2)+                  &
     &                                 dAdz(i,j+1,k1))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  cff=0.5_r8*(Kh(i,j-1)+Kh(i,j))
                  cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                  cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                  cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                  cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                  FZ(i,j,k2)=cff*                                       &
     &                       (cff1*(cff1*dAdz(i,j,k2)-dAdx(i  ,j,k1))+  &
     &                        cff2*(cff2*dAdz(i,j,k2)-dAdx(i+1,j,k2))+  &
     &                        cff3*(cff3*dAdz(i,j,k2)-dAdx(i  ,j,k2))+  &
     &                        cff4*(cff4*dAdz(i,j,k2)-dAdx(i+1,j,k1)))
                  cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                  cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                  cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                  cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                  FZ(i,j,k2)=FZ(i,j,k2)+                                &
     &                       cff*                                       &
     &                       (cff1*(cff1*dAdz(i,j,k2)-dAde(i,j-1,k1))+  &
     &                        cff2*(cff2*dAdz(i,j,k2)-dAde(i,j  ,k2))+  &
     &                        cff3*(cff3*dAdz(i,j,k2)-dAde(i,j-1,k2))+  &
     &                        cff4*(cff4*dAdz(i,j,k2)-dAde(i,j  ,k1)))
                END DO
              END DO
            END IF
!
!  Time-step harmonic, geopotential diffusion term (m Tunits).
!
            DO j=JstrV,Jend
              DO i=Istr,Iend
                Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                      &
     &                           Hfac(i,j)*                             &
     &                           (FX(i+1,j)-FX(i,j  )+                  &
     &                            FE(i  ,j)-FE(i,j-1))+                 &
     &                           DTsizeH*                               &
     &                           (FZ(i,j,k2)-FZ(i,j,k1))
              END DO
            END DO
          END IF
        END DO K_LOOP

# else

!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              FX(i,j)=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+       &
     &                                     Kh(i-1,j-1)+Kh(i,j-1))*      &
     &                (Awrk(i,j,k,Nold)-Awrk(i-1,j,k,Nold))
#  ifdef MASKING
              FX(i,j)=FX(i,j)*pmask(i,j)
#  endif
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
              FE(i,j)=pnom_r(i,j)*Kh(i,j)*                              &
     &                (Awrk(i,j+1,k,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
!
!  Time-step horizontal diffusion equation.
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         Hfac(i,j)*                               &
     &                         (FX(i+1,j)-FX(i,j)+                      &
     &                          FE(i,j)-FE(i,j-1))
            END DO
          END DO
        END DO
# endif
!
!  Apply boundary conditions.
!
        CALL dabc_v3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Awrk(:,:,:,Nnew))
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Nghost,                                     &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Awrk(:,:,:,Nnew))
# endif
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
#   ifdef SPLINES_VCONV
!
!-----------------------------------------------------------------------
!  Integrate vertical diffusion equation implicitly using parabolic
!  splines.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
            END DO
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Kv(i,j,k)
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
     &                         DTsizeV*oHz(i,j,k)*                      &
     &                         (DC(i,k)-DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   else
!
!-----------------------------------------------------------------------
!  Integerate vertical diffusion equation implicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the diffusion equation into DC.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
              BC(i,k)=cff-FC(i,j,k)-FC(i,j,k-1)
              DC(i,k)=Awrk(i,j,k,Nold)*cff
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=(DC(i,N(ng))-                                   &
     &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
     &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
#    ifdef MASKING
            Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*vmask(i,j)
#    endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              Awrk(i,j,k,Nnew)=DC(i,k)
#    ifdef MASKING
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*vmask(i,j)
#    endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   endif
#  else
!
!-----------------------------------------------------------------------
!  Integerate vertical diffusion equation explicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Compute vertical diffusive flux.  Notice that "FC" is assumed to
!  be time invariant.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FS(i,k)=FC(i,j,k)*(Awrk(i,j,k+1,Nold)-                    &
     &                           Awrk(i,j,k  ,Nold))
#  ifdef MASKING
              FS(i,k)=FS(i,k)*vmask(i,j)
#  endif
            END DO
          END DO
          DO i=Istr,Iend
            FS(i,0)=0.0_r8
            FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
               Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                       &
     &                          oHz(i,j,k)*(FS(i,k  )-                  &
     &                                      FS(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Load convolved solution.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
            A(i,j,k)=Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO
      CALL dabc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    A)
# endif

      RETURN
      END SUBROUTINE conv_v3d_tile
#endif
      END MODULE conv_3d_mod
