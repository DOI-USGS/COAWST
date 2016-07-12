#include "cppdefs.h"

      MODULE conv_3d_bry_mod

#if defined NONLINEAR && defined FOUR_DVAR && defined SOLVE3D && \
    defined ADJUST_BOUNDARY
!
!svn $Id: conv_bry3d.F 795 2016-05-11 01:42:43Z arango $
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
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number                                    !
!     tile       Tile partition                                        !
!     model      Calling model identifier                              !
!     boundary   Boundary edge to convolve                             !
!     edge       Boundary edges index                                  !
!     LBij       Lower bound MIN(I,J)-dimension                        !
!     LBij       Lower bound MAX(I,J)-dimension                        !
!     LBi        I-dimension lower bound                               !
!     UBi        I-dimension upper bound                               !
!     LBj        J-dimension lower bound                               !
!     UBj        J-dimension upper bound                               !
!     LBk        K-dimension lower bound                               !
!     UBk        K-dimension upper bound                               !
!     Nghost     Number of ghost points                                !
!     NHsteps    Number of horizontal diffusion integration steps      !
!     NVsteps    Number of vertical   diffusion integration steps      !
!     DTsizeH    Horizontal diffusion pseudo time-step size            !
!     DTsizeV    Vertical   diffusion pseudo time-step size            !
!     Kh         Horizontal diffusion coefficients                     !
!     Kv         Vertical   diffusion coefficients                     !
!     A          3D boundary state variable to diffuse                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Convolved 3D boundary state variable                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    conv_r3d_bry_tile  Nonlinear 3D boundary convolution at RHO-points!
!    conv_u3d_bry_tile  Nonlinear 3D boundary convolution at U-points  !
!    conv_v3d_bry_tile  Nonlinear 3D boundary convolution at V-points  !
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
      SUBROUTINE conv_r3d_bry_tile (ng, tile, model, boundary,          &
     &                              edge, LBij, UBij,                   &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Nghost, NHsteps, NVsteps,           &
     &                              DTsizeH, DTsizeV,                   &
     &                              Kh, Kv,                             &
     &                              pm, pn, pmon_u, pnom_v,             &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
     &                              Hz, z_r,                            &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_bry3d_mod, ONLY: bc_r3d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps, NVsteps

      real(r8), intent(in) :: DTsizeH, DTsizeV
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: A(LBij:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: cff, cff1

      real(r8), dimension(LBij:UBij,LBk:UBk,2) :: Awrk

      real(r8), dimension(JminS:JmaxS,LBk:UBk) :: FE
      real(r8), dimension(IminS:ImaxS,LBk:UBk) :: FX
      real(r8), dimension(LBij:UBij) :: Hfac
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng)) :: oHz
#  endif
#  if defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#   endif
#  else
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 2D state variable
!  at RHO-points.
!-----------------------------------------------------------------------
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Compute metrics factors.  Notice that "z_r" and "Hz" are assumed to
!  be time invariant in the vertical convolution.  Scratch array are
!  used for efficiency.
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          i=edge(boundary)
          DO j=Jstr-1,Jend+1
            Hfac(j)=DTsizeH*pm(i,j)*pn(i,j)
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
            FC(j,N(ng))=0.0_r8
            DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
              FC(j,k)=-DTsizeV*Kv(i,j,k)/                               &
     &                (z_r(i,j,k+1)-z_r(i,j,k))
#   else
              FC(j,k)=DTsizeV*Kv(i,j,k)/                                &
     &                (z_r(i,j,k+1)-z_r(i,j,k))
#   endif
            END DO
            FC(j,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
            DO k=1,N(ng)
              oHz(j,k)=1.0_r8/Hz(i,j,k)
            END DO
#  endif
# endif
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          j=edge(boundary)
          DO i=Istr-1,Iend+1
            Hfac(i)=DTsizeH*pm(i,j)*pn(i,j)
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
            FC(i,N(ng))=0.0_r8
            DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
              FC(i,k)=-DTsizeV*Kv(i,j,k)/                               &
     &                (z_r(i,j,k+1)-z_r(i,j,k))
#   else
              FC(i,k)=DTsizeV*Kv(i,j,k)/                                &
     &                (z_r(i,j,k+1)-z_r(i,j,k))
#   endif
            END DO
            FC(i,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
            DO k=1,N(ng)
              oHz(i,k)=1.0_r8/Hz(i,j,k)
            END DO
#  endif
# endif
          END DO
        END IF
      END IF
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2

      CALL bc_r3d_bry_tile (ng, tile, boundary,                         &
     &                      LBij, UBij, 1, N(ng),                       &
     &                      A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Nghost,                                   &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        A)
# endif
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              Awrk(j,k,Nold)=A(j,k)
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr-1,Iend+1
              Awrk(i,k,Nold)=A(i,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Compute XI- and ETA-components of diffusive flux.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend+1
                FE(j,k)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*         &
     &                  (Awrk(j  ,k,Nold)-                              &
     &                   Awrk(j-1,k,Nold))
# ifdef MASKING
                FE(j,k)=FE(j,k)*vmask(i,j)
# endif
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend+1
                FX(i,k)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*         &
     &                  (Awrk(i  ,k,Nold)-                              &
     &                   Awrk(i-1,k,Nold))
# ifdef MASKING
                FX(i,k)=FX(i,k)*umask(i,j)
# endif
              END DO
            END DO
          END IF
        END IF
!
!  Time-step horizontal diffusion equation.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         Hfac(j)*                                 &
     &                         (FE(j+1,k)-FE(j,k))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         Hfac(i)*                                 &
     &                         (FX(i+1,k)-FX(i,k))
              END DO
            END DO
          END IF
        END IF
!
!  Apply boundary conditions.
!
        CALL bc_r3d_bry_tile (ng, tile, boundary,                       &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Awrk(:,:,Nnew))
# ifdef DISTRIBUTE
        CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          Nghost,                                 &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          Awrk(:,:,Nnew))
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            cff1=1.0_r8/6.0_r8
            DO j=Jstr,Jend
              DO k=1,N(ng)-1
                FC(j,k)=cff1*Hz(i,j,k  )-                               &
     &                  DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                CF(j,k)=cff1*Hz(i,j,k+1)-                               &
     &                  DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
              END DO
              CF(j,0)=0.0_r8
              DC(j,0)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            cff1=1.0_r8/6.0_r8
            DO i=Istr,Iend
              DO k=1,N(ng)-1
                FC(i,k)=cff1*Hz(i,j,k  )-                               &
     &                  DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                CF(i,k)=cff1*Hz(i,j,k+1)-                               &
     &                  DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
              END DO
              CF(i,0)=0.0_r8
              DC(i,0)=0.0_r8
            END DO
          END IF
!
!  LU decomposition and forward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            cff1=1.0_r8/3.0_r8
            DO k=1,N(ng)-1
              DO j=Jstr,Jend
                BC(j,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                   &
     &                  DTsizeV*Kv(i,j,k)*                              &
     &                  (oHz(j,k)+oHz(j,k+1))
                cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                CF(j,k)=cff*CF(j,k)
                DC(j,k)=cff*(Awrk(j,k+1,Nold)-                          &
     &                       Awrk(j,k  ,Nold)-                          &
     &                       FC(j,k)*DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            cff1=1.0_r8/3.0_r8
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                   &
     &                  DTsizeV*Kv(i,j,k)*                              &
     &                  (oHz(i,k)+oHz(i,k+1))
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                CF(i,k)=cff*CF(i,k)
                DC(i,k)=cff*(Awrk(i,k+1,Nold)-                          &
     &                       Awrk(i,k  ,Nold)-                          &
     &                       FC(i,k)*DC(i,k-1))
              END DO
            END DO
          END IF
!
!  Backward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
              DC(j,N(ng))=0.0_r8
            END DO
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
              END DO
            END DO
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DC(j,k)=DC(j,k)*Kv(i,j,k)
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         DTsizeV*oHz(j,k)*                        &
     &                         (DC(j,k)-DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
              DC(i,N(ng))=0.0_r8
            END DO
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)*Kv(i,j,k)
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         DTsizeV*oHz(i,k)*                        &
     &                         (DC(i,k)-DC(i,k-1))
              END DO
            END DO
          END IF
        END IF
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
                BC(j,k)=Hz(i,j,k)-FC(j,k)-FC(j,k-1)
                DC(j,k)=Awrk(j,k,Nold)*Hz(i,j,k)
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
                BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
                DC(i,k)=Awrk(i,k,Nold)*Hz(i,j,k)
              END DO
            END DO
          END IF
!
!  Solve the tridiagonal system.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO j=Jstr,Jend
              cff=1.0_r8/BC(j,1)
              CF(j,1)=cff*FC(j,1)
              DC(j,1)=cff*DC(j,1)
            END DO
            DO k=2,N(ng)-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                CF(j,k)=cff*FC(j,k)
                DC(j,k)=cff*(DC(j,k)-FC(j,k-1)*DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO i=Istr,Iend
              cff=1.0_r8/BC(i,1)
              CF(i,1)=cff*FC(i,1)
              DC(i,1)=cff*DC(i,1)
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                CF(i,k)=cff*FC(i,k)
                DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
              END DO
            END DO
          END IF
!
!  Compute new solution by back substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
              DC(j,N(ng))=(DC(j,N(ng))-                                 &
     &                     FC(j,N(ng)-1)*DC(j,N(ng)-1))/                &
     &                    (BC(j,N(ng))-                                 &
     &                     FC(j,N(ng)-1)*CF(j,N(ng)-1))
              Awrk(j,N(ng),Nnew)=DC(j,N(ng))
#    ifdef MASKING
              Awrk(j,N(ng),Nnew)=Awrk(j,N(ng),Nnew)*                    &
     &                           rmask(i,j)
#    endif
            END DO
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
                Awrk(j,k,Nnew)=DC(j,k)
#    ifdef MASKING
                Awrk(j,k,Nnew)=Awrk(j,k,Nnew)*rmask(i,j)
#    endif
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
              DC(i,N(ng))=(DC(i,N(ng))-                                 &
     &                     FC(i,N(ng)-1)*DC(i,N(ng)-1))/                &
     &                    (BC(i,N(ng))-                                 &
     &                     FC(i,N(ng)-1)*CF(i,N(ng)-1))
              Awrk(i,N(ng),Nnew)=DC(i,N(ng))
#    ifdef MASKING
              Awrk(i,N(ng),Nnew)=Awrk(i,N(ng),Nnew)*rmask(i,j)
#    endif
            END DO
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
                Awrk(i,k,Nnew)=DC(i,k)
#    ifdef MASKING
                Awrk(i,k,Nnew)=Awrk(i,k,Nnew)*rmask(i,j)
#    endif
              END DO
            END DO
          END IF
        END IF
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
              DO k=1,N(ng)-1
                FS(j,k)=FC(j,k)*(Awrk(j,k+1,Nold)-                      &
     &                           Awrk(j,k  ,Nold))
#   ifdef MASKING
                FS(j,k)=FS(j,k)*rmask(i,j)
#   endif
              END DO
              FS(j,0)=0.0_r8
              FS(j,N(ng))=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
              DO k=1,N(ng)-1
                FS(i,k)=FC(i,k)*(Awrk(i,k+1,Nold)-                      &
     &                           Awrk(i,k  ,Nold))
#   ifdef MASKING
                FS(i,k)=FS(i,k)*rmask(i,j)
#   endif
              END DO
              FS(i,0)=0.0_r8
              FS(i,N(ng))=0.0_r8
            END DO
          END IF
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         oHz(j,k)*(FS(j,k  )-                     &
     &                                   FS(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         oHz(i,k)*(FS(i,k  )-                     &
     &                                   FS(i,k-1))
              END DO
            END DO
          END IF
        END IF
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
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              A(j,k)=Awrk(j,k,Nold)
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              A(i,k)=Awrk(i,k,Nold)
            END DO
          END DO
        END IF
      END IF
      CALL bc_r3d_bry_tile (ng, tile, boundary,                         &
     &                      LBij, UBij, 1, N(ng),                       &
     &                      A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Nghost,                                   &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        A)
# endif

      RETURN
      END SUBROUTINE conv_r3d_bry_tile

!
!***********************************************************************
      SUBROUTINE conv_u3d_bry_tile (ng, tile, model, boundary,          &
     &                              edge, LBij, UBij,                   &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Nghost, NHsteps, NVsteps,           &
     &                              DTsizeH, DTsizeV,                   &
     &                              Kh, Kv,                             &
     &                              pm, pn, pmon_r, pnom_p,             &
# ifdef MASKING
     &                              umask, pmask,                       &
# endif
     &                              Hz, z_r,                            &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_bry3d_mod, ONLY: bc_u3d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps, NVsteps

      real(r8), intent(in) :: DTsizeH, DTsizeV
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: A(LBij:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: cff, cff1

      real(r8), dimension(LBij:UBij,LBk:UBk,2) :: Awrk

      real(r8), dimension(JminS:JmaxS,LBk:UBk) :: FE
      real(r8), dimension(IminS:ImaxS,LBk:UBk) :: FX
      real(r8), dimension(LBij:UBij) :: Hfac
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng)) :: oHz
#  endif
#  if defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#   endif
#  else
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 2D state variable
!  at RHO-points.
!-----------------------------------------------------------------------
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Compute metrics factors.  Notice that "z_r" and "Hz" are assumed to
!  be time invariant in the vertical convolution.  Scratch array are
!  used for efficiency.
!
      IF (Lconvolve(boundary)) THEN
        cff=DTsizeH*0.25_r8
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          i=edge(boundary)
          DO j=Jstr-1,Jend+1
            Hfac(j)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
            FC(j,N(ng))=0.0_r8
            DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
              FC(j,k)=-DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                 &
     &                (z_r(i-1,j,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i-1,j,k  )-z_r(i,j,k  ))
#   else
              FC(j,k)=DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                  &
     &                (z_r(i-1,j,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i-1,j,k  )-z_r(i,j,k  ))
#   endif
            END DO
            FC(j,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
            DO k=1,N(ng)
              oHz(j,k)=2.0_r8/(Hz(i-1,j,k)+Hz(i,j,k))
            END DO
#  endif
# endif
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          j=edge(boundary)
          DO i=IstrU-1,Iend+1
            Hfac(i)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
            FC(i,N(ng))=0.0_r8
            DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
              FC(i,k)=-DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                 &
     &                (z_r(i-1,j,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i-1,j,k  )-z_r(i,j,k  ))
#   else
              FC(i,k)=DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                  &
     &                (z_r(i-1,j,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i-1,j,k  )-z_r(i,j,k  ))
#   endif
            END DO
            FC(i,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
            DO k=1,N(ng)
              oHz(i,k)=2.0_r8/(Hz(i-1,j,k)+Hz(i,j,k))
            END DO
#  endif
# endif
          END DO
        END IF
      END IF
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2

      CALL bc_u3d_bry_tile (ng, tile, boundary,                         &
     &                      LBij, UBij, 1, N(ng),                       &
     &                      A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Nghost,                                   &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        A)
# endif
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              Awrk(j,k,Nold)=A(j,k)
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend+1
              Awrk(i,k,Nold)=A(i,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Compute XI- and ETA-components of diffusive flux.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend+1
                FE(j,k)=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+     &
     &                                       Kh(i-1,j-1)+Kh(i,j-1))*    &
     &                  (Awrk(j  ,k,Nold)-                              &
     &                   Awrk(j-1,k,Nold))
# ifdef MASKING
                FE(j,k)=FE(j,k)*pmask(i,j)
# endif
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=IstrU-1,Iend
                FX(i,k)=pmon_r(i,j)*Kh(i,j)*                            &
     &                  (Awrk(i+1,k,Nold)-                              &
     &                   Awrk(i  ,k,Nold))
              END DO
            END DO
          END IF
        END IF
!
!  Time-step horizontal diffusion equation.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         Hfac(j)*                                 &
     &                         (FE(j+1,k)-FE(j,k))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         Hfac(i)*                                 &
     &                         (FX(i,k)-FX(i-1,k))
              END DO
            END DO
          END IF
        END IF
!
!  Apply boundary conditions.
!
        CALL bc_u3d_bry_tile (ng, tile, boundary,                       &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Awrk(:,:,Nnew))
# ifdef DISTRIBUTE
        CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          Nghost,                                 &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          Awrk(:,:,Nnew))
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            cff1=0.5_r8*(1.0_r8/6.0_r8)
            DO j=Jstr,Jend
              DO k=1,N(ng)-1
                FC(j,k)=cff1*(Hz(i-1,j,k  )+Hz(i,j,k  ))-               &
     &                  DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                CF(j,k)=cff1*(Hz(i-1,j,k+1)+Hz(i,j,k+1))-               &
     &                  DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
              END DO
              CF(j,0)=0.0_r8
              DC(j,0)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            cff1=0.5_r8*(1.0_r8/6.0_r8)
            DO i=IstrU,Iend
              DO k=1,N(ng)-1
                FC(i,k)=cff1*(Hz(i-1,j,k  )+Hz(i,j,k  ))-               &
     &                  DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                CF(i,k)=cff1*(Hz(i-1,j,k+1)+Hz(i,j,k+1))-               &
     &                  DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
              END DO
              CF(i,0)=0.0_r8
              DC(i,0)=0.0_r8
            END DO
          END IF
!
!  LU decomposition and forward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            cff1=0.5_r8*(1.0_r8/3.0_r8)
            DO k=1,N(ng)-1
              DO j=Jstr,Jend
                BC(j,k)=cff1*(Hz(i-1,j,k  )+Hz(i,j,k  )+                &
     &                        Hz(i-1,j,k+1)+Hz(i,j,k+1))+               &
     &                  DTsizeV*Kv(i,j,k)*                              &
     &                  (oHz(j,k)+oHz(j,k+1))
                cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                CF(j,k)=cff*CF(j,k)
                DC(j,k)=cff*(Awrk(j,k+1,Nold)-                          &
     &                       Awrk(j,k  ,Nold)-                          &
     &                       FC(j,k)*DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            cff1=0.5_r8*(1.0_r8/3.0_r8)
            DO k=1,N(ng)-1
              DO i=IstrU,Iend
                BC(i,k)=cff1*(Hz(i-1,j,k  )+Hz(i,j,k  )+                &
     &                        Hz(i-1,j,k+1)+Hz(i,j,k+1))+               &
     &                  DTsizeV*Kv(i,j,k)*                              &
     &                  (oHz(i,k)+oHz(i,k+1))
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                CF(i,k)=cff*CF(i,k)
                DC(i,k)=cff*(Awrk(i,k+1,Nold)-                          &
     &                       Awrk(i,k  ,Nold)-                          &
     &                       FC(i,k)*DC(i,k-1))
              END DO
            END DO
          END IF
!
!  Backward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
              DC(j,N(ng))=0.0_r8
            END DO
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
              END DO
            END DO
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DC(j,k)=DC(j,k)*Kv(i,j,k)
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         DTsizeV*oHz(j,k)*                        &
     &                         (DC(j,k)-DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=IstrU,Iend
              DC(i,N(ng))=0.0_r8
            END DO
            DO k=N(ng)-1,1,-1
              DO i=IstrU,Iend
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              END DO
            END DO
            DO k=1,N(ng)
              DO i=IstrU,Iend
                DC(i,k)=DC(i,k)*Kv(i,j,k)
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         DTsizeV*oHz(i,k)*                        &
     &                         (DC(i,k)-DC(i,k-1))
              END DO
            END DO
          END IF
        END IF
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                BC(j,k)=cff-FC(j,k)-FC(j,k-1)
                DC(j,k)=Awrk(j,k,Nold)*cff
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                BC(i,k)=cff-FC(i,k)-FC(i,k-1)
                DC(i,k)=Awrk(i,k,Nold)*cff
              END DO
            END DO
          END IF
!
!  Solve the tridiagonal system.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO j=Jstr,Jend
              cff=1.0_r8/BC(j,1)
              CF(j,1)=cff*FC(j,1)
              DC(j,1)=cff*DC(j,1)
            END DO
            DO k=2,N(ng)-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                CF(j,k)=cff*FC(j,k)
                DC(j,k)=cff*(DC(j,k)-FC(j,k-1)*DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO i=IstrU,Iend
              cff=1.0_r8/BC(i,1)
              CF(i,1)=cff*FC(i,1)
              DC(i,1)=cff*DC(i,1)
            END DO
            DO k=2,N(ng)-1
              DO i=IstrU,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                CF(i,k)=cff*FC(i,k)
                DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
              END DO
            END DO
          END IF
!
!  Compute new solution by back substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
              DC(j,N(ng))=(DC(j,N(ng))-                                 &
     &                     FC(j,N(ng)-1)*DC(j,N(ng)-1))/                &
     &                    (BC(j,N(ng))-                                 &
     &                     FC(j,N(ng)-1)*CF(j,N(ng)-1))
              Awrk(j,N(ng),Nnew)=DC(j,N(ng))
#    ifdef MASKING
              Awrk(j,N(ng),Nnew)=Awrk(j,N(ng),Nnew)*umask(i,j)
#    endif
            END DO
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
                Awrk(j,k,Nnew)=DC(j,k)
#    ifdef MASKING
                Awrk(j,k,Nnew)=Awrk(j,k,Nnew)*umask(i,j)
#    endif
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=IstrU,Iend
              DC(i,N(ng))=(DC(i,N(ng))-                                 &
     &                     FC(i,N(ng)-1)*DC(i,N(ng)-1))/                &
     &                    (BC(i,N(ng))-                                 &
     &                     FC(i,N(ng)-1)*CF(i,N(ng)-1))
              Awrk(i,N(ng),Nnew)=DC(i,N(ng))
#    ifdef MASKING
              Awrk(i,N(ng),Nnew)=Awrk(i,N(ng),Nnew)*umask(i,j)
#    endif
            END DO
            DO k=N(ng)-1,1,-1
              DO i=IstrU,Iend
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
                Awrk(i,k,Nnew)=DC(i,k)
#    ifdef MASKING
                Awrk(i,k,Nnew)=Awrk(i,k,Nnew)*umask(i,j)
#    endif
              END DO
            END DO
          END IF
        END IF
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
              DO k=1,N(ng)-1
                FS(j,k)=FC(j,k)*(Awrk(j,k+1,Nold)-                      &
     &                           Awrk(j,k  ,Nold))
#   ifdef MASKING
                FS(j,k)=FS(j,k)*umask(i,j)
#   endif
              END DO
              FS(j,0)=0.0_r8
              FS(j,N(ng))=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=IstrU,Iend
              DO k=1,N(ng)-1
                FS(i,k)=FC(i,k)*(Awrk(i,k+1,Nold)-                      &
     &                           Awrk(i,k  ,Nold))
#   ifdef MASKING
                FS(i,k)=FS(i,k)*umask(i,j)
#   endif
              END DO
              FS(i,0)=0.0_r8
              FS(i,N(ng))=0.0_r8
            END DO
          END IF
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         oHz(j,k)*(FS(j,k  )-                     &
     &                                   FS(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         oHz(i,k)*(FS(i,k  )-                     &
     &                                   FS(i,k-1))
              END DO
            END DO
          END IF
        END IF
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
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              A(j,k)=Awrk(j,k,Nold)
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              A(i,k)=Awrk(i,k,Nold)
            END DO
          END DO
        END IF
      END IF
      CALL bc_u3d_bry_tile (ng, tile, boundary,                         &
     &                      LBij, UBij, 1, N(ng),                       &
     &                      A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Nghost,                                   &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        A)
# endif

      RETURN
      END SUBROUTINE conv_u3d_bry_tile

!
!***********************************************************************
      SUBROUTINE conv_v3d_bry_tile (ng, tile, model, boundary,          &
     &                              edge, LBij, UBij,                   &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Nghost, NHsteps, NVsteps,           &
     &                              DTsizeH, DTsizeV,                   &
     &                              Kh, Kv,                             &
     &                              pm, pn, pmon_p, pnom_r,             &
# ifdef MASKING
     &                              vmask, pmask,                       &
# endif
     &                              Hz, z_r,                            &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_bry3d_mod, ONLY: bc_v3d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps, NVsteps

      real(r8), intent(in) :: DTsizeH, DTsizeV
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)

      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: A(LBij:,LBk:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: cff, cff1

      real(r8), dimension(LBij:UBij,LBk:UBk,2) :: Awrk

      real(r8), dimension(JminS:JmaxS,LBk:UBk) :: FE
      real(r8), dimension(IminS:ImaxS,LBk:UBk) :: FX
      real(r8), dimension(LBij:UBij) :: Hfac
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng)) :: oHz
#  endif
#  if defined IMPLICIT_VCONV || defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#   endif
#  else
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 2D state variable
!  at RHO-points.
!-----------------------------------------------------------------------
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Compute metrics factors.  Notice that "z_r" and "Hz" are assumed to
!  be time invariant in the vertical convolution.  Scratch array are
!  used for efficiency.
!
      IF (Lconvolve(boundary)) THEN
        cff=DTsizeH*0.25_r8
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          i=edge(boundary)
          DO j=JstrV-1,Jend+1
            Hfac(j)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
            FC(j,N(ng))=0.0_r8
            DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
              FC(j,k)=-DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                 &
     &                (z_r(i,j-1,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i,j-1,k  )-z_r(i,j,k  ))
#   else
              FC(j,k)=DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                  &
     &                (z_r(i,j-1,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i,j-1,k  )-z_r(i,j,k  ))
#   endif
            END DO
            FC(j,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
            DO k=1,N(ng)
              oHz(j,k)=2.0_r8/(Hz(i,j-1,k)+Hz(i,j,k))
            END DO
#  endif
# endif
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          j=edge(boundary)
          DO i=Istr-1,Iend+1
            Hfac(i)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
# ifdef VCONVOLUTION
#  ifndef SPLINES_VCONV
            FC(i,N(ng))=0.0_r8
            DO k=1,N(ng)-1
#   ifdef IMPLICIT_VCONV
              FC(i,k)=-DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                 &
     &                (z_r(i,j-1,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i,j-1,k  )-z_r(i,j,k  ))
#   else
              FC(i,k)=DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                  &
     &                (z_r(i,j-1,k+1)+z_r(i,j,k+1)-                     &
     &                 z_r(i,j-1,k  )-z_r(i,j,k  ))
#   endif
            END DO
            FC(i,0)=0.0_r8
#  endif
#  if !defined IMPLICIT_VCONV || defined SPLINES_VCONV
            DO k=1,N(ng)
              oHz(i,k)=2.0_r8/(Hz(i,j-1,k)+Hz(i,j,k))
            END DO
#  endif
# endif
          END DO
        END IF
      END IF
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2

      CALL bc_v3d_bry_tile (ng, tile, boundary,                         &
     &                      LBij, UBij, 1, N(ng),                       &
     &                      A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Nghost,                                   &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        A)
# endif
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend+1
              Awrk(j,k,Nold)=A(j,k)
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr-1,Iend+1
              Awrk(i,k,Nold)=A(i,k)
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Compute XI- and ETA-components of diffusive flux.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=JstrV-1,Jend
                FE(j,k)=pnom_r(i,j)*Kh(i,j)*                            &
     &                  (Awrk(j+1,k,Nold)-                              &
     &                   Awrk(j  ,k,Nold))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend+1
                FX(i,k)=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+     &
     &                                       Kh(i-1,j-1)+Kh(i,j-1))*    &
     &                  (Awrk(i  ,k,Nold)-                              &
     &                   Awrk(i-1,k,Nold))
# ifdef MASKING
                FX(i,k)=FX(i,k)*pmask(i,j)
# endif
              END DO
            END DO
          END IF
        END IF
!
!  Time-step horizontal diffusion equation.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         Hfac(j)*                                 &
     &                         (FE(j,k)-FE(j-1,k))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         Hfac(i)*                                 &
     &                         (FX(i+1,k)-FX(i,k))
              END DO
            END DO
          END IF
        END IF
!
!  Apply boundary conditions.
!
        CALL bc_v3d_bry_tile (ng, tile, boundary,                       &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Awrk(:,:,Nnew))
# ifdef DISTRIBUTE
        CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          Nghost,                                 &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          Awrk(:,:,Nnew))
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            cff1=0.5_r8*(1.0_r8/6.0_r8)
            DO j=JstrV,Jend
              DO k=1,N(ng)-1
                FC(j,k)=cff1*(Hz(i,j-1,k  )+Hz(i,j,k  ))-               &
     &                  DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                CF(j,k)=cff1*(Hz(i,j-1,k+1)+Hz(i,j,k+1))-               &
     &                  DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
              END DO
              CF(j,0)=0.0_r8
              DC(j,0)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            cff1=0.5_r8*(1.0_r8/6.0_r8)
            DO i=Istr,Iend
              DO k=1,N(ng)-1
                FC(i,k)=cff1*(Hz(i,j-1,k  )+Hz(i,j,k  ))-               &
     &                  DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                CF(i,k)=cff1*(Hz(i,j-1,k+1)+Hz(i,j,k+1))-               &
     &                  DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
              END DO
              CF(i,0)=0.0_r8
              DC(i,0)=0.0_r8
            END DO
          END IF
!
!  LU decomposition and forward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            cff1=0.5_r8*(1.0_r8/3.0_r8)
            DO k=1,N(ng)-1
              DO j=JstrV,Jend
                BC(j,k)=cff1*(Hz(i,j-1,k  )+Hz(i,j,k  )+                &
     &                        Hz(i,j-1,k+1)+Hz(i,j,k+1))+               &
     &                  DTsizeV*Kv(i,j,k)*                              &
     &                  (oHz(j,k)+oHz(j,k+1))
                cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                CF(j,k)=cff*CF(j,k)
                DC(j,k)=cff*(Awrk(j,k+1,Nold)-                          &
     &                       Awrk(j,k  ,Nold)-                          &
     &                       FC(j,k)*DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            cff1=0.5_r8*(1.0_r8/3.0_r8)
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                BC(i,k)=cff1*(Hz(i,j-1,k  )+Hz(i,j,k  )+                &
     &                        Hz(i,j-1,k+1)+Hz(i,j,k+1))+               &
     &                  DTsizeV*Kv(i,j,k)*                              &
     &                  (oHz(i,k)+oHz(i,k+1))
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                CF(i,k)=cff*CF(i,k)
                DC(i,k)=cff*(Awrk(i,k+1,Nold)-                          &
     &                       Awrk(i,k  ,Nold)-                          &
     &                       FC(i,k)*DC(i,k-1))
              END DO
            END DO
          END IF
!
!  Backward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=JstrV,Jend
              DC(j,N(ng))=0.0_r8
            END DO
            DO k=N(ng)-1,1,-1
              DO j=JstrV,Jend
                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
              END DO
            END DO
            DO k=1,N(ng)
              DO j=JstrV,Jend
                DC(j,k)=DC(j,k)*Kv(i,j,k)
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         DTsizeV*oHz(j,k)*                        &
     &                         (DC(j,k)-DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
              DC(i,N(ng))=0.0_r8
            END DO
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)*Kv(i,j,k)
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         DTsizeV*oHz(i,k)*                        &
     &                         (DC(i,k)-DC(i,k-1))
              END DO
            END DO
          END IF
        END IF
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                BC(j,k)=cff-FC(j,k)-FC(j,k-1)
                DC(j,k)=Awrk(j,k,Nold)*cff
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                BC(i,k)=cff-FC(i,k)-FC(i,k-1)
                DC(i,k)=Awrk(i,k,Nold)*cff
              END DO
            END DO
          END IF
!
!  Solve the tridiagonal system.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO j=JstrV,Jend
              cff=1.0_r8/BC(j,1)
              CF(j,1)=cff*FC(j,1)
              DC(j,1)=cff*DC(j,1)
            END DO
            DO k=2,N(ng)-1
              DO j=JstrV,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                CF(j,k)=cff*FC(j,k)
                DC(j,k)=cff*(DC(j,k)-FC(j,k-1)*DC(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO i=Istr,Iend
              cff=1.0_r8/BC(i,1)
              CF(i,1)=cff*FC(i,1)
              DC(i,1)=cff*DC(i,1)
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                CF(i,k)=cff*FC(i,k)
                DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
              END DO
            END DO
          END IF
!
!  Compute new solution by back substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=JstrV,Jend
              DC(j,N(ng))=(DC(j,N(ng))-                                 &
     &                     FC(j,N(ng)-1)*DC(j,N(ng)-1))/                &
     &                    (BC(j,N(ng))-                                 &
     &                     FC(j,N(ng)-1)*CF(j,N(ng)-1))
              Awrk(j,N(ng),Nnew)=DC(j,N(ng))
#    ifdef MASKING
              Awrk(j,N(ng),Nnew)=Awrk(j,N(ng),Nnew)*vmask(i,j)
#    endif
            END DO
            DO k=N(ng)-1,1,-1
              DO j=JstrV,Jend
                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
                Awrk(j,k,Nnew)=DC(j,k)
#    ifdef MASKING
                Awrk(j,k,Nnew)=Awrk(j,k,Nnew)*vmask(i,j)
#    endif
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
              DC(i,N(ng))=(DC(i,N(ng))-                                 &
     &                     FC(i,N(ng)-1)*DC(i,N(ng)-1))/                &
     &                    (BC(i,N(ng))-                                 &
     &                     FC(i,N(ng)-1)*CF(i,N(ng)-1))
              Awrk(i,N(ng),Nnew)=DC(i,N(ng))
#    ifdef MASKING
              Awrk(i,N(ng),Nnew)=Awrk(i,N(ng),Nnew)*vmask(i,j)
#    endif
            END DO
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
                Awrk(i,k,Nnew)=DC(i,k)
#    ifdef MASKING
                Awrk(i,k,Nnew)=Awrk(i,k,Nnew)*vmask(i,j)
#    endif
              END DO
            END DO
          END IF
        END IF
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
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=JstrV,Jend
              DO k=1,N(ng)-1
                FS(j,k)=FC(j,k)*(Awrk(j,k+1,Nold)-                      &
     &                           Awrk(j,k  ,Nold))
#   ifdef MASKING
                FS(j,k)=FS(j,k)*vmask(i,j)
#   endif
              END DO
              FS(j,0)=0.0_r8
              FS(j,N(ng))=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
              DO k=1,N(ng)-1
                FS(i,k)=FC(i,k)*(Awrk(i,k+1,Nold)-                      &
     &                           Awrk(i,k  ,Nold))
#   ifdef MASKING
                FS(i,k)=FS(i,k)*vmask(i,j)
#   endif
              END DO
              FS(i,0)=0.0_r8
              FS(i,N(ng))=0.0_r8
            END DO
          END IF
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                          &
     &                         oHz(j,k)*(FS(j,k  )-                     &
     &                                   FS(j,k-1))
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                          &
     &                         oHz(i,k)*(FS(i,k  )-                     &
     &                                   FS(i,k-1))
              END DO
            END DO
          END IF
        END IF
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
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              A(j,k)=Awrk(j,k,Nold)
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              A(i,k)=Awrk(i,k,Nold)
            END DO
          END DO
        END IF
      END IF
      CALL bc_v3d_bry_tile (ng, tile, boundary,                         &
     &                      LBij, UBij, 1, N(ng),                       &
     &                      A)
# ifdef DISTRIBUTE
      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        Nghost,                                   &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        A)
# endif

      RETURN
      END SUBROUTINE conv_v3d_bry_tile

#endif
      END MODULE conv_3d_bry_mod
