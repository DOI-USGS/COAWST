#include "cppdefs.h"
      MODULE ad_conv_bry2d_mod

#if defined ADJOINT && defined FOUR_DVAR && defined ADJUST_BOUNDARY
!
!svn $Id: ad_conv_bry2d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines applies the background error covariance to data      !
!  assimilation fields via the space convolution of the diffusion      !
!  equation (filter) for 3D state variables. The diffusion filter      !
!  is solved using an explicit (inefficient) algorithm.                !
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
!     LBi        I-dimension Lower bound                               !
!     UBi        I-dimension Upper bound                               !
!     LBj        J-dimension Lower bound                               !
!     UBj        J-dimension Upper bound                               !
!     Nghost     Number of ghost points                                !
!     NHsteps    Number of horizontal diffusion integration steps      !
!     DTsizeH    Horizontal diffusion pseudo time-step size            !
!     Kh         Horizontal diffusion coefficients                     !
!     ad_A       2D boundary state variable to diffuse                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A       Convolved 2D boundary state variable                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_conv_r2d_bry_tile  Tangent linear 2D boundary convolution at   !
!                            RHO-points                                !
!    ad_conv_u2d_bry_tile  Tangent linear 2D boundary convolution at   !
!                            U-points                                  !
!    ad_conv_v2d_bry_tile  Tangent linear 2D boundary convolution at   !
!                            V-points                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_conv_r2d_bry_tile (ng, tile, model, boundary,       &
     &                                 edge, LBij, UBij,                &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Nghost, NHsteps, DTsizeH,        &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_u, pnom_v,          &
# ifdef MASKING
     &                                 rmask, umask, vmask,             &
# endif
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_bry2d_mod, ONLY: ad_bc_r2d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps

      real(r8), intent(in) :: DTsizeH
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
      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(inout) :: ad_A(LBij:)
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
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_A(LBij:UBij)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, step

      real(r8) :: adfac

      real(r8), dimension(LBij:UBij,2) :: ad_Awrk

      real(r8), dimension(JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS) :: ad_FX
      real(r8), dimension(LBij:UBij) :: Hfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBij:UBij,1:2)=0.0_r8

      ad_FE(JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS)=0.0_r8
!
!-----------------------------------------------------------------------
!  Adjoint space convolution of the diffusion equation for a 1D state
!  variable at RHO-points.
!-----------------------------------------------------------------------
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set integration indices and initial conditions.
!
        Nold=1
        Nnew=2
!
!  Compute metrics factor.
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          i=edge(boundary)
          DO j=Jstr,Jend
            Hfac(j)=DTsizeH*pm(i,j)*pn(i,j)
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          j=edge(boundary)
          DO i=Istr,Iend
            Hfac(i)=DTsizeH*pm(i,j)*pn(i,j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of load convolved solution.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij,                               &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij,                            &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_r2d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij,                                 &
!>   &                      tl_A)
!>
      CALL ad_bc_r2d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij,                              &
     &                         ad_A)
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO j=Jstr,Jend
!>          tl_A(j)=tl_Awrk(j,Nold)
!>
            ad_Awrk(j,Nold)=ad_Awrk(j,Nold)+ad_A(j)
            ad_A(j)=0.0_r8
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO i=Istr,Iend
!>          tl_A(i)=tl_Awrk(i,Nold)
!>
            ad_Awrk(i,Nold)=ad_Awrk(i,Nold)+ad_A(i)
            ad_A(i)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Integrate adjoint horizontal diffusion terms.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Apply adjoint boundary conditions. If applicable, exchange boundary
!  data.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,           &
!>   &                          LBij, UBij,                             &
!>   &                          Nghost,                                 &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Awrk(:,Nnew))
!>
        CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,        &
     &                             LBij, UBij,                          &
     &                             Nghost,                              &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Awrk(:,Nnew))
# endif
!>      CALL bc_r2d_bry_tile (ng, tile, boundary,                       &
!>   &                        LBij, UBij,                               &
!>   &                        tl_Awrk(:,Nnew))
!>
        CALL ad_bc_r2d_bry_tile (ng, tile, boundary,                    &
     &                           LBij, UBij,                            &
     &                           ad_Awrk(:,Nnew))
!
!  Time-step adjoint horizontal diffusion terms.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO j=Jstr,Jend
!>            tl_Awrk(j,Nnew)=tl_Awrk(j,Nold)+                          &
!>   &                        Hfac(j)*                                  &
!>   &                        (tl_FE(j+1)-tl_FE(j))
!>
              adfac=Hfac(j)*ad_Awrk(j,Nnew)
              ad_FE(j  )=ad_FE(j  )-adfac
              ad_FE(j+1)=ad_FE(j+1)+adfac
              ad_Awrk(j,Nold)=ad_Awrk(j,Nold)+                          &
     &                        ad_Awrk(j,Nnew)
              ad_Awrk(j,Nnew)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO i=Istr,Iend
!>            tl_Awrk(i,Nnew)=tl_Awrk(i,Nold)+                          &
!>   &                        Hfac(i)*                                  &
!>   &                        (tl_FX(i+1)-tl_FX(i))
!>
              adfac=Hfac(i)*ad_Awrk(i,Nnew)
              ad_FX(i  )=ad_FX(i  )-adfac
              ad_FX(i+1)=ad_FX(i+1)+adfac
              ad_Awrk(i,Nold)=ad_Awrk(i,Nold)+                          &
     &                        ad_Awrk(i,Nnew)
              ad_Awrk(i,Nnew)=0.0_r8
            END DO
          END IF
        END IF
!
!  Compute XI- and ETA-components of adjoint diffusive flux.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend+1
# ifdef MASKING
!>            tl_FE(j)=tl_FE(j)*vmask(i,j)
!>
              ad_FE(j)=ad_FE(j)*vmask(i,j)
# endif
!>            tl_FE(j)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*          &
!>   &                 (tl_Awrk(j  ,Nold)-                              &
!>   &                  tl_Awrk(j-1,Nold))
!>
              adfac=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*ad_FE(j)
              ad_Awrk(j-1,Nold)=ad_Awrk(j-1,Nold)-adfac
              ad_Awrk(j  ,Nold)=ad_Awrk(j  ,Nold)+adfac
              ad_FE(j)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend+1
# ifdef MASKING
!>            tl_FX(i)=tl_FX(i)*umask(i,j)
!>
              ad_FX(i)=ad_FX(i)*umask(i,j)
# endif
!>            tl_FX(i)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*          &
!>   &                 (tl_Awrk(i  ,Nold)-                              &
!>   &                  tl_Awrk(i-1,Nold))
!>
              adfac=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*ad_FX(i)
              ad_Awrk(i-1,Nold)=ad_Awrk(i-1,Nold)-adfac
              ad_Awrk(i  ,Nold)=ad_Awrk(i  ,Nold)+adfac
              ad_FX(i)=0.0_r8
            END DO
          END IF
        END IF
      END DO
!
!  Set adjoint initial conditions.
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO j=Jstr-1,Jend+1
!>          tl_Awrk(j,Nold)=tl_A(j)
!>
            ad_A(j)=ad_A(j)+ad_Awrk(j,Nold)
            ad_Awrk(j,Nold)=0.0_r8
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO i=Istr-1,Iend+1
!>          tl_Awrk(i,Nold)=tl_A(i)
!>
            ad_A(i)=ad_A(i)+ad_Awrk(i,Nold)
            ad_Awrk(i,Nold)=0.0_r8
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij,                               &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij,                            &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_r2d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij,                                 &
!>   &                      tl_A)
!>
      CALL ad_bc_r2d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij,                              &
     &                         ad_A)

      RETURN
      END SUBROUTINE ad_conv_r2d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_conv_u2d_bry_tile (ng, tile, model, boundary,       &
     &                                 edge, LBij, UBij,                &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Nghost, NHsteps, DTsizeH,        &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_r, pnom_p,          &
# ifdef MASKING
     &                                 umask, pmask,                    &
# endif
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_bry2d_mod, ONLY: ad_bc_u2d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps

      real(r8), intent(in) :: DTsizeH
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
      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(inout) :: ad_A(LBij:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_A(LBij:UBij)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, step

      real(r8) :: adfac, cff

      real(r8), dimension(LBij:UBij,2) :: ad_Awrk

      real(r8), dimension(JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS) :: ad_FX
      real(r8), dimension(LBij:UBij) :: Hfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBij:UBij,1:2)=0.0_r8

      ad_FE(JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS)=0.0_r8
!
!-----------------------------------------------------------------------
!  Adjoint space convolution of the diffusion equation for a 1D state
!  variable at U-points.
!-----------------------------------------------------------------------
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
!
!  Compute metrics factor.
!
      cff=DTsizeH*0.25_r8
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          i=edge(boundary)
          DO j=Jstr,Jend
            Hfac(j)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          j=edge(boundary)
          DO i=IstrU,Iend
            Hfac(i)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of load convolved solution.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij,                               &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij,                            &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_u2d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij,                                 &
!>   &                      tl_A)
!>
      CALL ad_bc_u2d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij,                              &
     &                         ad_A)
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO j=Jstr,Jend
!>          tl_A(j)=tl_Awrk(j,Nold)
!>
            ad_Awrk(j,Nold)=ad_Awrk(j,Nold)+ad_A(j)
            ad_A(j)=0.0_r8
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO i=IstrU,Iend
!>          tl_A(i)=tl_Awrk(i,Nold)
!>
            ad_Awrk(i,Nold)=ad_Awrk(i,Nold)+ad_A(i)
            ad_A(i)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Integrate adjoint horizontal diffusion terms.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Apply adjoint boundary conditions. If applicable, exchange boundary
!  data.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,           &
!>   &                          LBij, UBij,                             &
!>   &                          Nghost,                                 &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Awrk(:,Nnew))
!>
        CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,        &
     &                             LBij, UBij,                          &
     &                             Nghost,                              &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Awrk(:,Nnew))
# endif
!>      CALL bc_u2d_bry_tile (ng, tile, boundary,                       &
!>   &                        LBij, UBij,                               &
!>   &                        tl_Awrk(:,Nnew))
!>
        CALL ad_bc_u2d_bry_tile (ng, tile, boundary,                    &
     &                           LBij, UBij,                            &
     &                           ad_Awrk(:,Nnew))
!
!  Time-step adjoint horizontal diffusion terms.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO j=Jstr,Jend
!>            tl_Awrk(j,Nnew)=tl_Awrk(j,Nold)+                          &
!>   &                        Hfac(j)*                                  &
!>   &                        (tl_FE(j+1)-tl_FE(j))
!>
              adfac=Hfac(j)*ad_Awrk(j,Nnew)
              ad_FE(j  )=ad_FE(j  )-adfac
              ad_FE(j+1)=ad_FE(j+1)+adfac
              ad_Awrk(j,Nold)=ad_Awrk(j,Nold)+                          &
     &                        ad_Awrk(j,Nnew)
              ad_Awrk(j,Nnew)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO i=IstrU,Iend
!>            tl_Awrk(i,Nnew)=tl_Awrk(i,Nold)+                          &
!>   &                        Hfac(i)*                                  &
!>   &                        (tl_FX(i)-tl_FX(i-1))
!>
              adfac=Hfac(i)*ad_Awrk(i,Nnew)
              ad_FX(i-1)=ad_FX(i-1)-adfac
              ad_FX(i  )=ad_FX(i  )+adfac
              ad_Awrk(i,Nold)=ad_Awrk(i,Nold)+                          &
     &                        ad_Awrk(i,Nnew)
              ad_Awrk(i,Nnew)=0.0_r8
            END DO
          END IF
        END IF
!
!  Compute XI- and ETA-components of diffusive flux.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend+1
# ifdef MASKING
!>            tl_FE(j)=tl_FE(j)*pmask(i,j)
!>
              ad_FE(j)=ad_FE(j)*pmask(i,j)
# endif
!>            tl_FE(j)=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+      &
!>   &                                      Kh(i-1,j-1)+Kh(i,j-1))*     &
!>   &                 (tl_Awrk(j  ,Nold)-                              &
!>   &                  tl_Awrk(j-1,Nold))
!>
              adfac=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+         &
     &                                   Kh(i-1,j-1)+Kh(i,j-1))*        &
     &              ad_FE(j)
              ad_Awrk(j-1,Nold)=ad_Awrk(j-1,Nold)-adfac
              ad_Awrk(j  ,Nold)=ad_Awrk(j  ,Nold)+adfac
              ad_FE(j)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=IstrU-1,Iend
!>            tl_FX(i)=pmon_r(i,j)*Kh(i,j)*                             &
!>   &                 (tl_Awrk(i+1,Nold)-                              &
!>   &                  tl_Awrk(i  ,Nold))
!>
              adfac=pmon_r(i,j)*Kh(i,j)*ad_FX(i)
              ad_Awrk(i  ,Nold)=ad_Awrk(i  ,Nold)-adfac
              ad_Awrk(i+1,Nold)=ad_Awrk(i+1,Nold)+adfac
              ad_FX(i)=0.0_r8
            END DO
          END IF
        END IF
      END DO
!
!  Set adjoint initial conditions.
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO j=Jstr-1,Jend+1
!>          tl_Awrk(j,Nold)=tl_A(j)
!>
            ad_A(j)=ad_A(j)+ad_Awrk(j,Nold)
            ad_Awrk(j,Nold)=0.0_r8
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO i=IstrU-1,Iend+1
!>          tl_Awrk(i,Nold)=tl_A(i)
!>
            ad_A(i)=ad_A(i)+ad_Awrk(i,Nold)
            ad_Awrk(i,Nold)=0.0_r8
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij,                               &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij,                            &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_u2d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij,                                 &
!>   &                      tl_A)
!>
      CALL ad_bc_u2d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij,                              &
     &                         ad_A)

      RETURN
      END SUBROUTINE ad_conv_u2d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_conv_v2d_bry_tile (ng, tile, model, boundary,       &
     &                                 edge, LBij, UBij,                &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Nghost, NHsteps, DTsizeH,        &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_p, pnom_r,          &
# ifdef MASKING
     &                                 vmask, pmask,                    &
# endif
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_bry2d_mod, ONLY: ad_bc_v2d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost, NHsteps

      real(r8), intent(in) :: DTsizeH
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
      real(r8), intent(in) :: Kh(LBi:,LBj:)
      real(r8), intent(inout) :: ad_A(LBij:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in)  :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_A(LBij:UBij)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, step

      real(r8) :: adfac, cff

      real(r8), dimension(LBij:UBij,2) :: ad_Awrk

      real(r8), dimension(JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS) :: ad_FX
      real(r8), dimension(LBij:UBij) :: Hfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBij:UBij,1:2)=0.0_r8

      ad_FE(JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS)=0.0_r8

!-----------------------------------------------------------------------
!  Adjoint space convolution of the diffusion equation for a 2D state
!  variable at RHO-points.
!-----------------------------------------------------------------------
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
!
!  Compute metrics factor.
!
      cff=DTsizeH*0.25_r8
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          i=edge(boundary)
          DO j=JstrV,Jend
            Hfac(j)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          j=edge(boundary)
          DO i=Istr,Iend
            Hfac(i)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint of load convolved solution.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij,                               &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij,                            &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_v2d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij,                                 &
!>   &                      tl_A)
!>
      CALL ad_bc_v2d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij,                              &
     &                         ad_A)
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO j=JstrV,Jend
!>          tl_A(j)=tl_Awrk(j,Nold)
!>
            ad_Awrk(j,Nold)=ad_Awrk(j,Nold)+ad_A(j)
            ad_A(j)=0.0_r8
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO i=Istr,Iend
!>          tl_A(i)=tl_Awrk(i,Nold)
!>
            ad_Awrk(i,Nold)=ad_Awrk(i,Nold)+ad_A(i)
            ad_A(i)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Integrate adjoint horizontal diffusion terms.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Apply boundary conditions. If applicable, exchange boundary data.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,           &
!>   &                          LBij, UBij,                             &
!>   &                          Nghost,                                 &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Awrk(:,Nnew))
!>
        CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,        &
     &                             LBij, UBij,                          &
     &                             Nghost,                              &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Awrk(:,Nnew))
# endif
!>      CALL bc_v2d_bry_tile (ng, tile, boundary,                       &
!>   &                        LBij, UBij,                               &
!>   &                        tl_Awrk(:,Nnew))
!>
        CALL ad_bc_v2d_bry_tile (ng, tile, boundary,                    &
     &                           LBij, UBij,                            &
     &                           ad_Awrk(:,Nnew))
!
!  Time-step adjoint horizontal diffusion terms.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO j=JstrV,Jend
!>            tl_Awrk(j,Nnew)=tl_Awrk(j,Nold)+                          &
!>   &                        Hfac(j)*                                  &
!>   &                        (tl_FE(j)-tl_FE(j-1))
!>
              adfac=Hfac(j)*ad_Awrk(j,Nnew)
              ad_FE(j-1)=ad_FE(j-1)-adfac
              ad_FE(j  )=ad_FE(j  )+adfac
              ad_Awrk(j,Nold)=ad_Awrk(j,Nold)+                          &
     &                        ad_Awrk(j,Nnew)
              ad_Awrk(j,Nnew)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO i=Istr,Iend
!>            tl_Awrk(i,Nnew)=tl_Awrk(i,Nold)+                          &
!>   &                        Hfac(i)*                                  &
!>   &                        (tl_FX(i+1)-tl_FX(i))
!>
              adfac=Hfac(i)*ad_Awrk(i,Nnew)
              ad_FX(i  )=ad_FX(i  )-adfac
              ad_FX(i+1)=ad_FX(i+1)+adfac
              ad_Awrk(i,Nold)=ad_Awrk(i,Nold)+                          &
     &                        ad_Awrk(i,Nnew)
              ad_Awrk(i,Nnew)=0.0_r8
            END DO
          END IF
        END IF
!
!  Compute XI- and ETA-components of adjoint diffusive flux.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=JstrV-1,Jend
!>            tl_FE(j)=pnom_r(i,j)*Kh(i,j)*                             &
!>   &                 (tl_Awrk(j+1,Nold)-                              &
!>   &                  tl_Awrk(j  ,Nold))
!>
              adfac=pnom_r(i,j)*Kh(i,j)*ad_FE(j)
              ad_Awrk(j  ,Nold)=ad_Awrk(j  ,Nold)-adfac
              ad_Awrk(j+1,Nold)=ad_Awrk(j+1,Nold)+adfac
              ad_FE(j)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend+1
# ifdef MASKING
!>            tl_FX(i)=tl_FX(i)*pmask(i,j)
!>
              ad_FX(i)=ad_FX(i)*pmask(i,j)
# endif
!>            tl_FX(i)=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+      &
!>   &                                      Kh(i-1,j-1)+Kh(i,j-1))*     &
!>   &                 (tl_Awrk(i  ,Nold)-                              &
!>   &                  tl_Awrk(i-1,Nold))
!>
              adfac=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+         &
     &                                   Kh(i-1,j-1)+Kh(i,j-1))*        &
     &              ad_FX(i)
              ad_Awrk(i-1,Nold)=ad_Awrk(i-1,Nold)-adfac
              ad_Awrk(i  ,Nold)=ad_Awrk(i  ,Nold)+adfac
              ad_FX(i)=0.0_r8
            END DO
          END IF
        END IF
      END DO
!
!  Set adjoint initial conditions.
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO j=JstrV-1,Jend+1
!>          tl_Awrk(j,Nold)=tl_A(j)
!>
            ad_A(j)=ad_A(j)+ad_Awrk(j,Nold)
            ad_Awrk(j,Nold)=0.0_r8
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO i=Istr-1,Iend+1
!>          tl_Awrk(i,Nold)=tl_A(i)
!>
            ad_A(i)=ad_A(i)+ad_Awrk(i,Nold)
            ad_Awrk(i,Nold)=0.0_r8
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij,                               &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij,                            &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_v2d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij,                                 &
!>   &                      tl_A)
!>
      CALL ad_bc_v2d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij,                              &
     &                         ad_A)

      RETURN
      END SUBROUTINE ad_conv_v2d_bry_tile
#endif
      END MODULE ad_conv_bry2d_mod
