#include "cppdefs.h"

      MODULE ad_conv_2d_mod

#if defined ADJOINT && defined FOUR_DVAR
!
!svn $Id: ad_conv_2d.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines applies the background error covariance to data      !
!  assimilation fields via the  adjoint space convolution  of the      !
!  diffusion equation (filter) for 3D state variables. The filter      !
!  is solved using an explicit (inefficient) algorithm.                !
!                                                                      !
!  For Gaussian (bell-shaped) correlations, the space convolution      !
!  of the diffusion operator is an efficient way  to estimate the      !
!  finite domain error covariances.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     Nghost     Number of ghost points.                               !
!     NHsteps    Number of horizontal diffusion integration steps.     !
!     DTsizeH    Horizontal diffusion pseudo time-step size.           !
!     Kh         Horizontal diffusion coefficients.                    !
!     ad_A       2D adjoint state variable to diffuse.                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A       Diffused 2D adjoint state variable.                   !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_conv_r2d_tile  Adjoint 2D convolution at RHO-points            !
!    ad_conv_u2d_tile  Adjoint 2D convolution at U-points              !
!    ad_conv_v2d_tile  Adjoint 2D convolution at V-points              !
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
      SUBROUTINE ad_conv_r2d_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Nghost, NHsteps, DTsizeH,            &
     &                             Kh,                                  &
     &                             pm, pn, pmon_u, pnom_v,              &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
     &                             ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_2d_mod, ONLY: ad_dabc_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
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
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, step

      real(r8) :: adfac, cff

      real(r8), dimension(LBi:UBi,LBj:UBj,2) :: ad_Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBi:UBi,LBj:UBj,1:2)=0.0_r8

      ad_FE(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS,JminS:JmaxS)=0.0_r8
!
!-----------------------------------------------------------------------
!  Adjoint space convolution of the diffusion equation for a 2D state
!  variable at RHO-points.
!-----------------------------------------------------------------------
!
!  Compute metrics factor.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Hfac(i,j)=DTsizeH*pm(i,j)*pn(i,j)
        END DO
      END DO
      Nold=1
      Nnew=2
!
!------------------------------------------------------------------------
!  Adjoint of load convolved solution.
!------------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_A)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Nghost,                                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
# endif
!>    CALL dabc_r2d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    tl_A)
!>
      CALL ad_dabc_r2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)
      DO j=Jstr,Jend
        DO i=Istr,Iend
!>        tl_A(i,j)=tl_Awrk(i,j,Nold)
!>
          ad_Awrk(i,j,Nold)=ad_Awrk(i,j,Nold)+ad_A(i,j)
          ad_A(i,j)=0.0_r8
        END DO
      END DO
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
!>      CALL mp_exchange2d (ng, tile, model, 1,                         &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      Nghost,                                     &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_Awrk(:,:,Nnew))
!>
        CALL ad_mp_exchange2d (ng, tile, model, 1,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Nghost,                                  &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Awrk(:,:,Nnew))
# endif
!>      CALL dabc_r2d_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      tl_Awrk(:,:,Nnew))
!>
        CALL ad_dabc_r2d_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         ad_Awrk(:,:,Nnew))
!
!  Time-step adjoint horizontal diffusion terms.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_Awrk(i,j,Nnew)=tl_Awrk(i,j,Nold)+                        &
!>   &                        Hfac(i,j)*                                &
!>   &                        (tl_FX(i+1,j)-tl_FX(i,j)+                 &
!>   &                         tl_FE(i,j+1)-tl_FE(i,j))
!>
            adfac=Hfac(i,j)*ad_Awrk(i,j,Nnew)
            ad_FE(i,j  )=ad_FE(i,j  )-adfac
            ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
            ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
            ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
            ad_Awrk(i,j,Nold)=ad_Awrk(i,j,Nold)+ad_Awrk(i,j,Nnew)
            ad_Awrk(i,j,Nnew)=0.0_r8
          END DO
        END DO
!
!  Compute XI- and ETA-components of the adjoint diffusive flux.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
# ifdef MASKING
!>          tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
!>
            ad_FE(i,j)=ad_FE(i,j)*vmask(i,j)
# endif
!>          tl_FE(i,j)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*          &
!>   &                 (tl_Awrk(i,j,Nold)-tl_Awrk(i,j-1,Nold))
!>
            adfac=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*ad_FE(i,j)
            ad_Awrk(i,j-1,Nold)=ad_Awrk(i,j-1,Nold)-adfac
            ad_Awrk(i,j  ,Nold)=ad_Awrk(i,j  ,Nold)+adfac
            ad_FE(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
# ifdef MASKING
!>          tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
!>
            ad_FX(i,j)=ad_FX(i,j)*umask(i,j)
# endif
!>          tl_FX(i,j)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*          &
!>   &                 (tl_Awrk(i,j,Nold)-tl_Awrk(i-1,j,Nold))
!>
            adfac=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*ad_FX(i,j)
            ad_Awrk(i-1,j,Nold)=ad_Awrk(i-1,j,Nold)-adfac
            ad_Awrk(i  ,j,Nold)=ad_Awrk(i  ,j,Nold)+adfac
            ad_FX(i,j)=0.0_r8
          END DO
        END DO

      END DO
!
!  Set adjoint initial conditions.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!>        tl_Awrk(i,j,Nold)=tl_A(i,j)
!>
          ad_A(i,j)=ad_A(i,j)+ad_Awrk(i,j,Nold)
          ad_Awrk(i,j,Nold)=0.0_r8
        END DO
      END DO
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    Nghost,                                       &
!>                        EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_A)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Nghost,                                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
# endif
!>    CALL dabc_r2d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    tl_A)
!>
      CALL ad_dabc_r2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)

      RETURN
      END SUBROUTINE ad_conv_r2d_tile
!
!***********************************************************************
      SUBROUTINE ad_conv_u2d_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Nghost, NHsteps, DTsizeH,            &
     &                             Kh,                                  &
     &                             pm, pn, pmon_r, pnom_p,              &
# ifdef MASKING
     &                             umask, pmask,                        &
# endif
     &                             ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_2d_mod, ONLY: ad_dabc_u2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
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
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, step

      real(r8) :: adfac, cff

      real(r8), dimension(LBi:UBi,LBj:UBj,2) :: ad_Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBi:UBi,LBj:UBj,1:2)=0.0_r8

      ad_FE(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS,JminS:JmaxS)=0.0_r8
!
!-----------------------------------------------------------------------
!  Adjoint space convolution of the diffusion equation for a 2D state
!  variable at U-points.
!-----------------------------------------------------------------------
!
!  Compute metrics factor.
!
      cff=DTsizeH*0.25_r8
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          Hfac(i,j)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
        END DO
      END DO
      Nold=1
      Nnew=2
!
!------------------------------------------------------------------------
!  Adjoint of load convolved solution.
!------------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_A)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Nghost,                                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
# endif
!>    CALL dabc_u2d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    tl_A)
!>
      CALL ad_dabc_u2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        tl_A(i,j)=tl_Awrk(i,j,Nold)
!>
          ad_Awrk(i,j,Nold)=ad_Awrk(i,j,Nold)+ad_A(i,j)
          ad_A(i,j)=0.0_r8
        END DO
      END DO
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
!>      CALL mp_exchange2d (ng, tile, model, 1,                         &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      Nghost,                                     &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_Awrk(:,:,Nnew))
!>
        CALL ad_mp_exchange2d (ng, tile, model, 1,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Nghost,                                  &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Awrk(:,:,Nnew))
# endif
!>      CALL dabc_u2d_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      tl_Awrk(:,:,Nnew))
!>
        CALL ad_dabc_u2d_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         ad_Awrk(:,:,Nnew))
!
!  Time-step adjoint horizontal diffusion terms.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!>          tl_Awrk(i,j,Nnew)=tl_Awrk(i,j,Nold)+                        &
!>   &                        Hfac(i,j)*                                &
!>   &                        (tl_FX(i,j)-tl_FX(i-1,j)+                 &
!>   &                         tl_FE(i,j+1)-tl_FE(i,j))
!>
            adfac=Hfac(i,j)*ad_Awrk(i,j,Nnew)
            ad_FE(i,j  )=ad_FE(i,j  )-adfac
            ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
            ad_FX(i-1,j)=ad_FX(i-1,j)-adfac
            ad_FX(i  ,j)=ad_FX(i  ,j)+adfac
            ad_Awrk(i,j,Nold)=ad_Awrk(i,j,Nold)+ad_Awrk(i,j,Nnew)
            ad_Awrk(i,j,Nnew)=0.0_r8
          END DO
        END DO
!
!  Compute XI- and ETA-components of the adjoint diffusive flux.
!
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
# ifdef MASKING
!>          tl_FE(i,j)=tl_FE(i,j)*pmask(i,j)
!>
            ad_FE(i,j)=ad_FE(i,j)*pmask(i,j)
# endif
!>          tl_FE(i,j)=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+      &
!>   &                                      Kh(i-1,j-1)+Kh(i,j-1))*     &
!>   &                 (tl_Awrk(i,j,Nold)-tl_Awrk(i,j-1,Nold))
!>
            adfac=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+           &
     &                                 Kh(i-1,j-1)+Kh(i,j-1))*          &
     &            ad_FE(i,j)
            ad_Awrk(i,j-1,Nold)=ad_Awrk(i,j-1,Nold)-adfac
            ad_Awrk(i,j  ,Nold)=ad_Awrk(i,j  ,Nold)+adfac
            ad_FE(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!>          tl_FX(i,j)=pmon_r(i,j)*Kh(i,j)*                             &
!>   &                 (tl_Awrk(i+1,j,Nold)-tl_Awrk(i,j,Nold))
!>
            adfac=pmon_r(i,j)*Kh(i,j)*ad_FX(i,j)
            ad_Awrk(i  ,j,Nold)=ad_Awrk(i  ,j,Nold)-adfac
            ad_Awrk(i+1,j,Nold)=ad_Awrk(i+1,j,Nold)+adfac
            ad_FX(i,j)=0.0_r8
          END DO
        END DO

      END DO
!
!  Set adjoint initial conditions.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
!>        tl_Awrk(i,j,Nold)=tl_A(i,j)
!>
          ad_A(i,j)=ad_A(i,j)+ad_Awrk(i,j,Nold)
          ad_Awrk(i,j,Nold)=0.0_r8
        END DO
      END DO
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_A)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Nghost,                                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
# endif
!>    CALL dabc_u2d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    tl_A)
!>
      CALL ad_dabc_u2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)

      RETURN
      END SUBROUTINE ad_conv_u2d_tile
!
!***********************************************************************
      SUBROUTINE ad_conv_v2d_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Nghost, NHsteps, DTsizeH,            &
     &                             Kh,                                  &
     &                             pm, pn, pmon_p, pnom_r,              &
# ifdef MASKING
     &                             vmask, pmask,                        &
# endif
     &                             ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_2d_mod, ONLY: ad_dabc_v2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, step

      real(r8) :: adfac, cff

      real(r8), dimension(LBi:UBi,LBj:UBj,2) :: ad_Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBi:UBi,LBj:UBj,1:2)=0.0_r8

      ad_FE(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_FX(IminS:ImaxS,JminS:JmaxS)=0.0_r8
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 2D state variable
!  at V-points.
!-----------------------------------------------------------------------
!
!  Compute metrics factor.
!
      cff=DTsizeH*0.25_r8
      DO j=JstrV,Jend
        DO i=Istr,Iend
          Hfac(i,j)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
        END DO
      END DO
      Nold=1
      Nnew=2
!
!------------------------------------------------------------------------
!  Adjoint of load convolved solution.
!------------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_A)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Nghost,                                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
# endif
!>    CALL dabc_v2d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    tl_A)
!>
      CALL ad_dabc_v2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        tl_A(i,j)=tl_Awrk(i,j,Nold)
!>
          ad_Awrk(i,j,Nold)=ad_Awrk(i,j,Nold)+ad_A(i,j)
          ad_A(i,j)=0.0_r8
        END DO
      END DO
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
!>      CALL mp_exchange2d (ng, tile, model, 1,                         &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      Nghost,                                     &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_Awrk(:,:,Nnew))
!>
        CALL ad_mp_exchange2d (ng, tile, model, 1,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Nghost,                                  &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Awrk(:,:,Nnew))
# endif
!>      CALL dabc_v2d_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      tl_Awrk(:,:,Nnew))
!>
        CALL ad_dabc_v2d_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         ad_Awrk(:,:,Nnew))
!
!  Time-step adjoint horizontal diffusion terms.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
!>          tl_Awrk(i,j,Nnew)=tl_Awrk(i,j,Nold)+                        &
!>   &                        Hfac(i,j)*                                &
!>   &                        (tl_FX(i+1,j)-tl_FX(i,j)+                 &
!>   &                         tl_FE(i,j)-tl_FE(i,j-1))
!>
            adfac=Hfac(i,j)*ad_Awrk(i,j,Nnew)
            ad_FE(i,j-1)=ad_FE(i,j-1)-adfac
            ad_FE(i,j  )=ad_FE(i,j  )+adfac
            ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
            ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
            ad_Awrk(i,j,Nold)=ad_Awrk(i,j,Nold)+ad_Awrk(i,j,Nnew)
            ad_Awrk(i,j,Nnew)=0.0_r8
          END DO
        END DO
!
!  Compute XI- and ETA-components of the adjoint diffusive flux.
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
!>          tl_FE(i,j)=pnom_r(i,j)*Kh(i,j)*                             &
!>   &                 (tl_Awrk(i,j+1,Nold)-tl_Awrk(i,j,Nold))
!>
            adfac=pnom_r(i,j)*Kh(i,j)*ad_FE(i,j)
            ad_Awrk(i,j  ,Nold)=ad_Awrk(i,j  ,Nold)-adfac
            ad_Awrk(i,j+1,Nold)=ad_Awrk(i,j+1,Nold)+adfac
            ad_FE(i,j)=0.0_r8
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
# ifdef MASKING
!>          tl_FX(i,j)=tl_FX(i,j)*pmask(i,j)
!>
            ad_FX(i,j)=ad_FX(i,j)*pmask(i,j)
# endif
!>          tl_FX(i,j)=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+      &
!>   &                                      Kh(i-1,j-1)+Kh(i,j-1))*     &
!>   &                 (tl_Awrk(i,j,Nold)-tl_Awrk(i-1,j,Nold))
!>
            adfac=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+           &
     &                                 Kh(i-1,j-1)+Kh(i,j-1))*          &
     &            ad_FX(i,j)
            ad_Awrk(i-1,j,Nold)=ad_Awrk(i-1,j,Nold)-adfac
            ad_Awrk(i  ,j,Nold)=ad_Awrk(i  ,j,Nold)+adfac
            ad_FX(i,j)=0.0_r8
          END DO
        END DO

      END DO
!
!  Set adjoint initial conditions.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
!>        tl_Awrk(i,j,Nold)=tl_A(i,j)
!>
          ad_A(i,j)=ad_A(i,j)+ad_Awrk(i,j,Nold)
          ad_Awrk(i,j,Nold)=0.0_r8
        END DO
      END DO
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_A)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Nghost,                                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
# endif
!>    CALL dabc_v2d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    tl_A)
!>
      CALL ad_dabc_v2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)

      RETURN
      END SUBROUTINE ad_conv_v2d_tile
#endif
      END MODULE ad_conv_2d_mod
