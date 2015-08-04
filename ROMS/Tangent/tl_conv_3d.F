#include "cppdefs.h"

      MODULE tl_conv_3d_mod

#if defined TANGENT && defined FOUR_DVAR && defined SOLVE3D
!
!svn $Id: tl_conv_3d.F 751 2015-01-07 22:56:36Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group       Andrew M. Moore   !
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
!     tl_A       3D tangent linear state variable to diffuse.          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     tl_A       Convolved 3D tangent linear state variable.           !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    tl_conv_r3d_tile  Tangent linear 3D convolution at RHO-points     !
!    tl_conv_u3d_tile  Tangent linear 3D convolution at U-points       !
!    tl_conv_v3d_tile  Tangent linear 3D convolution at V-points       !
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
      SUBROUTINE tl_conv_r3d_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Nghost, NHsteps, NVsteps,            &
     &                             DTsizeH, DTsizeV,                    &
     &                             Kh, Kv,                              &
     &                             pm, pn, pmon_u, pnom_v,              &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
     &                             Hz, z_r,                             &
     &                             tl_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY: bc_r3d_tile
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

      real(r8), intent(inout) :: tl_A(LBi:,LBj:,LBk:)
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
      real(r8), intent(inout) :: tl_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: cff, cff1

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk,2) :: tl_Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
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
#   ifdef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
#   endif
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DC
#  else
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 2D state variable
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
!>    CALL bc_r3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
!>   &                  A)
!>
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  tl_A)
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    A)
!>
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
# endif
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!>          Awrk(i,j,k,Nold)=A(i,j,k)
!>
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
!>            FX(i,j)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*           &
!>   &                (Awrk(i,j,k,Nold)-Awrk(i-1,j,k,Nold))
!>
              tl_FX(i,j)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*        &
     &                   (tl_Awrk(i,j,k,Nold)-tl_Awrk(i-1,j,k,Nold))
# ifdef MASKING
!>            FX(i,j)=FX(i,j)*umask(i,j)
!>
              tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
# endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
!>            FE(i,j)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*           &
!>   &                (Awrk(i,j,k,Nold)-Awrk(i,j-1,k,Nold))
!>
              tl_FE(i,j)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*        &
     &                   (tl_Awrk(i,j,k,Nold)-tl_Awrk(i,j-1,k,Nold))
# ifdef MASKING
!>            FE(i,j)=FE(i,j)*vmask(i,j)
!>
              tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
# endif
            END DO
          END DO
!
!  Time-step horizontal diffusion equation.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         Hfac(i,j)*                               &
!>   &                         (FX(i+1,j)-FX(i,j)+                      &
!>   &                          FE(i,j+1)-FE(i,j))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            Hfac(i,j)*                            &
     &                            (tl_FX(i+1,j)-tl_FX(i,j)+             &
     &                             tl_FE(i,j+1)-tl_FE(i,j))
            END DO
          END DO
        END DO
!
!  Apply boundary conditions.
!
!>      CALL bc_r3d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Awrk(:,:,:,Nnew))
!>
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    tl_Awrk(:,:,:,Nnew))
# ifdef DISTRIBUTE
!>      CALL mp_exchange3d (ng, tile, model, 1,                         &
!>   &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
!>   &                      Nghost,                                     &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      Awrk(:,:,:,Nnew))
!>
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Nghost,                                     &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_Awrk(:,:,:,Nnew))
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
!>          DC(i,0)=0.0_r8
!>
            tl_DC(i,0)=0.0_r8
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
!>            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
!>   &                     FC(i,k)*DC(i,k-1))
!>
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
!>          DC(i,N(ng))=0.0_r8
!>
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            DC(i,k)=DC(i,k)*Kv(i,j,k)
!>
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         DTsizeV*oHz(i,j,k)*                      &
!>   &                         (DC(i,k)-DC(i,k-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
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
!  terms for the tangent diffusion equation into tl_DC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,j,k)-FC(i,j,k-1)
!>            DC(i,k)=Awrk(i,j,k,Nold)*Hz(i,j,k)
!>
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*Hz(i,j,k)
            END DO
          END DO
!
!  Solve the tangent linear tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!>          DC(i,1)=cff*DC(i,1)
!>
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!>            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!>
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
!>          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!>   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!>   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!>
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!>          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!>
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
#    ifdef MASKING
!>          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*rmask(i,j)
!>
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*rmask(i,j)
#    endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>            Awrk(i,j,k,Nnew)=DC(i,k)
!>
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
#    ifdef MASKING
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*rmask(i,j)
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*rmask(i,j)
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
!>            FS(i,k)=FC(i,j,k)*(Awrk(i,j,k+1,Nold)-                    &
!>   &                           Awrk(i,j,k  ,Nold))
!>
              tl_FS(i,k)=FC(i,j,k)*(tl_Awrk(i,j,k+1,Nold)-              &
     &                              tl_Awrk(i,j,k  ,Nold))
#   ifdef MASKING
!>            FS(i,k)=FS(i,k)*rmask(i,j)
!>
              tl_FS(i,k)=tl_FS(i,k)*rmask(i,j)
#   endif
            END DO
          END DO
          DO i=Istr,Iend
!>          FS(i,0)=0.0_r8
!>
            tl_FS(i,0)=0.0_r8
!>          FS(i,N(ng))=0.0_r8
!>
            tl_FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         oHz(i,j,k)*(FS(i,k  )-                   &
!>   &                                     FS(i,k-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            oHz(i,j,k)*(tl_FS(i,k  )-             &
     &                                        tl_FS(i,k-1))
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
!>          A(i,j,k)=Awrk(i,j,k,Nold)
!>
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO
!>    CALL bc_r3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
!>   &                  A)
!>
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  tl_A)
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    A)
!>
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
# endif

      RETURN
      END SUBROUTINE tl_conv_r3d_tile
!
!***********************************************************************
      SUBROUTINE tl_conv_u3d_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Nghost, NHsteps, NVsteps,            &
     &                             DTsizeH, DTsizeV,                    &
     &                             Kh, Kv,                              &
     &                             pm, pn, pmon_r, pnom_p,              &
# ifdef MASKING
     &                             umask, pmask,                        &
# endif
     &                             Hz, z_r,                             &
     &                             tl_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY: bc_u3d_tile
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

      real(r8), intent(inout) :: tl_A(LBi:,LBj:,LBk:)
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
      real(r8), intent(inout) :: tl_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: cff, cff1

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk,2) :: tl_Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
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
#   ifdef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
#   endif
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DC
#  else
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FS
#  endif
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
!>    CALL bc_u3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
!>   &                  A)
!>
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  tl_A)
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    A)
!>
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
# endif
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
!>          Awrk(i,j,k,Nold)=A(i,j,k)
!>
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
!>            FX(i,j)=pmon_r(i,j)*Kh(i,j)*                              &
!>   &                (Awrk(i+1,j,k,Nold)-Awrk(i,j,k,Nold))
!>
              tl_FX(i,j)=pmon_r(i,j)*Kh(i,j)*                           &
     &                   (tl_Awrk(i+1,j,k,Nold)-tl_Awrk(i,j,k,Nold))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
!>            FE(i,j)=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+       &
!>   &                                     Kh(i-1,j-1)+Kh(i,j-1))*      &
!>   &                (Awrk(i,j,k,Nold)-Awrk(i,j-1,k,Nold))
!>
              tl_FE(i,j)=pnom_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+    &
     &                                        Kh(i-1,j-1)+Kh(i,j-1))*   &
     &                   (tl_Awrk(i,j,k,Nold)-tl_Awrk(i,j-1,k,Nold))
# ifdef MASKING
!>            FE(i,j)=FE(i,j)*pmask(i,j)
!>
              tl_FE(i,j)=tl_FE(i,j)*pmask(i,j)
# endif
            END DO
          END DO
!
!  Time-step horizontal diffusion equation.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         Hfac(i,j)*                               &
!>   &                         (FX(i,j)-FX(i-1,j)+                      &
!>   &                          FE(i,j+1)-FE(i,j))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            Hfac(i,j)*                            &
     &                            (tl_FX(i,j)-tl_FX(i-1,j)+             &
     &                             tl_FE(i,j+1)-tl_FE(i,j))
            END DO
          END DO
        END DO
!
!  Apply boundary conditions.
!
!>      CALL bc_u3d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Awrk(:,:,:,Nnew))
!>
        CALL bc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    tl_Awrk(:,:,:,Nnew))
# ifdef DISTRIBUTE
!>      CALL mp_exchange3d (ng, tile, model, 1,                         &
!>   &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
!>   &                      Nghost,                                     &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      Awrk(:,:,:,Nnew))
!>
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Nghost,                                     &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_Awrk(:,:,:,Nnew))
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
!>          DC(i,0)=0.0_r8
!>
            tl_DC(i,0)=0.0_r8
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
!>            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
!>   &                     FC(i,k)*DC(i,k-1))
!>
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=IstrU,Iend
!>          DC(i,N(ng))=0.0_r8
!>
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
!>            DC(i,k)=DC(i,k)*Kv(i,j,k)
!>
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         DTsizeV*oHz(i,j,k)*                      &
!>   &                         (DC(i,k)-DC(i,k-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
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
!  terms for the tangent linear diffusion equation into tl_DC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
              BC(i,k)=cff-FC(i,j,k)-FC(i,j,k-1)
!>            DC(i,k)=Awrk(i,j,k,Nold)*cff
!>
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*cff
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=IstrU,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!>          DC(i,1)=cff*DC(i,1)
!>
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=IstrU,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!>            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!>
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=IstrU,Iend
!>          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!>   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!>   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!>
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!>          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!>
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
#    ifdef MASKING
!>          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*umask(i,j)
!>
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*umask(i,j)
#    endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>            Awrk(i,j,k,Nnew)=DC(i,k)
!>
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
#    ifdef MASKING
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*umask(i,j)
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*umask(i,j)
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
!>            FS(i,k)=FC(i,j,k)*(Awrk(i,j,k+1,Nold)-                    &
!>   &                           Awrk(i,j,k  ,Nold))
!>
              tl_FS(i,k)=FC(i,j,k)*(tl_Awrk(i,j,k+1,Nold)-              &
     &                              tl_Awrk(i,j,k  ,Nold))
#   ifdef MASKING
!>            FS(i,k)=FS(i,k)*umask(i,j)
!>
              tl_FS(i,k)=tl_FS(i,k)*umask(i,j)
#   endif
            END DO
          END DO
          DO i=IstrU,Iend
!>          FS(i,0)=0.0_r8
!>
            tl_FS(i,0)=0.0_r8
!>          FS(i,N(ng))=0.0_r8
!>
            tl_FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         oHz(i,j,k)*(FS(i,k  )-                   &
!>   &                                     FS(i,k-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            oHz(i,j,k)*(tl_FS(i,k  )-             &
     &                                        tl_FS(i,k-1))
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
!>          A(i,j,k)=Awrk(i,j,k,Nold)
!>
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO
!>    CALL bc_u3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
!>   &                  A)
!>
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  tl_A)
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    A)
!>
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
# endif

      RETURN
      END SUBROUTINE tl_conv_u3d_tile
!
!***********************************************************************
      SUBROUTINE tl_conv_v3d_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Nghost, NHsteps, NVsteps,            &
     &                             DTsizeH, DTsizeV,                    &
     &                             Kh, Kv,                              &
     &                             pm, pn, pmon_p, pnom_r,              &
# ifdef MASKING
     &                             vmask, pmask,                        &
# endif
     &                             Hz, z_r,                             &
     &                             tl_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY: bc_v3d_tile
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

      real(r8), intent(inout) :: tl_A(LBi:,LBj:,LBk:)
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

      real(r8), intent(inout) :: tl_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: cff, cff1

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk,2) :: tl_Awrk

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
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
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DC
#   ifdef SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
#   endif
#  else
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FS
#  endif
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
!>    CALL bc_v3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
!>   &                  A)
!>
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  tl_A)
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    A)
!>
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
# endif
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
!>          Awrk(i,j,k,Nold)=A(i,j,k)
!>
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Integrate horizontal diffusion equation.
!-----------------------------------------------------------------------
!
      DO step=1,NHsteps
!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
!>            FX(i,j)=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+       &
!>   &                                     Kh(i-1,j-1)+Kh(i,j-1))*      &
!>   &                (Awrk(i,j,k,Nold)-Awrk(i-1,j,k,Nold))
!>
              tl_FX(i,j)=pmon_p(i,j)*0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+    &
     &                                        Kh(i-1,j-1)+Kh(i,j-1))*   &
     &                   (tl_Awrk(i,j,k,Nold)-tl_Awrk(i-1,j,k,Nold))
# ifdef MASKING
!>            FX(i,j)=FX(i,j)*pmask(i,j)
!>
              tl_FX(i,j)=tl_FX(i,j)*pmask(i,j)
# endif
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
!>            FE(i,j)=pnom_r(i,j)*Kh(i,j)*                              &
!>   &                (Awrk(i,j+1,k,Nold)-Awrk(i,j,k,Nold))
!>
              tl_FE(i,j)=pnom_r(i,j)*Kh(i,j)*                           &
     &                   (tl_Awrk(i,j+1,k,Nold)-tl_Awrk(i,j,k,Nold))
            END DO
          END DO
!
!  Time-step horizontal diffusion equation.
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         Hfac(i,j)*                               &
!>   &                         (FX(i+1,j)-FX(i,j)+                      &
!>   &                          FE(i,j)-FE(i,j-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            Hfac(i,j)*                            &
     &                            (tl_FX(i+1,j)-tl_FX(i,j)+             &
     &                             tl_FE(i,j)-tl_FE(i,j-1))
            END DO
          END DO
        END DO
!
!  Apply boundary conditions.
!
!>      CALL bc_v3d_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Awrk(:,:,:,Nnew))
!>
        CALL bc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    tl_Awrk(:,:,:,Nnew))
# ifdef DISTRIBUTE
!>      CALL mp_exchange3d (ng, tile, model, 1,                         &
!>   &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
!>   &                      Nghost,                                     &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      Awrk(:,:,:,Nnew))
!>
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
     &                      Nghost,                                     &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_Awrk(:,:,:,Nnew))
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
!>          DC(i,0)=0.0_r8
!>
            tl_DC(i,0)=0.0_r8
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
!>            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-                          &
!>   &                     Awrk(i,j,k  ,Nold)-                          &
!>   &                     FC(i,k)*DC(i,k-1))
!>
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
!>          DC(i,N(ng))=0.0_r8
!>
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            DC(i,k)=DC(i,k)*Kv(i,j,k)
!>
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         DTsizeV*oHz(i,j,k)*                      &
!>   &                         (DC(i,k)-DC(i,k-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
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
!  terms for the tangent linear diffusion equation into tl_DC.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
              BC(i,k)=cff-FC(i,j,k)-FC(i,j,k-1)
!>            DC(i,k)=Awrk(i,j,k,Nold)*cff
!>
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*cff
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!>          DC(i,1)=cff*DC(i,1)
!>
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!>            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!>
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
!>          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!>   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!>   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!>
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!>          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!>
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
#    ifdef MASKING
!>          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*vmask(i,j)
!>
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*vmask(i,j)
#    endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>            Awrk(i,j,k,Nnew)=DC(i,k)
!>
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
#    ifdef MASKING
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*vmask(i,j)
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*vmask(i,j)
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
!>            FS(i,k)=FC(i,j,k)*(Awrk(i,j,k+1,Nold)-                    &
!>   &                           Awrk(i,j,k  ,Nold))
!>
              tl_FS(i,k)=FC(i,j,k)*(tl_Awrk(i,j,k+1,Nold)-              &
     &                              tl_Awrk(i,j,k  ,Nold))
#   ifdef MASKING
!>            FS(i,k)=FS(i,k)*vmask(i,j)
!>
              tl_FS(i,k)=tl_FS(i,k)*vmask(i,j)
#   endif
            END DO
          END DO
          DO i=Istr,Iend
!>          FS(i,0)=0.0_r8
!>
            tl_FS(i,0)=0.0_r8
!>          FS(i,N(ng))=0.0_r8
!>
            tl_FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term. Notice that "oHz" is assumed to
!  be time invariant.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!>   &                         oHz(i,j,k)*(FS(i,k  )-                   &
!>   &                                     FS(i,k-1))
!>
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            oHz(i,j,k)*(tl_FS(i,k  )-             &
     &                                        tl_FS(i,k-1))
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
!>          A(i,j,k)=Awrk(i,j,k,Nold)
!>
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO
!>    CALL bc_v3d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
!>   &                  A)
!>
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  tl_A)
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
!>   &                    Nghost,                                       &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    A)
!>
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
# endif

      RETURN
      END SUBROUTINE tl_conv_v3d_tile
#endif
      END MODULE tl_conv_3d_mod
