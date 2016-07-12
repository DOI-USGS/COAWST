#include "cppdefs.h"
      MODULE ad_conv_bry3d_mod

#if defined ADJOINT && defined FOUR_DVAR && defined SOLVE3D && \
    defined ADJUST_BOUNDARY
!
!svn $Id: ad_conv_bry3d.F 795 2016-05-11 01:42:43Z arango $
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
!     ad_A       3D boundary state variable to diffuse                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A       Convolved 3D boundary state variable                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_conv_r3d_bry_tile  Tangent linear 3D boundary convolution at   !
!                            RHO-points                                !
!    ad_conv_u3d_bry_tile  Tangent linear 3D boundary convolution at   !
!                            U-points                                  !
!    ad_conv_v3d_bry_tile  Tangent linear 3D boundary convolution at   !
!                            V-points                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_conv_r3d_bry_tile (ng, tile, model, boundary,       &
     &                                 edge, LBij, UBij,                &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Nghost, NHsteps, NVsteps,        &
     &                                 DTsizeH, DTsizeV,                &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_u, pnom_v,          &
# ifdef MASKING
     &                                 rmask, umask, vmask,             &
# endif
     &                                 Hz, z_r,                         &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_bry3d_mod, ONLY: ad_bc_r3d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d_bry
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
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)
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
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: adfac, cff, cff1

      real(r8), dimension(LBij:UBij,LBk:UBk,2) :: ad_Awrk

      real(r8), dimension(JminS:JmaxS,LBk:UBk) :: ad_FE
      real(r8), dimension(IminS:ImaxS,LBk:UBk) :: ad_FX
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
#   ifdef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#   endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_DC
#  else
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBij:UBij,LBk:UBk,1:2)=0.0_r8
# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
      ad_DC(LBij:UBij,0:N(ng))=0.0_r8
#  else
      ad_FS(LBij:UBij,0:N(ng))=0.0_r8
#  endif
# endif
      ad_FE(JminS:JmaxS,LBk:UBk)=0.0_r8
      ad_FX(IminS:ImaxS,LBk:UBk)=0.0_r8
!
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
!-----------------------------------------------------------------------
!  Adjoint of load convolved solution.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_r3d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij, 1, N(ng),                       &
!>   &                      tl_A)
!>
      CALL ad_bc_r3d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij, 1, N(ng),                    &
     &                         ad_A)
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
!>            tl_A(j,k)=tl_Awrk(j,k,Nold)
!>
              ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                      &
     &                          ad_A(j,k)
              ad_A(j,k)=0.0_r8
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_A(i,k)=tl_Awrk(i,k,Nold)
!>
              ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                      &
     &                          ad_A(i,k)
              ad_A(i,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF

# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
#   ifdef SPLINES_VCONV
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation implicitly using
!  parabolic splines.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute basic state time-invariant coefficients.
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
              END DO
            END DO
          END IF
!
!  Adjoint backward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            DTsizeV*oHz(j,k)*                     &
!>   &                            (tl_DC(j,k)-tl_DC(j,k-1))
!>
                adfac=DTsizeV*oHz(j,k)*ad_Awrk(j,k,Nnew)
                ad_DC(j,k-1)=ad_DC(j,k-1)-adfac
                ad_DC(j,k  )=ad_DC(j,k  )+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
!>              tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!>
                ad_DC(j,k)=ad_DC(j,k)*Kv(i,j,k)
              END DO
            END DO
            DO k=1,N(ng)-1
              DO j=Jstr,Jend
!>              tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!>
                ad_DC(j,k+1)=ad_DC(j,k+1)-CF(j,k)*ad_DC(j,k)
              END DO
            END DO
            DO j=Jstr,Jend
!>            tl_DC(j,N(ng))=0.0_r8
!>
              ad_DC(j,N(ng))=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            DTsizeV*oHz(i,k)*                     &
!>   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!>
                adfac=DTsizeV*oHz(i,k)*ad_Awrk(i,k,Nnew)
                ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
                ad_DC(i,k  )=ad_DC(i,k  )+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
!>              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!>
                ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
              END DO
            END DO
            DO k=1,N(ng)-1
              DO i=Istr,Iend
!>              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
                ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
              END DO
            END DO
            DO i=Istr,Iend
!>            tl_DC(i,N(ng))=0.0_r8
!>
              ad_DC(i,N(ng))=0.0_r8
            END DO
          END IF
!
!  Adjoint of LU decomposition and forward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
!>              tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                    &
!>   &                          tl_Awrk(j,k  ,Nold)-                    &
!>   &                          FC(j,k)*tl_DC(j,k-1))
!>
                adfac=cff*ad_DC(j,k)
                ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k)*adfac
                ad_DC(j,k)=0.0_r8
              END DO
            END DO
            DO j=Jstr,Jend
!>            tl_DC(j,0)=0.0_r8
!>
              ad_DC(j,0)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!>              tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                    &
!>   &                          tl_Awrk(i,k  ,Nold)-                    &
!>   &                          FC(i,k)*tl_DC(i,k-1))
!>
                adfac=cff*ad_DC(i,k)
                ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
                ad_DC(i,k)=0.0_r8
              END DO
            END DO
            DO i=Istr,Iend
!>            tl_DC(i,0)=0.0_r8
!>
              ad_DC(i,0)=0.0_r8
            END DO
          END IF
        END IF
      END DO

#   else
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation implicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Compute diagonal matrix coefficients BC and CF.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
                BC(j,k)=Hz(i,j,k)-FC(j,k)-FC(j,k-1)
              END DO
            END DO
            DO j=Jstr,Jend
              cff=1.0_r8/BC(j,1)
              CF(j,1)=cff*FC(j,1)
            END DO
            DO k=2,N(ng)-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                CF(j,k)=cff*FC(j,k)
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
                BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              END DO
            END DO
            DO i=Istr,Iend
              cff=1.0_r8/BC(i,1)
              CF(i,1)=cff*FC(i,1)
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                CF(i,k)=cff*FC(i,k)
              END DO
            END DO
          END IF
!
!  Adjoint of compute new solution by back substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
!>          DO k=N(ng)-1,1,-1
!>
            DO k=1,N(ng)-1
              DO j=Jstr,Jend
#    ifdef MASKING
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*rmask(i,j)
!>
                ad_Awrk(j,k,Nnew)=ad_Awrk(j,k,Nnew)*rmask(i,j)
#    endif
!>              tl_Awrk(j,k,Nnew)=tl_DC(j,k)
!>
                ad_DC(j,k)=ad_DC(j,k)+                                  &
     &                     ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
!>              tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!>
                ad_DC(j,k+1)=-CF(j,k)*ad_DC(j,k)
              END DO
            END DO
            DO j=Jstr,Jend
#    ifdef MASKING
!>            tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*rmask(i,j)
!>
              ad_Awrk(j,N(ng),Nnew)=ad_Awrk(j,N(ng),Nnew)*rmask(i,j)

#    endif
!>            tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
!>
              ad_DC(j,N(ng))=ad_DC(j,N(ng))+                            &
     &                       ad_Awrk(j,N(ng),Nnew)
              ad_Awrk(j,N(ng),Nnew)=0.0_r8
!>            tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                           &
!>   &                        FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/          &
!>   &                       (BC(j,N(ng))-                              &
!>   &                        FC(j,N(ng)-1)*CF(j,N(ng)-1))
!>
              adfac=ad_DC(j,N(ng))/                                     &
     &              (BC(j,N(ng))-                                       &
     &               FC(j,N(ng)-1)*CF(j,N(ng)-1))
              ad_DC(j,N(ng)-1)=ad_DC(j,N(ng)-1)-                        &
     &                         FC(j,N(ng)-1)*adfac
              ad_DC(j,N(ng)  )=adfac
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
!>          DO k=N(ng)-1,1,-1
!>
            DO k=1,N(ng)-1
              DO i=Istr,Iend
#    ifdef MASKING
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*rmask(i,j)
!>
                ad_Awrk(i,k,Nnew)=ad_Awrk(i,k,Nnew)*rmask(i,j)
#    endif
!>              tl_Awrk(i,k,Nnew)=tl_DC(i,k)
!>
                ad_DC(i,k)=ad_DC(i,k)+                                  &
     &                     ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
!>              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
                ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
              END DO
            END DO
            DO i=Istr,Iend
#    ifdef MASKING
!>            tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*rmask(i,j)
!>
              ad_Awrk(i,N(ng),Nnew)=ad_Awrk(i,N(ng),Nnew)*rmask(i,j)
#    endif
!>            tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
!>
              ad_DC(i,N(ng))=ad_DC(i,N(ng))+                            &
     &                       ad_Awrk(i,N(ng),Nnew)
              ad_Awrk(i,N(ng),Nnew)=0.0_r8
!>            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                           &
!>   &                        FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!>   &                       (BC(i,N(ng))-                              &
!>   &                        FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>
              adfac=ad_DC(i,N(ng))/                                     &
     &              (BC(i,N(ng))-                                       &
     &               FC(i,N(ng)-1)*CF(i,N(ng)-1))
              ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-                        &
     &                         FC(i,N(ng)-1)*adfac
              ad_DC(i,N(ng)  )=adfac
            END DO
          END IF
!
!  Solve the adjoint tridiagonal system.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
!>          DO k=2,N(ng)-1
!>
            DO k=N(ng)-1,2,-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
!>              tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
!>
                adfac=cff*ad_DC(j,k)
                ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k-1)*adfac
                ad_DC(j,k  )=adfac
              END DO
            END DO
            DO j=Jstr,Jend
              cff=1.0_r8/BC(j,1)
!>            tl_DC(j,1)=cff*tl_DC(j,1)
!>
              ad_DC(j,1)=cff*ad_DC(j,1)
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
!>          DO k=2,N(ng)-1
!>
            DO k=N(ng)-1,2,-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!>              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
!>
                adfac=cff*ad_DC(i,k)
                ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
                ad_DC(i,k  )=adfac
              END DO
            END DO
            DO i=Istr,Iend
              cff=1.0_r8/BC(i,1)
!>            tl_DC(i,1)=cff*tl_DC(i,1)
!>
              ad_DC(i,1)=cff*ad_DC(i,1)
            END DO
          END IF
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_DC(j,k)=tl_Awrk(j,k,Nold)*Hz(i,j,k)
!>
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            Hz(i,j,k)*ad_DC(j,k)
                ad_DC(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_DC(i,k)=tl_Awrk(i,k,Nold)*Hz(i,j,k)
!>
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            Hz(i,j,k)*ad_DC(i,k)
                ad_DC(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
#   endif

#  else
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation explicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Time-step adjoint vertical diffusive term. Notice that "oHz" is
!  assumed to be time invariant.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            oHz(j,k)*(tl_FS(j,k  )-               &
!>   &                                      tl_FS(j,k-1))
!>
                adfac=oHz(j,k)*ad_Awrk(j,k,Nnew)
                ad_FS(j,k-1)=ad_FS(j,k-1)-adfac
                ad_FS(j,k  )=ad_FS(j,k  )+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            oHz(i,k)*(tl_FS(i,k  )-               &
!>   &                                      tl_FS(i,k-1))
!>
                adfac=oHz(i,k)*ad_Awrk(i,k,Nnew)
                ad_FS(i,k-1)=ad_FS(i,k-1)-adfac
                ad_FS(i,k  )=ad_FS(i,k  )+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute adjoint vertical diffusive flux.  Notice that "FC" is
!  assumed to be time invariant.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
!>            tl_FS(j,N(ng))=0.0_r8
!>
              ad_FS(j,N(ng))=0.0_r8
!>            tl_FS(j,0)=0.0_r8
!>
              ad_FS(j,0)=0.0_r8
              DO k=1,N(ng)-1
#   ifdef MASKING
!>              tl_FS(j,k)=tl_FS(j,k)*rmask(i,j)
!>
                ad_FS(j,k)=ad_FS(j,k)*rmask(i,j)
#   endif
!>              tl_FS(j,k)=FC(j,k)*(tl_Awrk(j,k+1,Nold)-                &
!>   &                              tl_Awrk(j,k  ,Nold))
!>
                adfac=FC(j,k)*ad_FS(j,k)
                ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                ad_FS(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
!>            tl_FS(i,N(ng))=0.0_r8
!>
              ad_FS(i,N(ng))=0.0_r8
!>            tl_FS(i,0)=0.0_r8
!>
              ad_FS(i,0)=0.0_r8
              DO k=1,N(ng)-1
#   ifdef MASKING
!>              tl_FS(i,k)=tl_FS(i,k)*rmask(i,j)
!>
                ad_FS(i,k)=ad_FS(i,k)*rmask(i,j)
#   endif
!>              tl_FS(i,k)=FC(i,k)*(tl_Awrk(i,k+1,Nold)-                &
!>   &                              tl_Awrk(i,k  ,Nold))
!>
                adfac=FC(i,k)*ad_FS(i,k)
                ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                ad_FS(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Integrate adjoint horizontal diffusion equation.
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
!  Apply adjoint boundary conditions.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
!>   &                          LBij, UBij, 1, N(ng),                   &
!>   &                          Nghost,                                 &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Awrk(:,:,Nnew))
!>
        CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,        &
     &                             LBij, UBij, 1, N(ng),                &
     &                             Nghost,                              &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Awrk(:,:,Nnew))
# endif
!>      CALL bc_r3d_bry_tile (ng, tile, boundary,                       &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        tl_Awrk(:,:,Nnew))
!>
        CALL ad_bc_r3d_bry_tile (ng, tile, boundary,                    &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           ad_Awrk(:,:,Nnew))
!
!  Time-step adjoint horizontal diffusion equation.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            Hfac(j)*                              &
!>   &                            (tl_FE(j+1,k)-tl_FE(j,k))
!>
                adfac=Hfac(j)*ad_Awrk(j,k,Nnew)
                ad_FE(j  ,k)=ad_FE(j  ,k)-adfac
                ad_FE(j+1,k)=ad_FE(j+1,k)+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            Hfac(i)*                              &
!>   &                            (tl_FX(i+1,k)-tl_FX(i,k))
!>
                adfac=Hfac(i)*ad_Awrk(i,k,Nnew)
                ad_FX(i  ,k)=ad_FX(i  ,k)-adfac
                ad_FX(i+1,k)=ad_FX(i+1,k)+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute XI- and ETA-components of the adjoint diffusive flux.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend+1
# ifdef MASKING
!>              tl_FE(j,k)=tl_FE(j,k)*vmask(i,j)
!>
                ad_FE(j,k)=ad_FE(j,k)*vmask(i,j)
# endif
!>              tl_FE(j,k)=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*      &
!>   &                     (tl_Awrk(j  ,k,Nold)-                        &
!>   &                      tl_Awrk(j-1,k,Nold))
!>
                adfac=pnom_v(i,j)*0.5_r8*(Kh(i,j-1)+Kh(i,j))*ad_FE(j,k)
                ad_Awrk(j-1,k,Nold)=ad_Awrk(j-1,k,Nold)-adfac
                ad_Awrk(j  ,k,Nold)=ad_Awrk(j  ,k,Nold)+adfac
                ad_FE(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend+1
# ifdef MASKING
!>              tl_FX(i,k)=tl_FX(i,k)*umask(i,j)
!>
                ad_FX(i,k)=ad_FX(i,k)*umask(i,j)
# endif
!>              tl_FX(i,k)=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*      &
!>   &                     (tl_Awrk(i  ,k,Nold)-                        &
!>   &                      tl_Awrk(i-1,k,Nold))
!>
                adfac=pmon_u(i,j)*0.5_r8*(Kh(i-1,j)+Kh(i,j))*ad_FX(i,k)
                ad_Awrk(i-1,k,Nold)=ad_Awrk(i-1,k,Nold)-adfac
                ad_Awrk(i  ,k,Nold)=ad_Awrk(i  ,k,Nold)+adfac
                ad_FX(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set adjoint initial conditions.
!-----------------------------------------------------------------------
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
!>            tl_Awrk(j,k,Nold)=tl_A(j,k)
!>
              ad_A(j,k)=ad_A(j,k)+ad_Awrk(j,k,Nold)
              ad_Awrk(j,k,Nold)=0.0_r8
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr-1,Iend+1
!>            tl_Awrk(i,k,Nold)=tl_A(i,k)
!>
              ad_A(i,k)=ad_A(i,k)+ad_Awrk(i,k,Nold)
              ad_Awrk(i,k,Nold)=0.0_r8
            END DO
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_r3d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij, 1, N(ng),                       &
!>   &                      tl_A)
!>
      CALL ad_bc_r3d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij, 1, N(ng),                    &
     &                         ad_A)

      RETURN
      END SUBROUTINE ad_conv_r3d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_conv_u3d_bry_tile (ng, tile, model, boundary,       &
     &                                 edge, LBij, UBij,                &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Nghost, NHsteps, NVsteps,        &
     &                                 DTsizeH, DTsizeV,                &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_r, pnom_p,          &
# ifdef MASKING
     &                                 umask, pmask,                    &
# endif
     &                                 Hz, z_r,                         &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_bry3d_mod, ONLY: ad_bc_u3d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d_bry
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary
      integer, intent(in) :: edge(4)
      integer, intent(in) :: LBij, UBij
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Nghost
      integer, intent(in) :: NHsteps, NVsteps

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
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)
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
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: adfac, cff, cff1

      real(r8), dimension(LBij:UBij,LBk:UBk,2) :: ad_Awrk

      real(r8), dimension(JminS:JmaxS,LBk:UBk) :: ad_FE
      real(r8), dimension(IminS:ImaxS,LBk:UBk) :: ad_FX
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
#   ifdef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#   endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_DC
#  else
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBij:UBij,LBk:UBk,1:2)=0.0_r8
# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
      ad_DC(LBij:UBij,0:N(ng))=0.0_r8
#  else
      ad_FS(LBij:UBij,0:N(ng))=0.0_r8
#  endif
# endif
      ad_FE(JminS:JmaxS,LBk:UBk)=0.0_r8
      ad_FX(IminS:ImaxS,LBk:UBk)=0.0_r8
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
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
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
!-----------------------------------------------------------------------
!  Adjoint of load convolved solution.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_u3d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij, 1, N(ng),                       &
!>   &                      tl_A)
!>
      CALL ad_bc_u3d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij, 1, N(ng),                    &
     &                         ad_A)
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
!>            tl_A(j,k)=tl_Awrk(j,k,Nold)
!>
              ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                      &
     &                          ad_A(j,k)
              ad_A(j,k)=0.0_r8
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_A(i,k)=tl_Awrk(i,k,Nold)
!>
              ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                      &
     &                          ad_A(i,k)
              ad_A(i,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF

# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
#   ifdef SPLINES_VCONV
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation implicitly using
!  parabolic splines.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute basic state time-invariant coefficients.
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
              END DO
            END DO
          END IF
!
!  Adjoint backward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            DTsizeV*oHz(j,k)*                     &
!>   &                            (tl_DC(j,k)-tl_DC(j,k-1))
!>
                adfac=DTsizeV*oHz(j,k)*ad_Awrk(j,k,Nnew)
                ad_DC(j,k-1)=ad_DC(j,k-1)-adfac
                ad_DC(j,k  )=ad_DC(j,k  )+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
!>              tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!>
                ad_DC(j,k)=ad_DC(j,k)*Kv(i,j,k)
              END DO
            END DO
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
!>              tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!>
                ad_DC(j,k+1)=ad_DC(j,k+1)-CF(j,k)*ad_DC(j,k)
              END DO
            END DO
            DO j=Jstr,Jend
!>            tl_DC(j,N(ng))=0.0_r8
!>
              ad_DC(j,N(ng))=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=IstrU,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            DTsizeV*oHz(i,k)*                     &
!>   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!>
                adfac=DTsizeV*oHz(i,k)*ad_Awrk(i,k,Nnew)
                ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
                ad_DC(i,k  )=ad_DC(i,k  )+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
!>              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!>
                ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
              END DO
            END DO
            DO k=N(ng)-1,1,-1
              DO i=IstrU,Iend
!>              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
                ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
              END DO
            END DO
            DO i=IstrU,Iend
!>            tl_DC(i,N(ng))=0.0_r8
!>
              ad_DC(i,N(ng))=0.0_r8
            END DO
          END IF
!
!  Adjoint of LU decomposition and forward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=N(ng)-1,1,-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
!>              tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                    &
!>   &                          tl_Awrk(j,k  ,Nold)-                    &
!>   &                          FC(j,k)*tl_DC(j,k-1))
!>
                adfac=cff*ad_DC(j,k)
                ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k)*adfac
                ad_DC(j,k)=0.0_r8
              END DO
            END DO
            DO j=Jstr,Jend
!>            tl_DC(j,0)=0.0_r8
!>
              ad_DC(j,0)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=N(ng)-1,1,-1
              DO i=IstrU,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!>              tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                    &
!>   &                          tl_Awrk(i,k  ,Nold)-                    &
!>   &                          FC(i,k)*tl_DC(i,k-1))
!>
                adfac=cff*ad_DC(i,k)
                ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
                ad_DC(i,k)=0.0_r8
              END DO
            END DO
            DO i=IstrU,Iend
!>            tl_DC(i,0)=0.0_r8
!>
              ad_DC(i,0)=0.0_r8
            END DO
          END IF
        END IF
      END DO

#   else
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation implicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Compute diagonal matrix coefficients BC and CF.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                BC(j,k)=cff-FC(j,k)-FC(j,k-1)
              END DO
            END DO
            DO j=Jstr,Jend
              cff=1.0_r8/BC(j,1)
              CF(j,1)=cff*FC(j,1)
            END DO
            DO k=2,N(ng)-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                CF(j,k)=cff*FC(j,k)
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                BC(i,k)=cff-FC(i,k)-FC(i,k-1)
              END DO
            END DO
            DO i=IstrU,Iend
              cff=1.0_r8/BC(i,1)
              CF(i,1)=cff*FC(i,1)
            END DO
            DO k=2,N(ng)-1
              DO i=IstrU,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                CF(i,k)=cff*FC(i,k)
              END DO
            END DO
          END IF
!
!  Compute new solution by back substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
!>          DO k=N(ng)-1,1,-1
!>
            DO k=1,N(ng)-1
              DO j=Jstr,Jend
#    ifdef MASKING
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*umask(i,j)
!>
                ad_Awrk(j,k,Nnew)=ad_Awrk(j,k,Nnew)*umask(i,j)
#    endif
!>              tl_Awrk(j,k,Nnew)=tl_DC(j,k)
!>
                ad_DC(j,k)=ad_DC(j,k)+                                  &
     &                     ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
!>              tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!>
                ad_DC(j,k+1)=-CF(j,k)*ad_DC(j,k)
              END DO
            END DO
            DO j=Jstr,Jend
#    ifdef MASKING
!>            tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*umask(i,j)
!>
              ad_Awrk(j,N(ng),Nnew)=ad_Awrk(j,N(ng),Nnew)*umask(i,j)
#    endif
!>            tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
!>
              ad_DC(j,N(ng))=ad_DC(j,N(ng))+                            &
     &                       ad_Awrk(j,N(ng),Nnew)
              ad_Awrk(j,N(ng),Nnew)=0.0_r8
!>            tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                           &
!>   &                        FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/          &
!>   &                       (BC(j,N(ng))-                              &
!>   &                        FC(j,N(ng)-1)*CF(j,N(ng)-1))
!>
              adfac=ad_DC(j,N(ng))/                                     &
     &              (BC(j,N(ng))-                                       &
     &               FC(j,N(ng)-1)*CF(j,N(ng)-1))
              ad_DC(j,N(ng)-1)=ad_DC(j,N(ng)-1)-                        &
     &                         FC(j,N(ng)-1)*adfac
              ad_DC(j,N(ng)  )=adfac
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
!>          DO k=N(ng)-1,1,-1
!>
            DO k=1,N(ng)-1
              DO i=IstrU,Iend
#    ifdef MASKING
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*umask(i,j)
!>
                ad_Awrk(i,k,Nnew)=ad_Awrk(i,k,Nnew)*umask(i,j)
#    endif
!>              tl_Awrk(i,k,Nnew)=tl_DC(i,k)
!>
                ad_DC(i,k)=ad_DC(i,k)+                                  &
     &                     ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
!>              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
                ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
              END DO
            END DO
            DO i=IstrU,Iend
#    ifdef MASKING
!>            tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*umask(i,j)
!>
              ad_Awrk(i,N(ng),Nnew)=ad_Awrk(i,N(ng),Nnew)*umask(i,j)
#    endif
!>            tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
!>
              ad_DC(i,N(ng))=ad_DC(i,N(ng))+                            &
     &                       ad_Awrk(i,N(ng),Nnew)
              ad_Awrk(i,N(ng),Nnew)=0.0_r8
!>            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                           &
!>   &                        FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!>   &                       (BC(i,N(ng))-                              &
!>   &                        FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>
              adfac=ad_DC(i,N(ng))/                                     &
     &              (BC(i,N(ng))-                                       &
     &               FC(i,N(ng)-1)*CF(i,N(ng)-1))
              ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-                        &
     &                         FC(i,N(ng)-1)*adfac
              ad_DC(i,N(ng)  )=adfac
            END DO
          END IF
!
!  Adjoint of solve the tridiagonal system.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
!>          DO k=2,N(ng)-1
!>
            DO k=N(ng)-1,2,-1
              DO j=Jstr,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
!>              tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
!>
                adfac=cff*ad_DC(j,k)
                ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k-1)*adfac
                ad_DC(j,k  )=adfac
              END DO
            END DO
            DO j=Jstr,Jend
              cff=1.0_r8/BC(j,1)
!>            tl_DC(j,1)=cff*tl_DC(j,1)
!>
              ad_DC(j,1)=cff*ad_DC(j,1)
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
!>          DO k=2,N(ng)-1
!>
            DO k=N(ng)-1,2,-1
              DO i=IstrU,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!>              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
!>
                adfac=cff*ad_DC(i,k)
                ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
                ad_DC(i,k  )=adfac
              END DO
            END DO
            DO i=IstrU,Iend
              cff=1.0_r8/BC(i,1)
!>            tl_DC(i,1)=cff*tl_DC(i,1)
!>
              ad_DC(i,1)=cff*ad_DC(i,1)
            END DO
          END IF
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
!>              tl_DC(j,k)=tl_Awrk(j,k,Nold)*cff
!>
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            cff*ad_DC(j,k)
                ad_DC(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
!>              tl_DC(i,k)=tl_Awrk(i,k,Nold)*cff
!>
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            cff*ad_DC(i,k)
                ad_DC(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
#   endif

#  else
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation explicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Time-step adjoint vertical diffusive term. Notice that "oHz" is
!  assumed to be time invariant.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            oHz(j,k)*(tl_FS(j,k  )-               &
!>   &                                      tl_FS(j,k-1))
!>
                adfac=oHz(j,k)*ad_Awrk(j,k,Nnew)
                ad_FS(j,k-1)=ad_FS(j,k-1)-adfac
                ad_FS(j,k  )=ad_FS(j,k  )+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            oHz(i,k)*(tl_FS(i,k  )-               &
!>   &                                      tl_FS(i,k-1))
!>
                adfac=oHz(i,k)*ad_Awrk(i,k,Nnew)
                ad_FS(i,k-1)=ad_FS(i,k-1)-adfac
                ad_FS(i,k  )=ad_FS(i,k  )+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute adjoint vertical diffusive flux.  Notice that "FC" is
!  assumed to be time invariant.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=Jstr,Jend
!>            tl_FS(j,N(ng))=0.0_r8
!>
              ad_FS(j,N(ng))=0.0_r8
!>            tl_FS(j,0)=0.0_r8
!>
              ad_FS(j,0)=0.0_r8
              DO k=1,N(ng)-1
#   ifdef MASKING
!>              tl_FS(j,k)=tl_FS(j,k)*umask(i,j)
!>
                ad_FS(j,k)=ad_FS(j,k)*umask(i,j)
#   endif
!>              tl_FS(j,k)=FC(j,k)*(tl_Awrk(j,k+1,Nold)-                &
!>   &                              tl_Awrk(j,k  ,Nold))
!>
                adfac=FC(j,k)*ad_FS(j,k)
                ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                ad_FS(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=IstrU,Iend
!>            tl_FS(i,N(ng))=0.0_r8
!>
              ad_FS(i,N(ng))=0.0_r8
!>            tl_FS(i,0)=0.0_r8
!>
              ad_FS(i,0)=0.0_r8
              DO k=1,N(ng)-1
#   ifdef MASKING
!>              tl_FS(i,k)=tl_FS(i,k)*umask(i,j)
!>
                ad_FS(i,k)=ad_FS(i,k)*umask(i,j)
#   endif
!>              tl_FS(i,k)=FC(i,k)*(tl_Awrk(i,k+1,Nold)-                &
!>   &                              tl_Awrk(i,k  ,Nold))
!>
                adfac=FC(i,k)*ad_FS(i,k)
                ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                ad_FS(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Integrate adjoint horizontal diffusion equation.
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
!  Apply boundary conditions.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
!>   &                          LBij, UBij, 1, N(ng),                   &
!>   &                          Nghost,                                 &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Awrk(:,:,Nnew))
!>
        CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,        &
     &                             LBij, UBij, 1, N(ng),                &
     &                             Nghost,                              &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Awrk(:,:,Nnew))
# endif
!>      CALL bc_u3d_bry_tile (ng, tile, boundary,                       &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        tl_Awrk(:,:,Nnew))
!>
        CALL ad_bc_u3d_bry_tile (ng, tile, boundary,                    &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           ad_Awrk(:,:,Nnew))
!
!  Time-step adjoint horizontal diffusion equation.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            Hfac(j)*                              &
!>   &                            (tl_FE(j+1,k)-tl_FE(j,k))
!>
                adfac=Hfac(j)*ad_Awrk(j,k,Nnew)
                ad_FE(j  ,k)=ad_FE(j  ,k)-adfac
                ad_FE(j+1,k)=ad_FE(j+1,k)+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            Hfac(i)*                              &
!>   &                            (tl_FX(i,k)-tl_FX(i-1,k))
!>
                adfac=Hfac(i)*ad_Awrk(i,k,Nnew)
                ad_FX(i-1,k)=ad_FX(i-1,k)-adfac
                ad_FX(i  ,k)=ad_FX(i  ,k)+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute XI- and ETA-components of adjoint diffusive flux.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=Jstr,Jend+1
# ifdef MASKING
!>              tl_FE(j,k)=tl_FE(j,k)*pmask(i,j)
!>
                ad_FE(j,k)=ad_FE(j,k)*pmask(i,j)
# endif
!>              tl_FE(j,k)=pnom_p(i,j)*                                 &
!>   &                     0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+              &
!>   &                              Kh(i-1,j-1)+Kh(i,j-1))*             &
!>   &                     (tl_Awrk(j  ,k,Nold)-                        &
!>   &                      tl_Awrk(j-1,k,Nold))
!>
                adfac=pnom_p(i,j)*                                      &
     &                0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+                   &
     &                         Kh(i-1,j-1)+Kh(i,j-1))*ad_FE(j,k)
                ad_Awrk(j-1,k,Nold)=ad_Awrk(j-1,k,Nold)-adfac
                ad_Awrk(j  ,k,Nold)=ad_Awrk(j  ,k,Nold)+adfac
                ad_FE(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=IstrU-1,Iend
!>              tl_FX(i,k)=pmon_r(i,j)*Kh(i,j)*                         &
!>   &                     (tl_Awrk(i+1,k,Nold)-                        &
!>   &                      tl_Awrk(i  ,k,Nold))
!>
                adfac=pmon_r(i,j)*Kh(i,j)*ad_FX(i,k)
                ad_Awrk(i  ,k,Nold)=ad_Awrk(i  ,k,Nold)-adfac
                ad_Awrk(i+1,k,Nold)=ad_Awrk(i+1,k,Nold)+adfac
                ad_FX(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set adjoint initial conditions.
!-----------------------------------------------------------------------
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
!>            tl_Awrk(j,k,Nold)=tl_A(j,k)
!>
              ad_A(j,k)=ad_A(j,k)+ad_Awrk(j,k,Nold)
              ad_Awrk(j,k,Nold)=0.0_r8
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend+1
!>            tl_Awrk(i,k,Nold)=tl_A(i,k)
!>
              ad_A(i,k)=ad_A(i,k)+ad_Awrk(i,k,Nold)
              ad_Awrk(i,k,Nold)=0.0_r8
            END DO
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_u3d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij, 1, N(ng),                       &
!>   &                      tl_A)
!>
      CALL ad_bc_u3d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij, 1, N(ng),                    &
     &                         ad_A)

      RETURN
      END SUBROUTINE ad_conv_u3d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_conv_v3d_bry_tile (ng, tile, model, boundary,       &
     &                                 edge, LBij, UBij,                &
     &                                 LBi, UBi, LBj, UBj, LBk, UBk,    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Nghost, NHsteps, NVsteps,        &
     &                                 DTsizeH, DTsizeV,                &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_p, pnom_r,          &
# ifdef MASKING
     &                                 vmask, pmask,                    &
# endif
     &                                 Hz, z_r,                         &
     &                                 ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_bry3d_mod, ONLY: ad_bc_v3d_bry_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d_bry
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
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)
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
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: Lconvolve

      integer :: Nnew, Nold, Nsav, i, j, k, step

      real(r8) :: adfac, cff, cff1

      real(r8), dimension(LBij:UBij,LBk:UBk,2) :: ad_Awrk

      real(r8), dimension(JminS:JmaxS,LBk:UBk) :: ad_FE
      real(r8), dimension(IminS:ImaxS,LBk:UBk) :: ad_FX
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
#   ifdef SPLINES_VCONV
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
#   endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_DC
#  else
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_FS
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Awrk(LBij:UBij,LBk:UBk,1:2)=0.0_r8
# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
      ad_DC(LBij:UBij,0:N(ng))=0.0_r8
#  else
      ad_FS(LBij:UBij,0:N(ng))=0.0_r8
#  endif
# endif
      ad_FE(JminS:JmaxS,LBk:UBk)=0.0_r8
      ad_FX(IminS:ImaxS,LBk:UBk)=0.0_r8
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
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
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
!-----------------------------------------------------------------------
!  Adjoint of load convolved solution.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_v3d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij, 1, N(ng),                       &
!>   &                      tl_A)
!>
      CALL ad_bc_v3d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij, 1, N(ng),                    &
     &                         ad_A)
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
!>            tl_A(j,k)=tl_Awrk(j,k,Nold)
!>
              ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                      &
     &                          ad_A(j,k)
              ad_A(j,k)=0.0_r8
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_A(i,k)=tl_Awrk(i,k,Nold)
!>
              ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                      &
     &                          ad_A(i,k)
              ad_A(i,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF

# ifdef VCONVOLUTION
#  ifdef IMPLICIT_VCONV
#   ifdef SPLINES_VCONV
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation implicitly using
!  parabolic splines.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute basic state time-invariant coefficients.
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
              END DO
            END DO
          END IF
!
!  Adjoint backward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=JstrV,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            DTsizeV*oHz(j,k)*                     &
!>   &                            (tl_DC(j,k)-tl_DC(j,k-1))
!>
                adfac=DTsizeV*oHz(j,k)*ad_Awrk(j,k,Nnew)
                ad_DC(j,k-1)=ad_DC(j,k-1)-adfac
                ad_DC(j,k  )=ad_DC(j,k  )+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
!>              tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!>
                ad_DC(j,k)=ad_DC(j,k)*Kv(i,j,k)
              END DO
            END DO
            DO k=1,N(ng)-1
              DO j=JstrV,Jend
!>              tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!>
                ad_DC(j,k+1)=ad_DC(j,k+1)-CF(j,k)*ad_DC(j,k)
              END DO
            END DO
            DO j=JstrV,Jend
!>            tl_DC(j,N(ng))=0.0_r8
!>
              ad_DC(j,N(ng))=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            DTsizeV*oHz(i,k)*                     &
!>   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!>
                adfac=DTsizeV*oHz(i,k)*ad_Awrk(i,k,Nnew)
                ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
                ad_DC(i,k  )=ad_DC(i,k  )+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
!>              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!>
                ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
              END DO
            END DO
            DO k=1,N(ng)-1
              DO i=Istr,Iend
!>              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
                ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
              END DO
            END DO
            DO i=Istr,Iend
!>            tl_DC(i,N(ng))=0.0_r8
!>
              ad_DC(i,N(ng))=0.0_r8
            END DO
          END IF
!
!  Adjoint of LU decomposition and forward substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=N(ng)-1,1,-1
              DO j=JstrV,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
!>              tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                    &
!>   &                          tl_Awrk(j,k  ,Nold)-                    &
!>   &                          FC(j,k)*tl_DC(j,k-1))
!>
                adfac=cff*ad_DC(j,k)
                ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k)*adfac
                ad_DC(j,k)=0.0_r8
              END DO
            END DO
            DO j=JstrV,Jend
!>            tl_DC(j,0)=0.0_r8
!>
              ad_DC(j,0)=0.0_r8
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!>              tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                    &
!>   &                          tl_Awrk(i,k  ,Nold)-                    &
!>   &                          FC(i,k)*tl_DC(i,k-1))
!>
                adfac=cff*ad_DC(i,k)
                ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
                ad_DC(i,k)=0.0_r8
              END DO
            END DO
            DO i=Istr,Iend
!>            tl_DC(i,0)=0.0_r8
!>
              ad_DC(i,0)=0.0_r8
            END DO
          END IF
        END IF
      END DO

#   else
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation implicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Compute diagonal matrix coefficients BC and CF.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                BC(j,k)=cff-FC(j,k)-FC(j,k-1)
              END DO
            END DO
            DO j=JstrV,Jend
              cff=1.0_r8/BC(j,1)
              CF(j,1)=cff*FC(j,1)
            END DO
            DO k=2,N(ng)-1
              DO j=JstrV,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                CF(j,k)=cff*FC(j,k)
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                BC(i,k)=cff-FC(i,k)-FC(i,k-1)
              END DO
            END DO
            DO i=Istr,Iend
              cff=1.0_r8/BC(i,1)
              CF(i,1)=cff*FC(i,1)
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                CF(i,k)=cff*FC(i,k)
              END DO
            END DO
          END IF
!
!  Compute new solution by back substitution.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
!>          DO k=N(ng)-1,1,-1
!>
            DO k=1,N(ng)-1
              DO j=JstrV,Jend
#    ifdef MASKING
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*vmask(i,j)
!>
                ad_Awrk(j,k,Nnew)=ad_Awrk(j,k,Nnew)*vmask(i,j)
#    endif
!>              tl_Awrk(j,k,Nnew)=tl_DC(j,k)
!>
                ad_DC(j,k)=ad_DC(j,k)+                                  &
     &                     ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
!>              tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!>
                ad_DC(j,k+1)=-CF(j,k)*ad_DC(j,k)
              END DO
            END DO
            DO j=JstrV,Jend
#    ifdef MASKING
!>            tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*vmask(i,j)
!>
              ad_Awrk(j,N(ng),Nnew)=ad_Awrk(j,N(ng),Nnew)*vmask(i,j)
#    endif
!>            tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
              ad_DC(j,N(ng))=ad_DC(j,N(ng))+ad_Awrk(j,N(ng),Nnew)
              ad_Awrk(j,N(ng),Nnew)=0.0_r8

!>            tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                           &
!>   &                        FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/          &
!>   &                       (BC(j,N(ng))-                              &
!>   &                        FC(j,N(ng)-1)*CF(j,N(ng)-1))
!>
              adfac=ad_DC(j,N(ng))/                                     &
     &                       (BC(j,N(ng))-                              &
     &                        FC(j,N(ng)-1)*CF(j,N(ng)-1))
              ad_DC(j,N(ng)-1)=ad_DC(j,N(ng)-1)-FC(j,N(ng)-1)*adfac
              ad_DC(j,N(ng))=adfac
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
!>          DO k=N(ng)-1,1,-1
!>
            DO k=1,N(ng)-1
              DO i=Istr,Iend
#    ifdef MASKING
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*vmask(i,j)
!>
                ad_Awrk(i,k,Nnew)=ad_Awrk(i,k,Nnew)*vmask(i,j)
#    endif
!>              tl_Awrk(i,k,Nnew)=tl_DC(i,k)
!>
                ad_DC(i,k)=ad_DC(i,k)+                                  &
     &                     ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
!>              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
                ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
              END DO
            END DO
            DO i=Istr,Iend
#    ifdef MASKING
!>            tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*vmask(i,j)
!>
              ad_Awrk(i,N(ng),Nnew)=ad_Awrk(i,N(ng),Nnew)*vmask(i,j)
#    endif
!>            tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
!>
              ad_DC(i,N(ng))=ad_DC(i,N(ng))+                            &
     &                       ad_Awrk(i,N(ng),Nnew)
              ad_Awrk(i,N(ng),Nnew)=0.0_r8
!>            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                           &
!>   &                        FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!>   &                       (BC(i,N(ng))-                              &
!>   &                        FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>
              adfac=ad_DC(i,N(ng))/                                     &
     &              (BC(i,N(ng))-                                       &
     &               FC(i,N(ng)-1)*CF(i,N(ng)-1))
              ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-                        &
     &                         FC(i,N(ng)-1)*adfac
              ad_DC(i,N(ng)  )=adfac
            END DO
          END IF
!
!  Solve the adjoint tridiagonal system.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
!>          DO k=2,N(ng)-1
!>
            DO k=N(ng)-1,2,-1
              DO j=JstrV,Jend
                cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
!>              tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
!>
                adfac=cff*ad_DC(j,k)
                ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k-1)*adfac
                ad_DC(j,k  )=adfac
              END DO
            END DO
            DO j=JstrV,Jend
              cff=1.0_r8/BC(j,1)
!>            tl_DC(j,1)=cff*tl_DC(j,1)
!>
              ad_DC(j,1)=cff*ad_DC(j,1)
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
!>          DO k=2,N(ng)-1
!>
            DO k=N(ng)-1,2,-1
              DO i=Istr,Iend
                cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!>              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
!>
                adfac=cff*ad_DC(i,k)
                ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
                ad_DC(i,k  )=adfac
              END DO
            END DO
            DO i=Istr,Iend
              cff=1.0_r8/BC(i,1)
!>            tl_DC(i,1)=cff*tl_DC(i,1)
!>
              ad_DC(i,1)=cff*ad_DC(i,1)
            END DO
          END IF
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
!>              tl_DC(j,k)=tl_Awrk(j,k,Nold)*cff
!>
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            cff*ad_DC(j,k)
                ad_DC(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
!>              tl_DC(i,k)=tl_Awrk(i,k,Nold)*cff
!>
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            cff*ad_DC(i,k)
                ad_DC(i,k)=0.0
              END DO
            END DO
          END IF
        END IF
      END DO
#   endif

#  else
!
!-----------------------------------------------------------------------
!  Integrate adjoint vertical diffusion equation explicitly.
!-----------------------------------------------------------------------
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Time-step adjoint vertical diffusive term. Notice that "oHz" is
!  assumed to be time invariant.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            oHz(j,k)*(tl_FS(j,k  )-               &
!>   &                                      tl_FS(j,k-1))
!>
                adfac=oHz(j,k)*ad_Awrk(j,k,Nnew)
                ad_FS(j,k-1)=ad_FS(j,k-1)-adfac
                ad_FS(j,k  )=ad_FS(j,k  )+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            oHz(i,k)*(tl_FS(i,k  )-               &
!>   &                                      tl_FS(i,k-1))
!>
                adfac=oHz(i,k)*ad_Awrk(i,k,Nnew)
                ad_FS(i,k-1)=ad_FS(i,k-1)-adfac
                ad_FS(i,k  )=ad_FS(i,k  )+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute adjoint vertical diffusive flux.  Notice that "FC" is
!  assumed to be time invariant.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO j=JstrV,Jend
!>            tl_FS(j,N(ng))=0.0_r8
!>
              ad_FS(j,N(ng))=0.0_r8
!>            tl_FS(j,0)=0.0_r8
!>
              ad_FS(j,0)=0.0_r8
              DO k=1,N(ng)-1
#   ifdef MASKING
!>              tl_FS(j,k)=tl_FS(j,k)*vmask(i,j)
!>
                ad_FS(j,k)=ad_FS(j,k)*vmask(i,j)
#   endif
!>              tl_FS(j,k)=FC(j,k)*(tl_Awrk(j,k+1,Nold)-                &
!>   &                              tl_Awrk(j,k  ,Nold))
!>
                adfac=FC(j,k)*ad_FS(j,k)
                ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                ad_FS(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO i=Istr,Iend
!>            tl_FS(i,N(ng))=0.0_r8
!>
              ad_FS(i,N(ng))=0.0_r8
!>            tl_FS(i,0)=0.0_r8
!>
              ad_FS(i,0)=0.0_r8
              DO k=1,N(ng)-1
#   ifdef MASKING
!>              tl_FS(i,k)=tl_FS(i,k)*vmask(i,j)
!>
                ad_FS(i,k)=ad_FS(i,k)*vmask(i,j)
#   endif
!>              tl_FS(i,k)=FC(i,k)*(tl_Awrk(i,k+1,Nold)-                &
!>   &                              tl_Awrk(i,k  ,Nold))
!>
                adfac=FC(i,k)*ad_FS(i,k)
                ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                ad_FS(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Integrate adjoint horizontal diffusion equation.
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
!  Apply boundary conditions.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
!>   &                          LBij, UBij, 1, N(ng),                   &
!>   &                          Nghost,                                 &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Awrk(:,:,Nnew))
!>
        CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,        &
     &                             LBij, UBij, 1, N(ng),                &
     &                             Nghost,                              &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Awrk(:,:,Nnew))
# endif
!>      CALL bc_v3d_bry_tile (ng, tile, boundary,                       &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        tl_Awrk(:,:,Nnew))
!>
        CALL ad_bc_v3d_bry_tile (ng, tile, boundary,                    &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           ad_Awrk(:,:,Nnew))
!
!  Time-step adjoint horizontal diffusion equation.
!
        IF (Lconvolve(boundary)) THEN
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
!>              tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                    &
!>   &                            Hfac(j)*                              &
!>   &                            (tl_FE(j,k)-tl_FE(j-1,k))
!>
                adfac=Hfac(j)*ad_Awrk(j,k,Nnew)
                ad_FE(j-1,k)=ad_FE(j-1,k)-adfac
                ad_FE(j  ,k)=ad_FE(j  ,k)+adfac
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_Awrk(j,k,Nnew)
                ad_Awrk(j,k,Nnew)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                    &
!>   &                            Hfac(i)*                              &
!>   &                            (tl_FX(i+1,k)-tl_FX(i,k))
!>
                adfac=Hfac(i)*ad_Awrk(i,k,Nnew)
                ad_FX(i  ,k)=ad_FX(i  ,k)-adfac
                ad_FX(i+1,k)=ad_FX(i+1,k)+adfac
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_Awrk(i,k,Nnew)
                ad_Awrk(i,k,Nnew)=0.0_r8
              END DO
            END DO
          END IF
!
!  Compute XI- and ETA-components of adjoint diffusive flux.
!
          IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
            i=edge(boundary)
            DO k=1,N(ng)
              DO j=JstrV-1,Jend
!>              tl_FE(j,k)=pnom_r(i,j)*Kh(i,j)*                         &
!>   &                     (tl_Awrk(j+1,k,Nold)-                        &
!>   &                      tl_Awrk(j  ,k,Nold))
!>
                adfac=pnom_r(i,j)*Kh(i,j)*ad_FE(j,k)
                ad_Awrk(j  ,k,Nold)=ad_Awrk(j  ,k,Nold)-adfac
                ad_Awrk(j+1,k,Nold)=ad_Awrk(j+1,k,Nold)+adfac
                ad_FE(j,k)=0.0_r8
              END DO
            END DO
          ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
            j=edge(boundary)
            DO k=1,N(ng)
              DO i=Istr,Iend+1
# ifdef MASKING
!>              tl_FX(i,k)=tl_FX(i,k)*pmask(i,j)
!>
                ad_FX(i,k)=ad_FX(i,k)*pmask(i,j)
# endif
!>              tl_FX(i,k)=pmon_p(i,j)*                                 &
!>   &                     0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+              &
!>   &                              Kh(i-1,j-1)+Kh(i,j-1))*             &
!>   &                     (tl_Awrk(i  ,k,Nold)-                        &
!>   &                      tl_Awrk(i-1,k,Nold))
!>
                adfac=pmon_p(i,j)*                                      &
     &                0.25_r8*(Kh(i-1,j  )+Kh(i,j  )+                   &
     &                         Kh(i-1,j-1)+Kh(i,j-1))*ad_FX(i,k)
                ad_Awrk(i-1,k,Nold)=ad_Awrk(i-1,k,Nold)-adfac
                ad_Awrk(i  ,k,Nold)=ad_Awrk(i  ,k,Nold)+adfac
                ad_FX(i,k)=0.0_r8
              END DO
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set adjoint initial conditions.
!-----------------------------------------------------------------------
!
      IF (Lconvolve(boundary)) THEN
        IF ((boundary.eq.iwest).or.(boundary.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend+1
!>            tl_Awrk(j,k,Nold)=tl_A(j,k)
!>
              ad_A(j,k)=ad_A(j,k)+ad_Awrk(j,k,Nold)
              ad_Awrk(j,k,Nold)=0.0_r8
            END DO
          END DO
        ELSE IF ((boundary.eq.isouth).or.(boundary.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Istr-1,Iend+1
!>            tl_Awrk(i,k,Nold)=tl_A(i,k)
!>
              ad_A(i,k)=ad_A(i,k)+ad_Awrk(i,k,Nold)
              ad_Awrk(i,k,Nold)=0.0_r8
            END DO
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,             &
!>   &                        LBij, UBij, 1, N(ng),                     &
!>   &                        Nghost,                                   &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_A)
!>
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, boundary,          &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           Nghost,                                &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!>    CALL bc_v3d_bry_tile (ng, tile, boundary,                         &
!>   &                      LBij, UBij, 1, N(ng),                       &
!>   &                      tl_A)
!>
      CALL ad_bc_v3d_bry_tile (ng, tile, boundary,                      &
     &                         LBij, UBij, 1, N(ng),                    &
     &                         ad_A)

      RETURN
      END SUBROUTINE ad_conv_v3d_bry_tile

#endif
      END MODULE ad_conv_bry3d_mod
