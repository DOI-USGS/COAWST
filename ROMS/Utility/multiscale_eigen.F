#include "cppdefs.h"
      MODULE multiscale_eigen_mod

#if defined FOUR_DVAR && defined MULTI_SCALE_B
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group           Andrew M. Moore    !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Multi-scale, Background-Error Covariace modeling:                   !
!                                                                      !
!  This module computes the extrema Ritz eigenvalues required by the   !
!  Chebyshev Iterations (CI) solver to invert the implicit horizontal  !
!  pseudo-diffusion operator, which models the spreading and smoothing !
!  effects of the multiscale background-error covariance matrix B.     !
!                                                                      !
!  For each variable in the control vector and for predetermined       !
!  correlation scales, the minimum and maximum eigenvalues are         !
!  estimated by iteratively running the Conjugate Gradient (CG)        !
!  solver and are stored in the "multiscale" CLASS structure.          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation  !
!      operators based on an implicitly formulated diffusion equation  !
!      solved with the Chebyshev iteration, Q.J.R. Meteorol. Soc.,     !
!      142, 455-471, doi:10.1002/qj.2664.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_fourdvar
      USE mod_grid
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
      USE mod_scalars
# if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
# endif
!
# ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod,        ONLY : bc_r2d_bry_tile,                  &
     &                                bc_u2d_bry_tile,                  &
     &                                bc_v2d_bry_tile
#  ifdef SOLVE3D
      USE bc_bry3d_mod,        ONLY : bc_r3d_bry_tile,                  &
     &                                bc_u3d_bry_tile,                  &
     &                                bc_v3d_bry_tile
#  endif
# endif
      USE bc_2d_mod,           ONLY : dabc_r2d_tile,                    &
     &                                dabc_u2d_tile,                    &
     &                                dabc_v2d_tile
# ifdef SOLVE3D
      USE bc_3d_mod,           ONLY : dabc_r3d_tile,                    &
     &                                dabc_u3d_tile,                    &
     &                                dabc_v3d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod,     ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod,     ONLY : mp_exchange3d
      USE mp_exchange_mod,     ONLY : mp_exchange4d
#  endif
# endif
      USE set_depth_mod,       ONLY : set_depth_tile       
      USE strings_mod,         ONLY : FoundError
      USE white_noise_mod,     ONLY : white_noise2d
# ifdef ADJUST_BOUNDARY
      USE white_noise_mod,     ONLY : white_noise2d_bry
# endif
# ifdef SOLVE3D
      USE white_noise_mod,     ONLY : white_noise3d
#  ifdef ADJUST_BOUNDARY
      USE white_noise_mod,     ONLY : white_noise3d_bry
#  endif
# endif
!
      USE roms_multiscale_mod, ONLY : multiscale         ! CLASS object
!
      implicit none
!
      PUBLIC  :: multiscale_eigen
      PUBLIC  :: multiscale_eigen_read
      PUBLIC  :: multiscale_eigen_write
      PRIVATE :: multiscale_eigen_tile
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE multiscale_eigen (ng, tile, Linp, ifac)
!***********************************************************************
!
      USE mod_stepping,        ONLY : nnew, nstp
      USE roms_multiscale_mod, ONLY : MSB
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, ifac
!
!  Local variable declarations.
!
      logical :: Lweak
!
# include "tile.h"
!
!   Set Lweak to FALSE here for now. Need to figure out the general
!   case for weak constrain model error (Lweak=TRUE).
!
      Lweak=.FALSE.
!
      CALL multiscale_eigen_tile (MSB(ng), ng, tile, iTLM,              &
     &                            ifac, Lweak,                          &
     &                            nstp(ng), nnew(ng), Linp,             &
     &                            LBi, UBi, LBj, UBj,                   &
# ifdef ADJUST_BOUNDARY
     &                            LBij, UBij,                           &
# endif
     &                            IminS, ImaxS, JminS, JmaxS,           &
# ifdef ADJUST_BOUNDARY
     &                            BOUNDARY(ng) % tl_ubar_obc,           &
     &                            BOUNDARY(ng) % tl_vbar_obc,           &
     &                            BOUNDARY(ng) % tl_zeta_obc,           &
#  ifdef SOLVE3D
     &                            BOUNDARY(ng) % tl_t_obc,              &
     &                            BOUNDARY(ng) % tl_u_obc,              &
     &                            BOUNDARY(ng) % tl_v_obc,              &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                            FORCES(ng) % tl_ustr,                 &
     &                            FORCES(ng) % tl_vstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                            FORCES(ng) % tl_tflux,                &
# endif
# ifdef SOLVE3D
     &                            OCEAN(ng) % tl_t,                     &
     &                            OCEAN(ng) % tl_u,                     &
     &                            OCEAN(ng) % tl_v,                     &
# endif
     &                            OCEAN(ng) % tl_ubar,                  &
     &                            OCEAN(ng) % tl_vbar,                  &
     &                            OCEAN(ng) % tl_zeta)
!
      RETURN
      END SUBROUTINE multiscale_eigen
!
!***********************************************************************
      SUBROUTINE multiscale_eigen_tile (self, ng, tile, model,          &
     &                                  ifac, Lweak,                    &
     &                                  nstp, nnew, Linp,               &
     &                                  LBi, UBi, LBj, UBj,             &
# ifdef ADJUST_BOUNDARY
     &                                  LBij, UBij,                     &
# endif
     &                                  IminS, ImaxS, JminS, JmaxS,     &
# ifdef ADJUST_BOUNDARY
     &                                  tl_ubar_obc, tl_vbar_obc,       &
     &                                  tl_zeta_obc,                    &
#  ifdef SOLVE3D
     &                                  tl_t_obc, tl_u_obc, tl_v_obc,   &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                                  tl_ustr, tl_vstr,               &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                  tl_tflux,                       &
# endif
# ifdef SOLVE3D
     &                                  tl_t, tl_u, tl_v,               &
# endif
     &                                  tl_ubar, tl_vbar,               &
     &                                  tl_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      logical, intent(in) :: Lweak
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: nstp, nnew, Linp, ifac
      integer, intent(in) :: LBi, UBi, LBj, UBj
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  endif

# else

#  ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
# endif
!
!  Local variable declarations.
!
      integer, parameter :: ir = 1
!
      integer :: i, ib,  j, k, ns
# ifdef SOLVE3D
      integer :: ifield, itrc
# endif
      real(r8) :: Amax, Amin
# ifdef ADJUST_BOUNDARY
      real(r8) :: Bmin, Bmax
      real(r8), parameter :: Aspv = 0.0_r8
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: A2d
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
# endif
# ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij)               :: B2d
#  ifdef SOLVE3D
      real(r8), dimension(LBij:UBij,1:N(ng))       :: B3d
#  endif
# endif

# include "set_bounds.h"

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8                   ! free surface
        END DO
      END DO
!
      CALL set_depth_tile (ng, tile, iTLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     GRID(ng) % h,                                &
#  ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     SEDBED(ng) % bed_thick,                      &
#  endif
     &                     A2d,                                         &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
# endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Initialize the right-hand-size of the control vector for the
!  Conjugate Gradient (CG) solver with random values.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  
      IF (Master) WRITE (stdout,10)
!
!  Tangent linear free-surface.
!
      CALL white_noise2d (ng, model, r2dvar, Rscheme(ng),               &
     &                    IstrR, IendR, JstrR, JendR,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          tl_zeta(i,j,Linp)=A2d(i,j)
        END DO
      END DO
!
      CALL dabc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_zeta(:,:,Linp))
!
!  Tangent linear 2D momentum.
!
      CALL white_noise2d (ng, model, u2dvar, Rscheme(ng),               &
     &                    Istr, IendR, JstrR, JendR,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          tl_ubar(i,j,Linp)=A2d(i,j)
        END DO
      END DO
!
      CALL white_noise2d (ng, model, v2dvar, Rscheme(ng),               &
     &                    IstrR, IendR, Jstr, JendR,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          tl_vbar(i,j,Linp)=A2d(i,j)
        END DO
      END DO
!
      CALL dabc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_ubar(:,:,Linp))
      CALL dabc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_vbar(:,:,Linp))

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,Linp),                            &
     &                    tl_ubar(:,:,Linp),                            &
     &                    tl_vbar(:,:,Linp))
# endif

# ifdef SOLVE3D
!
!  Tangent linear 3D momentum.
!
      CALL white_noise3d (ng, model, u3dvar, Rscheme(ng),               &
     &                    Istr, IendR, JstrR, JendR,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Amin, Amax, A3d)
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          DO k=1,N(ng)
            tl_u(i,j,k,Linp)=A3d(i,j,k)
          END DO
        END DO
      END DO
!
      CALL white_noise3d (ng, model, v3dvar, Rscheme(ng),               &
     &                    IstrR, IendR, Jstr, JendR,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Amin, Amax, A3d)
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          DO k=1,N(ng)
            tl_v(i,j,k,Linp)=A3d(i,j,k)
          END DO
        END DO
      END DO
!
      CALL dabc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    tl_u(:,:,:,Linp))
      CALL dabc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    tl_v(:,:,:,Linp))

#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,Linp),                             &
     &                    tl_v(:,:,:,Linp))
#  endif
!
!  Tangent linear tracer variables.
!
      CALL white_noise3d (ng, model, r3dvar, Rscheme(ng),               &
     &                    IstrR, IendR, JstrR, JendR,                   &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Amin, Amax, A3d)
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              tl_t(i,j,k,Linp,itrc)=A3d(i,j,k)
            END DO
          END DO
        END DO
!
        CALL dabc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      tl_t(:,:,:,Linp,itrc))
      END DO

#  ifdef DISTRIBUTE
!
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,Linp,:))
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!  Tangent linear free-surface open boundaries.
!
      tl_zeta_obc(:,:,ir,Linp)=Aspv
      DO ib=1,4
        SELECT CASE (ib)
          CASE (iwest, ieast)
            CALL white_noise2d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              JstrT, JendT,                       &
     &                              LBij, UBij,                         &
     &                              Bmin, Bmax, B2d)
            DO j=JstrT,JendT
              tl_zeta_obc(j,ib,ir,Linp)=B2d(j)
            END DO
          CASE (isouth, inorth)
            CALL white_noise2d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              IstrT, IendT,                       &
     &                              LBij, UBij,                         &
     &                              Bmin, Bmax, B2d)
            DO i=IstrT,IendT
              tl_zeta_obc(i,ib,ir,Linp)=B2d(i)
            END DO
        END SELECT
        CALL bc_r2d_bry_tile (ng, tile, ib,                             &
     &                        LBij, UBij,                               &
     &                        tl_zeta_obc(:,ib,ir,Linp))
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      tl_ubar_obc(:,:,ir,Linp)=Aspv
      DO ib=1,4
        SELECT CASE (ib)
          CASE (iwest, ieast)
            CALL white_noise2d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              JstrT, JendT,                       &
     &                              LBij, UBij,                         &
     &                              Bmin, Bmax, B2d)
            DO j=JstrT,JendT
              tl_ubar_obc(j,ib,ir,Linp)=B2d(j)
            END DO
          CASE (isouth, inorth)
            CALL white_noise2d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              IstrP, IendT,                       &
     &                              LBij, UBij,                         &
     &                              Bmin, Bmax, B2d)
            DO i=IstrP,IendT
              tl_ubar_obc(i,ib,ir,Linp)=B2d(i)
            END DO
        END SELECT
        CALL bc_u2d_bry_tile (ng, tile, ib,                             &
     &                        LBij, UBij,                               &
     &                        tl_ubar_obc(:,ib,ir,Linp))
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      tl_vbar_obc(:,:,ir,Linp)=Aspv
      DO ib=1,4
        SELECT CASE (ib)
          CASE (iwest, ieast)
            CALL white_noise2d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              JstrP, JendT,                       &
     &                              LBij, UBij,                         &
     &                              Bmin, Bmax, B2d)
            DO j=JstrP,JendT
              tl_vbar_obc(j,ib,ir,Linp)=B2d(j)
            END DO
          CASE (isouth, inorth)
            CALL white_noise2d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              IstrT, IendT,                       &
     &                              LBij, UBij,                         &
     &                              Bmin, Bmax, B2d)
            DO i=IstrT,IendT
              tl_vbar_obc(i,ib,ir,Linp)=B2d(i)
            END DO
        END SELECT
        CALL bc_v2d_bry_tile (ng, tile, ib,                             &
     &                        LBij, UBij,                               &
     &                        tl_vbar_obc(:,ib,ir,Linp))
      END DO

#  ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      tl_u_obc(:,:,:,ir,Linp)=Aspv
      DO ib=1,4
        SELECT CASE (ib)
          CASE (iwest, ieast)
            CALL white_noise3d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              JstrT, JendT,                       &
     &                              LBij, UBij, 1, N(ng),               &
     &                              Bmin, Bmax, B3d)
            DO k=1,N(ng)
              DO j=JstrT,JendT
                tl_u_obc(j,k,ib,ir,Linp)=B3d(j,k)
              END DO
            END DO
          CASE (isouth, inorth)
            CALL white_noise3d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              IstrP, IendT,                       &
     &                              LBij, UBij, 1, N(ng),               &
     &                              Bmin, Bmax, B3d)
            DO k=1,N(ng)
              DO i=IstrP,IendT
                tl_u_obc(i,k,ib,ir,Linp)=B3d(i,k)
              END DO
            END DO
        END SELECT
        CALL bc_u3d_bry_tile (ng, tile, ib,                             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        tl_u_obc(:,:,ib,ir,Linp))
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      tl_v_obc(:,:,:,ir,Linp)=Aspv
      DO ib=1,4
        SELECT CASE (ib)
          CASE (iwest, ieast)
            CALL white_noise3d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              JstrP, JendT,                       &
     &                              LBij, UBij, 1, N(ng),               &
     &                              Bmin, Bmax, B3d)
            DO k=1,N(ng)
              DO j=JstrP,JendT
                tl_v_obc(j,k,ib,ir,Linp)=B3d(j,k)
              END DO
            END DO
          CASE (isouth, inorth)
            CALL white_noise3d_bry (ng, tile, model, ib,                &
     &                              Rscheme(ng),                        &
     &                              IstrT, IendT,                       &
     &                              LBij, UBij, 1, N(ng),               &
     &                              Bmin, Bmax, B3d)
            DO k=1,N(ng)
              DO i=IstrT,IendT
                tl_v_obc(i,k,ib,ir,Linp)=B3d(i,k)
              END DO
            END DO
        END SELECT
        CALL bc_v3d_bry_tile (ng, tile, ib,                             &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        tl_v_obc(:,:,ib,ir,Linp))
      END DO
!
!  Tangent linear tracer variables open boundaries.
!
      DO itrc=1,NT(ng)
        tl_t_obc(:,:,:,ir,Linp,itrc)=Aspv
        DO ib=1,4
          SELECT CASE (ib)
            CASE (iwest, ieast)
              CALL white_noise3d_bry (ng, tile, model, ib,              &
     &                                Rscheme(ng),                      &
     &                                JstrT, JendT,                     &
     &                                LBij, UBij, 1, N(ng),             &
     &                                Bmin, Bmax, B3d)
              DO k=1,N(ng)
                DO j=JstrT,JendT
                  tl_t_obc(j,k,ib,ir,Linp,itrc)=B3d(j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              CALL white_noise3d_bry (ng, tile, model, ib,              &
     &                                Rscheme(ng),                      &
     &                                IstrR, IendR,                     &
     &                                LBij, UBij, 1, N(ng),             &
     &                                Bmin, Bmax, B3d)
              DO k=1,N(ng)
                DO i=IstrT,IendT
                  tl_t_obc(i,k,ib,ir,Linp,itrc)=B3d(i,k)
                END DO
              END DO
          END SELECT
          CALL bc_r3d_bry_tile (ng, tile, ib,                           &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          tl_t_obc(:,:,ib,ir,Linp,itrc))
        END DO
      END DO
#  endif
# endif

# ifdef ADJUST_WSTRESS
!
!  Tangent linear surface momentum stress.
!
      CALL white_noise2d (ng, model, u2dvar, Rscheme(ng),               &
     &                    IstrP, IendT, JstrT, JendT,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          tl_ustr(i,j,1,Linp)=A2d(i,j)
        END DO
      END DO
!
      CALL white_noise2d (ng, model, v2dvar, Rscheme(ng),               &
     &                    IstrT, IendT, JstrP, JendT,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          tl_vstr(i,j,1,Linp)=A2d(i,j)
        END DO
      END DO
!
      CALL dabc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_ustr(:,:,ir,Linp))
      CALL dabc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_vstr(:,:,ir,Linp))

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ustr(:,:,ir,Linp),                         &
     &                    tl_vstr(:,:,ir,Linp))
#  endif
# endif

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Tangent linear surface tracer variables flux.
!
      DO itrc=1,NT(ng)
        CALL white_noise2d (ng, model, r2dvar, Rscheme(ng),             &
     &                      IstrT, IendT, JstrT, JendT,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Amin, Amax, A2d)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            tl_tflux(i,j,ir,Linp,itrc)=A2d(i,j)
          END DO
        END DO
!
        CALL dabc_r2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      tl_tflux(:,:,ir,Linp,itrc))

#  ifdef DISTRIBUTE
!
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_tflux(:,:,ir,Linp,itrc))
#  endif
      END DO
# endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Prior/model multiscale error convariance: Convolve tangent linear
!  state vector with a generalized implicit diffusion operators to
!  compute the minimum and maximum Ritz eigenvalues required by the
!  Chebyshev Iterations (CI) solver.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Loop over the different length scales.
!
      MS_LOOP : DO ns=1,Nscale(ng)
!
!  Tangent linear free-surface, inplicit CG convolution.
!
        IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isFsur))
!
        CALL self%tl_CG_2d (ng, tile, model, isFsur, r2dvar,            &
     &                      ns, NiterCG(ns,ng), ifac, Lweak,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      tl_zeta(:,:,Linp))
!
!  Tangent linear 2D U-momentum, implicit CG convolution.
!
        IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isUbar))
!
        CALL self%tl_CG_2d (ng, tile, model, isUbar, u2dvar,            &
     &                      ns, NiterCG(ns,ng), ifac, Lweak,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      tl_ubar(:,:,Linp))
!
!  Tangent linear 2D V-momentum, implicit CG convolution.
!
        IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isVbar))
!
        CALL self%tl_CG_2d (ng, tile, model, isVbar, v2dvar,            &
     &                      ns, NiterCG(ns,ng), ifac, Lweak,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      tl_vbar(:,:,Linp))

# ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum, implicit CG convolution.
!
        IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isUvel))
!
        CALL self%tl_CG_3d (ng, tile, model, isUvel, u3dvar,            &
     &                      ns, NiterCG(ns,ng), ifac, Lweak,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      tl_u(:,:,:,Linp))
!
!  Tangent linear 3D V-momentum, implicit CG convolution.
!
        IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isVvel))
!
        CALL self%tl_CG_3d (ng, tile, model, isVvel, v3dvar,            &
     &                      ns, NiterCG(ns,ng), ifac, Lweak,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      tl_v(:,:,:,Linp))
!
!  Tangent linear tracer variables, implicit CG convolution.
!
        DO itrc=1,NT(ng)
          ifield=isTvar(itrc)
          IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(ifield))
!
          CALL self%tl_CG_3d (ng, tile, model, ifield, r3dvar,          &
     &                        ns, NiterCG(ns,ng), ifac, Lweak,          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        tl_t(:,:,:,Linp,itrc))
        END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  Tangent linear free-surface open boundaries, implicit GC onvolution.
!
        IF (Master) WRITE (stdout,30) ns, TRIM(StateVarName(isFsur))
!
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
            CALL self%tl_CG_b1d (ng, tile, model, isFsur, ib,           &
     &                           r2dvar, ns, NiterCG(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_zeta_obc(:,ib,ir,Linp))
          END IF
        END DO
!
!  Tangent linear 2D U-momentum open boundaries, implicit CG convolution.
!
        IF (Master) WRITE (stdout,30) ns, TRIM(StateVarName(isUbar))
!
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
            CALL self%tl_CG_b1d (ng, tile, model, isUbar, ib,           &
     &                           u2dvar, ns, NiterCG(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_ubar_obc(:,ib,ir,Linp))
          END IF
        END DO
!
!  Tangent linear 2D V-momentum open boundaries, implicit CG convolution.
!
        IF (Master) WRITE (stdout,30) ns, TRIM(StateVarName(isVbar))
!
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
            CALL self%tl_CG_b1d (ng, tile, model, isVbar, ib,           &
     &                           v2dvar, ns, NiterCG(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_vbar_obc(:,ib,ir,Linp))
          END IF
        END DO

#  ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries, implicit CG convolution.
!
        IF (Master) WRITE (stdout,30) ns, TRIM(StateVarName(isUvel))
!
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
            CALL self%tl_CG_b2d (ng, tile, model, isUvel, ib,           &
     &                           u3dvar, ns, NiterCG(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_u_obc(:,:,ib,1,Linp))
          END IF
        END DO
!
!  Tangent linear 3D V-momentum open boundaries, implicit CG convolution.
!
        IF (Master) WRITE (stdout,30) ns, TRIM(StateVarName(isVvel))
!
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
            CALL self%tl_CG_b2d (ng, tile, model, isVvel, ib,           &
     &                           v3dvar, ns, NiterCG(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_v_obc(:,:,ib,1,Linp))
          END IF
        END DO
!
!  Tangent linear tracers open boundaries, implicit CG convolution.
!
        DO itrc=1,NT(ng)
          ifield=isTvar(itrc)
        IF (Master) WRITE (stdout,30) ns, TRIM(StateVarName(ifield))
!
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,ifield,ng)) THEN
              CALL self%tl_CG_b2d (ng, tile, model, ifield, ib,         &
     &                             r3dvar, ns, NiterCG(ns,ng), ifac,    &
     &                             LBij, UBij,                          &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             tl_t_obc(:,:,ib,ir,Linp,itrc))

             END IF
          END DO
        END DO
#  endif
# endif

# ifdef ADJUST_WSTRESS
!
!  Tangent linear surface U-momentum stress, implicit CG convolution.
!
        IF (.not.Lweak) THEN
          IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isUstr))
!
          CALL self%tl_CG_2d (ng, tile, model, isUstr, u2dvar,          &
     &                        ns, NiterCG(ns,ng), ifac, Lweak,          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        tl_ustr(:,:,ir,Linp))
        END IF
!
!  Tangent linear surface V-momentum stress, implicit CG convolution.
!
        IF (.not.Lweak) THEN
          IF (Master) WRITE (stdout,20) ns, TRIM(StateVarName(isVstr))
!
          CALL self%tl_CG_2d (ng, tile, model, isVstr, v2dvar,          &
     &                        ns, NiterCG(ns,ng), ifac, Lweak,          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        tl_vstr(:,:,ir,Linp))
        END IF
# endif

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Tangent linear surface tracers flux, implicit CG convolution.
!
        IF (.not.Lweak) THEN
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              ifield=isTsur(itrc)
              IF (Master) WRITE (stdout,20) ns,                         &
     &                                      TRIM(StateVarName(ifield))
!
              CALL self%tl_CG_2d (ng, tile, model, ifield, r2dvar,      &
     &                            ns, NiterCG(ns,ng), ifac, Lweak,      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            tl_tflux(:,:,ir,Linp,itrc))
            END IF
          END DO
        END IF
# endif
      END DO MS_LOOP
!
  10  FORMAT (/,' Ritz Extrema Eigenvalues Required by Implicit ',      &
     &        'Chebyshev Iterations (CI) Solver:',/) 
  20  FORMAT (4x,'B-scale = ',i2.2,', computing eigenspectrum for: ',a)
# ifdef ADJUST_BOUNDARY
  30  FORMAT (4x,'B-scale = ',i2.2,', computing eigenspectrum for: ',a, &
     &        ' lateral boundaries')
# endif
!
      RETURN
      END SUBROUTINE multiscale_eigen_tile
!
!***********************************************************************
      SUBROUTINE multiscale_eigen_read (self, ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: i, ifld, itrc, nl, ns
      integer :: omode
!
      real (r8) :: Fmin, Fmax
!
      character (len=80 ) :: Text
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_eigen_read"
!
      SourceFile=MyFile
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Read in Extrema K-Laplacian eigenvalues for Prior/Model background-
!  error background error convariance into output Normalization NetCDF
!  file(s).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      nl=MAXVAL(self%Mlap)  ! nunber of k-Laplacian applications
      ns=self%Nscale        ! number of correlation scales
!
      FILE_LOOP : DO i=1,NSA

        LREAD_EIGEN : IF (.not.LwrtNRM(i,ng)) THEN
!
          ncname=NRM(i,ng)%name
          IF (i.eq.1) THEN
            Text='prior B-normalization file:'
          ELSE IF (i.eq.2) THEN
            Text='model Q-normalization file:'
          END IF
          IF (Master) WRITE (stdout,10) TRIM(Text), TRIM(ncname)
!
!  If appropriate, open input NetCDF file.
!
          IF (NRM(i,ng)%IOtype.eq.io_nf90) THEN
            IF (NRM(i,ng)%ncid.eq.-1) THEN
              omode=1
              CALL netcdf_open (ng, model, ncname, omode,               &
     &                          NRM(i,ng)%ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
# if defined PIO_LIB && defined DISTRIBUTE
          ELSE IF (NRM(i,ng)%IOtype.eq.io_pio) THEN
            IF (NRM(i,ng)%pioFile%fh.eq.-1) THEN
              omode=1
              CALL pio_netcdf_open (ng, model, ncname, omode,           &
     &                              NRM(i,ng)%pioFile)
            END IF
# endif
          END IF
!
!  Read in minimum and maximum Ritz eigenvalues.
!
          SELECT CASE (NRM(i,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieFsur),                    &
     &                              self%zeta_eigen,                    &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              start = (/1,1,1/),                  &
     &                              total = (/nl,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_eigen,                            &
     &                      MASK = self%zeta_eigen > 0.0)
                Fmax=MAXVAL(self%zeta_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieFsur)), Fmin, Fmax
              END IF
!
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieUbar),                    &
     &                              self%ubar_eigen,                    &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              start = (/1,1,1/),                  &
     &                              total = (/nl,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_eigen,                            &
     &                      MASK = self%ubar_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUbar)), Fmin, Fmax
              END IF
!
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieVbar),                    &
     &                              self%vbar_eigen,                    &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              start = (/1,1,1/),                  &
     &                              total = (/nl,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_eigen,                            &
     &                      MASK = self%vbar_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVbar)), Fmin, Fmax
              END IF

# ifdef SOLVE3D
!
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieUvel),                    &
     &                              self%u_eigen,                       &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,N(ng),ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_eigen,                               &
     &                      MASK = self%u_eigen > 0.0)
                Fmax=MAXVAL(self%u_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUvel)), Fmin, Fmax
              END IF
!
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieVvel),                    &
     &                              self%v_eigen,                       &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,N(ng),ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_eigen,                               &
     &                      MASK = self%v_eigen > 0.0)
                Fmax=MAXVAL(self%v_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVvel)), Fmin, Fmax
              END IF
!
              DO itrc=1,NT(ng)
                ifld=ieTvar(itrc)
                CALL netcdf_get_fvar (ng, model, ncname,                &
     &                                Vname(1,ifld),                    &
     &                                self%t_eigen(:,:,:,:,itrc),       &
     &                                ncid = NRM(i,ng)%ncid,            &
     &                                start = (/1,1,1,1/),              &
     &                                total = (/nl,N(ng),ns,2/))    
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_eigen(:,:,:,:,itrc),               &
     &                        MASK = self%t_eigen(:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_eigen(:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END DO
# endif

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieFsur),                &
     &                                  self%zeta_eigen,                &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                                  start = (/1,1,1/),              &
     &                                  total = (/nl,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_eigen,                            &
     &                      MASK = self%zeta_eigen > 0.0)
                Fmax=MAXVAL(self%zeta_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieFsur)), Fmin, Fmax
              END IF
!
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieUbar),                &
     &                                  self%ubar_eigen,                &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                                  start = (/1,1,1/),              &
     &                                  total = (/nl,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_eigen,                            &
     &                      MASK = self%ubar_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUbar)), Fmin, Fmax
              END IF
!
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieVbar),                &
     &                                  self%vbar_eigen,                &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                                  start = (/1,1,1/),              &
     &                                  total = (/nl,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_eigen,                            &
     &                      MASK = self%vbar_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVbar)), Fmin, Fmax
              END IF

#  ifdef SOLVE3D
!
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieUvel),                &
     &                                  self%u_eigen,                   &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,N(ng),ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_eigen,                               &
     &                      MASK = self%u_eigen > 0.0)
                Fmax=MAXVAL(self%u_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUvel)), Fmin, Fmax
              END IF
!
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieVvel),                &
     &                                  self%v_eigen,                   &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,N(ng),ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_eigen,                               &
     &                      MASK = self%v_eigen > 0.0)
                Fmax=MAXVAL(self%v_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVvel)), Fmin, Fmax
              END IF
!
              DO itrc=1,NT(ng)
                ifld=ieTvar(itrc)
                CALL pio_netcdf_get_fvar (ng, model, ncname,            &
     &                                    Vname(1,ifld),                &
     &                                    self%t_eigen(:,:,:,:,itrc),   &
     &                                    pioFile = NRM(i,ng)%pioFile,  &
     &                                    start = (/1,1,1,1/),          &
     &                                    total = (/nl,N(ng),ns,2/))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_eigen(:,:,:,:,itrc),               &
     &                        MASK = self%t_eigen(:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_eigen(:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END DO
#  endif
# endif
          END SELECT
!
        END IF LREAD_EIGEN
      END DO FILE_LOOP

# ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Read Extrema K-Laplacian eigenvalues for open boundaries background-
!  error background error convariance into output Normalization NetCDF
!  file.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      BRY_LREAD_EIGEN : IF (.not.LwrtNRM(3,ng)) THEN
!
        ncname=NRM(3,ng)%name
        Text='boundary B-normalization file:'
        IF (Master) WRITE (stdout,10) TRIM(Text), TRIM(ncname)
!
!  If appropriate, open input NetCDF file.
!
        IF (NRM(3,ng)%IOtype.eq.io_nf90) THEN
          IF (NRM(3,ng)%ncid.eq.-1) THEN
            omode=1
            CALL netcdf_open (ng, model, ncname, omode,                 &
     &                        NRM(3,ng)%ncid)
            IF (FoundError(exit_flag, NoError,                          &
     &                      __LINE__, MyFile)) RETURN
          END IF
#  if defined PIO_LIB && defined DISTRIBUTE
        ELSE IF (NRM(3,ng)%IOtype.eq.io_pio) THEN
          IF (NRM(3,ng)%pioFile%fh.eq.-1) THEN
            omode=1
            CALL pio_netcdf_open (ng, model, ncname, omode,             &
     &                            NRM(3,ng)%pioFile)
          END IF
#  endif
        END IF
!
!  Read in minimum and maximum Ritz eigenvalues.
!
        SELECT CASE (NRM(3,ng)%IOtype)
!
          CASE (io_nf90)
!
            IF (ANY(Lobc(:,isFsur,ng))) THEN
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieZbry),                    &
     &                              self%zeta_obc_eigen,                &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_obc_eigen,                        &
     &                      MASK = self%zeta_obc_eigen > 0.0)
                Fmax=MAxVAL(self%zeta_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieZbry)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isUbar,ng))) THEN
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieU2bc),                    &
     &                              self%ubar_obc_eigen,                &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_obc_eigen,                        &
     &                      MASK = self%ubar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU2bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVbar,ng))) THEN
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieV2bc),                    &
     &                              self%vbar_obc_eigen,                &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_obc_eigen,                        &
     &                      MASK = self%vbar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV2bc)), Fmin, Fmax
              END IF
            END IF

#  ifdef SOLVE3D
!
            IF (ANY(Lobc(:,isUvel,ng))) THEN
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieU3bc),                    &
     &                              self%u_obc_eigen,                   &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              start = (/1,1,1,1,1/),              &
     &                              total = (/nl,N(ng),4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_obc_eigen,                           &
     &                      MASK = self%u_obc_eigen > 0.0)
                Fmax=MAXVAL(self%u_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU3bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVvel,ng))) THEN
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieV3bc),                    &
     &                              self%v_obc_eigen,                   &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              start = (/1,1,1,1,1/),              &
     &                              total = (/nl,N(ng),4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_obc_eigen,                           &
     &                      MASK = self%v_obc_eigen > 0.0)
                Fmax=MAXVAL(self%v_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV3bc)), Fmin, Fmax
              END IF
            END IF
!
            DO itrc=1,NT(ng)
              IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
                ifld=ieTbry(itrc)
                CALL netcdf_get_fvar (ng, model, ncname,                &
     &                                Vname(1,ifld),                    &
     &                                self%t_obc_eigen(:,:,:,:,:,itrc), &
     &                                ncid = NRM(3,ng)%ncid,            &
     &                                start = (/1,1,1,1,1/),            &
     &                                total = (/nl,N(ng),4,ns,2/))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_obc_eigen(:,:,:,:,:,itrc),         &
     &                    MASK = self%t_obc_eigen(:,:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_obc_eigen(:,:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#  endif

#  if defined PIO_LIB && defined DISTRIBUTE
!
          CASE (io_pio)
!
            IF (ANY(Lobc(:,isFsur,ng))) THEN
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieZbry),                &
     &                                  self%zeta_obc_eigen,            &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_obc_eigen,                        &
     &                      MASK = self%zeta_obc_eigen > 0.0)
                Fmax=MAXVAL(self%zeta_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieZbry)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isUbar,ng))) THEN
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieU2bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_obc_eigen,                        &
     &                      MASK = self%ubar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU2bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVbar,ng))) THEN
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieV2bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_obc_eigen,                        &
     &                      MASK = self%vbar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV2bc)), Fmin, Fmax
              END IF
            END IF

#   ifdef SOLVE3D
!
            IF (ANY(Lobc(:,isUvel,ng))) THEN
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieU3bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                                  start = (/1,1,1,1,1/),          &
     &                                  total = (/nl,N(ng),4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_obc_eigen,                           &
     &                      MASK = self%u_obc_eigen > 0.0)
                Fmax=MAXVAL(self%u_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU3bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVvel,ng))) THEN
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieV3bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                                  start = (/1,1,1,1,1/),          &
     &                                  total = (/nl,N(ng),4,ns,2/))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_obc_eigen,                           &
     &                      MASK = self%v_obc_eigen > 0.0)
                Fmax=MAXVAL(self%v_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV3bc)), Fmin, Fmax
              END IF
            END IF
!
            DO itrc=1,NT(ng)
              IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
                ifld=ieTbry(itrc)
                CALL pio_netcdf_get_fvar (ng, model, ncname,            &
     &                                    Vname(1,ifld),                &
     &                                    self%t_obc_eigen(:,:,:,:,:,   &
     &                                                     itrc),       &
     &                                    pioFile = NRM(3,ng)%pioFile,  &
     &                                    start = (/1,1,1,1,1/),        &
     &                                    total = (/nl,N(ng),4,ns,2/))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_obc_eigen(:,:,:,:,:,itrc),         &
     &                    MASK = self%t_obc_eigen(:,:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_obc_eigen(:,:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#   endif
#  endif
        END SELECT
!
      END IF BRY_LREAD_EIGEN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Read Extrema K-Laplacian eigenvalues for open boundaries background-
!  error background error convariance into output Normalization NetCDF
!  file.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      FRC_LREAD_EIGEN : IF (.not.LwrtNRM(4,ng)) THEN
!
        ncname=NRM(4,ng)%name
        Text='forcing B-normalization file:'
        IF (Master) WRITE (stdout,10) TRIM(Text), TRIM(ncname)
!
!  If appropriate, open input NetCDF file.
!
        IF (NRM(4,ng)%IOtype.eq.io_nf90) THEN
          IF (NRM(4,ng)%ncid.eq.-1) THEN
            omode=1
            CALL netcdf_open (ng, model, ncname, omode,                 &
     &                        NRM(4,ng)%ncid)
            IF (FoundError(exit_flag, NoError,                          &
     &                      __LINE__, MyFile)) RETURN
          END IF
#  if defined PIO_LIB && defined DISTRIBUTE
        ELSE IF (NRM(4,ng)%IOtype.eq.io_pio) THEN
          IF (NRM(4,ng)%pioFile%fh.eq.-1) THEN
            omode=1
            CALL pio_netcdf_open (ng, model, ncname, omode,             &
     &                            NRM(4,ng)%pioFile)
          END IF
#  endif
        END IF
!
!  Read in minimum and maximum Ritz eigenvalues.
!
        SELECT CASE (NRM(4,ng)%IOtype)
!
          CASE (io_nf90)

#  ifdef ADJUST_WSTRESS
!
            CALL netcdf_get_fvar (ng, model, ncname,                    &
     &                            Vname(1,ieUsms),                      &
     &                            self%sustr_eigen,                     &
     &                            ncid = NRM(4,ng)%ncid,                &
     &                            start = (/1,1,1/),                    &
     &                            total = (/nl,ns,2/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%sustr_eigen,                             &
     &                    MASK = self%sustr_eigen > 0.0)
              Fmax=MAXVAL(self%sustr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieUsms)), Fmin, Fmax
            END IF
!
            CALL netcdf_get_fvar (ng, model, ncname,                    &
     &                            Vname(1,ieVsms),                      &
     &                            self%svstr_eigen,                     &
     &                            ncid = NRM(4,ng)%ncid,                &
     &                            start = (/1,1,1/),                    &
     &                            total = (/nl,ns,2/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%svstr_eigen,                             &
     &                    MASK = self%svstr_eigen > 0.0)
              Fmax=MAXVAL(self%svstr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieVsms)), Fmin, Fmax
            END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
            DO itrc=1,NT(ng)
              IF (Lstflux(itrc,ng)) THEN
                ifld=ieTsur(itrc)
                CALL netcdf_get_fvar (ng, model, ncname,                &
     &                                Vname(1,ifld),                    &
     &                                self%stflux_eigen(:,:,:,itrc),    &
     &                                ncid = NRM(4,ng)%ncid,            &
     &                                start = (/1,1,1/),                &
     &                                total = (/nl,ns,2/))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%stflux_eigen(:,:,:,itrc),            &
     &                    MASK = self%stflux_eigen(:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%stflux_eigen(:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#  endif
#  if defined PIO_LIB && defined DISTRIBUTE
!
          CASE (io_pio)

#   ifdef ADJUST_WSTRESS
!
            CALL pio_netcdf_get_fvar (ng, model, ncname,                &
     &                                Vname(1,ieUsms),                  &
     &                                self%sustr_eigen,                 &
     &                                pioFile = NRM(4,ng)%pioFile,      &
     &                                start = (/1,1,1/),                &
     &                                total = (/nl,ns,2/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%sustr_eigen,                             &
     &                    MASK = self%sustr_eigen > 0.0)
              Fmax=MAXVAL(self%sustr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieUsms)), Fmin, Fmax
            END IF
!
            CALL pio_netcdf_get_fvar (ng, model, ncname,                &
     &                                Vname(1,ieVsms),                  &
     &                                self%svstr_eigen,                 &
     &                                pioFile = NRM(4,ng)%pioFile,      &
     &                                start = (/1,1,1/),                &
     &                                total = (/nl,ns,2/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%svstr_eigen,                             &
     &                    MASK = self%svstr_eigen > 0.0)
              Fmax=MAXVAL(self%svstr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieVsms)), Fmin, Fmax
            END IF
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
!
            DO itrc=1,NT(ng)
              IF (Lstflux(itrc,ng)) THEN
                ifld=ieTsur(itrc)
                CALL pio_netcdf_get_fvar (ng, model, ncname,            &
     &                                    Vname(1,ifld),                &
     &                                    self%stflux_eigen(:,:,:,itrc),&
     &                                    pioFile = NRM(4,ng)%pioFile,  &
     &                                    start = (/1,1,1/),            &
     &                                    total = (/nl,ns,2/))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%stflux_eigen(:,:,:,itrc),            &
     &                    MASK = self%stflux_eigen(:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%stflux_eigen(:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#   endif
#  endif
        END SELECT
!
      END IF FRC_LREAD_EIGEN
# endif
!
  10  FORMAT (/,1x,'Reading extrema Ritz eigenvalues from ',a,2x,a/,3x, &
     &        '(Required by Implicit Chebyshev Iterations, CI, '        &
     &        'Solver)',/) 
  20  FORMAT (4x,'- ',a,':',t35,'Min = ',1p,e15.8,',  Max = ',1p,e15.8)
!
      RETURN
      END SUBROUTINE multiscale_eigen_read
!
!***********************************************************************
      SUBROUTINE multiscale_eigen_write (self, ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: i, ifld, itrc, nl, ns
!
      real (r8) :: Fmin, Fmax
!
      character (len=80 ) :: Text
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_eigen_write"
!
      SourceFile=MyFile
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Write Extrema K-Laplacian eigenvalues for Prior/Model background-
!  error background error convariance into output Normalization NetCDF
!  file(s).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      nl=MAXVAL(self%Mlap)  ! nunber of k-Laplacian applications
      ns=self%Nscale        ! number of correlation scales
!
      FILE_LOOP : DO i=1,NSA

        LWRITE_EIGEN : IF (LwrtNRM(i,ng)) THEN
!
          ncname=NRM(i,ng)%name
          IF (i.eq.1) THEN
            Text='prior B-normalization file:'
          ELSE IF (i.eq.2) THEN
            Text='model Q-normalization file:'
          END IF
          IF (Master) WRITE (stdout,10) TRIM(Text), TRIM(ncname)
!
          SELECT CASE (NRM(i,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieFsur),                    &
     &                              self%zeta_eigen,                    &
     &                              start = (/1,1,1/),                  &
     &                              total = (/nl,ns,2/),                &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              varid = NRM(i,ng)%Vid(ieFsur))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_eigen,                            &
     &                      MASK = self%zeta_eigen > 0.0)
                Fmax=MAXVAL(self%zeta_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieFsur)), Fmin, Fmax
              END IF
!
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieUbar),                    &
     &                              self%ubar_eigen,                    &
     &                              start = (/1,1,1/),                  &
     &                              total = (/nl,ns,2/),                &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              varid = NRM(i,ng)%Vid(ieUbar))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_eigen,                            &
     &                      MASK = self%ubar_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUbar)), Fmin, Fmax
              END IF
!
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieVbar),                    &
     &                              self%vbar_eigen,                    &
     &                              start = (/1,1,1/),                  &
     &                              total = (/nl,ns,2/),                &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              varid = NRM(i,ng)%Vid(ieVbar))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_eigen,                            &
     &                      MASK = self%vbar_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVbar)), Fmin, Fmax
              END IF

# ifdef SOLVE3D
!
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieUvel),                    &
     &                              self%u_eigen,                       &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,N(ng),ns,2/),          &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              varid = NRM(i,ng)%Vid(ieUvel))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_eigen,                               &
     &                      MASK = self%u_eigen > 0.0)
                Fmax=MAXVAL(self%u_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUvel)), Fmin, Fmax
              END IF
!
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieVvel),                    &
     &                              self%v_eigen,                       &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,N(ng),ns,2/),          &
     &                              ncid = NRM(i,ng)%ncid,              &
     &                              varid = NRM(i,ng)%Vid(ieVvel))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_eigen,                               &
     &                      MASK = self%v_eigen > 0.0)
                Fmax=MAXVAL(self%v_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVvel)), Fmin, Fmax
              END IF
!
              DO itrc=1,NT(ng)
                ifld=ieTvar(itrc)
                CALL netcdf_put_fvar (ng, model, ncname,                &
     &                                Vname(1,ifld),                    &
     &                                self%t_eigen(:,:,:,:,itrc),       &
     &                                start = (/1,1,1,1/),              &
     &                                total = (/nl,N(ng),ns,2/),        &
     &                                ncid = NRM(i,ng)%ncid,            &
     &                                varid = NRM(i,ng)%Vid(ifld))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_eigen(:,:,:,:,itrc),               &
     &                        MASK = self%t_eigen(:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_eigen(:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END DO
# endif

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieFsur),                &
     &                                  self%zeta_eigen,                &
     &                                  start = (/1,1,1/),              &
     &                                  total = (/nl,ns,2/),            &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                             pioVar = NRM(i,ng)%pioVar(ieFsur)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_eigen,                            &
     &                      MASK = self%zeta_eigen > 0.0)
                Fmax=MAXVAL(self%zeta_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieFsur)), Fmin, Fmax
              END IF
!
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieUbar),                &
     &                                  self%ubar_eigen,                &
     &                                  start = (/1,1,1/),              &
     &                                  total = (/nl,ns,2/),            &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                             pioVar = NRM(i,ng)%pioVar(ieUbar)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_eigen,                            &
     &                      MASK = self%ubar_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUbar)), Fmin, Fmax
              END IF
!
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieVbar),                &
     &                                  self%vbar_eigen,                &
     &                                  start = (/1,1,1/),              &
     &                                  total = (/nl,ns,2/),            &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                             pioVar = NRM(i,ng)%pioVar(ieVbar)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_eigen,                            &
     &                      MASK = self%vbar_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVbar)), Fmin, Fmax
              END IF

#  ifdef SOLVE3D
!
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieUvel),                &
     &                                  self%u_eigen,                   &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,N(ng),ns,2/),      &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                             pioVar = NRM(i,ng)%pioVar(ieUvel)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_eigen,                               &
     &                      MASK = self%u_eigen > 0.0)
                Fmax=MAXVAL(self%u_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieUvel)), Fmin, Fmax
              END IF
!
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieVvel),                &
     &                                  self%v_eigen,                   &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,N(ng),ns,2/),      &
     &                                  pioFile = NRM(i,ng)%pioFile,    &
     &                             pioVar = NRM(i,ng)%pioVar(ieVvel)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_eigen,                               &
     &                      MASK = self%v_eigen > 0.0)
                Fmax=MAXVAL(self%v_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieVvel)), Fmin, Fmax
              END IF
!
              DO itrc=1,NT(ng)
                ifld=ieTvar(itrc)
                CALL pio_netcdf_put_fvar (ng, model, ncname,            &
     &                                    Vname(1,ifld),                &
     &                                    self%t_eigen(:,:,:,:,itrc),   &
     &                                    start = (/1,1,1,1/),          &
     &                                    total = (/nl,N(ng),ns,2/),    &
     &                                    pioFile = NRM(i,ng)%pioFile,  &
     &                             pioVar = NRM(i,ng)%pioVar(ifld)%vd)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_eigen(:,:,:,:,itrc),               &
     &                        MASK = self%t_eigen(:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_eigen(:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END DO
#  endif
# endif
          END SELECT
!
        END IF LWRITE_EIGEN
      END DO FILE_LOOP

# ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Write Extrema K-Laplacian eigenvalues for open boundaries background-
!  error background error convariance into output Normalization NetCDF
!  file.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      BRY_LWRITE_EIGEN : IF (LwrtNRM(3,ng)) THEN
!
        ncname=NRM(3,ng)%name
        Text='boundary B-normalization file:'
        IF (Master) WRITE (stdout,10) TRIM(Text), TRIM(ncname)
!
        SELECT CASE (NRM(3,ng)%IOtype)
!
          CASE (io_nf90)
!
            IF (ANY(Lobc(:,isFsur,ng))) THEN
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieZbry),                    &
     &                              self%zeta_obc_eigen,                &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,4,ns,2/),              &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              varid = NRM(3,ng)%Vid(ieZbry))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_obc_eigen,                        &
     &                      MASK = self%zeta_obc_eigen > 0.0)
                Fmax=MAxVAL(self%zeta_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieZbry)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isUbar,ng))) THEN
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieU2bc),                    &
     &                              self%ubar_obc_eigen,                &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,4,ns,2/),              &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              varid = NRM(3,ng)%Vid(ieU2bc))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_obc_eigen,                        &
     &                      MASK = self%ubar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU2bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVbar,ng))) THEN
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieV2bc),                    &
     &                              self%vbar_obc_eigen,                &
     &                              start = (/1,1,1,1/),                &
     &                              total = (/nl,4,ns,2/),              &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              varid = NRM(3,ng)%Vid(ieV2bc))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_obc_eigen,                        &
     &                      MASK = self%vbar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV2bc)), Fmin, Fmax
              END IF
            END IF

#  ifdef SOLVE3D
!
            IF (ANY(Lobc(:,isUvel,ng))) THEN
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieU3bc),                    &
     &                              self%u_obc_eigen,                   &
     &                              start = (/1,1,1,1,1/),              &
     &                              total = (/nl,N(ng),4,ns,2/),        &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              varid = NRM(3,ng)%Vid(ieU3bc))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_obc_eigen,                           &
     &                      MASK = self%u_obc_eigen > 0.0)
                Fmax=MAXVAL(self%u_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU3bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVvel,ng))) THEN
              CALL netcdf_put_fvar (ng, model, ncname,                  &
     &                              Vname(1,ieV3bc),                    &
     &                              self%v_obc_eigen,                   &
     &                              start = (/1,1,1,1,1/),              &
     &                              total = (/nl,N(ng),4,ns,2/),        &
     &                              ncid = NRM(3,ng)%ncid,              &
     &                              varid = NRM(3,ng)%Vid(ieV3bc))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_obc_eigen,                           &
     &                      MASK = self%v_obc_eigen > 0.0)
                Fmax=MAXVAL(self%v_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV3bc)), Fmin, Fmax
              END IF
            END IF
!
            DO itrc=1,NT(ng)
              IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
                ifld=ieTbry(itrc)
                CALL netcdf_put_fvar (ng, model, ncname,                &
     &                                Vname(1,ifld),                    &
     &                                self%t_obc_eigen(:,:,:,:,:,itrc), &
     &                                start = (/1,1,1,1,1/),            &
     &                                total = (/nl,N(ng),4,ns,2/),      &
     &                                ncid = NRM(3,ng)%ncid,            &
     &                                varid = NRM(3,ng)%Vid(ifld))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_obc_eigen(:,:,:,:,:,itrc),         &
     &                    MASK = self%t_obc_eigen(:,:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_obc_eigen(:,:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#  endif

#  if defined PIO_LIB && defined DISTRIBUTE
!
          CASE (io_pio)
!
            IF (ANY(Lobc(:,isFsur,ng))) THEN
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieZbry),                &
     &                                  self%zeta_obc_eigen,            &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,4,ns,2/),          &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                             pioVar = NRM(3,ng)%pioVar(ieZbry)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%zeta_obc_eigen,                        &
     &                      MASK = self%zeta_obc_eigen > 0.0)
                Fmax=MAXVAL(self%zeta_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieZbry)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isUbar,ng))) THEN
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieU2bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,4,ns,2/),          &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                             pioVar = NRM(3,ng)%pioVar(ieU2bc)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%ubar_obc_eigen,                        &
     &                      MASK = self%ubar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%ubar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU2bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVbar,ng))) THEN
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieV2bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  start = (/1,1,1,1/),            &
     &                                  total = (/nl,4,ns,2/),          &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                             pioVar = NRM(3,ng)%pioVar(ieV2bc)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%vbar_obc_eigen,                        &
     &                      MASK = self%vbar_obc_eigen > 0.0)
                Fmax=MAXVAL(self%vbar_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV2bc)), Fmin, Fmax
              END IF
            END IF

#   ifdef SOLVE3D
!
            IF (ANY(Lobc(:,isUvel,ng))) THEN
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieU3bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  start = (/1,1,1,1,1/),          &
     &                                  total = (/nl,N(ng),4,ns,2/),    &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                             pioVar = NRM(3,ng)%pioVar(ieU3bc)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%u_obc_eigen,                           &
     &                      MASK = self%u_obc_eigen > 0.0)
                Fmax=MAXVAL(self%u_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieU3bc)), Fmin, Fmax
              END IF
            END IF
!
            IF (ANY(Lobc(:,isVvel,ng))) THEN
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,ieV3bc),                &
     &                                  self%ubar_obc_eigen,            &
     &                                  start = (/1,1,1,1,1/),          &
     &                                  total = (/nl,N(ng),4,ns,2/),    &
     &                                  pioFile = NRM(3,ng)%pioFile,    &
     &                             pioVar = NRM(3,ng)%pioVar(ieV3bc)%vd)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
              IF (Master) THEN 
                Fmin=MINVAL(self%v_obc_eigen,                           &
     &                      MASK = self%v_obc_eigen > 0.0)
                Fmax=MAXVAL(self%v_obc_eigen)
                WRITE (stdout,20) TRIM(Vname(1,ieV3bc)), Fmin, Fmax
              END IF
            END IF
!
            DO itrc=1,NT(ng)
              IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
                ifld=ieTbry(itrc)
                CALL pio_netcdf_put_fvar (ng, model, ncname,            &
     &                                    Vname(1,ifld),                &
     &                                    self%t_obc_eigen(:,:,:,:,:,   &
     &                                                     itrc),       &
     &                                    start = (/1,1,1,1,1/),        &
     &                                    total = (/nl,N(ng),4,ns,2/),  &
     &                                    pioFile = NRM(3,ng)%pioFile,  &
     &                               pioVar = NRM(3,ng)%pioVar(ifld)%vd)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%t_obc_eigen(:,:,:,:,:,itrc),         &
     &                    MASK = self%t_obc_eigen(:,:,:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%t_obc_eigen(:,:,:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#   endif
#  endif
        END SELECT
!
      END IF BRY_LWRITE_EIGEN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Write Extrema K-Laplacian eigenvalues for open boundaries background-
!  error background error convariance into output Normalization NetCDF
!  file.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      FRC_LWRITE_EIGEN : IF (LwrtNRM(4,ng)) THEN
!
        ncname=NRM(4,ng)%name
        Text='forcing B-normalization file:'
        IF (Master) WRITE (stdout,10) TRIM(Text), TRIM(ncname)
!
        SELECT CASE (NRM(4,ng)%IOtype)
!
          CASE (io_nf90)

#  ifdef ADJUST_WSTRESS
!
            CALL netcdf_put_fvar (ng, model, ncname,                    &
     &                            Vname(1,ieUsms),                      &
     &                            self%sustr_eigen,                     &
     &                            start = (/1,1,1/),                    &
     &                            total = (/nl,ns,2/),                  &
     &                            ncid = NRM(4,ng)%ncid,                &
     &                            varid = NRM(4,ng)%Vid(ieUsms))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%sustr_eigen,                             &
     &                    MASK = self%sustr_eigen > 0.0)
              Fmax=MAXVAL(self%sustr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieUsms)), Fmin, Fmax
            END IF
!
            CALL netcdf_put_fvar (ng, model, ncname,                    &
     &                            Vname(1,ieVsms),                      &
     &                            self%svstr_eigen,                     &
     &                            start = (/1,1,1/),                    &
     &                            total = (/nl,ns,2/),                  &
     &                            ncid = NRM(4,ng)%ncid,                &
     &                            varid = NRM(4,ng)%Vid(ieVsms))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%svstr_eigen,                             &
     &                    MASK = self%svstr_eigen > 0.0)
              Fmax=MAXVAL(self%svstr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieVsms)), Fmin, Fmax
            END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
            DO itrc=1,NT(ng)
              IF (Lstflux(itrc,ng)) THEN
                ifld=ieTsur(itrc)
                CALL netcdf_put_fvar (ng, model, ncname,                &
     &                                Vname(1,ifld),                    &
     &                                self%stflux_eigen(:,:,:,itrc),    &
     &                                start = (/1,1,1/),                &
     &                                total = (/nl,ns,2/),              &
     &                                ncid = NRM(4,ng)%ncid,            &
     &                                varid = NRM(4,ng)%Vid(ifld))
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%stflux_eigen(:,:,:,itrc),            &
     &                    MASK = self%stflux_eigen(:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%stflux_eigen(:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#  endif
#  if defined PIO_LIB && defined DISTRIBUTE
!
          CASE (io_pio)

#   ifdef ADJUST_WSTRESS
!
            CALL pio_netcdf_put_fvar (ng, model, ncname,                &
     &                                Vname(1,ieUsms),                  &
     &                                self%sustr_eigen,                 &
     &                                start = (/1,1,1/),                &
     &                                total = (/nl,ns,2/),              &
     &                                pioFile = NRM(4,ng)%pioFile,      &
     &                           pioVar = NRM(4,ng)%pioVar(ieUsms)%vd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%sustr_eigen,                             &
     &                    MASK = self%sustr_eigen > 0.0)
              Fmax=MAXVAL(self%sustr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieUsms)), Fmin, Fmax
            END IF
!
            CALL pio_netcdf_put_fvar (ng, model, ncname,                &
     &                                Vname(1,ieVsms),                  &
     &                                self%svstr_eigen,                 &
     &                                start = (/1,1,1/),                &
     &                                total = (/nl,ns,2/),              &
     &                                pioFile = NRM(4,ng)%pioFile,      &
     &                           pioVar = NRM(4,ng)%pioVar(ieVsms)%vd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (Master) THEN 
              Fmin=MINVAL(self%svstr_eigen,                             &
     &                    MASK = self%svstr_eigen > 0.0)
              Fmax=MAXVAL(self%svstr_eigen)
              WRITE (stdout,20) TRIM(Vname(1,ieVsms)), Fmin, Fmax
            END IF
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
!
            DO itrc=1,NT(ng)
              IF (Lstflux(itrc,ng)) THEN
                ifld=ieTsur(itrc)
                CALL pio_netcdf_put_fvar (ng, model, ncname,            &
     &                                    Vname(1,ifld),                &
     &                                    self%stflux_eigen(:,:,:,itrc),&
     &                                    start = (/1,1,1/),            &
     &                                    total = (/nl,ns,2/),          &
     &                                    pioFile = NRM(4,ng)%pioFile,  &
     &                               pioVar = NRM(4,ng)%pioVar(ifld)%vd)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
                IF (Master) THEN 
                  Fmin=MINVAL(self%stflux_eigen(:,:,:,itrc),            &
     &                    MASK = self%stflux_eigen(:,:,:,itrc) > 0.0)
                  Fmax=MAXVAL(self%stflux_eigen(:,:,:,itrc))
                  WRITE (stdout,20) TRIM(Vname(1,ifld)), Fmin, Fmax
                END IF
              END IF
            END DO
#   endif
#  endif
        END SELECT
!
      END IF FRC_LWRITE_EIGEN
# endif
!
  10  FORMAT (/,1x,'Writing extrema Ritz eigenvalues into ',a,2x,a/,3x, &
     &        '(Required by Implicit Chebyshev Iterations, CI, '        &
     &        'Solver)',/) 
  20  FORMAT (4x,'- ',a,':',t35,'Min = ',1p,e15.8,',  Max = ',1p,e15.8)
!
      RETURN
      END SUBROUTINE multiscale_eigen_write
#endif
      END MODULE multiscale_eigen_mod
