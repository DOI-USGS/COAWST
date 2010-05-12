#include "cppdefs.h"
      MODULE random_ic_mod

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
    defined POSTERIOR_ERROR_F)

# ifdef EW_PERIODIC
#  define IR_RANGE Istr-1,Iend+1
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JR_RANGE Jstr-1,Jend+1
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif

!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module initializes the tangent linear state variables with     !
!  a random vector for use when computing the eigenvectors of the      !
!  posterior analysis error covariance matrix.                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: random_ic

      CONTAINS
!
!***********************************************************************
      SUBROUTINE random_ic (ng, tile, model, innLoop, outLoop, Lout,    &
     &                      Ltrace)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, innLoop, outLoop, Lout
      logical, intent(in) :: Ltrace
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
!
      CALL random_ic_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lout, outLoop, Ltrace,                       &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % tl_t_obc,                     &
     &                     BOUNDARY(ng) % tl_u_obc,                     &
     &                     BOUNDARY(ng) % tl_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % tl_ubar_obc,                  &
     &                     BOUNDARY(ng) % tl_vbar_obc,                  &
     &                     BOUNDARY(ng) % tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
#  endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE random_ic
!
!***********************************************************************
      SUBROUTINE random_ic_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lout, outLoop, Ltrace,                 &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           tl_t_obc, tl_u_obc, tl_v_obc,          &
#  endif
     &                           tl_ubar_obc, tl_vbar_obc,              &
     &                           tl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
#  endif
     &                           tl_t, tl_u, tl_v,                      &
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
!
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE white_noise_mod
# ifdef DISTRIBUTE
#  ifdef ADJUST_BOUNDARY
      USE distribute_mod, ONLY : mp_collect
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lout, outLoop
      logical, intent(in) :: Ltrace
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j, ir, Zscheme
# ifdef SOLVE3D
      integer :: itrc, k
# endif

      real(r8) :: Amax, Amin, Bmax, Bmin

      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
# ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
# endif
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
#  ifdef ADJUST_BOUNDARY
      integer :: ib
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
#  endif
# endif

# include "set_bounds.h"
!
      SourceFile='random_ic.F, random_ic_tile'
!
!-----------------------------------------------------------------------
!  Generates random initial vectors for computation of the posterior
!  analysis error covariance matrix EOFs.
!-----------------------------------------------------------------------
!
!  Always use Gaussian distribution between -1 and +1, by Zscheme=1.
!
      Zscheme=1
!
!  2D random initialization at RHO-points.
!
      CALL white_noise2d (ng, iTLM, r2dvar, Zscheme,                    &
     &                    IstrR, IendR, JstrR, JendR,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      IF (.not.Ltrace) THEN
        DO j=JR_RANGE
           DO i=IR_RANGE
             tl_zeta(i,j,Lout)=A2d(i,j)
# ifdef MASKING
             tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
# endif
           END DO
        END DO
      ELSE
        DO j=JR_RANGE
           DO i=IR_RANGE
             tl_zeta(i,j,Lout)=DSIGN(1.0_r8,A2d(i,j))
# ifdef MASKING
             tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
# endif
           END DO
        END DO
      ENDIF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_zeta(:,:,Lout))
# endif
# ifndef SOLVE3D
!
!  2D random initialization at U-points.
!
      CALL white_noise2d (ng, iTLM, u2dvar, Zscheme,                    &
     &                    Istr, IendR, JstrR, JendR,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      IF (.not.Ltrace) THEN
        DO j=JU_RANGE
          DO i=IU_RANGE
            tl_ubar(i,j,Lout)=A2d(i,j)
# ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
# endif
          END DO
        END DO
      ELSE
        DO j=JU_RANGE
          DO i=IU_RANGE
            tl_ubar(i,j,Lout)=DSIGN(1.0_r8,A2d(i,j))
# ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
# endif
          END DO
        END DO
      END IF
!
!  2D random initialization at V-points.
!
      CALL white_noise2d (ng, iTLM, v2dvar, Zscheme,                    &
     &                    IstrR, IendR, Jstr, JendR,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Amin, Amax, A2d)
      IF (.not.Ltrace) THEN
        DO j=JV_RANGE
          DO i=IV_RANGE
            tl_vbar(i,j,Lout)=A2d(i,j)
# ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
# endif
          END DO
        END DO
      ELSE
        DO j=JV_RANGE
          DO i=IV_RANGE
            tl_vbar(i,j,Lout)=DSIGN(1.0_r8,A2d(i,j))
# ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
# endif
          END DO
        END DO
      END IF
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_ubar(:,:,Lout),                            &
     &                    tl_vbar(:,:,Lout))
#  endif
# endif

# ifdef SOLVE3D
!
!  3D random initialization U-points.
!
      CALL white_noise3d (ng, iTLM, u3dvar, Zscheme,                    &
     &                    Istr, IendR, JstrR, JendR,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Amin, Amax, A3d)
      IF (.not.Ltrace) THEN
        DO k=1,N(ng)
          DO j=JU_RANGE
            DO i=IU_RANGE
              tl_u(i,j,k,Lout)=A3d(i,j,k)
#  ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#  endif
            END DO
          END DO
        END DO
      ELSE
        DO k=1,N(ng)
          DO j=JU_RANGE
            DO i=IU_RANGE
              tl_u(i,j,k,Lout)=DSIGN(1.0_r8,A3d(i,j,k))
#  ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#  endif
            END DO
          END DO
        END DO
      END IF
!
!  3D random initialization at V-points.
!
       CALL white_noise3d (ng, iTLM, v3dvar, Zscheme,                   &
     &                     IstrR, IendR, Jstr, JendR,                   &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     Amin, Amax, A3d)
      IF (.not.Ltrace) THEN
        DO k=1,N(ng)
          DO j=JV_RANGE
            DO i=IV_RANGE
              tl_v(i,j,k,Lout)=A3d(i,j,k)
#  ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#  endif
            END DO
          END DO
        END DO
      ELSE
        DO k=1,N(ng)
          DO j=JV_RANGE
            DO i=IV_RANGE
              tl_v(i,j,k,Lout)=DSIGN(1.0_r8,A3d(i,j,k))
#  ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#  endif
            END DO
          END DO
        END DO
      ENDIF
#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_u(:,:,:,Lout), tl_v(:,:,:,Lout))
#  endif
!
!  3D random initialization at RHO-points.
!
       DO itrc=1,NT(ng)
         CALL white_noise3d (ng, iTLM, r3dvar, Zscheme,                 &
     &                       IstrR, IendR, JstrR, JendR,                &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       Amin, Amax, A3d)
         IF (.not.Ltrace) THEN
           DO k=1,N(ng)
             DO j=JR_RANGE
               DO i=IR_RANGE
                 tl_t(i,j,k,Lout,itrc)=A3d(i,j,k)
#  ifdef MASKING
                 tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
#  endif
               END DO
             END DO
           END DO
         ELSE
           DO k=1,N(ng)
             DO j=JR_RANGE
               DO i=IR_RANGE
                 tl_t(i,j,k,Lout,itrc)=DSIGN(1.0_r8,A3d(i,j,k))
# ifdef MASKING
                 tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
# endif
               END DO
             END DO
           END DO
         END IF
       END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_t(:,:,:,Lout,:))
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D boundary random initialization at RHO-points.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isFsur,ng)) THEN
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              CALL white_noise2d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                JstrR, JendR,                     &
     &                                LBij, UBij,                       &
     &                                Bmin, Bmax, B2d)
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              CALL white_noise2d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                IstrR, IendR,                     &
     &                                LBij, UBij,                       &
     &                                Bmin, Bmax, B2d)
            END IF
            IF (((ib.eq.iwest).and.WESTERN_EDGE).or.                    &
     &          ((ib.eq.ieast).and.EASTERN_EDGE)) THEN
              i=BOUNDS(ng)%edge(ib,r2dvar)
              IF (.not.Ltrace) THEN
                DO j=JR_RANGE
                  tl_zeta_obc(j,ib,ir,Lout)=B2d(j)
#  ifdef MASKING
                  tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Lout)*  &
     &                                      rmask(i,j)
#  endif
                END DO
              ELSE
                DO j=JR_RANGE
                  tl_zeta_obc(j,ib,ir,Lout)=DSIGN(1.0_r8,B2d(j))
#  ifdef MASKING
                  tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Lout)*  &
     &                                      rmask(i,j)
#  endif
                END DO
              END IF
            ELSE IF (((ib.eq.isouth).and.SOUTHERN_EDGE).or.             &
     &               ((ib.eq.inorth).and.NORTHERN_EDGE)) THEN
              j=BOUNDS(ng)%edge(ib,r2dvar)
              IF (.not.Ltrace) THEN
                DO i=IR_RANGE
                  tl_zeta_obc(i,ib,ir,Lout)=B2d(i)
#  ifdef MASKING
                  tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Lout)*  &
     &                                      rmask(i,j)
#  endif
                END DO
              ELSE
                DO i=IR_RANGE
                  tl_zeta_obc(i,ib,ir,Lout)=DSIGN(1.0_r8,B2d(i))
#  ifdef MASKING
                  tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Lout)*  &
     &                                      rmask(i,j)
#  endif
                END DO
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_bry (ng, tile, iTLM, 1, ib,              &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic, NSperiodic,             &
     &                              tl_zeta_obc(:,ib,ir,Lout))
#  endif
          END IF
        END DO
      END DO
!
!  2D boundary random initialization at U-points.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isUbar,ng)) THEN
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              CALL white_noise2d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                JstrR, JendR,                     &
     &                                LBij, UBij,                       &
     &                                Bmin, Bmax, B2d)
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              CALL white_noise2d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                Istr, IendR,                      &
     &                                LBij, UBij,                       &
     &                                Bmin, Bmax, B2d)
            END IF
            IF (((ib.eq.iwest).and.WESTERN_EDGE).or.                    &
     &          ((ib.eq.ieast).and.EASTERN_EDGE)) THEN
              i=BOUNDS(ng)%edge(ib,u2dvar)
              IF (.not.Ltrace) THEN
                DO j=JU_RANGE
                  tl_ubar_obc(j,ib,ir,Lout)=B2d(j)
#  ifdef MASKING
                  tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Lout)*  &
     &                                      umask(i,j)
#  endif
                END DO
              ELSE
                DO j=JU_RANGE
                  tl_ubar_obc(j,ib,ir,Lout)=DSIGN(1.0_r8,B2d(j))
#  ifdef MASKING
                  tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Lout)*  &
     &                                      umask(i,j)
#  endif
                END DO
              END IF
            ELSE IF (((ib.eq.isouth).and.SOUTHERN_EDGE).or.             &
     &               ((ib.eq.inorth).and.NORTHERN_EDGE)) THEN
              j=BOUNDS(ng)%edge(ib,u2dvar)
              IF (.not.Ltrace) THEN
                DO i=IU_RANGE
                  tl_ubar_obc(i,ib,ir,Lout)=B2d(i)
#  ifdef MASKING
                  tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Lout)*  &
     &                                      umask(i,j)
#  endif
                END DO
              ELSE
                DO i=IU_RANGE
                  tl_ubar_obc(i,ib,ir,Lout)=DSIGN(1.0_r8,B2d(i))
#  ifdef MASKING
                  tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Lout)*  &
     &                                      umask(i,j)
#  endif
                END DO
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_bry (ng, tile, iTLM, 1, ib,              &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic, NSperiodic,             &
     &                              tl_ubar_obc(:,ib,ir,Lout))
#  endif
          END IF
        END DO
      END DO
!
!  2D boundary random initialization at V-points.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isVbar,ng)) THEN
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              CALL white_noise2d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                Jstr, JendR,                      &
     &                                LBij, UBij,                       &
     &                                Bmin, Bmax, B2d)
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              CALL white_noise2d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                IstrR, IendR,                     &
     &                                LBij, UBij,                       &
     &                                Bmin, Bmax, B2d)
            END IF
            IF (((ib.eq.iwest).and.WESTERN_EDGE).or.                    &
     &          ((ib.eq.ieast).and.EASTERN_EDGE)) THEN
              i=BOUNDS(ng)%edge(ib,v2dvar)
              IF (.not.Ltrace) THEN
                DO j=JV_RANGE
                  tl_vbar_obc(j,ib,ir,Lout)=B2d(j)
#  ifdef MASKING
                  tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Lout)*  &
     &                                      vmask(i,j)
#  endif
                END DO
              ELSE
                DO j=JV_RANGE
                  tl_vbar_obc(j,ib,ir,Lout)=DSIGN(1.0_r8,B2d(j))
#  ifdef MASKING
                  tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Lout)*  &
     &                                      vmask(i,j)
#  endif
                END DO
              END IF
            ELSE IF (((ib.eq.isouth).and.SOUTHERN_EDGE).or.             &
     &               ((ib.eq.inorth).and.NORTHERN_EDGE)) THEN
              j=BOUNDS(ng)%edge(ib,v2dvar)
              IF (.not.Ltrace) THEN
                DO i=IV_RANGE
                  tl_vbar_obc(i,ib,ir,Lout)=B2d(i)
#  ifdef MASKING
                  tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Lout)*  &
     &                                      vmask(i,j)
#  endif
                END DO
              ELSE
                DO i=IV_RANGE
                  tl_vbar_obc(i,ib,ir,Lout)=DSIGN(1.0_r8,B2d(i))
#  ifdef MASKING
                  tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Lout)*  &
     &                                      vmask(i,j)
#  endif
                END DO
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange2d_bry (ng, tile, iTLM, 1, ib,              &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic, NSperiodic,             &
     &                              tl_vbar_obc(:,ib,ir,Lout))
#  endif
          END IF
        END DO
      END DO

#  ifdef SOLVE3D
!
!  3D boundary norm at U-points.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isUvel,ng)) THEN
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              CALL white_noise3d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                JstrR, JendR,                     &
     &                                LBij, UBij, 1, N(ng),             &
     &                                Bmin, Bmax, B3d)
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              CALL white_noise3d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                Istr, IendR,                      &
     &                                LBij, UBij, 1, N(ng),             &
     &                                Bmin, Bmax, B3d)
            END IF
            IF (((ib.eq.iwest).and.WESTERN_EDGE).or.                    &
     &          ((ib.eq.ieast).and.EASTERN_EDGE)) THEN
              i=BOUNDS(ng)%edge(ib,u2dvar)
              IF (.not.Ltrace) THEN
                DO k=1,N(ng)
                  DO j=JU_RANGE
                    tl_u_obc(j,k,ib,ir,Lout)=B3d(j,k)
#   ifdef MASKING
                    tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Lout)*  &
     &                                       umask(i,j)
#   endif
                  END DO
                END DO
              ELSE
                DO k=1,N(ng)
                  DO j=JU_RANGE
                    tl_u_obc(j,k,ib,ir,Lout)=DSIGN(1.0_r8,B3d(j,k))
#   ifdef MASKING
                    tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Lout)*  &
     &                                       umask(i,j)
#   endif
                  END DO
                END DO
              END IF
            ELSE IF (((ib.eq.isouth).and.SOUTHERN_EDGE).or.             &
     &               ((ib.eq.inorth).and.NORTHERN_EDGE)) THEN
              j=BOUNDS(ng)%edge(ib,u2dvar)
              IF (.not.Ltrace) THEN
                DO k=1,N(ng)
                  DO i=IU_RANGE
                    tl_u_obc(i,k,ib,ir,Lout)=B3d(i,k)
#   ifdef MASKING
                    tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Lout)*  &
     &                                       umask(i,j)
#   endif
                  END DO
                END DO
              ELSE
                DO k=1,N(ng)
                  DO i=IU_RANGE
                    tl_u_obc(i,k,ib,ir,Lout)=DSIGN(1.0_r8,B3d(i,k))
#   ifdef MASKING
                    tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Lout)*  &
     &                                       umask(i,j)
#   endif
                  END DO
                END DO
              END IF
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange3d_bry (ng, tile, iTLM, 1, ib,              &
     &                              LBij, UBij, 1, N(ng),               &
     &                              NghostPoints,                       &
     &                              EWperiodic, NSperiodic,             &
     &                              tl_u_obc(:,:,ib,ir,Lout))
#   endif
          END IF
        END DO
      END DO
!
!  3D boundary random initialization at V-points.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isVvel,ng)) THEN
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              CALL white_noise3d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                Jstr, JendR,                      &
     &                                LBij, UBij, 1, N(ng),             &
     &                                Bmin, Bmax, B3d)
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              CALL white_noise3d_bry (ng, tile, iTLM, ib,               &
     &                                Zscheme,                          &
     &                                IstrR, IendR,                     &
     &                                LBij, UBij, 1, N(ng),             &
     &                                Bmin, Bmax, B3d)
            END IF
            IF (((ib.eq.iwest).and.WESTERN_EDGE).or.                    &
     &          ((ib.eq.ieast).and.EASTERN_EDGE)) THEN
              i=BOUNDS(ng)%edge(ib,v2dvar)
              IF (.not.Ltrace) THEN
                DO k=1,N(ng)
                  DO j=JV_RANGE
                    tl_v_obc(j,k,ib,ir,Lout)=B3d(j,k)
#   ifdef MASKING
                    tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Lout)*  &
     &                                       vmask(i,j)
#   endif
                  END DO
                END DO
              ELSE
                DO k=1,N(ng)
                  DO j=JV_RANGE
                    tl_v_obc(j,k,ib,ir,Lout)=DSIGN(1.0_r8,B3d(j,k))
#   ifdef MASKING
                    tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Lout)*  &
     &                                       vmask(i,j)
#   endif
                  END DO
                END DO
              END IF
            ELSE IF (((ib.eq.isouth).and.SOUTHERN_EDGE).or.             &
     &               ((ib.eq.inorth).and.NORTHERN_EDGE)) THEN
              j=BOUNDS(ng)%edge(ib,v2dvar)
              IF (.not.Ltrace) THEN
                DO k=1,N(ng)
                  DO i=IV_RANGE
                    tl_v_obc(i,k,ib,ir,Lout)=B3d(i,k)
#   ifdef MASKING
                    tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Lout)*  &
     &                                       vmask(i,j)
#   endif
                  END DO
                END DO
              ELSE
                DO k=1,N(ng)
                  DO i=IV_RANGE
                    tl_v_obc(i,k,ib,ir,Lout)=DSIGN(1.0_r8,B3d(i,k))
#   ifdef MASKING
                    tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Lout)*  &
     &                                       vmask(i,j)
#   endif
                  END DO
                END DO
              END IF
            END IF
#   ifdef DISTRIBUTE
            CALL mp_exchange3d_bry (ng, tile, iTLM, 1, ib,              &
     &                              LBij, UBij, 1, N(ng),               &
     &                              NghostPoints,                       &
     &                              EWperiodic, NSperiodic,             &
     &                              tl_v_obc(:,:,ib,ir,Lout))
#   endif
          END IF
        END DO
      END DO
!
!  3D boundary random initialization at RHO-points.
!
      DO itrc=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (Lobc(ib,isTvar(itrc),ng)) THEN
              IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ib,r2dvar)
                CALL white_noise3d_bry (ng, tile, iTLM, ib,             &
     &                                  Zscheme,                        &
     &                                  JstrR, JendR,                   &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  Bmin, Bmax, B3d)
              ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ib,r2dvar)
                CALL white_noise3d_bry (ng, tile, iTLM, ib,             &
     &                                  Zscheme,                        &
     &                                  IstrR, IendR,                   &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  Bmin, Bmax, B3d)
              END IF
              IF (((ib.eq.iwest).and.WESTERN_EDGE).or.                  &
     &            ((ib.eq.ieast).and.EASTERN_EDGE)) THEN
                i=BOUNDS(ng)%edge(ib,r2dvar)
                IF (.not.Ltrace) THEN
                  DO k=1,N(ng)
                    DO j=JR_RANGE
                      tl_t_obc(j,k,ib,ir,Lout,itrc)=B3d(j,k)
#   ifdef MASKING
                      tl_t_obc(j,k,ib,ir,Lout,itrc)=                    &
     &                         tl_t_obc(j,k,ib,ir,Lout,itrc)*rmask(i,j)
#   endif
                    END DO
                  END DO
                ELSE
                  DO k=1,N(ng)
                    DO j=JR_RANGE
                      tl_t_obc(j,k,ib,ir,Lout,itrc)=                    &
     &                                        DSIGN(1.0_r8,B3d(j,k))
#   ifdef MASKING
                      tl_t_obc(j,k,ib,ir,Lout,itrc)=                    &
     &                         tl_t_obc(j,k,ib,ir,Lout,itrc)*rmask(i,j)
#   endif
                    END DO
                  END DO
                END IF
              ELSE IF (((ib.eq.isouth).and.SOUTHERN_EDGE).or.           &
     &                 ((ib.eq.inorth).and.NORTHERN_EDGE)) THEN
                j=BOUNDS(ng)%edge(ib,r2dvar)
                IF (.not.Ltrace) THEN
                  DO k=1,N(ng)
                    DO i=IR_RANGE
                      tl_t_obc(i,k,ib,ir,Lout,itrc)=B3d(i,k)
#   ifdef MASKING
                      tl_t_obc(i,k,ib,ir,Lout,itrc)=                    &
     &                         tl_t_obc(i,k,ib,ir,Lout,itrc)*rmask(i,j)
#   endif
                    END DO
                  END DO
                ELSE
                  DO k=1,N(ng)
                    DO i=IR_RANGE
                      tl_t_obc(i,k,ib,ir,Lout,itrc)=                    &
     &                                        DSIGN(1.0_r8,B3d(i,k))
#   ifdef MASKING
                      tl_t_obc(i,k,ib,ir,Lout,itrc)=                    &
     &                         tl_t_obc(i,k,ib,ir,Lout,itrc)*rmask(i,j)
#   endif
                    END DO
                  END DO
                END IF
              END IF
#   ifdef DISTRIBUTE
              CALL mp_exchange3d_bry (ng, tile, iTLM, 1, ib,            &
     &                                LBij, UBij, 1, N(ng),             &
     &                                NghostPoints,                     &
     &                                EWperiodic, NSperiodic,           &
     &                                tl_t_obc(:,:,ib,ir,Lout,itrc))
#   endif
            END IF
          END DO
        END DO
      END DO
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
#  ifdef ADJUST_WSTRESS
!
!  2D random initialization at U-stress points.
!
      DO ir=1,Nfrec(ng)
        CALL white_noise2d (ng, iTLM, u2dvar, Zscheme,                  &
     &                      Istr, IendR, JstrR, JendR,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Amin, Amax, A2d)
        IF (.not.Ltrace) THEN
          DO j=JU_RANGE
            DO i=IU_RANGE
              tl_ustr(i,j,ir,Lout)=A2d(i,j)
#   ifdef MASKING
              tl_ustr(i,j,ir,Lout)=tl_ustr(i,j,ir,Lout)*umask(i,j)
#   endif
            END DO
          END DO
        ELSE
          DO j=JU_RANGE
            DO i=IU_RANGE
              tl_ustr(i,j,ir,Lout)=DSIGN(1.0_r8,A2d(i,j))
#   ifdef MASKING
              tl_ustr(i,j,ir,Lout)=tl_ustr(i,j,ir,Lout)*umask(i,j)
#   endif
            END DO
          END DO
        END IF
      END DO
!
!  2D random initialization at V-stress points.
!
      DO ir=1,Nfrec(ng)
        CALL white_noise2d (ng, iTLM, v2dvar, Zscheme,                  &
     &                      IstrR, IendR, Jstr, JendR,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Amin, Amax, A2d)
        IF (.not.Ltrace) THEN
          DO j=JV_RANGE
            DO i=IV_RANGE
              tl_vstr(i,j,ir,Lout)=A2d(i,j)
#   ifdef MASKING
              tl_vstr(i,j,ir,Lout)=tl_vstr(i,j,ir,Lout)*vmask(i,j)
#   endif
            END DO
          END DO
        ELSE
          DO j=JV_RANGE
            DO i=IV_RANGE
              tl_vstr(i,j,ir,Lout)=DSIGN(1.0_r8,A2d(i,j))
#   ifdef MASKING
              tl_vstr(i,j,ir,Lout)=tl_vstr(i,j,ir,Lout)*vmask(i,j)
#   endif
            END DO
          END DO
        END IF
      END DO
#   ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nfrec(ng),             &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_ustr(:,:,:,Lout),                          &
     &                    tl_vstr(:,:,:,Lout))
#   endif
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  2D random initialization at surface tracer flux points.
!
      DO itrc=1,NT(ng)
        DO ir=1,Nfrec(ng)
          CALL white_noise2d (ng, iTLM, r2dvar, Zscheme,                &
     &                        IstrR, IendR, JstrR, JendR,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Amin, Amax, A2d)
          IF (.not.Ltrace) THEN
            DO j=JR_RANGE
              DO i=IR_RANGE
                tl_tflux(i,j,ir,Lout,itrc)=A2d(i,j)
#   ifdef MASKING
                tl_tflux(i,j,ir,Lout,itrc)=tl_tflux(i,j,ir,Lout,itrc)*  &
     &                                     rmask(i,j)
#   endif
              END DO
            END DO
          ELSE
            DO j=JR_RANGE
              DO i=IR_RANGE
                tl_tflux(i,j,ir,Lout,itrc)=DSIGN(1.0_r8,A2d(i,j))
#   ifdef MASKING
                tl_tflux(i,j,ir,Lout,itrc)=tl_tflux(i,j,ir,Lout,itrc)*  &
     &                                     rmask(i,j)
#   endif
              END DO
            END DO
          END IF
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, Nfrec(ng),         &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        tl_tflux(:,:,:,Lout,itrc))
#   endif
        END DO
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE random_ic_tile
#endif

      END MODULE random_ic_mod
