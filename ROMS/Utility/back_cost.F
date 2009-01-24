#include "cppdefs.h"
      MODULE back_cost_mod
#ifdef BACKGROUND
!
!svn $Id: back_cost.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the background cost function (Jb) as the      !
!  misfit (squared difference) between model and background state:     !
!                                                                      !
!    Jb = 1/2 * transpose(X - Xb) * B^(-1) * (X - Xb)                  !
!                                                                      !
!  where                                                               !
!                                                                      !
!    Xb : background or reference state (first guess)                  !
!    X  : nonlinear model initial state                                !
!    B  : background-error covariance                                  !
!                                                                      !
!  In incremental 4DVAR, the initial conditions estimate is done in    !
!  in minimization space by defining the following transformation:     !
!                                                                      !
!               deltaV = B^(-1/2) deltaX      (v-space)                !
!  or                                                                  !
!               deltaX = B^(1/2) deltaV       (x-space)                !
!                                                                      !
!  where                                                               !
!                                                                      !
!                    B = tanspose{B^(1/2)} B^(1/2)                     !
!                                                                      !
!  Then, the background cost function becomes:                         !
!                                                                      !
!    Jb = 1/2 * transpose{(X - Xb)} * B^(-1) * (X - Xb)     x-space    !
!                                                                      !
!  or                                                                  !
!                                                                      !
!    Jb = 1/2 * transpose(deltaX) * B^(-1) * deltaX         x-space    !
!                                                                      !
!  or                                                                  !
!                                                                      !
!    Jb = 1/2 * transpose(deltaV) * deltaV                  v-space    !
!                                                                      !
!  Therefore, in v-space the background cost function gradient is:     !
!                                                                      !
!    GRADv(Jb) = deltaV                                                !
!                                                                      !
!  Notice that initially, Jb is zero since the model is initialized    !
!  with the background state (X=Xb) and the minimization increment,    !
!  deltaX is zero.                                                     !
!                                                                      !
!                                                                      !
!  WARNING:                                                            !
!  -------                                                             !
!                                                                      !
!  The background cost function term is computed in v-space. Recall    !
!  that in the inner loop the increment vector deltaX (x-space) and    !
!  deltaV (v-space) are written into ITLname NetCDF file at records    !
!  1 and 2, respectively.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE back_cost (ng, tile, Lbck)
!***********************************************************************
!
      USE mod_param
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL back_cost_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lbck,                                        &
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
#  endif
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta)
      RETURN
      END SUBROUTINE back_cost
!
!***********************************************************************
      SUBROUTINE back_cost_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lbck,                                  &
# ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
# endif
# ifdef SOLVE3D
     &                           tl_t,                                  &
     &                           tl_u,                                  &
     &                           tl_v,                                  &
#  ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
#  endif
# else
     &                           tl_ubar,                               &
     &                           tl_vbar,                               &
# endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lbck
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
#  else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#   endif
#  else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

      real(r8), dimension(0:NstateVar(ng)) :: my_BackCost

      real(r8) :: cff

# ifdef DISTRIBUTE
      character (len=3), dimension(0:NstateVar(ng)) :: op_handle
# endif

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute background cost function (Jb) as misfit between model and
!  background states at initial time of the assimilation window.
!  Initially, the misfit innovation matrix (X-Xb) is zero. As the
!  assimilation algorithm iterates, Jb becomes greater than zero.
!----------------------------------------------------------------------
!
      DO i=0,NstateVar(ng)
        my_BackCost(i)=0.0_r8
      END DO
!
!  Free-surface contribution.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=0.5_r8*tl_zeta(i,j,Lbck)*tl_zeta(i,j,Lbck)
          my_BackCost(0)=my_BackCost(0)+cff
          my_BackCost(isFsur)=my_BackCost(isFsur)+cff
        END DO
      END DO

# if !defined SOLVE3D
!
!  2D momentum contribution.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff=0.5_r8*tl_ubar(i,j,Lbck)*tl_ubar(i,j,Lbck)
          my_BackCost(0)=my_BackCost(0)+cff
          my_BackCost(isUbar)=my_BackCost(isUbar)+cff
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff=0.5_r8*tl_vbar(i,j,Lbck)*tl_vbar(i,j,Lbck)
          my_BackCost(0)=my_BackCost(0)+cff
          my_BackCost(isVbar)=my_BackCost(isVbar)+cff
        END DO
      END DO
# else
!
!  3D momentum contribution.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=0.5_r8*tl_u(i,j,k,Lbck)*tl_u(i,j,k,Lbck)
            my_BackCost(0)=my_BackCost(0)+cff
            my_BackCost(isUvel)=my_BackCost(isUvel)+cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=0.5_r8*tl_v(i,j,k,Lbck)*tl_v(i,j,k,Lbck)
            my_BackCost(0)=my_BackCost(0)+cff
            my_BackCost(isVvel)=my_BackCost(isVvel)+cff
          END DO
        END DO
      END DO
!
!  Tracers contribution.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=0.5_r8*tl_t(i,j,k,Lbck,itrc)*tl_t(i,j,k,Lbck,itrc)
              my_BackCost(0)=my_BackCost(0)+cff
              my_BackCost(isTvar(itrc))=my_BackCost(isTvar(itrc))+cff
            END DO
          END DO
        END DO
      END DO
#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux contribution.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=0.5_r8*tl_tflux(i,j,k,Lbck,itrc)*                     &
     &                   tl_tflux(i,j,k,Lbck,itrc)
              my_BackCost(0)=my_BackCost(0)+cff
              my_BackCost(isTsur(itrc))=my_BackCost(isTsur(itrc))+cff
            END DO
          END DO
        END DO
      END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress contribution.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=0.5_r8*tl_ustr(i,j,k,Lbck)*tl_ustr(i,j,k,Lbck)
            my_BackCost(0)=my_BackCost(0)+cff
            my_BackCost(isUstr)=my_BackCost(isUstr)+cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=0.5_r8*tl_vstr(i,j,k,Lbck)*tl_vstr(i,j,k,Lbck)
            my_BackCost(0)=my_BackCost(0)+cff
            my_BackCost(isVstr)=my_BackCost(isVstr)+cff
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
! Compute global background cost function.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (BACKCOST)
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%BackCost(i)=my_BackCost(i)
        END DO
      ELSE
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%BackCost(i)=FOURDVAR(ng)%BackCost(i)+            &
     &                             my_BackCost(i)
        END DO
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        DO i=0,NstateVar(ng)
          op_handle(i)='SUM'
        END DO 
        CALL mp_reduce (ng, iNLM, NstateVar(ng)+1,                      &
     &                  FOURDVAR(ng)%BackCost(0:),  op_handle(0:))
# endif
      END IF
!$OMP END CRITICAL (BACKCOST)
!
!  If start of minimization, add background contribution to total cost
!  function scales used to report normalized values.
!
      IF (Nrun.eq.1) THEN
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%CostNorm(i)=FOURDVAR(ng)%CostNorm(i)+            &
     &                             FOURDVAR(ng)%BackCost(i)
        END DO
      END IF

      RETURN
      END SUBROUTINE back_cost_tile
#endif
      END MODULE back_cost_mod
