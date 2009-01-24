#include "cppdefs.h"
      MODULE cost_norm_mod
#if defined IS4DVAR || defined IS4DVAR_OLD
!
!svn $Id: cost_norm.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the cost function gradient norm by taking     !
!  dot product of the adjoint solution:                                !
!                                                                      !
!      gradient norm = SQRT ( transpose{grad(J)} * grad(J) )           !
!                                                                      !
!  It is only used as a diagnostic.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC cost_norm

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cost_norm (ng, tile, Linp)
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
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL cost_norm_tile (ng, tile,                                    &  
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Linp,                                        &
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ad_ustr,                        &
     &                     FORCES(ng) % ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % ad_tflux,                       &
#  endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
# else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
# endif
     &                     OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE cost_norm
!
!***********************************************************************
      SUBROUTINE cost_norm_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Linp,                                  &
# ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           ad_tflux,                              &
#  endif
     &                           ad_t,                                  &
     &                           ad_u,                                  &
     &                           ad_v,                                  &
# else
     &                           ad_ubar,                               &
     &                           ad_vbar,                               &
# endif
     &                           ad_zeta)
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
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

      real(r8), dimension(0:NstateVar(ng)) :: my_CostGrad

      real(r8) :: cff

# ifdef DISTRIBUTE
      character (len=3), dimension(0:NstateVar(ng)) :: op_handle
# endif

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute cost function gradient norm.
!----------------------------------------------------------------------
!
      DO i=0,NstateVar(ng)
        my_CostGrad(i)=0.0_r8
      END DO
!
!  Free-surface gradient norm.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=ad_zeta(i,j,Linp)*ad_zeta(i,j,Linp)
          my_CostGrad(0)=my_CostGrad(0)+cff
          my_CostGrad(isFsur)=my_CostGrad(isFsur)+cff
        END DO
      END DO

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress gradient norm.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=ad_ustr(i,j,k,Linp)*ad_ustr(i,j,k,Linp)
            my_CostGrad(0)=my_CostGrad(0)+cff
            my_CostGrad(isUstr)=my_CostGrad(isUstr)+cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=ad_vstr(i,j,k,Linp)*ad_vstr(i,j,k,Linp)
            my_CostGrad(0)=my_CostGrad(0)+cff
            my_CostGrad(isVstr)=my_CostGrad(isVstr)+cff
          END DO
        END DO
      END DO
# endif
# if !defined SOLVE3D
!
!  2D momentum gradient norm.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff=ad_ubar(i,j,Linp)*ad_ubar(i,j,Linp)
          my_CostGrad(0)=my_CostGrad(0)+cff
          my_CostGrad(isUbar)=my_CostGrad(isUbar)+cff
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff=ad_vbar(i,j,Linp)*ad_vbar(i,j,Linp)
          my_CostGrad(0)=my_CostGrad(0)+cff
          my_CostGrad(isVbar)=my_CostGrad(isVbar)+cff
        END DO
      END DO
# else
!
!  3D momentum gradient norm.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=ad_u(i,j,k,Linp)*ad_u(i,j,k,Linp)
            my_CostGrad(0)=my_CostGrad(0)+cff
            my_CostGrad(isUvel)=my_CostGrad(isUvel)+cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=ad_v(i,j,k,Linp)*ad_v(i,j,k,Linp)
            my_CostGrad(0)=my_CostGrad(0)+cff
            my_CostGrad(isVvel)=my_CostGrad(isVvel)+cff
          END DO
        END DO
      END DO
!
!  Tracers gradient norm.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=ad_t(i,j,k,Linp,itrc)*ad_t(i,j,k,Linp,itrc)
              my_CostGrad(0)=my_CostGrad(0)+cff
              my_CostGrad(isTvar(itrc))=my_CostGrad(isTvar(itrc))+cff
            END DO
          END DO
        END DO
      END DO
#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux gradient norm.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=ad_tflux(i,j,k,Linp,itrc)*ad_tflux(i,j,k,Linp,itrc)
              my_CostGrad(0)=my_CostGrad(0)+cff
              my_CostGrad(isTsur(itrc))=my_CostGrad(isTsur(itrc))+cff
            END DO
          END DO
        END DO
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
! Compute cost function gradient norm.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (COST_GRAD)
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%CostGrad(i)=0.0_r8
        END DO  
      END IF
      DO i=0,NstateVar(ng)
        FOURDVAR(ng)%CostGrad(i)=FOURDVAR(ng)%CostGrad(i)+my_CostGrad(i)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        DO i=0,NstateVar(ng)
          op_handle(i)='SUM'
        END DO 
        CALL mp_reduce (ng, iADM, NstateVar(ng)+1,                      &
     &                  FOURDVAR(ng)%CostGrad(0:),  op_handle(0:))
# endif
      END IF
!$OMP END CRITICAL (COST_GRAD)
!
!  Take norm square root.
!
      DO i=0,NstateVar(ng)
        FOURDVAR(ng)%CostGrad(i)=SQRT(FOURDVAR(ng)%CostGrad(i))
      END DO

      RETURN
      END SUBROUTINE cost_norm_tile
#endif
      END MODULE cost_norm_mod
