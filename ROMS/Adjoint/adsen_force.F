#include "cppdefs.h"
      MODULE adsen_force_mod

#if defined ADJOINT          && \
   (defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS)
!
!svn $Id: adsen_force.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine forces the adjoint state with the functional whose     !
!  sensitivity is required.                                            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: adsen_force

      CONTAINS
!
!***********************************************************************
      SUBROUTINE adsen_force (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL adsen_force_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       knew(ng),                                  &
# ifdef SOLVE3D
     &                       nnew(ng),                                  &
# endif
     &                       GRID(ng) % Rscope,                         &
     &                       GRID(ng) % Uscope,                         &
     &                       GRID(ng) % Vscope,                         &
# ifdef SOLVE3D
     &                       CLIMA(ng) % u_ads,                         &
     &                       CLIMA(ng) % v_ads,                         &
     &                       CLIMA(ng) % t_ads,                         &
# endif
     &                       CLIMA(ng) % ubar_ads,                      &
     &                       CLIMA(ng) % vbar_ads,                      &
     &                       CLIMA(ng) % zeta_ads,                      &
# ifdef SOLVE3D
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
     &                       OCEAN(ng) % ad_t,                          &
     &                       COUPLING(ng) % ad_Zt_avg1,                 &
# else
     &                       OCEAN(ng) % ad_zeta,                       &
# endif
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar)

      RETURN
      END SUBROUTINE adsen_force
!
!***********************************************************************
      SUBROUTINE adsen_force_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             knew,                                &
# ifdef SOLVE3D
     &                             nnew,                                &
# endif
     &                             Rscope, Uscope, Vscope,              &
# ifdef SOLVE3D
     &                             u_ads, v_ads, t_ads,                 &
# endif
     &                             ubar_ads, vbar_ads, zeta_ads,        &
# ifdef SOLVE3D
     &                             ad_u, ad_v, ad_t,                    &
     &                             ad_Zt_avg1,                          &
# else
     &                             ad_zeta,                             &
# endif
     &                             ad_ubar, ad_vbar)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
# ifdef SOLVE3D
      integer, intent(in) :: nnew
# endif
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Rscope(LBi:,LBj:)
      real(r8), intent(in) :: Uscope(LBi:,LBj:)
      real(r8), intent(in) :: Vscope(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: u_ads(LBi:,LBj:,:)
      real(r8), intent(in) :: v_ads(LBi:,LBj:,:)
      real(r8), intent(in) :: t_ads(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ubar_ads(LBi:,LBj:)
      real(r8), intent(in) :: vbar_ads(LBi:,LBj:)
      real(r8), intent(in) :: zeta_ads(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Rscope(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uscope(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vscope(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: u_ads(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_ads(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t_ads(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
      real(r8), intent(in) :: ubar_ads(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_ads(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta_ads(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Kfrc, i, itrc, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint staye with the functional whose sensitivity is
!  required.  Use functional loaded into first record of climatological
!  arrays.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntend(ng)) THEN
        Kfrc=knew
      ELSE
        Kfrc=1
      END IF
!
!  Free-surface.
!
      IF (SCALARS(ng)%Lstate(isFsur)) THEN
# ifdef SOLVE3D
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                            &
     &                      zeta_ads(i,j)*Rscope(i,j)
          END DO
        END DO
# else
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Kfrc)=ad_zeta(i,j,Kfrc)+                        &
     &                        zeta_ads(i,j)*Rscope(i,j)
          END DO
        END DO
# endif
      END IF
!
!  2D Momentum.
!
      IF (SCALARS(ng)%Lstate(isUbar)) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,Kfrc)=ad_ubar(i,j,Kfrc)+                        &
     &                        ubar_ads(i,j)*Uscope(i,j)
          END DO
        END DO
      END IF
!
      IF (SCALARS(ng)%Lstate(isVbar)) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Kfrc)=ad_vbar(i,j,Kfrc)+                        &
     &                        vbar_ads(i,j)*Vscope(i,j)
          END DO
        END DO
      END IF
# ifdef SOLVE3D
!
!  3D Momentum.
!
      IF (SCALARS(ng)%Lstate(isUvel)) THEN
        DO k=KstrS(ng),KendS(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+                        &
     &                         u_ads(i,j,k)*Uscope(i,j)
            END DO
          END DO
        END DO
      END IF    
!
      IF (SCALARS(ng)%Lstate(isVvel)) THEN
        DO k=KstrS(ng),KendS(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+                        &
     &                         v_ads(i,j,k)*Vscope(i,j)
            END DO
          END DO
        END DO
      END IF
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Lstate(isTvar(itrc))) THEN
          DO k=KstrS(ng),KendS(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_t(i,j,k,nnew,itrc)=ad_t(i,j,k,nnew,itrc)+            &
     &                                t_ads(i,j,k,itrc)*Rscope(i,j)
              END DO
            END DO
          END DO
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE adsen_force_tile
#endif
      END MODULE adsen_force_mod
