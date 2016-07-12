#include "cppdefs.h"
      MODULE adsen_initial_mod

#if defined ADJOINT          && \
   (defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI)
!
!svn $Id: adsen_initial.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes the adjoint state with the functional      !
!  whose sensitivity is required.                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: adsen_initial

      CONTAINS
!
!***********************************************************************
      SUBROUTINE adsen_initial (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
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
      CALL adsen_initial_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         knew(ng),                                &
# ifdef SOLVE3D
     &                         nstp(ng),                                &
# endif
     &                         GRID(ng) % Rscope,                       &
     &                         GRID(ng) % Uscope,                       &
     &                         GRID(ng) % Vscope,                       &
# ifdef SOLVE3D
     &                         CLIMA(ng) % u_adsG,                      &
     &                         CLIMA(ng) % v_adsG,                      &
     &                         CLIMA(ng) % t_adsG,                      &
# endif
     &                         CLIMA(ng) % ubar_adsG,                   &
     &                         CLIMA(ng) % vbar_adsG,                   &
     &                         CLIMA(ng) % zeta_adsG,                   &
# ifdef SOLVE3D
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % ad_v,                        &
     &                         OCEAN(ng) % ad_t,                        &
# endif
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ad_vbar,                     &
     &                         OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE adsen_initial
!
!***********************************************************************
      SUBROUTINE adsen_initial_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               knew,                              &
# ifdef SOLVE3D
     &                               nstp,                              &
# endif
     &                               Rscope, Uscope, Vscope,            &
# ifdef SOLVE3D
     &                               u_adsG, v_adsG, t_adsG,            &
# endif
     &                               ubar_adsG, vbar_adsG, zeta_adsG,   &
# ifdef SOLVE3D
     &                               ad_u, ad_v, ad_t,                  &
# endif
     &                               ad_ubar, ad_vbar, ad_zeta)
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
      integer, intent(in) :: nstp
# endif
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Rscope(LBi:,LBj:)
      real(r8), intent(in) :: Uscope(LBi:,LBj:)
      real(r8), intent(in) :: Vscope(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: u_adsG(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v_adsG(LBi:,LBj:,:,:)
      real(r8), intent(in) :: t_adsG(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(in) :: ubar_adsG(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_adsG(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta_adsG(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Rscope(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uscope(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vscope(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: u_adsG(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v_adsG(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: t_adsG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))
#  endif
      real(r8), intent(in) :: ubar_adsG(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vbar_adsG(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: zeta_adsG(LBi:UBi,LBj:UBj,2)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Initialize adjoint staye with the functional whose sensitivity is
!  required.  Use functional loaded into first record of climatological
!  arrays.
!------------------------------------------------------------------------
!
!  Free-surface.
!
      IF (SCALARS(ng)%Lstate(isFsur)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,knew)=zeta_adsG(i,j,1)*Rscope(i,j)
          END DO
        END DO
      END IF
!
!  2D Momentum.
!
      IF (SCALARS(ng)%Lstate(isUbar)) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,knew)=ubar_adsG(i,j,1)*Uscope(i,j)
          END DO
        END DO
      END IF
!
      IF (SCALARS(ng)%Lstate(isVbar)) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,knew)=vbar_adsG(i,j,1)*Vscope(i,j)
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
              ad_u(i,j,k,nstp)=u_adsG(i,j,k,1)*Uscope(i,j)
            END DO
          END DO
        END DO
      END IF
!
      IF (SCALARS(ng)%Lstate(isVvel)) THEN
        DO k=KstrS(ng),KendS(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,nstp)=v_adsG(i,j,k,1)*Vscope(i,j)
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
                ad_t(i,j,k,nstp,itrc)=t_adsG(i,j,k,1,itrc)*Rscope(i,j)
              END DO
            END DO
          END DO
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE adsen_initial_tile
#endif
      END MODULE adsen_initial_mod
