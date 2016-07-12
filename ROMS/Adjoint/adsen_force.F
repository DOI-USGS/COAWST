#include "cppdefs.h"
      MODULE adsen_force_mod

#if defined ADJOINT          && \
   (defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR)
!
!svn $Id: adsen_force.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
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
     &                       nnew(ng), nstp(ng),                        &
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
     &                             nnew, nstp,                          &
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
      USE mod_parallel
      USE mod_iounits
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
      integer, intent(in) :: nnew, nstp
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
      integer :: Kfrc, Nfrc, i, itrc, j, k

      real(r8) :: adFac

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
        Nfrc=nstp
      ELSE
        Kfrc=1
        Nfrc=nnew
      END IF

# ifdef AD_IMPULSE
!
!  Impose the forcing in the adjoint model as impulses in a manner that
!  is consistent with the definition of the sensitivity functional.
!  Notice that nTLM is used to control the forcing since nADJ is
!  used to control strong versus weak constraint so cannot be
!  arbitrarily defined.
!
      adFac=0.0_r8
#  ifdef IS4DVAR_SENSITIVITY
      IF ((mod(iic(ng)-1,nHIS(ng)).eq.0).and.                           &
     &    (DendS(ng).ge.tdays(ng)).and.(tdays(ng).ge.DstrS(ng))) THEN
#  else
      IF ((mod(iic(ng)-1,nTLM(ng)).eq.0).and.                           &
     &    (DendS(ng).ge.tdays(ng)).and.(tdays(ng).ge.DstrS(ng))) THEN
#  endif
          adFac=1.0_r8
          IF (Master) THEN
            WRITE (stdout,10) iic(ng)-1
 10         FORMAT (2x,' ADSEN_FORCE  - forcing Adjoint model at',      &
     &                 ' TimeStep: ', i8.8)
          END IF
      END IF
# else
      adFac=1.0_r8
# endif
!
!  Free-surface.
!
      IF (SCALARS(ng)%Lstate(isFsur)) THEN
# ifdef SOLVE3D
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                            &
     &                      adFac*zeta_ads(i,j)*Rscope(i,j)
          END DO
        END DO
# else
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Kfrc)=ad_zeta(i,j,Kfrc)+                        &
     &                        adFac*zeta_ads(i,j)*Rscope(i,j)
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
     &                        adFac*ubar_ads(i,j)*Uscope(i,j)
          END DO
        END DO
      END IF
!
      IF (SCALARS(ng)%Lstate(isVbar)) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Kfrc)=ad_vbar(i,j,Kfrc)+                        &
     &                        adFac*vbar_ads(i,j)*Vscope(i,j)
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
              ad_u(i,j,k,Nfrc)=ad_u(i,j,k,Nfrc)+                        &
     &                         adFac*u_ads(i,j,k)*Uscope(i,j)
            END DO
          END DO
        END DO
      END IF
!
      IF (SCALARS(ng)%Lstate(isVvel)) THEN
        DO k=KstrS(ng),KendS(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,Nfrc)=ad_v(i,j,k,Nfrc)+                        &
     &                         adFac*v_ads(i,j,k)*Vscope(i,j)
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
                ad_t(i,j,k,Nfrc,itrc)=ad_t(i,j,k,Nfrc,itrc)+            &
     &                                adFac*t_ads(i,j,k,itrc)*          &
     &                                Rscope(i,j)
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
