#include "cppdefs.h"
      MODULE ad_force_dual_mod

#if defined ADJOINT && defined WEAK_CONSTRAINT
!
!svn $Id: ad_force_dual.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group   Emanuele Di Lorenzo   !
!    Licensed under a MIT/X style license            Andrew M. Moore   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is used to force the adjoint state equations.          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_force_dual

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_force_dual (ng, tile, Kfrc, Nfrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kfrc, Nfrc
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_force_dual_tile (ng, tile,                                &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Kfrc, Nfrc,                                 &
# ifdef SOLVE3D
     &                      OCEAN(ng) % f_t,                            &
     &                      OCEAN(ng) % f_u,                            &
     &                      OCEAN(ng) % f_v,                            &
# else
     &                      OCEAN(ng) % f_ubar,                         &
     &                      OCEAN(ng) % f_vbar,                         &
# endif
     &                      OCEAN(ng) % f_zeta,                         &
# ifdef SOLVE3D
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
# else
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
# endif
     &                      OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_force_dual
!
!***********************************************************************
      SUBROUTINE ad_force_dual_tile (ng, tile,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Kfrc, Nfrc,                           &
# ifdef SOLVE3D
     &                            f_t, f_u, f_v,                        &
# else
     &                            f_ubar, f_vbar,                       &
# endif
     &                            f_zeta,                               &
# ifdef SOLVE3D
     &                            ad_t, ad_u, ad_v,                     &
# else
     &                            ad_ubar, ad_vbar,                     &
# endif
     &                            ad_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kfrc
      integer, intent(in) :: Nfrc
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add given forcing to adjoint state.
!-----------------------------------------------------------------------
!
!  Adjoint free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Kfrc)=ad_zeta(i,j,Kfrc)+f_zeta(i,j)
          f_zeta(i,j)=0.0_r8
        END DO
      END DO

# ifndef SOLVE3D
!
!  Adjoint 2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Kfrc)=ad_ubar(i,j,Kfrc)+f_ubar(i,j)
          f_ubar(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Kfrc)=ad_vbar(i,j,Kfrc)+f_vbar(i,j)
          f_vbar(i,j)=0.0_r8
        END DO
      END DO
# else
!
!  Adjoint 3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Nfrc)=ad_u(i,j,k,Nfrc)+f_u(i,j,k)
            f_u(i,j,k)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Nfrc)=ad_v(i,j,k,Nfrc)+f_v(i,j,k)
            f_v(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
!
!  Adjoint tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Nfrc,itrc)=ad_t(i,j,k,Nfrc,itrc)+              &
     &                              f_t(i,j,k,itrc)
              f_t(i,j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE ad_force_dual_tile
#endif
      END MODULE ad_force_dual_mod
