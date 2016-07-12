#include "cppdefs.h"
      MODULE sum_imp_mod
#if defined WEAK_CONSTRAINT && defined RPCG
!
!svn $Id: sum_imp.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC sum_imp

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sum_imp (ng, tile, Lout)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL sum_imp_tile (ng, tile,                                      &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lout,                                         &
# ifdef SOLVE3D
     &                    OCEAN(ng) % tl_t,                             &
     &                    OCEAN(ng) % tl_u,                             &
     &                    OCEAN(ng) % tl_v,                             &
     &                    OCEAN(ng) % f_t,                              &
     &                    OCEAN(ng) % f_u,                              &
     &                    OCEAN(ng) % f_v,                              &
# else
     &                    OCEAN(ng) % tl_ubar,                          &
     &                    OCEAN(ng) % tl_vbar,                          &
     &                    OCEAN(ng) % f_ubar,                           &
     &                    OCEAN(ng) % f_vbar,                           &
# endif
     &                    OCEAN(ng) % tl_zeta,                          &
     &                    OCEAN(ng) % f_zeta)
      RETURN
      END SUBROUTINE sum_imp
!
!***********************************************************************
      SUBROUTINE sum_imp_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lout,                                   &
# ifdef SOLVE3D
     &                          tl_t, tl_u, tl_v,                       &
     &                          f_t, f_u, f_v,                          &
# else
     &                          tl_ubar, tl_vbar,                       &
     &                          f_ubar, f_vbar,                         &
# endif
     &                          tl_zeta, f_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
# else
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, ib, ir, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lout)=f_zeta(i,j)+tl_zeta(i,j,Lout)
        END DO
      END DO

# ifndef SOLVE3D
!
!  2D U-momentum component.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Lout)=f_ubar(i,j)+tl_ubar(i,j,Lout)
        END DO
      END DO
!
!  2D V-momentum.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lout)=f_vbar(i,j)+tl_vbar(i,j,Lout)
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_u(i,j,k,Lout)=f_u(i,j,k)+tl_u(i,j,k,Lout)
          END DO
        END DO
      END DO

!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lout)=f_v(i,j,k)+tl_v(i,j,k,Lout)
          END DO
        END DO
      END DO

!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
             tl_t(i,j,k,Lout,itrc)=f_t(i,j,k,itrc)+tl_t(i,j,k,Lout,itrc)
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE sum_imp_tile
#endif
      END MODULE sum_imp_mod
