#include "cppdefs.h"
      MODULE forcing_mod

#if defined NLM_OUTER  || defined SENSITIVITY_4DVAR || \
    defined TL_W4DPSAS || defined W4DPSAS
!
!svn $Id: forcing.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is used to force the nonlinear state equations in      !
!  weak constraint variational data assimilation.                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: forcing

      CONTAINS
!
!***********************************************************************
      SUBROUTINE forcing (ng, tile, Kfrc, Nfrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
# ifdef SOLVE3D
      USE mod_coupling
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kfrc, Nfrc
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL forcing_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   Kfrc, Nfrc,                                    &
# ifdef SOLVE3D
     &                   OCEAN(ng) % f_t,                               &
     &                   OCEAN(ng) % f_u,                               &
     &                   OCEAN(ng) % f_v,                               &
# endif
     &                   OCEAN(ng) % f_ubar,                            &
     &                   OCEAN(ng) % f_vbar,                            &
     &                   OCEAN(ng) % f_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
# endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
# ifdef SOLVE3D
     &                   COUPLING(ng) % Zt_avg1,                        &
# endif
     &                   OCEAN(ng) % zeta)

      RETURN
      END SUBROUTINE forcing
!
!***********************************************************************
      SUBROUTINE forcing_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Kfrc, Nfrc,                              &
# ifdef SOLVE3D
     &                         f_t, f_u, f_v,                           &
# endif
     &                         f_ubar, f_vbar,                          &
     &                         f_zeta,                                  &
# ifdef SOLVE3D
     &                         t, u, v,                                 &
# endif
     &                         ubar, vbar,                              &
# ifdef SOLVE3D
     &                         Zt_avg1,                                 &
# endif
     &                         zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
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
      real(r8), intent(in) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: f_u(LBi:,LBj:,:)
      real(r8), intent(in) :: f_v(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: f_ubar(LBi:,LBj:)
      real(r8), intent(in) :: f_vbar(LBi:,LBj:)
      real(r8), intent(in) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: f_v(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
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
!  Add adjoint impulse forcing to nonlinear linear state.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        IF (Master) WRITE (stdout,10) time_code(ng)
 10     FORMAT (3x,'NL_FORCING   - added convolved adjoint impulse,',   &
     &          t62,'t = ', a)
      END IF
!
!  Free-surface.
!
# ifdef SOLVE3D
      IF (iic(ng).eq.ntstart(ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zeta(i,j,Kfrc)=zeta(i,j,Kfrc)+f_zeta(i,j)
          END DO
        END DO
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+f_zeta(i,j)
          END DO
        END DO
      END IF
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,Kfrc)=zeta(i,j,Kfrc)+f_zeta(i,j)
        END DO
      END DO
# endif

# ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,Kfrc)=ubar(i,j,Kfrc)+f_ubar(i,j)
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,Kfrc)=vbar(i,j,Kfrc)+f_vbar(i,j)
        END DO
      END DO

# else
!
!  Tangent linear 2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,1)=ubar(i,j,1)+f_ubar(i,j)
          ubar(i,j,2)=ubar(i,j,2)+f_ubar(i,j)
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,1)=vbar(i,j,1)+f_vbar(i,j)
          vbar(i,j,2)=vbar(i,j,2)+f_vbar(i,j)
        END DO
      END DO
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            u(i,j,k,Nfrc)=u(i,j,k,Nfrc)+f_u(i,j,k)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,Nfrc)=v(i,j,k,Nfrc)+f_v(i,j,k)
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
              t(i,j,k,Nfrc,itrc)=t(i,j,k,Nfrc,itrc)+                    &
     &                           f_t(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE forcing_tile
#endif
      END MODULE forcing_mod
