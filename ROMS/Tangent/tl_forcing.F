#include "cppdefs.h"
      MODULE tl_forcing_mod

#if defined TANGENT && (defined WEAK_CONSTRAINT || defined FORCING_SV)
!
!svn $Id: tl_forcing.F 751 2015-01-07 22:56:36Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group   Emanuele Di Lorenzo   !
!    Licensed under a MIT/X style license            Andrew M. Moore   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is used to force the tangent linear state equations.   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tl_forcing

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_forcing (ng, tile, Kfrc, Nfrc)
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
      CALL tl_forcing_tile (ng, tile,                                   &
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
     &                      OCEAN(ng) % tl_t,                           &
     &                      OCEAN(ng) % tl_u,                           &
     &                      OCEAN(ng) % tl_v,                           &
# else
     &                      OCEAN(ng) % tl_ubar,                        &
     &                      OCEAN(ng) % tl_vbar,                        &
# endif
# ifdef SOLVE3D
     &                      COUPLING(ng) % tl_Zt_avg1,                  &
# endif
     &                      OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_forcing
!
!***********************************************************************
      SUBROUTINE tl_forcing_tile (ng, tile,                             &
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
     &                            tl_t, tl_u, tl_v,                     &
# else
     &                            tl_ubar, tl_vbar,                     &
# endif
# ifdef SOLVE3D
     &                            tl_Zt_avg1,                           &
# endif
     &                            tl_zeta)
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
#  else
      real(r8), intent(in) :: f_ubar(LBi:,LBj:)
      real(r8), intent(in) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: f_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: f_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
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
!  Add given forcing to tangent linear state.
!-----------------------------------------------------------------------
!
# ifdef WEAK_CONSTRAINT
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        IF (Master) THEN
          IF (FrequentImpulse(ng)) THEN
            IF (iic(ng).eq.ntstart(ng)) THEN
              WRITE (stdout,10)
 10           FORMAT (3x,'TL_FORCING   - adding weak-constraint ',      &
     &                'forcing at each time-step.')
            END IF
          ELSE
            WRITE (stdout,20) time_code(ng)
 20         FORMAT (3x,'TL_FORCING   - added weak-constraint forcing,', &
     &              t62,'t = ', a)
          END IF
        END IF
      END IF
# endif
!
!  Tangent linear free-surface. The two different cases in the
!  case of SOLVE3D are due to the fact that tl_ini_fields is
!  also called on the first timestep. tl_forcing MUST be called
!  before tl_ini_fields.
!
# ifdef SOLVE3D
      IF(iic(ng).eq.ntstart(ng)) THEN
       DO j=JstrR,JendR
         DO i=IstrR,IendR
            tl_zeta(i,j,Kfrc)=tl_zeta(i,j,Kfrc)+f_zeta(i,j)
         END DO
       END DO
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_Zt_avg1(i,j)=tl_Zt_avg1(i,j)+f_zeta(i,j)
          END DO
        END DO
      END IF
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Kfrc)=tl_zeta(i,j,Kfrc)+f_zeta(i,j)
        END DO
      END DO
# endif

# ifndef SOLVE3D
!
!  Tangent linear 2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Kfrc)=tl_ubar(i,j,Kfrc)+f_ubar(i,j)
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Kfrc)=tl_vbar(i,j,Kfrc)+f_vbar(i,j)
        END DO
      END DO

# else
!
!  Tangent linear 3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_u(i,j,k,Nfrc)=tl_u(i,j,k,Nfrc)+f_u(i,j,k)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Nfrc)=tl_v(i,j,k,Nfrc)+f_v(i,j,k)
          END DO
        END DO
      END DO
!
!  Tangent linear tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_t(i,j,k,Nfrc,itrc)=tl_t(i,j,k,Nfrc,itrc)+              &
     &                              f_t(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE tl_forcing_tile
#endif
      END MODULE tl_forcing_mod
