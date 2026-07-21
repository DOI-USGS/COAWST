#include "cppdefs.h"
      MODULE wec_wave_mix_mod
#if defined SOLVE3D && defined WEC && defined WAVE_MIXING
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.md                            Nirnimesh Kumar    !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the wave induced mixing contribution to       !
!  vertical viscosity.                                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Uchiyama, Y., McWilliams, J.C., and Shchepetkin, A.F. (2010).       !
!  Wave current interacation in an oceanic circulation model with a    !
!  vortex-force formalism: Applications to surf zone, Ocean Modeling,  !
!  34, 16-35.                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
# if defined DIAGNOSTICS_UV
      USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY :  mp_exchange3d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: wec_wave_mix
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_wave_mix (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 43, __LINE__, MyFile)
# endif
      CALL wec_wave_mix_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),      &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        FORCES(ng) % Dissip_break,                &
     &                        FORCES(ng) % Dissip_wcap,                 &
# ifdef WEC_ROLLER
     &                        FORCES(ng) % Dissip_roller,               &
# endif
     &                        FORCES(ng) % Hwave,                       &
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % z_w,                           &
     &                        MIXING(ng) % Akv)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 43, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE wec_wave_mix
!
!***********************************************************************
      SUBROUTINE wec_wave_mix_tile (ng, tile, LBi, UBi, LBj, UBj, UBk,  &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Dissip_break,                       &
     &                              Dissip_wcap,                        &
# ifdef WEC_ROLLER
     &                              Dissip_roller,                      &
# endif
     &                              Hwave,                              &
     &                              Hz, z_r, z_w,                       &
     &                              Akv)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_wcap(LBi:,LBj:)
#  ifdef WEC_ROLLER
      real(r8), intent(in) :: Dissip_roller(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dissip_wcap(LBi:UBi,LBj:UBj)
#  ifdef WEC_ROLLER
      real(r8), intent(in) :: Dissip_roller(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, numits, it
      real(r8) :: fac1, fac2, cff1, cff2, cff3, cff4
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: sinb=0.1_r8
      real(r8), parameter :: Cb=0.002_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oDstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oroller

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute wave induced mixing contribution to vertical viscosity.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Compute total depth
!
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
          oDstp(i,j)=1.0_r8/Dstp(i,j)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          oroller(i,j)=0.0_r8
          gamr(i,j)=MIN(0.707_r8*Dstp(i,j)/                             &
     &              (1.25_r8*Hwave(i,j)+eps),1.0_r8)
          DO k=1,N(ng)
            cff2=(z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)*gamr(i,j)
            oroller(i,j)=oroller(i,j)+Hz(i,j,k)*                        &
     &                   COSH(2.0_r8*pi*cff2)
          END DO
          oroller(i,j)=1.0_r8/(oroller(i,j)+eps)
        END DO
      END DO   
! 
!  Compute contribution of wave breaking induced mixing to vertical
!  viscosity  (Akv).
!
      DO k=1,N(ng)
        fac1=0.5_r8*SQRT(2.0_r8)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            fac2=(z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)
            cff2=fac2*gamr(i,j)
            cff3=COSH(2.0_r8*pi*cff2)
!           cff1=Dissip_wcap(i,j)+Dissip_break(i,j)
            cff1=Dissip_break(i,j)
            cff4=((1.0_r8-wec_alpha(ng))*cff1)**(1.0_r8/3.0_r8)
# ifdef WEC_ROLLER
            cff4=cff4+(Dissip_roller(i,j))**(1.0_r8/3.0_r8)
# endif
            Akv(i,j,k)=Akv(i,j,k)+                                      &
     &                 Cb*cff4*fac1*Hwave(i,j)*                         &
     &                 Dstp(i,j)*cff3*oroller(i,j)
          END DO
        END DO
      END DO   
!
!  Exchange lateral boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
      END IF
# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akv)
# endif
!
      RETURN
      END SUBROUTINE wec_wave_mix_tile
#endif
      END MODULE wec_wave_mix_mod
