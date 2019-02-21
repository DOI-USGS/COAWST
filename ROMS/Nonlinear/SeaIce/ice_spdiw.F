#include "cppdefs.h"
      MODULE ice_spdiw_mod
#ifdef ICE_MODEL
!
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module computes the magnitude of the shear between the ice
!  and the surface water. In this case, the surface water is defined
!  as the water in a surface mixed layer, so that velocity must be
!  computed first.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_spdiw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_spdiw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_ice
      USE mod_coupling
#  ifdef LMD_SKPP
      USE mod_mixing
#  endif
      USE mod_stepping
!
      implicit none
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
# endif
      CALL ice_spdiw_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef ICE_MODEL
     &                     liuol(ng),                                   &
# endif
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
#ifdef LMD_SKPP
     &                     MIXING(ng) % hsbl,                           &
#endif
     &                     ICE(ng) % ui,                                &
     &                     ICE(ng) % vi,                                &
     &                     ICE(ng) % uwater,                            &
     &                     ICE(ng) % vwater,                            &
     &                     ICE(ng) % spd_iw                             &
     &                     )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
# endif
      RETURN
      END SUBROUTINE ice_spdiw
!
!***********************************************************************
      SUBROUTINE ice_spdiw_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef ICE_MODEL
     &                           liuol,                                 &
# endif
     &                           z_r, z_w,                              &
     &                           u, v,                                  &
# ifdef LMD_SKPP
     &                           hsbl,                                  &
# endif
     &                           ui, vi,                                &
     &                           uwater, vwater,                        &
     &                           spd_iw)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(in) :: liuol

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(out) :: uwater(LBi:,LBj:)
      real(r8), intent(out) :: vwater(LBi:,LBj:)
      real(r8), intent(out) :: spd_iw(LBi:,LBj:)
# else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: uwater(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vwater(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: spd_iw(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      integer :: nlio, nbotu, nbotv, k
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: nbot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vw

      real(r8) :: mlio
      real(r8) :: dml
      real(r8) :: totml

# include "set_bounds.h"

#  define I_RANGE MAX(Istr-2,0),MIN(Iend+2,Lm(ng)+1)
#  define J_RANGE MAX(Jstr-2,0),MIN(Jend+2,Mm(ng)+1)
      do j=J_RANGE
        do i=I_RANGE
!         sl_dpth = lmd_epsilon*(z_w(i,j,N(ng))-hsbl(i,j))
#  ifdef LMD_SKPP
! hsbl is now a positive quantity
          mlio = min(-hsbl(i,j),-10._r8)
#  else
          mlio = -10._r8
#  endif
          nbot(i,j) = 1
          do k=N(ng),1,-1
            if(z_r(i,j,k).lt.mlio) then
              nbot(i,j) = min(k,N(ng))
              nbot(i,j) = max(nbot(i,j),1)
              goto 1111
            endif
          enddo
 1111   continue
        enddo
      enddo
#undef I_RANGE
#undef J_RANGE
#  define I_RANGE MAX(Istr-1,1),Iend+1             
      do j=Jstr,Jend
        do i=I_RANGE
          nlio = 0
          nbotu = NINT(0.5_r8*(nbot(i-1,j)+nbot(i,j)))
          nbotu = max(min(nbotu,N(ng)),1)
          uw(i,j) = 0._r8
          totml = 0._r8
          do k=N(ng),nbotu,-1
            nlio = nlio + 1
            dml = 0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)               &
     &                      + z_w(i,j,k)-z_w(i,j,k-1))
            uw(i,j) = uw(i,j) + u(i,j,k,nrhs)*dml
            totml = totml + dml
          enddo
          uw(i,j) = uw(i,j)/totml
!         uw(i,j) =  u(i,j,N,nrhs)
        enddo
      enddo

#  define J_RANGE MAX(Jstr-1,1),Jend+1 
      do j=J_RANGE
        do i=Istr,Iend
          nlio = 0
          nbotv = NINT(0.5_r8*(nbot(i,j-1)+nbot(i,j)))
          nbotv = max(min(nbotv,N(ng)),1)
          vw(i,j) = 0._r8
          totml = 0._r8
          do k=N(ng),nbotv,-1
            nlio = nlio + 1
            dml = 0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)               &
     &                      + z_w(i,j,k)-z_w(i,j,k-1))
            vw(i,j) = vw(i,j) + v(i,j,k,nrhs)*dml
            totml = totml + dml
          enddo
          vw(i,j) = vw(i,j)/totml
!         vw(i,j) =  v(i,j,N,nrhs)
        enddo
      enddo
#undef I_RANGE
#undef J_RANGE
      do j=Jstr,Jend
        do i=Istr,Iend
          spd_iw(i,j) = 0.5*sqrt((uw(i,j)-ui(i,j,liuol)                 &
     &                 +  uw(i+1,j)-ui(i+1,j,liuol))**2                 &
     &                  +(vw(i,j)-vi(i,j,liuol)                         &
     &                 +  vw(i,j+1)-vi(i,j+1,liuol))**2)
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrP,Iend
           uwater(i,j) = uw(i,j)
        enddo
      enddo
      do j=JstrP,Jend
        do i=Istr,Iend
           vwater(i,j) = vw(i,j)
        enddo
      enddo
!
!  Apply boundary conditions.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          spd_iw)
        CALL bc_u2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          uwater)
        CALL bc_v2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vwater)
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    spd_iw, uwater, vwater)
#endif

      RETURN
      END SUBROUTINE ice_spdiw_tile
#endif
      END MODULE ice_spdiw_mod
