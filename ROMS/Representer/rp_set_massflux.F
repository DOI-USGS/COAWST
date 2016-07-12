#include "cppdefs.h"
      MODULE rp_set_massflux_mod

#if defined TL_IOMS && defined SOLVE3D
!
!svn $Id: rp_set_massflux.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes representers tangent linear horizontal mass   !
!  fluxes, Hz*u/n and Hz*v/m.                                          !
!                                                                      !
!   BASIC STATE variables required:  Hz, u, v                          !
!   Dependend variables:  tl_Huon, tl_Hvom                             !
!   Independend variables:  tl_Hz, tl_u, tl_v                          !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rp_set_massflux
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_set_massflux (ng, tile)
!***********************************************************************
!
      USE mod_param
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
# ifdef PROFILE
      CALL wclock_on (ng, iRPM, 12)
# endif
      CALL rp_set_massflux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs(ng),                              &
     &                           OCEAN(ng) % u,                         &
     &                           OCEAN(ng) % v,                         &
     &                           OCEAN(ng) % tl_u,                      &
     &                           OCEAN(ng) % tl_v,                      &
# ifdef NEARSHORE_MELLOR
     &                           OCEAN(ng) % u_stokes,                  &
     &                           OCEAN(ng) % v_stokes,                  &
     &                           OCEAN(ng) % tl_u_stokes,               &
     &                           OCEAN(ng) % tl_v_stokes,               &
# endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % tl_Hz,                      &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % tl_Huon,                    &
     &                           GRID(ng) % tl_Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 12)
# endif

      RETURN
      END SUBROUTINE rp_set_massflux
!
!***********************************************************************
      SUBROUTINE rp_set_massflux_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 nrhs,                            &
     &                                 u, v,                            &
     &                                 tl_u, tl_v,                      &
# ifdef NEARSHORE_MELLOR
     &                                 u_stokes, v_stokes,              &
     &                                 tl_u_stokes, tl_v_stokes,        &
# endif
     &                                 Hz, tl_Hz,                       &
     &                                 om_v, on_u,                      &
     &                                 tl_Huon, tl_Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(out) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_Hvom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
!  Compute horizontal mass fluxes.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
!>          Huon(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   &
!>   &                  on_u(i,j)
!>
            tl_Huon(i,j,k)=0.5_r8*on_u(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      tl_u(i,j,k,nrhs)+                           &
     &                      (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      u(i,j,k,nrhs))-                             &
# ifdef TL_IOMS
     &                     0.5_r8*on_u(i,j)*                            &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)
# endif
# ifdef NEARSHORE_MELLOR
!>          Huon(i,j,k)=Huon(i,j,k)+                                    &
!>   &                  0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                 &
!>   &                  u_stokes(i,j,k)*on_u(i,j)
!>
            tl_Huon(i,j,k)=tl_Huon(i,j,k)+                              &
     &                     0.5_r8*on_u(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      tl_u_stokes(i,j,k)+                         &
     &                      (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      u_stokes(i,j,k))-                           &
#  ifdef TL_IOMS
     &                     0.5_r8*on_u(i,j)*                            &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*u_stokes(i,j,k)
#  endif
# endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
!>          Hvom(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   &
!>   &                  om_v(i,j)
!>
            tl_Hvom(i,j,k)=0.5_r8*om_v(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      tl_v(i,j,k,nrhs)+                           &
     &                      (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      v(i,j,k,nrhs))-                             &
#  ifdef TL_IOMS
     &                     0.5_r8*om_v(i,j)*                            &
     &                     (Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)
#  endif
# ifdef NEARSHORE_MELLOR
!>          Hvom(i,j,k)=Hvom(i,j,k)+                                    &
!>   &                  0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                 &
!>   &                  v_stokes(i,j,k)*om_v(i,j)
!>
            tl_Hvom(i,j,k)=tl_Hvom(i,j,k)+                              &
     &                     0.5_r8*om_v(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      tl_v_stokes(i,j,k)+                         &
     &                      (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      v_stokes(i,j,k))-                           &
#  ifdef TL_IOMS
     &                     0.5_r8*om_v(i,j)*                            &
     &                     (Hz(i,j,k)+Hz(i,j-1,k))*v_stokes(i,j,k)
#  endif
# endif
          END DO
        END DO
      END DO
!
!  Exchange boundary information.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          Huon)
!>
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_Huon)
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          Hvom)
!>
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_Hvom)
      END IF

# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    Huon, Hvom)
!>
      CALL mp_exchange3d (ng, tile, iRPM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_Huon, tl_Hvom)
# endif

      RETURN
      END SUBROUTINE rp_set_massflux_tile
#endif
      END MODULE rp_set_massflux_mod
