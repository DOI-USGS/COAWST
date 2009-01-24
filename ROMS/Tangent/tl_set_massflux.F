#include "cppdefs.h"
      MODULE tl_set_massflux_mod

#if defined TANGENT && defined SOLVE3D

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: tl_set_massflux.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine computes tangent linear horizontal mass fluxes,        !
!  Hz*u/n and Hz*v/m.                                                  !
!                                                                      !
!   BASIC STATE variables required:  Hz, u, v                          !
!   Dependend variables:  tl_Huon, tl_Hvom                             !
!   Independend variables:  tl_Hz, tl_u, tl_v                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tl_set_massflux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_set_massflux (ng, tile)
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
      CALL wclock_on (ng, iTLM, 12)
# endif
      CALL tl_set_massflux_tile (ng, tile,                              &
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
      CALL wclock_off (ng, iTLM, 12)
# endif
      RETURN
      END SUBROUTINE tl_set_massflux
!
!***********************************************************************
      SUBROUTINE tl_set_massflux_tile (ng, tile,                        &
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

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod
# endif
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
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
!>          Huon(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   &
!>   &                  on_u(i,j)
!>
            tl_Huon(i,j,k)=0.5_r8*on_u(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      tl_u(i,j,k,nrhs)+                           &
     &                      (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      u(i,j,k,nrhs))
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
     &                      u_stokes(i,j,k))
# endif
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
!>          Hvom(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   &
!>   &                  om_v(i,j)
!>
            tl_Hvom(i,j,k)=0.5_r8*om_v(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      tl_v(i,j,k,nrhs)+                           &
     &                      (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      v(i,j,k,nrhs))
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
     &                      v_stokes(i,j,k))
# endif
          END DO
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_u3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        Huon)
!>
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        tl_Huon)
!>    CALL exchange_v3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        Hvom)
!>
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        tl_Hvom)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_Huon, tl_Hvom)
#  endif
# endif
      RETURN
      END SUBROUTINE tl_set_massflux_tile
#endif
      END MODULE tl_set_massflux_mod

