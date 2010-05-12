#include "cppdefs.h"
      MODULE ad_obc_adjust_mod

#ifdef ADJUST_BOUNDARY
!
!svn $Id: ad_obc_adjust.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the adjoint of time-interpolation of 4DVar open     !
!  boundary increments.  The increments can be constant  (Nbrec=1)     !
!  or time interpolated between snapshots (Nbrec>1).                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      4DVar increment time index to process.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ad_obc_adjust
# ifdef SOLVE3D
      PUBLIC :: ad_obc2d_adjust
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_obc_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_forces
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
# endif
      CALL ad_obc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp,                                    &
# ifdef MASKING
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef SOLVE3D
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % Hz_bry,                       &
     &                         GRID(ng) % ad_Hz,                        &
     &                         GRID(ng) % ad_Hz_bry,                    &
# endif
# ifdef WEST_FSOBC
     &                         BOUNDARY(ng) % zeta_west,                &
     &                         BOUNDARY(ng) % ad_zeta_west,             &
# endif
# ifdef EAST_FSOBC
     &                         BOUNDARY(ng) % zeta_east,                &
     &                         BOUNDARY(ng) % ad_zeta_east,             &
# endif
# ifdef SOUTH_FSOBC
     &                         BOUNDARY(ng) % zeta_south,               &
     &                         BOUNDARY(ng) % ad_zeta_south,            &
# endif
# ifdef NORTH_FSOBC
     &                         BOUNDARY(ng) % zeta_north,               &
     &                         BOUNDARY(ng) % ad_zeta_north,            &
# endif
# ifdef WEST_M2OBC
     &                         BOUNDARY(ng) % ubar_west,                &
     &                         BOUNDARY(ng) % vbar_west,                &
     &                         BOUNDARY(ng) % ad_ubar_west,             &
     &                         BOUNDARY(ng) % ad_vbar_west,             &
# endif
# ifdef EAST_M2OBC
     &                         BOUNDARY(ng) % ubar_east,                &
     &                         BOUNDARY(ng) % vbar_east,                &
     &                         BOUNDARY(ng) % ad_ubar_east,             &
     &                         BOUNDARY(ng) % ad_vbar_east,             &
# endif
# ifdef SOUTH_M2OBC
     &                         BOUNDARY(ng) % ubar_south,               &
     &                         BOUNDARY(ng) % vbar_south,               &
     &                         BOUNDARY(ng) % ad_ubar_south,            &
     &                         BOUNDARY(ng) % ad_vbar_south,            &
# endif
# ifdef NORTH_M2OBC
     &                         BOUNDARY(ng) % ubar_north,               &
     &                         BOUNDARY(ng) % vbar_north,               &
     &                         BOUNDARY(ng) % ad_ubar_north,            &
     &                         BOUNDARY(ng) % ad_vbar_north,            &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                         BOUNDARY(ng) % u_west,                   &
     &                         BOUNDARY(ng) % v_west,                   &
     &                         BOUNDARY(ng) % ad_u_west,                &
     &                         BOUNDARY(ng) % ad_v_west,                &
#  endif
#  ifdef EAST_M3OBC
     &                         BOUNDARY(ng) % u_east,                   &
     &                         BOUNDARY(ng) % v_east,                   &
     &                         BOUNDARY(ng) % ad_u_east,                &
     &                         BOUNDARY(ng) % ad_v_east,                &
#  endif
#  ifdef SOUTH_M3OBC
     &                         BOUNDARY(ng) % u_south,                  &
     &                         BOUNDARY(ng) % v_south,                  &
     &                         BOUNDARY(ng) % ad_u_south,               &
     &                         BOUNDARY(ng) % ad_v_south,               &
#  endif
#  ifdef NORTH_M3OBC
     &                         BOUNDARY(ng) % u_north,                  &
     &                         BOUNDARY(ng) % v_north,                  &
     &                         BOUNDARY(ng) % ad_u_north,               &
     &                         BOUNDARY(ng) % ad_v_north,               &
#  endif
#  ifdef WEST_TOBC
     &                         BOUNDARY(ng) % ad_t_west,                &
#  endif
#  ifdef EAST_TOBC
     &                         BOUNDARY(ng) % ad_t_east,                &
#  endif
#  ifdef SOUTH_TOBC
     &                         BOUNDARY(ng) % ad_t_south,               &
#  endif
#  ifdef NORTH_TOBC
     &                         BOUNDARY(ng) % ad_t_north,               &
#  endif
# endif
# ifdef SOLVE3D
     &                         BOUNDARY(ng) % ad_t_obc,                 &
     &                         BOUNDARY(ng) % ad_u_obc,                 &
     &                         BOUNDARY(ng) % ad_v_obc,                 &
# endif
     &                         BOUNDARY(ng) % ad_ubar_obc,              &
     &                         BOUNDARY(ng) % ad_vbar_obc,              &
     &                         BOUNDARY(ng) % ad_zeta_obc)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
# endif

      RETURN
      END SUBROUTINE ad_obc_adjust
!
!***********************************************************************
      SUBROUTINE ad_obc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp,                              &
# ifdef MASKING
     &                               umask, vmask,                      &
# endif
# ifdef SOLVE3D
     &                               Hz, Hz_bry,                        &
     &                               ad_Hz, ad_Hz_bry,                  &
# endif
# ifdef WEST_FSOBC
     &                               zeta_west, ad_zeta_west,           &
# endif
# ifdef EAST_FSOBC
     &                               zeta_east, ad_zeta_east,           &
# endif
# ifdef SOUTH_FSOBC
     &                               zeta_south, ad_zeta_south,         &
# endif
# ifdef NORTH_FSOBC
     &                               zeta_north, ad_zeta_north,         &
# endif
# ifdef WEST_M2OBC
     &                               ubar_west, vbar_west,              &
     &                               ad_ubar_west, ad_vbar_west,        &
# endif
# ifdef EAST_M2OBC
     &                               ubar_east, vbar_east,              &
     &                               ad_ubar_east, ad_vbar_east,        &
# endif
# ifdef SOUTH_M2OBC
     &                               ubar_south, vbar_south,            &
     &                               ad_ubar_south, ad_vbar_south,      &
# endif
# ifdef NORTH_M2OBC
     &                               ubar_north, vbar_north,            &
     &                               ad_ubar_north, ad_vbar_north,      &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                               u_west, v_west,                    &
     &                               ad_u_west, ad_v_west,              &
#  endif
#  ifdef EAST_M3OBC
     &                               u_east, v_east,                    &
     &                               ad_u_east, ad_v_east,              &
#  endif
#  ifdef SOUTH_M3OBC
     &                               u_south, v_south,                  &
     &                               ad_u_south, ad_v_south,            &
#  endif
#  ifdef NORTH_M3OBC
     &                               u_north, v_north,                  &
     &                               ad_u_north, ad_v_north,            &
#  endif
#  ifdef WEST_TOBC
     &                               ad_t_west,                         &
#  endif
#  ifdef EAST_TOBC
     &                               ad_t_east,                         &
#  endif
#  ifdef SOUTH_TOBC
     &                               ad_t_south,                        &
#  endif
#  ifdef NORTH_TOBC
     &                               ad_t_north,                        &
#  endif
# endif
# ifdef SOLVE3D
     &                               ad_t_obc, ad_u_obc, ad_v_obc,      &
# endif
     &                               ad_ubar_obc, ad_vbar_obc,          &
     &                               ad_zeta_obc)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz_bry(LBij:,:,:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hz_bry(LBij:,:,:)
#  endif
#  ifdef WEST_FSOBC
      real(r8), intent(in) :: zeta_west(0:)
      real(r8), intent(inout) :: ad_zeta_west(0:)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(in) :: zeta_east(0:)
      real(r8), intent(inout) :: ad_zeta_east(0:)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(in) :: zeta_south(0:)
      real(r8), intent(inout) :: ad_zeta_south(0:)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(in) :: zeta_north(0:)
      real(r8), intent(inout) :: ad_zeta_north(0:)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(in) :: ubar_west(0:)
      real(r8), intent(in) :: vbar_west(0:)
      real(r8), intent(inout) :: ad_ubar_west(0:)
      real(r8), intent(inout) :: ad_vbar_west(0:)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(in) :: ubar_east(0:)
      real(r8), intent(in) :: vbar_east(0:)
      real(r8), intent(inout) :: ad_ubar_east(0:)
      real(r8), intent(inout) :: ad_vbar_east(0:)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ubar_south(0:)
      real(r8), intent(in) :: vbar_south(0:)
      real(r8), intent(inout) :: ad_ubar_south(0:)
      real(r8), intent(inout) :: ad_vbar_south(0:)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(in) :: ubar_north(0:)
      real(r8), intent(in) :: vbar_north(0:)
      real(r8), intent(inout) :: ad_ubar_north(0:)
      real(r8), intent(inout) :: ad_vbar_north(0:)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(inout) :: ad_u_west(0:,:)
      real(r8), intent(inout) :: ad_v_west(0:,:)
      real(r8), intent(in) :: u_west(0:,:)
      real(r8), intent(in) :: v_west(0:,:)
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: u_east(0:,:)
      real(r8), intent(in) :: v_east(0:,:)
      real(r8), intent(inout) :: ad_u_east(0:,:)
      real(r8), intent(inout) :: ad_v_east(0:,:)
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: u_south(0:,:)
      real(r8), intent(in) :: v_south(0:,:)
      real(r8), intent(inout) :: ad_u_south(0:,:)
      real(r8), intent(inout) :: ad_v_south(0:,:)
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: u_north(0:,:)
      real(r8), intent(in) :: v_north(0:,:)
      real(r8), intent(inout) :: ad_u_north(0:,:)
      real(r8), intent(inout) :: ad_v_north(0:,:)
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(inout) :: ad_t_west(0:,:,:)
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(inout) :: ad_t_east(0:,:,:)
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(inout) :: ad_t_south(0:,:,:)
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(inout) :: ad_t_north(0:,:,:)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz_bry(LBij:UBij,N(ng),4)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hz_bry(LBij:UBij,N(ng),4)
#  endif
#  ifdef WEST_FSOBC
      real(r8), intent(in) :: zeta_west(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_zeta_west(0:Jm(ng)+1)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(in) :: zeta_east(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_zeta_east(0:Jm(ng)+1)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(in) :: zeta_south(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_zeta_south(0:Im(ng)+1)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(inout) :: ad_zeta_north(0:Im(ng)+1)
      real(r8), intent(in) :: zeta_north(0:Im(ng)+1)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(in) :: ubar_west(0:Jm(ng)+1)
      real(r8), intent(in) :: vbar_west(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_ubar_west(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_vbar_west(0:Jm(ng)+1)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(in) :: ubar_east(0:Jm(ng)+1)
      real(r8), intent(in) :: vbar_east(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_ubar_east(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_vbar_east(0:Jm(ng)+1)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ubar_south(0:Im(ng)+1)
      real(r8), intent(in) :: vbar_south(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_ubar_south(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_vbar_south(0:Im(ng)+1)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(in) :: ubar_north(0:Im(ng)+1)
      real(r8), intent(in) :: vbar_north(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_ubar_north(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_vbar_north(0:Im(ng)+1)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(in) :: u_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(in) :: v_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_west(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: u_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(in) :: v_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_east(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: u_south(0:Im(ng)+1,N(ng))
      real(r8), intent(in) :: v_south(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_south(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_south(0:Im(ng)+1,N(ng))
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: u_north(0:Im(ng)+1,N(ng))
      real(r8), intent(in) :: v_north(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_north(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_north(0:Im(ng)+1,N(ng))
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(inout) :: ad_t_west(0:Jm(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(inout) :: ad_t_east(0:Jm(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(inout) :: ad_t_south(0:Im(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(inout) :: ad_t_north(0:Im(ng)+1,N(ng),NT(ng))
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, ib, it1, it2, j
# ifdef SOLVE3D
      integer :: it, k
# endif
      real(r8) :: adfac, fac, fac1, fac2
      real(r8) :: cff1, ad_cff1, ad_cff2

# ifdef SOLVE3D
      real(r8), dimension(0:N(ng)) :: CF
      real(r8), dimension(0:N(ng)) :: DC

      real(r8), dimension(0:N(ng)) :: ad_CF
      real(r8), dimension(0:N(ng)) :: ad_DC
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint of adjust open boundary fields with 4DVAR increments.
!-----------------------------------------------------------------------
!
!  Set time records and interpolation factor, if any.
!
      IF (Nbrec(ng).eq.1) THEN
        it1=1
        it2=1
        fac1=1.0_r8
        fac2=0.0_r8
      ELSE
        it1=MAX(0,(iic(ng)-1)/nOBC(ng))+1
        it2=MIN(it1+1,Nbrec(ng))
        fac1=OBC_time(it2,ng)-(time(ng)+dt(ng))
        fac2=(time(ng)+dt(ng))-OBC_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF
!
!   Clear arrays and constants
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_CF=0.0_r8
      ad_DC=0.0_r8

# ifndef SOLVE3D
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
#  ifdef WEST_M2OBC
        IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
          DO j=Jstr,Jend
!>          tl_ubar_west(j)=fac1*tl_ubar_obc(j,iwest,it1,Linp)+         &
!>   &                      fac2*tl_ubar_obc(j,iwest,it2,Linp)
!>
            ad_ubar_obc(j,iwest,it1,Linp)=ad_ubar_obc(j,iwest,it1,Linp)+&
     &                                    fac1*ad_ubar_west(j)
            ad_ubar_obc(j,iwest,it2,Linp)=ad_ubar_obc(j,iwest,it2,Linp)+&
     &                                    fac2*ad_ubar_west(j)
            ad_ubar_west(j)=0.0_r8
          END DO
        END IF
#  endif
#  ifdef EAST_M2OBC
        IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
          DO j=Jstr,Jend
!>          tl_ubar_east(j)=fac1*tl_ubar_obc(j,ieast,it1,Linp)+         &
!>   &                      fac2*tl_ubar_obc(j,ieast,it2,Linp)
!>
            ad_ubar_obc(j,ieast,it1,Linp)=ad_ubar_obc(j,ieast,it1,Linp)+&
     &                       fac1*ad_ubar_east(j)
            ad_ubar_obc(j,ieast,it2,Linp)=ad_ubar_obc(j,ieast,it2,Linp)+&
     &                       fac2*ad_ubar_east(j)
            ad_ubar_east(j)=0.0_r8
          END DO
        END IF
#  endif
#  ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
          DO i=IstrU,Iend
!>          tl_ubar_south(i)=fac1*tl_ubar_obc(i,isouth,it1,Linp)+       &
!>   &                       fac2*tl_ubar_obc(i,isouth,it2,Linp)
!>
            ad_ubar_obc(i,isouth,it1,Linp)=                             &
     &          ad_ubar_obc(i,isouth,it1,Linp)+fac1*ad_ubar_south(i)
            ad_ubar_obc(i,isouth,it2,Linp)=                             &
     &          ad_ubar_obc(i,isouth,it2,Linp)+fac2*ad_ubar_south(i)
            ad_ubar_south(i)=0.0_r8
          END DO
        END IF
#  endif
#  ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
          DO i=IstrU,Iend
!>          tl_ubar_north(i)=fac1*tl_ubar_obc(i,inorth,it1,Linp)+       &
!>   &                       fac2*tl_ubar_obc(i,inorth,it2,Linp)
!>
            ad_ubar_obc(i,inorth,it1,Linp)=                             &
     &            ad_ubar_obc(i,inorth,it1,Linp)+fac1*ad_ubar_north(i)
            ad_ubar_obc(i,inorth,it2,Linp)=                             &
     &            ad_ubar_obc(i,inorth,it2,Linp)+fac2*ad_ubar_north(i)
            ad_ubar_north(i)=0.0_r8
          END DO
        END IF
#  endif
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
#  ifdef WEST_M2OBC
        IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
          DO j=JstrV,Jend
!>          tl_vbar_west(j)=fac1*tl_vbar_obc(j,iwest,it1,Linp)+         &
!>   &                      fac2*tl_vbar_obc(j,iwest,it2,Linp)
!>
            ad_vbar_obc(j,iwest,it1,Linp)=                              &
     &         ad_vbar_obc(j,iwest,it1,Linp)+fac1*ad_vbar_west(j)
            ad_vbar_obc(j,iwest,it2,Linp)=                              &
     &         ad_vbar_obc(j,iwest,it2,Linp)+fac2*ad_vbar_west(j)
            ad_vbar_west(j)=0.0_r8
          END DO
        END IF
#  endif
#  ifdef EAST_M2OBC
        IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
          DO j=JstrV,Jend
!>          tl_vbar_east(j)=fac1*tl_vbar_obc(j,ieast,it1,Linp)+         &
!>   &                      fac2*tl_vbar_obc(j,ieast,it2,Linp)
!>
            ad_vbar_obc(j,ieast,it1,Linp)=                              &
     &               ad_vbar_obc(j,ieast,it1,Linp)+fac1*ad_vbar_east(j)
            ad_vbar_obc(j,ieast,it2,Linp)=                              &
     &               ad_vbar_obc(j,ieast,it2,Linp)+fac2*ad_vbar_east(j)
            ad_vbar_east(j)=0.0_r8
          END DO
        END IF
#  endif
#  ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
!>          tl_vbar_south(i)=fac1*tl_vbar_obc(i,isouth,it1,Linp)+       &
!>   &                       fac2*tl_vbar_obc(i,isouth,it2,Linp)
!>
            ad_vbar_obc(i,isouth,it1,Linp)=                             &
     &            ad_vbar_obc(i,isouth,it1,Linp)+fac1*ad_vbar_south(i)
            ad_vbar_obc(i,isouth,it2,Linp)=                             &
     &            ad_vbar_obc(i,isouth,it2,Linp)+fac2*ad_vbar_south(i)
            ad_vbar_south(i)=0.0_r8
          END DO
        END IF
#  endif
#  ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
          DO i=Istr,Iend
!>          tl_vbar_north(i)=fac1*tl_vbar_obc(i,inorth,it1,Linp)+       &
!>   &                       fac2*tl_vbar_obc(i,inorth,it2,Linp)
!>
            ad_vbar_obc(i,inorth,it1,Linp)=                             &
     &             ad_vbar_obc(i,inorth,it1,Linp)+fac1*ad_vbar_north(i)
            ad_vbar_obc(i,inorth,it2,Linp)=                             &
     &             ad_vbar_obc(i,inorth,it2,Linp)+fac2*ad_vbar_north(i)
            ad_vbar_north(i)=0.0_r8
          END DO
        END IF
#  endif
      END IF
# endif
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
# ifdef WEST_FSOBC
        IF ((Lobc(iwest,isFsur,ng)).and.WESTERN_EDGE) THEN
          ib=iwest
          DO j=Jstr,Jend
!>          tl_zeta_west(j)=fac1*tl_zeta_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_zeta_obc(j,ib,it2,Linp)
!>
            ad_zeta_obc(j,ib,it1,Linp)=ad_zeta_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_west(j)
            ad_zeta_obc(j,ib,it2,Linp)=ad_zeta_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_west(j)
            ad_zeta_west(j)=0.0_r8
          END DO
        END IF
# endif
# ifdef EAST_FSOBC
        IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO j=Jstr,Jend
!>          tl_zeta_east(j)=fac1*tl_zeta_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_zeta_obc(j,ib,it2,Linp)
!>
            ad_zeta_obc(j,ib,it1,Linp)=ad_zeta_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_east(j)
            ad_zeta_obc(j,ib,it2,Linp)=ad_zeta_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_east(j)
            ad_zeta_east(j)=0.0_r8
          END DO
        END IF
# endif
# ifdef SOUTH_FSOBC
        IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO i=Istr,Iend
!>          tl_zeta_south(i)=fac1*tl_zeta_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_zeta_obc(i,ib,it2,Linp)
!>
            ad_zeta_obc(i,ib,it1,Linp)=ad_zeta_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_south(i)
            ad_zeta_obc(i,ib,it2,Linp)=ad_zeta_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_south(i)
            ad_zeta_south(i)=0.0_r8
          END DO
        END IF
# endif
# ifdef NORTH_FSOBC
        IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO i=Istr,Iend
!>          tl_zeta_north(i)=fac1*tl_zeta_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_zeta_obc(i,ib,it2,Linp)
!>
            ad_zeta_obc(i,ib,it1,Linp)=ad_zeta_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_north(i)
            ad_zeta_obc(i,ib,it2,Linp)=ad_zeta_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_north(i)
            ad_zeta_north(i)=0.0_r8
          END DO
        END IF
# endif
      END IF

# ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
#  ifdef WEST_M3OBC
        IF ((Lobc(iwest,isUvel,ng)).and.WESTERN_EDGE) THEN
          ib=iwest
          DO k=1,N(ng)
            DO j=Jstr,Jend
!>            tl_u_west(j,k)=fac1*tl_u_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_u_obc(j,k,ib,it2,Linp)
!>
              ad_u_obc(j,k,ib,it1,Linp)=ad_u_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_west(j,k)
              ad_u_obc(j,k,ib,it2,Linp)=ad_u_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_west(j,k)
              ad_u_west(j,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO k=1,N(ng)
            DO j=Jstr,Jend
!>            tl_u_east(j,k)=fac1*tl_u_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_u_obc(j,k,ib,it2,Linp)
!>
              ad_u_obc(j,k,ib,it1,Linp)=ad_u_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_east(j,k)
              ad_u_obc(j,k,ib,it2,Linp)=ad_u_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_east(j,k)
              ad_u_east(j,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO k=1,N(ng)
            DO i=IstrU,Iend
!>            tl_u_south(i,k)=fac1*tl_u_obc(i,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_u_obc(i,k,ib,it2,Linp)
!>
              ad_u_obc(i,k,ib,it1,Linp)=ad_u_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_south(i,k)
              ad_u_obc(i,k,ib,it2,Linp)=ad_u_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_south(i,k)
              ad_u_south(i,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO k=1,N(ng)
            DO i=IstrU,Iend
!>            tl_u_north(i,k)=fac1*tl_u_obc(i,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_u_obc(i,k,ib,it2,Linp)
!>
              ad_u_obc(i,k,ib,it1,Linp)=ad_u_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_north(i,k)
              ad_u_obc(i,k,ib,it2,Linp)=ad_u_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_north(i,k)
              ad_u_north(i,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
      END IF
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
#  ifdef WEST_M3OBC
        IF ((Lobc(iwest,isVvel,ng)).and.WESTERN_EDGE) THEN
          ib=iwest
          DO k=1,N(ng)
            DO j=JstrV,Jend
!>            tl_v_west(j,k)=fac1*tl_v_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_v_obc(j,k,ib,it2,Linp)
!>
              ad_v_obc(j,k,ib,it1,Linp)=ad_v_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_west(j,k)
              ad_v_obc(j,k,ib,it2,Linp)=ad_v_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_west(j,k)
              ad_v_west(j,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO k=1,N(ng)
            DO j=JstrV,Jend
!>            tl_v_east(j,k)=fac1*tl_v_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_v_obc(j,k,ib,it2,Linp)
!>
              ad_v_obc(j,k,ib,it1,Linp)=ad_v_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_east(j,k)
              ad_v_obc(j,k,ib,it2,Linp)=ad_v_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_east(j,k)
              ad_v_east(j,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_v_south(i,k)=fac1*tl_v_obc(i,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_v_obc(i,k,ib,it2,Linp)
!>
              ad_v_obc(i,k,ib,it1,Linp)=ad_v_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_south(i,k)
              ad_v_obc(i,k,ib,it2,Linp)=ad_v_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_south(i,k)
              ad_v_south(i,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_v_north(i,k)=fac1*tl_v_obc(i,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_v_obc(i,k,ib,it2,Linp)
!>
              ad_v_obc(i,k,ib,it1,Linp)=ad_v_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_north(i,k)
              ad_v_obc(i,k,ib,it2,Linp)=ad_v_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_north(i,k)
              ad_v_north(i,k)=0.0_r8
            END DO
          END DO
        END IF
#  endif
      END IF
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
#  ifdef WEST_TOBC
          IF ((Lobc(iwest,isTvar(it),ng)).and.WESTERN_EDGE) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_t_west(j,k,it)=fac1*tl_t_obc(j,k,ib,it1,Linp,it)+    &
!>   &                            fac2*tl_t_obc(j,k,ib,it2,Linp,it)
!>
                ad_t_obc(j,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_west(j,k,it)
                ad_t_obc(j,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_west(j,k,it)
                ad_t_west(j,k,it)=0.0_r8
              END DO
            END DO
          END IF
#  endif
#  ifdef EAST_TOBC
          IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_t_east(j,k,it)=fac1*tl_t_obc(j,k,ib,it1,Linp,it)+    &
!>   &                            fac2*tl_t_obc(j,k,ib,it2,Linp,it)
!>
                ad_t_obc(j,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_east(j,k,it)
                ad_t_obc(j,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_east(j,k,it)
                ad_t_east(j,k,it)=0.0_r8
              END DO
            END DO
          END IF
#  endif
#  ifdef SOUTH_TOBC
          IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_t_south(i,k,it)=fac1*tl_t_obc(i,k,ib,it1,Linp,it)+   &
!>   &                             fac2*tl_t_obc(i,k,ib,it2,Linp,it)
!>
                ad_t_obc(i,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(i,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_south(i,k,it)
                ad_t_obc(i,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(i,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_south(i,k,it)
                ad_t_south(i,k,it)=0.0_r8
              END DO
            END DO
          END IF
#  endif
#  ifdef NORTH_TOBC
          IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_t_north(i,k,it)=fac1*tl_t_obc(i,k,ib,it1,Linp,it)+   &
!>   &                             fac2*tl_t_obc(i,k,ib,it2,Linp,it)
!>
                ad_t_obc(i,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(i,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_north(i,k,it)
                ad_t_obc(i,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(i,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_north(i,k,it)
                ad_t_north(i,k,it)=0.0_r8
              END DO
            END DO
          END IF
#  endif
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE ad_obc_adjust_tile

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE ad_obc2d_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_forces
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      CALL ad_obc2d_adjust_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Linp,                                  &
#  ifdef MASKING
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % vmask,                      &
#  endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % Hz_bry,                     &
     &                           GRID(ng) % ad_Hz,                      &
     &                           GRID(ng) % ad_Hz_bry,                  &
#  ifdef WEST_M2OBC
     &                           BOUNDARY(ng) % ubar_west,              &
     &                           BOUNDARY(ng) % vbar_west,              &
     &                           BOUNDARY(ng) % ad_ubar_west,           &
     &                           BOUNDARY(ng) % ad_vbar_west,           &
#  endif
#  ifdef EAST_M2OBC
     &                           BOUNDARY(ng) % ubar_east,              &
     &                           BOUNDARY(ng) % vbar_east,              &
     &                           BOUNDARY(ng) % ad_ubar_east,           &
     &                           BOUNDARY(ng) % ad_vbar_east,           &
#  endif
#  ifdef SOUTH_M2OBC
     &                           BOUNDARY(ng) % ubar_south,             &
     &                           BOUNDARY(ng) % vbar_south,             &
     &                           BOUNDARY(ng) % ad_ubar_south,          &
     &                           BOUNDARY(ng) % ad_vbar_south,          &
#  endif
#  ifdef NORTH_M2OBC
     &                           BOUNDARY(ng) % ubar_north,             &
     &                           BOUNDARY(ng) % vbar_north,             &
     &                           BOUNDARY(ng) % ad_ubar_north,          &
     &                           BOUNDARY(ng) % ad_vbar_north,          &
#  endif
#  ifdef WEST_M3OBC
     &                           BOUNDARY(ng) % u_west,                 &
     &                           BOUNDARY(ng) % v_west,                 &
     &                           BOUNDARY(ng) % ad_u_west,              &
     &                           BOUNDARY(ng) % ad_v_west,              &
#  endif
#  ifdef EAST_M3OBC
     &                           BOUNDARY(ng) % u_east,                 &
     &                           BOUNDARY(ng) % v_east,                 &
     &                           BOUNDARY(ng) % ad_u_east,              &
     &                           BOUNDARY(ng) % ad_v_east,              &
#  endif
#  ifdef SOUTH_M3OBC
     &                           BOUNDARY(ng) % u_south,                &
     &                           BOUNDARY(ng) % v_south,                &
     &                           BOUNDARY(ng) % ad_u_south,             &
     &                           BOUNDARY(ng) % ad_v_south,             &
#  endif
#  ifdef NORTH_M3OBC
     &                           BOUNDARY(ng) % u_north,                &
     &                           BOUNDARY(ng) % v_north,                &
     &                           BOUNDARY(ng) % ad_u_north,             &
     &                           BOUNDARY(ng) % ad_v_north,             &
#  endif
     &                           BOUNDARY(ng) % ad_u_obc,               &
     &                           BOUNDARY(ng) % ad_v_obc)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif

      RETURN
      END SUBROUTINE ad_obc2d_adjust
!
!***********************************************************************
      SUBROUTINE ad_obc2d_adjust_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBij, UBij,  &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Linp,                            &
#  ifdef MASKING
     &                                 umask, vmask,                    &
#  endif
     &                                 Hz, Hz_bry,                      &
     &                                 ad_Hz, ad_Hz_bry,                &
#  ifdef WEST_M2OBC
     &                                 ubar_west, vbar_west,            &
     &                                 ad_ubar_west, ad_vbar_west,      &
#  endif
#  ifdef EAST_M2OBC
     &                                 ubar_east, vbar_east,            &
     &                                 ad_ubar_east, ad_vbar_east,      &
#  endif
#  ifdef SOUTH_M2OBC
     &                                 ubar_south, vbar_south,          &
     &                                 ad_ubar_south, ad_vbar_south,    &
#  endif
#  ifdef NORTH_M2OBC
     &                                 ubar_north, vbar_north,          &
     &                                 ad_ubar_north, ad_vbar_north,    &
#  endif
#  ifdef WEST_M3OBC
     &                                 u_west, v_west,                  &
     &                                 ad_u_west, ad_v_west,            &
#  endif
#  ifdef EAST_M3OBC
     &                                 u_east, v_east,                  &
     &                                 ad_u_east, ad_v_east,            &
#  endif
#  ifdef SOUTH_M3OBC
     &                                 u_south, v_south,                &
     &                                 ad_u_south, ad_v_south,          &
#  endif
#  ifdef NORTH_M3OBC
     &                                 u_north, v_north,                &
     &                                 ad_u_north, ad_v_north,          &
#  endif
     &                                 ad_u_obc, ad_v_obc)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz_bry(LBij:,:,:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hz_bry(LBij:,:,:)
#   ifdef WEST_M2OBC
      real(r8), intent(in) :: ubar_west(0:)
      real(r8), intent(in) :: vbar_west(0:)
      real(r8), intent(inout) :: ad_ubar_west(0:)
      real(r8), intent(inout) :: ad_vbar_west(0:)
#   endif
#   ifdef EAST_M2OBC
      real(r8), intent(in) :: ubar_east(0:)
      real(r8), intent(in) :: vbar_east(0:)
      real(r8), intent(inout) :: ad_ubar_east(0:)
      real(r8), intent(inout) :: ad_vbar_east(0:)
#   endif
#   ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ubar_south(0:)
      real(r8), intent(in) :: vbar_south(0:)
      real(r8), intent(inout) :: ad_ubar_south(0:)
      real(r8), intent(inout) :: ad_vbar_south(0:)
#   endif
#   ifdef NORTH_M2OBC
      real(r8), intent(in) :: ubar_north(0:)
      real(r8), intent(in) :: vbar_north(0:)
      real(r8), intent(inout) :: ad_ubar_north(0:)
      real(r8), intent(inout) :: ad_vbar_north(0:)
#   endif
#   ifdef WEST_M3OBC
      real(r8), intent(inout) :: ad_u_west(0:,:)
      real(r8), intent(inout) :: ad_v_west(0:,:)
      real(r8), intent(in) :: u_west(0:,:)
      real(r8), intent(in) :: v_west(0:,:)
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: u_east(0:,:)
      real(r8), intent(in) :: v_east(0:,:)
      real(r8), intent(inout) :: ad_u_east(0:,:)
      real(r8), intent(inout) :: ad_v_east(0:,:)
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: u_south(0:,:)
      real(r8), intent(in) :: v_south(0:,:)
      real(r8), intent(inout) :: ad_u_south(0:,:)
      real(r8), intent(inout) :: ad_v_south(0:,:)
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: u_north(0:,:)
      real(r8), intent(in) :: v_north(0:,:)
      real(r8), intent(inout) :: ad_u_north(0:,:)
      real(r8), intent(inout) :: ad_v_north(0:,:)
#   endif
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)

#  else

#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz_bry(LBij:UBij,N(ng),4)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hz_bry(LBij:UBij,N(ng),4)
#   ifdef WEST_M2OBC
      real(r8), intent(in) :: ubar_west(0:Jm(ng)+1)
      real(r8), intent(in) :: vbar_west(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_ubar_west(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_vbar_west(0:Jm(ng)+1)
#   endif
#   ifdef EAST_M2OBC
      real(r8), intent(in) :: ubar_east(0:Jm(ng)+1)
      real(r8), intent(in) :: vbar_east(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_ubar_east(0:Jm(ng)+1)
      real(r8), intent(inout) :: ad_vbar_east(0:Jm(ng)+1)
#   endif
#   ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ubar_south(0:Im(ng)+1)
      real(r8), intent(in) :: vbar_south(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_ubar_south(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_vbar_south(0:Im(ng)+1)
#   endif
#   ifdef NORTH_M2OBC
      real(r8), intent(in) :: ubar_north(0:Im(ng)+1)
      real(r8), intent(in) :: vbar_north(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_ubar_north(0:Im(ng)+1)
      real(r8), intent(inout) :: ad_vbar_north(0:Im(ng)+1)
#   endif
#   ifdef WEST_M3OBC
      real(r8), intent(in) :: u_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(in) :: v_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_west(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: u_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(in) :: v_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_east(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: u_south(0:Im(ng)+1,N(ng))
      real(r8), intent(in) :: v_south(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_south(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_south(0:Im(ng)+1,N(ng))
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: u_north(0:Im(ng)+1,N(ng))
      real(r8), intent(in) :: v_north(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_u_north(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: ad_v_north(0:Im(ng)+1,N(ng))
#   endif
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ib, it1, it2, j, k
      real(r8) :: adfac, fac, fac1, fac2
      real(r8) :: cff1, ad_cff1, ad_cff2

      real(r8), dimension(0:N(ng)) :: CF
      real(r8), dimension(0:N(ng)) :: DC

      real(r8), dimension(0:N(ng)) :: ad_CF
      real(r8), dimension(0:N(ng)) :: ad_DC

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint of adjust open boundary fields with 4DVAR increments.
!-----------------------------------------------------------------------
!
!  Set time records and interpolation factor, if any.
!
      IF (Nbrec(ng).eq.1) THEN
        it1=1
        it2=1
        fac1=1.0_r8
        fac2=0.0_r8
      ELSE
        it1=MAX(0,(iic(ng)-1)/nOBC(ng))+1
        it2=MIN(it1+1,Nbrec(ng))
        fac1=OBC_time(it2,ng)-(time(ng)+dt(ng))
        fac2=(time(ng)+dt(ng))-OBC_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF
!
!   Clear arrays and constants
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_CF=0.0_r8
      ad_DC=0.0_r8

!
!  2D U-momentum open boundaries: integrate 3D U-momentum at the open
!  boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
#  ifdef WEST_M2OBC
        IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=Jstr,Jend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz_bry(j,k,iwest)+                          &
     &                      Hz(i+1,j,k))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*u_west(j,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_ubar_west(j)=tl_cff2
!>
            ad_cff2=ad_ubar_west(j)
            ad_ubar_west(j)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*umask(i,j)
!>
            ad_cff2=ad_cff2*umask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_CF(0)=ad_CF(0)+ad_cff2*cff1
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*u_west(j,k)+                            &
!>   &                 DC(k)*tl_u_west(j,k)
!>
              ad_DC(k)=ad_DC(k)+u_west(j,k)*ad_CF(0)
              ad_u_west(j,k)=ad_u_west(j,k)+DC(k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz_bry(j,k,iwest)+                    &
!>   &                         tl_Hz(i+1,j,k))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz(i+1,j,k)=ad_Hz(i+1,j,k)+adfac
              ad_Hz_bry(j,k,iwest)=ad_Hz_bry(j,k,iwest)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif

#  ifdef EAST_M2OBC
        IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=Jstr,Jend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz(i-1,j,k)+                                &
     &                      Hz_bry(j,k,ieast))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*u_east(j,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_ubar_east(j)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_ubar_east(j)
            ad_ubar_east(j)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*umask(i,j)
!>
            ad_cff2=ad_cff2*umask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*u_east(j,k)+                            &
!>   &                 DC(k)*tl_u_east(j,k)
!>
              ad_DC(k)=ad_DC(k)+u_east(j,k)*ad_CF(0)
              ad_u_east(j,k)=ad_u_east(j,k)+DC(k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz(i-1,j,k)+                          &
!>   &                         tl_Hz_bry(j,k,ieast))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
              ad_Hz_bry(j,k,ieast)=ad_Hz_bry(j,k,ieast)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif

#  ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=Istr,Iend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz_bry(i-1,k,isouth)+                       &
     &                      Hz_bry(i  ,k,isouth))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*u_south(i,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_ubar_south(i)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_ubar_south(i)
            ad_ubar_south(i)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*umask(i,j)
!>
            ad_cff2=ad_cff2*umask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*u_south(i,k)+                           &
!>   &                 DC(k)*tl_u_south(i,k)
!>
              ad_DC(k)=ad_DC(k)+u_south(i,k)*ad_CF(0)
              ad_u_south(i,k)=ad_u_south(i,k)+DC(k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz_bry(i-1,k,isouth)+                 &
!>   &                         tl_Hz_bry(i  ,k,isouth))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz_bry(i-1,k,isouth)=ad_Hz_bry(i-1,k,isouth)+adfac
              ad_Hz_bry(i  ,k,isouth)=ad_Hz_bry(i  ,k,isouth)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif

#  ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=Istr,Iend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz_bry(i-1,k,inorth)+                       &
     &                      Hz_bry(i  ,k,inorth))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*u_north(i,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_ubar_north(i)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_ubar_north(i)
            ad_ubar_north(i)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*umask(i,j)
!>
            ad_cff2=ad_cff2*umask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*u_north(i,k)+                           &
!>   &                 DC(k)*tl_u_north(i,k)
!>
              ad_DC(k)=ad_DC(k)+u_north(i,k)*ad_CF(0)
              ad_u_north(i,k)=ad_u_north(i,k)+DC(k)*ad_CF(0)
              ad_CF(0)=0.0_r8
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz_bry(i-1,k,inorth)+                 &
!>   &                         tl_Hz_bry(i  ,k,inorth))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz_bry(i-1,k,inorth)=ad_Hz_bry(i-1,k,inorth)+adfac
              ad_Hz_bry(i  ,k,inorth)=ad_Hz_bry(i  ,k,inorth)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif
      END IF
!
!  3D V-momentum open boundaries: integrate 3D V-momentum at the open
!  boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
#  ifdef WEST_M2OBC
        IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrV,Jend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz_bry(j-1,k,iwest)+                        &
     &                      Hz_bry(j  ,k,iwest))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*v_west(j,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_vbar_west(j)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar_west(j)
            ad_vbar_west(j)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*vmask(i,j)
!>
            ad_cff2=ad_cff2*vmask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*v_west(j,k)+                            &
!>   &                 DC(k)*tl_v_west(j,k)
!>
              ad_v_west(j,k)=ad_v_west(j,k)+DC(k)*ad_CF(0)
              ad_DC(k)=ad_DC(k)+v_west(j,k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz_bry(j-1,k,iwest)+                  &
!>   &                         tl_Hz_bry(j  ,k,iwest))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz_bry(j-1,k,iwest)=ad_Hz_bry(j-1,k,iwest)+adfac
              ad_Hz_bry(j  ,k,iwest)=ad_Hz_bry(j  ,k,iwest)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
          END DO
        END IF
#  endif

#  ifdef EAST_M2OBC
        IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrV,Jend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz_bry(j-1,k,ieast)+                        &
     &                      Hz_bry(j  ,k,ieast))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*v_east(j,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_vbar_east(j)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar_east(j)
            ad_vbar_east(j)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*vmask(i,j)
!>
            ad_cff2=ad_cff2*vmask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*v_east(j,k)+                            &
!>   &                 DC(k)*tl_v_east(j,k)
!>
              ad_v_east(j,k)=ad_v_east(j,k)+DC(k)*ad_CF(0)
              ad_DC(k)=ad_DC(k)+v_east(j,k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
              ad_DC(0)=0.0_r8
!>            tl_DC(k)=0.5_r8*(tl_Hz_bry(j-1,k,ieast)+                  &
!>   &                         tl_Hz_bry(j  ,k,ieast))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz_bry(j-1,k,ieast)=ad_Hz_bry(j-1,k,ieast)+adfac
              ad_Hz_bry(j  ,k,ieast)=ad_Hz_bry(j  ,k,ieast)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif

#  ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=Istr,Iend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz_bry(i,k,isouth)+                         &
     &                      Hz(i+1,j,k))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*v_south(i,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_vbar_south(i)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar_south(i)
            ad_vbar_south(i)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*vmask(i,j)
!>
            ad_cff2=ad_cff2*vmask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*v_south(i,k)+                           &
!>   &                 DC(k)*tl_v_south(i,k)
!>
              ad_v_south(i,k)=ad_v_south(i,k)+DC(k)*ad_CF(0)
              ad_DC(k)=ad_DC(k)+v_south(i,k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz_bry(i,k,isouth)+                   &
!>   &                         tl_Hz(i+1,j,k))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz_bry(i,k,isouth)=ad_Hz_bry(i,k,isouth)+adfac
              ad_Hz(i+1,j,k)=ad_Hz(i+1,j,k)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif

#  ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=Istr,Iend
            DC(0)=0.0_r8
            CF(0)=0.0_r8
            DO k=1,N(ng)
              DC(k)=0.5_r8*(Hz(i,j-1,k)+                                &
     &                      Hz_bry(i,k,inorth))
              DC(0)=DC(0)+DC(k)
              CF(0)=CF(0)+DC(k)*v_north(i,k)
            END DO
            cff1=1.0_r8/DC(0)
!>          tl_vbar_north(i)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar_north(i)
            ad_vbar_north(i)=0.0_r8
#   ifdef MASKING
!>          tl_cff2=tl_cff2*vmask(i,j)
!>
            ad_cff2=ad_cff2*vmask(i,j)
#   endif
!>          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
!>
            ad_CF(0)=ad_CF(0)+cff1*ad_cff2
            ad_cff1=ad_cff1+CF(0)*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(0)
!>
            ad_DC(0)=ad_DC(0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
            DO k=N(ng),1,-1
!>            tl_CF(0)=tl_CF(0)+                                        &
!>   &                 tl_DC(k)*v_north(i,k)+                           &
!>   &                 DC(k)*tl_v_north(i,k)
!>
              ad_v_north(i,k)=ad_v_north(i,k)+DC(k)*ad_CF(0)
              ad_DC(k)=ad_DC(k)+v_north(i,k)*ad_CF(0)
!>            tl_DC(0)=tl_DC(0)+tl_DC(k)
!>
              ad_DC(k)=ad_DC(k)+ad_DC(0)
!>            tl_DC(k)=0.5_r8*(tl_Hz(i,j-1,k)+                          &
!>   &                         tl_Hz_bry(i,k,inorth))
!>
              adfac=0.5_r8*ad_DC(k)
              ad_Hz_bry(i,k,inorth)=ad_Hz_bry(i,k,inorth)+adfac
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
              ad_DC(k)=0.0_r8
            END DO
!>          tl_DC(0)=0.0_r8
!>
            ad_DC(0)=0.0_r8
!>          tl_CF(0)=0.0_r8
!>
            ad_CF(0)=0.0_r8
          END DO
        END IF
#  endif
      END IF

      RETURN
      END SUBROUTINE ad_obc2d_adjust_tile
# endif
#endif
      END MODULE ad_obc_adjust_mod
