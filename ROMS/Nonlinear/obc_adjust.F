#include "cppdefs.h"
      MODULE obc_adjust_mod

#ifdef ADJUST_BOUNDARY
!
!svn $Id: obc_adjust.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates 4DVar open boundary increments       !
!  and adds them to  nonlinear model  open boundary arrays.  The       !
!  increments can be constant  (Nbrec=1)  or  time  interpolated       !
!  between snapshots (Nbrec>1).                                        !
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
      PUBLIC :: obc_adjust
      PUBLIC :: load_obc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE obc_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
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
      CALL wclock_on (ng, iNLM, 7)
# endif
      CALL obc_adjust_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Linp,                                       &
# ifdef WEST_FSOBC
     &                      BOUNDARY(ng) % zeta_west,                   &
# endif
# ifdef EAST_FSOBC
     &                      BOUNDARY(ng) % zeta_east,                   &
# endif
# ifdef SOUTH_FSOBC
     &                      BOUNDARY(ng) % zeta_south,                  &
# endif
# ifdef NORTH_FSOBC
     &                      BOUNDARY(ng) % zeta_north,                  &
# endif
# ifdef WEST_M2OBC
     &                      BOUNDARY(ng) % ubar_west,                   &
     &                      BOUNDARY(ng) % vbar_west,                   &
# endif
# ifdef EAST_M2OBC
     &                      BOUNDARY(ng) % ubar_east,                   &
     &                      BOUNDARY(ng) % vbar_east,                   &
# endif
# ifdef SOUTH_M2OBC
     &                      BOUNDARY(ng) % ubar_south,                  &
     &                      BOUNDARY(ng) % vbar_south,                  &
# endif
# ifdef NORTH_M2OBC
     &                      BOUNDARY(ng) % ubar_north,                  &
     &                      BOUNDARY(ng) % vbar_north,                  &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                      BOUNDARY(ng) % u_west,                      &
     &                      BOUNDARY(ng) % v_west,                      &
#  endif
#  ifdef EAST_M3OBC
     &                      BOUNDARY(ng) % u_east,                      &
     &                      BOUNDARY(ng) % v_east,                      &
#  endif
#  ifdef SOUTH_M3OBC
     &                      BOUNDARY(ng) % u_south,                     &
     &                      BOUNDARY(ng) % v_south,                     &
#  endif
#  ifdef NORTH_M3OBC
     &                      BOUNDARY(ng) % u_north,                     &
     &                      BOUNDARY(ng) % v_north,                     &
#  endif
#  ifdef WEST_TOBC
     &                      BOUNDARY(ng) % t_west,                      &
#  endif
#  ifdef EAST_TOBC
     &                      BOUNDARY(ng) % t_east,                      &
#  endif
#  ifdef SOUTH_TOBC
     &                      BOUNDARY(ng) % t_south,                     &
#  endif
#  ifdef NORTH_TOBC
     &                      BOUNDARY(ng) % t_north,                     &
#  endif
# endif
# ifdef SOLVE3D
     &                      BOUNDARY(ng) % tl_t_obc,                    &
     &                      BOUNDARY(ng) % tl_u_obc,                    &
     &                      BOUNDARY(ng) % tl_v_obc,                    &
# endif
     &                      BOUNDARY(ng) % tl_ubar_obc,                 &
     &                      BOUNDARY(ng) % tl_vbar_obc,                 &
     &                      BOUNDARY(ng) % tl_zeta_obc)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
# endif
      RETURN
      END SUBROUTINE obc_adjust
!
!***********************************************************************
      SUBROUTINE obc_adjust_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Linp,                                 &
# ifdef WEST_FSOBC
     &                            zeta_west,                            &
# endif
# ifdef EAST_FSOBC
     &                            zeta_east,                            &
# endif
# ifdef SOUTH_FSOBC
     &                            zeta_south,                           &
# endif
# ifdef NORTH_FSOBC
     &                            zeta_north,                           &
# endif
# ifdef WEST_M2OBC
     &                            ubar_west, vbar_west,                 &
# endif
# ifdef EAST_M2OBC
     &                            ubar_east, vbar_east,                 &
# endif
# ifdef SOUTH_M2OBC
     &                            ubar_south, vbar_south,               &
# endif
# ifdef NORTH_M2OBC
     &                            ubar_north, vbar_north,               &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                            u_west, v_west,                       &
#  endif
#  ifdef EAST_M3OBC
     &                            u_east, v_east,                       &
#  endif
#  ifdef SOUTH_M3OBC
     &                            u_south, v_south,                     &
#  endif
#  ifdef NORTH_M3OBC
     &                            u_north, v_north,                     &
#  endif
#  ifdef WEST_TOBC
     &                            t_west,                               &
#  endif
#  ifdef EAST_TOBC
     &                            t_east,                               &
#  endif
#  ifdef SOUTH_TOBC
     &                            t_south,                              &
#  endif
#  ifdef NORTH_TOBC
     &                            t_north,                              &
#  endif
# endif
# ifdef SOLVE3D
     &                            tl_t_obc, tl_u_obc, tl_v_obc,         &
# endif
     &                            tl_ubar_obc, tl_vbar_obc,             &
     &                            tl_zeta_obc)
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
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: tl_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_zeta_obc(LBij:,:,:,:)
#  ifdef WEST_FSOBC
      real(r8), intent(inout) :: zeta_west(0:)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(inout) :: zeta_east(0:)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(inout) :: zeta_south(0:)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(inout) :: zeta_north(0:)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(inout) :: ubar_west(0:)
      real(r8), intent(inout) :: vbar_west(0:)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(inout) :: ubar_east(0:)
      real(r8), intent(inout) :: vbar_east(0:)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(inout) :: ubar_south(0:)
      real(r8), intent(inout) :: vbar_south(0:)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(inout) :: ubar_north(0:)
      real(r8), intent(inout) :: vbar_north(0:)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(inout) :: u_west(0:,:)
      real(r8), intent(inout) :: v_west(0:,:)
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(inout) :: u_east(0:,:)
      real(r8), intent(inout) :: v_east(0:,:)
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(inout) :: u_south(0:,:)
      real(r8), intent(inout) :: v_south(0:,:)
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(inout) :: u_north(0:,:)
      real(r8), intent(inout) :: v_north(0:,:)
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(inout) :: t_west(0:,:,:)
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(inout) :: t_east(0:,:,:)
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(inout) :: t_south(0:,:,:)
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(inout) :: t_north(0:,:,:)
#   endif
#  endif

# else

#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  ifdef WEST_FSOBC
      real(r8), intent(inout) :: zeta_west(0:Jm(ng)+1)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(inout) :: zeta_east(0:Jm(ng)+1)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(inout) :: zeta_south(0:Im(ng)+1)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(inout) :: zeta_north(0:Im(ng)+1)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(inout) :: ubar_west(0:Jm(ng)+1)
      real(r8), intent(inout) :: vbar_west(0:Jm(ng)+1)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(inout) :: ubar_east(0:Jm(ng)+1)
      real(r8), intent(inout) :: vbar_east(0:Jm(ng)+1)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(inout) :: ubar_south(0:Im(ng)+1)
      real(r8), intent(inout) :: vbar_south(0:Im(ng)+1)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(inout) :: ubar_north(0:Im(ng)+1)
      real(r8), intent(inout) :: vbar_north(0:Im(ng)+1)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(inout) :: u_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: v_west(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(inout) :: u_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(inout) :: v_east(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(inout) :: u_south(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: v_south(0:Im(ng)+1,N(ng))
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(inout) :: u_north(0:Im(ng)+1,N(ng))
      real(r8), intent(inout) :: v_north(0:Im(ng)+1,N(ng))
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(inout) :: t_west(0:Jm(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(inout) :: t_east(0:Jm(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(inout) :: t_south(0:Im(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(inout) :: t_north(0:Im(ng)+1,N(ng),NT(ng))
#   endif
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j
# ifdef SOLVE3D
      integer :: k, it
# endif
      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust nonlinear open boundaries with 4DVar increments.
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
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
# ifdef WEST_FSOBC
        IF ((Lobc(iwest,isFsur,ng)).and.WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            zeta_west(j)=zeta_west(j)+                                  &
     &                   fac1*tl_zeta_obc(j,iwest,it1,Linp)+            &
     &                   fac2*tl_zeta_obc(j,iwest,it2,Linp)
          END DO
        END IF
# endif
# ifdef EAST_FSOBC
        IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            zeta_east(j)=zeta_east(j)+                                  &
     &                   fac1*tl_zeta_obc(j,ieast,it1,Linp)+            &
     &                   fac2*tl_zeta_obc(j,ieast,it2,Linp)
          END DO
        END IF
# endif
# ifdef SOUTH_FSOBC
        IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            zeta_south(i)=zeta_south(i)+                                &
     &                    fac1*tl_zeta_obc(i,isouth,it1,Linp)+          &
     &                    fac2*tl_zeta_obc(i,isouth,it2,Linp)
          END DO
        END IF
# endif
# ifdef NORTH_FSOBC
        IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            zeta_north(i)=zeta_north(i)+                                &
     &                    fac1*tl_zeta_obc(i,inorth,it1,Linp)+          &
     &                    fac2*tl_zeta_obc(i,inorth,it2,Linp)
          END DO
        END IF
# endif
      END IF
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
# ifdef WEST_M2OBC
        IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            ubar_west(j)=ubar_west(j)+                                  &
     &                   fac1*tl_ubar_obc(j,iwest,it1,Linp)+            &
     &                   fac2*tl_ubar_obc(j,iwest,it2,Linp)
          END DO
        END IF
# endif
# ifdef EAST_M2OBC
        IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            ubar_east(j)=ubar_east(j)+                                  &
     &                   fac1*tl_ubar_obc(j,ieast,it1,Linp)+            &
     &                   fac2*tl_ubar_obc(j,ieast,it2,Linp)
          END DO
        END IF
# endif
# ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
          DO i=IstrU,Iend
            ubar_south(i)=ubar_south(i)+                                &
     &                    fac1*tl_ubar_obc(i,isouth,it1,Linp)+          &
     &                    fac2*tl_ubar_obc(i,isouth,it2,Linp)
          END DO
        END IF
# endif
# ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
          DO i=IstrU,Iend
            ubar_north(i)=ubar_north(i)+                                &
     &                    fac1*tl_ubar_obc(i,inorth,it1,Linp)+          &
     &                    fac2*tl_ubar_obc(i,inorth,it2,Linp)
          END DO
        END IF
# endif
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
# ifdef WEST_M2OBC
        IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
          DO j=JstrV,Jend
            vbar_west(j)=vbar_west(j)+                                  &
     &                   fac1*tl_vbar_obc(j,iwest,it1,Linp)+            &
     &                   fac2*tl_vbar_obc(j,iwest,it2,Linp)
          END DO
        END IF
# endif
# ifdef EAST_M2OBC
        IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
          DO j=JstrV,Jend
            vbar_east(j)=vbar_east(j)+                                  &
     &                   fac1*tl_vbar_obc(j,ieast,it1,Linp)+            &
     &                   fac2*tl_vbar_obc(j,ieast,it2,Linp)
          END DO
        END IF
# endif
# ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            vbar_south(i)=vbar_south(i)+                                &
     &                    fac1*tl_vbar_obc(i,isouth,it1,Linp)+          &
     &                    fac2*tl_vbar_obc(i,isouth,it2,Linp)
          END DO
        END IF
# endif
# ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            vbar_north(i)=vbar_north(i)+                                &
     &                    fac1*tl_vbar_obc(i,inorth,it1,Linp)+          &
     &                    fac2*tl_vbar_obc(i,inorth,it2,Linp)
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
          DO k=1,N(ng)
            DO j=Jstr,Jend
              u_west(j,k)=u_west(j,k)+                                  &
     &                    fac1*tl_u_obc(j,k,iwest,it1,Linp)+            &
     &                    fac2*tl_u_obc(j,k,iwest,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              u_east(j,k)=u_east(j,k)+                                  &
     &                    fac1*tl_u_obc(j,k,ieast,it1,Linp)+            &
     &                    fac2*tl_u_obc(j,k,ieast,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              u_south(i,k)=u_south(i,k)+                                &
     &                     fac1*tl_u_obc(i,k,isouth,it1,Linp)+          &
     &                     fac2*tl_u_obc(i,k,isouth,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              u_north(i,k)=u_north(i,k)+                                &
     &                     fac1*tl_u_obc(i,k,inorth,it1,Linp)+          &
     &                     fac2*tl_u_obc(i,k,inorth,it2,Linp)
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
          DO k=1,N(ng)
            DO j=JstrV,Jend
              v_west(j,k)=v_west(j,k)+                                  &
     &                    fac1*tl_v_obc(j,k,iwest,it1,Linp)+            &
     &                    fac2*tl_v_obc(j,k,iwest,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              v_east(j,k)=v_east(j,k)+                                  &
     &                    fac1*tl_v_obc(j,k,ieast,it1,Linp)+            &
     &                    fac2*tl_v_obc(j,k,ieast,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              v_south(i,k)=v_south(i,k)+                                &
     &                     fac1*tl_v_obc(i,k,isouth,it1,Linp)+          &
     &                     fac2*tl_v_obc(i,k,isouth,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              v_north(i,k)=v_north(i,k)+                                &
     &                     fac1*tl_v_obc(i,k,inorth,it1,Linp)+          &
     &                     fac2*tl_v_obc(i,k,inorth,it2,Linp)
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
            DO k=1,N(ng)
              DO j=Jstr,Jend
                t_west(j,k,it)=t_west(j,k,it)+                          &
     &                         fac1*tl_t_obc(j,k,iwest,it1,Linp,it)+    &
     &                         fac2*tl_t_obc(j,k,iwest,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef EAST_TOBC
          IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                t_east(j,k,it)=t_east(j,k,it)+                          &
     &                         fac1*tl_t_obc(j,k,ieast,it1,Linp,it)+    &
     &                         fac2*tl_t_obc(j,k,ieast,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef SOUTH_TOBC
          IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                t_south(i,k,it)=t_south(i,k,it)+                        &
     &                          fac1*tl_t_obc(i,k,isouth,it1,Linp,it)+  &
     &                          fac2*tl_t_obc(i,k,isouth,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef NORTH_TOBC
          IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                t_north(i,k,it)=t_north(i,k,it)+                        &
     &                          fac1*tl_t_obc(i,k,inorth,it1,Linp,it)+  &
     &                          fac2*tl_t_obc(i,k,inorth,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE obc_adjust_tile

      SUBROUTINE load_obc (ng, tile, Lout)
!
!=======================================================================
!                                                                      !
!  This routine loads open boundaries into nonlinear storage arrays.   !
!  In 4DVAR open boundary adjustment, the boundary values are stored   !
!  in arrays with extra dimensions to aid minimization at other times  !
!  in addition to initialization time.                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Lout      Time index to process in storage arrays.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 8)
# endif
      CALL load_obc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lout,                                         &
# ifdef WEST_FSOBC
     &                    BOUNDARY(ng) % zeta_west,                     &
# endif
# ifdef EAST_FSOBC
     &                    BOUNDARY(ng) % zeta_east,                     &
# endif
# ifdef SOUTH_FSOBC
     &                    BOUNDARY(ng) % zeta_south,                    &
# endif
# ifdef NORTH_FSOBC
     &                    BOUNDARY(ng) % zeta_north,                    &
# endif
# ifdef WEST_M2OBC
     &                    BOUNDARY(ng) % ubar_west,                     &
     &                    BOUNDARY(ng) % vbar_west,                     &
# endif
# ifdef EAST_M2OBC
     &                    BOUNDARY(ng) % ubar_east,                     &
     &                    BOUNDARY(ng) % vbar_east,                     &
# endif
# ifdef SOUTH_M2OBC
     &                    BOUNDARY(ng) % ubar_south,                    &
     &                    BOUNDARY(ng) % vbar_south,                    &
# endif
# ifdef NORTH_M2OBC
     &                    BOUNDARY(ng) % ubar_north,                    &
     &                    BOUNDARY(ng) % vbar_north,                    &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                    BOUNDARY(ng) % u_west,                        &
     &                    BOUNDARY(ng) % v_west,                        &
#  endif
#  ifdef EAST_M3OBC
     &                    BOUNDARY(ng) % u_east,                        &
     &                    BOUNDARY(ng) % v_east,                        &
#  endif
#  ifdef SOUTH_M3OBC
     &                    BOUNDARY(ng) % u_south,                       &
     &                    BOUNDARY(ng) % v_south,                       &
#  endif
#  ifdef NORTH_M3OBC
     &                    BOUNDARY(ng) % u_north,                       &
     &                    BOUNDARY(ng) % v_north,                       &
#  endif
#  ifdef WEST_TOBC
     &                    BOUNDARY(ng) % t_west,                        &
#  endif
#  ifdef EAST_TOBC
     &                    BOUNDARY(ng) % t_east,                        &
#  endif
#  ifdef SOUTH_TOBC
     &                    BOUNDARY(ng) % t_south,                       &
#  endif
#  ifdef NORTH_TOBC
     &                    BOUNDARY(ng) % t_north,                       &
#  endif
# endif
# ifdef SOLVE3D
     &                    BOUNDARY(ng) % t_obc,                         &
     &                    BOUNDARY(ng) % u_obc,                         &
     &                    BOUNDARY(ng) % v_obc,                         &
# endif
     &                    BOUNDARY(ng) % ubar_obc,                      &
     &                    BOUNDARY(ng) % vbar_obc,                      &
     &                    BOUNDARY(ng) % zeta_obc)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 8)
# endif

      RETURN
      END SUBROUTINE load_obc
!
!***********************************************************************
      SUBROUTINE load_obc_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lout,                                   &
# ifdef WEST_FSOBC
     &                          zeta_west,                              &
# endif
# ifdef EAST_FSOBC
     &                          zeta_east,                              &
# endif
# ifdef SOUTH_FSOBC
     &                          zeta_south,                             &
# endif
# ifdef NORTH_FSOBC
     &                          zeta_north,                             &
# endif
# ifdef WEST_M2OBC
     &                          ubar_west, vbar_west,                   &
# endif
# ifdef EAST_M2OBC
     &                          ubar_east, vbar_east,                   &
# endif
# ifdef SOUTH_M2OBC
     &                          ubar_south, vbar_south,                 &
# endif
# ifdef NORTH_M2OBC
     &                          ubar_north, vbar_north,                 &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                          u_west, v_west,                         &
#  endif
#  ifdef EAST_M3OBC
     &                          u_east, v_east,                         &
#  endif
#  ifdef SOUTH_M3OBC
     &                          u_south, v_south,                       &
#  endif
#  ifdef NORTH_M3OBC
     &                          u_north, v_north,                       &
#  endif
#  ifdef WEST_TOBC
     &                          t_west,                                 &
#  endif
#  ifdef EAST_TOBC
     &                          t_east,                                 &
#  endif
#  ifdef SOUTH_TOBC
     &                          t_south,                                &
#  endif
#  ifdef NORTH_TOBC
     &                          t_north,                                &
#  endif
# endif
# ifdef SOLVE3D
     &                          t_obc, u_obc, v_obc,                    &
# endif
     &                          ubar_obc, vbar_obc,                     &
     &                          zeta_obc)
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
      integer, intent(in) :: Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef WEST_FSOBC
      real(r8), intent(in) :: zeta_west(0:)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(in) :: zeta_east(0:)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(in) :: zeta_south(0:)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(in) :: zeta_north(0:)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(in) :: ubar_west(0:)
      real(r8), intent(in) :: vbar_west(0:)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(in) :: ubar_east(0:)
      real(r8), intent(in) :: vbar_east(0:)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ubar_south(0:)
      real(r8), intent(in) :: vbar_south(0:)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(in) :: ubar_north(0:)
      real(r8), intent(in) :: vbar_north(0:)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(in) :: u_west(0:,:)
      real(r8), intent(in) :: v_west(0:,:)
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: u_east(0:,:)
      real(r8), intent(in) :: v_east(0:,:)
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: u_south(0:,:)
      real(r8), intent(in) :: v_south(0:,:)
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: u_north(0:,:)
      real(r8), intent(in) :: v_north(0:,:)
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(in) :: t_west(0:,:,:)
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(in) :: t_east(0:,:,:)
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(in) :: t_south(0:,:,:)
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(in) :: t_north(0:,:,:)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: zeta_obc(LBij:,:,:,:)

# else

#  ifdef WEST_FSOBC
      real(r8), intent(in) :: zeta_west(0:Jm(ng)+1)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(in) :: zeta_east(0:Jm(ng)+1)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(in) :: zeta_south(0:Im(ng)+1)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(in) :: zeta_north(0:Im(ng)+1)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(in) :: ubar_west(0:Jm(ng)+1)
      real(r8), intent(in) :: vbar_west(0:Jm(ng)+1)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(in) :: ubar_east(0:Jm(ng)+1)
      real(r8), intent(in) :: vbar_east(0:Jm(ng)+1)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ubar_south(0:Im(ng)+1)
      real(r8), intent(in) :: vbar_south(0:Im(ng)+1)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(in) :: ubar_north(0:Im(ng)+1)
      real(r8), intent(in) :: vbar_north(0:Im(ng)+1)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(in) :: u_west(0:Jm(ng)+1,N(ng))
      real(r8), intent(in) :: v_west(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: u_east(0:Jm(ng)+1,N(ng))
      real(r8), intent(in) :: v_east(0:Jm(ng)+1,N(ng))
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: u_south(0:Im(ng)+1,N(ng))
      real(r8), intent(in) :: v_south(0:Im(ng)+1,N(ng))
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: u_north(0:Im(ng)+1,N(ng))
      real(r8), intent(in) :: v_north(0:Im(ng)+1,N(ng))
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(in) :: t_west(0:Jm(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(in) :: t_east(0:Jm(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(in) :: t_south(0:Im(ng)+1,N(ng),NT(ng))
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(in) :: t_north(0:Im(ng)+1,N(ng),NT(ng))
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef SOLVE3D
      integer :: it, k
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Load nonlinear open boundary fields into storage arrays.
!-----------------------------------------------------------------------
!
      LOAD : IF (MOD(iic(ng)-1,nOBC(ng)).eq.0) THEN
        OBCcount(ng)=OBCcount(ng)+1
!
!  Free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
# ifdef WEST_FSOBC
          IF ((Lobc(iwest,isFsur,ng)).and.WESTERN_EDGE) THEN
            DO j=Jstr,Jend
              zeta_obc(j,iwest,OBCcount(ng),Lout)=zeta_west(j)
            END DO
          END IF
# endif
# ifdef EAST_FSOBC
          IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
            DO j=Jstr,Jend
              zeta_obc(j,ieast,OBCcount(ng),Lout)=zeta_east(j)
            END DO
          END IF
# endif
# ifdef SOUTH_FSOBC
          IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
            DO i=Istr,Iend
              zeta_obc(i,isouth,OBCcount(ng),Lout)=zeta_south(i)
            END DO
          END IF
# endif
# ifdef NORTH_FSOBC
          IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
            DO i=Istr,Iend
              zeta_obc(i,inorth,OBCcount(ng),Lout)=zeta_north(i)
            END DO
          END IF
# endif
        END IF
!
!  2D U-momentum open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
# ifdef WEST_M2OBC
          IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
            DO j=Jstr,Jend
              ubar_obc(j,iwest,OBCcount(ng),Lout)=ubar_west(j)
            END DO
          END IF
# endif
# ifdef EAST_M2OBC
          IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
            DO j=Jstr,Jend
              ubar_obc(j,ieast,OBCcount(ng),Lout)=ubar_east(j)
            END DO
          END IF
# endif
# ifdef SOUTH_M2OBC
          IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
            DO i=IstrU,Iend
              ubar_obc(i,isouth,OBCcount(ng),Lout)=ubar_south(i)
            END DO
          END IF
# endif
# ifdef NORTH_M2OBC
          IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
            DO i=IstrU,Iend
              ubar_obc(i,inorth,OBCcount(ng),Lout)=ubar_north(i)
            END DO
          END IF
# endif
        END IF
!
!  2D V-momentum open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
# ifdef WEST_M2OBC
          IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
            DO j=JstrV,Jend
              vbar_obc(j,iwest,OBCcount(ng),Lout)=vbar_west(j)
            END DO
          END IF
# endif
# ifdef EAST_M2OBC
          IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
            DO j=JstrV,Jend
              vbar_obc(j,ieast,OBCcount(ng),Lout)=vbar_east(j)
            END DO
          END IF
# endif
# ifdef SOUTH_M2OBC
          IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
            DO i=Istr,Iend
              vbar_obc(i,isouth,OBCcount(ng),Lout)=vbar_south(i)
            END DO
          END IF
# endif
# ifdef NORTH_M2OBC
          IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
            DO i=Istr,Iend
              vbar_obc(i,inorth,OBCcount(ng),Lout)=vbar_north(i)
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
            DO k=1,N(ng)
              DO j=Jstr,Jend
                u_obc(j,k,iwest,OBCcount(ng),Lout)=u_west(j,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef EAST_M3OBC
          IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                u_obc(j,k,ieast,OBCcount(ng),Lout)=u_east(j,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef SOUTH_M3OBC
          IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                u_obc(i,k,isouth,OBCcount(ng),Lout)=u_south(i,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef NORTH_M3OBC
          IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                u_obc(i,k,inorth,OBCcount(ng),Lout)=u_north(i,k)
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
            DO k=1,N(ng)
              DO j=JstrV,Jend
                v_obc(j,k,iwest,OBCcount(ng),Lout)=v_west(j,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef EAST_M3OBC
          IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                v_obc(j,k,ieast,OBCcount(ng),Lout)=v_east(j,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef SOUTH_M3OBC
          IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                v_obc(i,k,isouth,OBCcount(ng),Lout)=v_south(i,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef NORTH_M3OBC
          IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                v_obc(i,k,inorth,OBCcount(ng),Lout)=v_north(i,k)
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
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  t_obc(j,k,iwest,OBCcount(ng),Lout,it)=t_west(j,k,it)
                END DO
              END DO
            END IF
#  endif
#  ifdef EAST_TOBC
            IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  t_obc(j,k,ieast,OBCcount(ng),Lout,it)=t_east(j,k,it)
                END DO
              END DO
            END IF
#  endif
#  ifdef SOUTH_TOBC
            IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  t_obc(i,k,isouth,OBCcount(ng),Lout,it)=t_south(i,k,it)
                END DO
              END DO
            END IF
#  endif
#  ifdef NORTH_TOBC
            IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  t_obc(i,k,inorth,OBCcount(ng),Lout,it)=t_north(i,k,it)
                END DO
              END DO
            END IF
#  endif
          END IF
        END DO
# endif
      END IF LOAD

      RETURN
      END SUBROUTINE load_obc_tile
#endif
      END MODULE obc_adjust_mod
