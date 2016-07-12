#include "cppdefs.h"
      MODULE rp_obc_adjust_mod

#ifdef ADJUST_BOUNDARY
!
!svn $Id: rp_obc_adjust.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates 4DVar tangent linear model open      !
!  boundary increments.  The increments can be constant (Nbrec=1)      !
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
!
      PRIVATE
      PUBLIC :: rp_obc_adjust
# ifdef SOLVE3D
      PUBLIC :: rp_obc2d_adjust
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_obc_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
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
      CALL wclock_on (ng, iRPM, 7)
# endif
      CALL rp_obc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 7)
# endif

      RETURN
      END SUBROUTINE rp_obc_adjust
!
!***********************************************************************
      SUBROUTINE rp_obc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
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
!  Local variable declarations.
!
      integer :: i, it1, it2, j
# ifdef SOLVE3D
      integer :: it, k
# endif
      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust tangent linear open boundary fields with 4DVar increments.
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

# ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
      IF (tl_LBC(iwest,isUvel,ng)%acquire.and.                          &
     &      Lobc(iwest,isUvel,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            BOUNDARY(ng)%tl_u_west(j,k)=BOUNDARY(ng)%tl_u_west(j,k)+    &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_u_obc(j,k,      &
     &                                    iwest,it1,Linp)+              &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_u_obc(j,k,      &
     &                                    iwest,it2,Linp)
          END DO
        END DO
      END IF

      IF (tl_LBC(ieast,isUvel,ng)%acquire.and.                          &
     &      Lobc(ieast,isUvel,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            BOUNDARY(ng)%tl_u_east(j,k)=BOUNDARY(ng)%tl_u_east(j,k)+    &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_u_obc(j,k,      &
     &                                    ieast,it1,Linp)+              &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_u_obc(j,k,      &
     &                                    ieast,it2,Linp)
          END DO
        END DO
      END IF

      IF (tl_LBC(isouth,isUvel,ng)%acquire.and.                         &
     &      Lobc(isouth,isUvel,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BOUNDARY(ng)%tl_u_south(i,k)=BOUNDARY(ng)%tl_u_south(i,k)+  &
     &                                   fac1*                          &
     &                                   BOUNDARY(ng)%tl_u_obc(i,k,     &
     &                                     isouth,it1,Linp)+            &
     &                                   fac2*                          &
     &                                   BOUNDARY(ng)%tl_u_obc(i,k,     &
     &                                     isouth,it2,Linp)
          END DO
        END DO
      END IF

      IF (tl_LBC(inorth,isUvel,ng)%acquire.and.                         &
     &      Lobc(inorth,isUvel,ng).and.                                 &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BOUNDARY(ng)%tl_u_north(i,k)=BOUNDARY(ng)%tl_u_north(i,k)+  &
     &                                   fac1*                          &
     &                                   BOUNDARY(ng)%tl_u_obc(i,k,     &
     &                                     inorth,it1,Linp)+            &
     &                                   fac2*                          &
     &                                   BOUNDARY(ng)%tl_u_obc(i,k,     &
     &                                     inorth,it2,Linp)
          END DO
        END DO
      END IF
!
!  3D V-momentum open boundaries.
!
      IF (tl_LBC(iwest,isVvel,ng)%acquire.and.                          &
     &      Lobc(iwest,isVvel,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            BOUNDARY(ng)%tl_v_west(j,k)=BOUNDARY(ng)%tl_v_west(j,k)+    &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_v_obc(j,k,      &
     &                                    iwest,it1,Linp)+              &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_v_obc(j,k,      &
     &                                    iwest,it2,Linp)
          END DO
        END DO
      END IF

      IF (tl_LBC(ieast,isVvel,ng)%acquire.and.                          &
     &      Lobc(ieast,isVvel,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            BOUNDARY(ng)%tl_v_east(j,k)=BOUNDARY(ng)%tl_v_east(j,k)+    &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_v_obc(j,k,      &
     &                                    ieast,it1,Linp)+              &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_v_obc(j,k,      &
     &                                    ieast,it2,Linp)
          END DO
        END DO
      END IF

      IF (tl_LBC(isouth,isVvel,ng)%acquire.and.                         &
     &      Lobc(isouth,isVvel,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            BOUNDARY(ng)%tl_v_south(i,k)=BOUNDARY(ng)%tl_v_south(i,k)+  &
     &                                   fac1*                          &
     &                                   BOUNDARY(ng)%tl_v_obc(i,k,     &
     &                                     isouth,it1,Linp)+            &
     &                                   fac2*                          &
     &                                   BOUNDARY(ng)%tl_v_obc(i,k,     &
     &                                     isouth,it2,Linp)
          END DO
        END DO
      END IF

      IF (tl_LBC(inorth,isVvel,ng)%acquire.and.                         &
     &      Lobc(inorth,isVvel,ng).and.                                 &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            BOUNDARY(ng)%tl_v_north(i,k)=BOUNDARY(ng)%tl_v_north(i,k)+  &
     &                                   fac1*                          &
     &                                   BOUNDARY(ng)%tl_v_obc(i,k,     &
     &                                     inorth,it1,Linp)+            &
     &                                   fac2*                          &
     &                                   BOUNDARY(ng)%tl_v_obc(i,k,     &
     &                                     inorth,it2,Linp)
          END DO
        END DO
      END IF
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (tl_LBC(iwest,isTvar(it),ng)%acquire.and.                    &
     &        Lobc(iwest,isTvar(it),ng).and.                            &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              BOUNDARY(ng)%tl_t_west(j,k,it)=BOUNDARY(ng)%tl_t_west(j,  &
     &                                         k,it)+                   &
     &                                       fac1*                      &
     &                                       BOUNDARY(ng)%tl_t_obc(j,   &
     &                                         k,iwest,it1,Linp,it)+    &
     &                                       fac2*                      &
     &                                       BOUNDARY(ng)%tl_t_obc(j,   &
     &                                         k,iwest,it2,Linp,it)
            END DO
          END DO
        END IF

        IF (tl_LBC(ieast,isTvar(it),ng)%acquire.and.                    &
     &        Lobc(ieast,isTvar(it),ng).and.                            &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              BOUNDARY(ng)%tl_t_east(j,k,it)=BOUNDARY(ng)%tl_t_east(j,  &
     &                                         k,it)+                   &
     &                                       fac1*                      &
     &                                       BOUNDARY(ng)%tl_t_obc(j,   &
     &                                         k,ieast,it1,Linp,it)+    &
     &                                       fac2*                      &
     &                                       BOUNDARY(ng)%tl_t_obc(j,   &
     &                                         k,ieast,it2,Linp,it)
            END DO
          END DO
        END IF

        IF (tl_LBC(isouth,isTvar(it),ng)%acquire.and.                   &
     &        Lobc(isouth,isTvar(it),ng).and.                           &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              BOUNDARY(ng)%tl_t_south(i,k,it)=BOUNDARY(ng)%tl_t_south(i,&
     &                                          k,it)+                  &
     &                                        fac1*                     &
     &                                        BOUNDARY(ng)%tl_t_obc(i,  &
     &                                          k,isouth,it1,Linp,it)+  &
     &                                        fac2*                     &
     &                                        BOUNDARY(ng)%tl_t_obc(i,  &
     &                                          k,isouth,it2,Linp,it)
            END DO
          END DO
        END IF

        IF (tl_LBC(inorth,isTvar(it),ng)%acquire.and.                   &
     &        Lobc(inorth,isTvar(it),ng).and.                           &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              BOUNDARY(ng)%tl_t_north(i,k,it)=BOUNDARY(ng)%tl_t_north(i,&
     &                                          k,it)+                  &
     &                                        fac1*                     &
     &                                        BOUNDARY(ng)%tl_t_obc(i,  &
     &                                          k,inorth,it1,Linp,it)+  &
     &                                        fac2*                     &
     &                                        BOUNDARY(ng)%tl_t_obc(i,  &
     &                                          k,inorth,it2,Linp,it)
            END DO
          END DO
        END IF
      END DO
# endif
!
!  Free-surface open boundaries.
!
      IF (tl_LBC(iwest,isFsur,ng)%acquire.and.                          &
     &      Lobc(iwest,isFsur,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%tl_zeta_west(j)=BOUNDARY(ng)%tl_zeta_west(j)+    &
     &                                 fac1*                            &
     &                                 BOUNDARY(ng)%tl_zeta_obc(j,      &
     &                                   iwest,it1,Linp)+               &
     &                                 fac2*                            &
     &                                 BOUNDARY(ng)%tl_zeta_obc(j,      &
     &                                   iwest,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(ieast,isFsur,ng)%acquire.and.                          &
     &      Lobc(ieast,isFsur,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%tl_zeta_east(j)=BOUNDARY(ng)%tl_zeta_east(j)+    &
     &                                 fac1*                            &
     &                                 BOUNDARY(ng)%tl_zeta_obc(j,      &
     &                                   ieast,it1,Linp)+               &
     &                                 fac2*                            &
     &                                 BOUNDARY(ng)%tl_zeta_obc(j,      &
     &                                   ieast,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(isouth,isFsur,ng)%acquire.and.                         &
     &      Lobc(isouth,isFsur,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%tl_zeta_south(i)=BOUNDARY(ng)%tl_zeta_south(i)+  &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_zeta_obc(i,     &
     &                                    isouth,it1,Linp)+             &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_zeta_obc(i,     &
     &                                    isouth,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(inorth,isFsur,ng)%acquire.and.                         &
     &      Lobc(inorth,isFsur,ng).and.                                 &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%tl_zeta_north(i)=BOUNDARY(ng)%tl_zeta_north(i)+  &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_zeta_obc(i,     &
     &                                    inorth,it1,Linp)+             &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_zeta_obc(i,     &
     &                                    inorth,it2,Linp)
        END DO
      END IF

# ifndef SOLVE3D
!
!  2D U-momentum open boundaries.
!
      IF (tl_LBC(iwest,isUbar,ng)%acquire.and.                          &
     &      Lobc(iwest,isUbar,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%tl_ubar_west(j)=BOUNDARY(ng)%tl_ubar_west(j)+    &
     &                                 fac1*                            &
     &                                 BOUNDARY(ng)%tl_ubar_obc(j,      &
     &                                   iwest,it1,Linp)+               &
     &                                 fac2*                            &
     &                                 BOUNDARY(ng)%tl_ubar_obc(j,      &
     &                                   iwest,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(ieast,isUbar,ng)%acquire.and.                          &
     &      Lobc(ieast,isUbar,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%tl_ubar_east(j)=BOUNDARY(ng)%tl_ubar_east(j)+    &
     &                                 fac1*                            &
     &                                 BOUNDARY(ng)%tl_ubar_obc(j,      &
     &                                   ieast,it1,Linp)+               &
     &                                 fac2*                            &
     &                                 BOUNDARY(ng)%tl_ubar_obc(j,      &
     &                                   ieast,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(isouth,isUbar,ng)%acquire.and.                         &
     &      Lobc(isouth,isUbar,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrU,Iend
          BOUNDARY(ng)%tl_ubar_south(i)=BOUNDARY(ng)%tl_ubar_south(i)+  &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_ubar_obc(i,     &
     &                                    isouth,it1,Linp)+             &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_ubar_obc(i,     &
     &                                    isouth,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(inorth,isUbar,ng)%acquire.and.                         &
     &      Lobc(inorth,isUbar,ng).and.                                 &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrU,Iend
          BOUNDARY(ng)%tl_ubar_north(i)=BOUNDARY(ng)%tl_ubar_north(i)+  &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_ubar_obc(i,     &
     &                                    inorth,it1,Linp)+             &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_ubar_obc(i,     &
     &                                    inorth,it2,Linp)
        END DO
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (tl_LBC(iwest,isVbar,ng)%acquire.and.                          &
     &      Lobc(iwest,isVbar,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrV,Jend
          BOUNDARY(ng)%tl_vbar_west(j)=BOUNDARY(ng)%tl_vbar_west(j)+    &
     &                                 fac1*                            &
     &                                 BOUNDARY(ng)%tl_vbar_obc(j,      &
     &                                   iwest,it1,Linp)+               &
     &                                 fac2*                            &
     &                                 BOUNDARY(ng)%tl_vbar_obc(j,      &
     &                                   iwest,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(ieast,isVbar,ng)%acquire.and.                          &
     &      Lobc(ieast,isVbar,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrV,Jend
          BOUNDARY(ng)%tl_vbar_east(j)=BOUNDARY(ng)%tl_vbar_east(j)+    &
     &                                 fac1*                            &
     &                                 BOUNDARY(ng)%tl_vbar_obc(j,      &
     &                                   ieast,it1,Linp)+               &
     &                                 fac2*                            &
     &                                 BOUNDARY(ng)%tl_vbar_obc(j,      &
     &                                   ieast,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(isouth,isVbar,ng)%acquire.and.                         &
     &      Lobc(isouth,isVbar,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%tl_vbar_south(i)=BOUNDARY(ng)%tl_vbar_south(i)+  &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_vbar_obc(i,     &
     &                                    isouth,it1,Linp)+             &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_vbar_obc(i,     &
     &                                    isouth,it2,Linp)
        END DO
      END IF

      IF (tl_LBC(inorth,isVbar,ng)%acquire.and.                         &
     &      Lobc(inorth,isVbar,ng).and.                                 &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%tl_vbar_north(i)=BOUNDARY(ng)%tl_vbar_north(i)+  &
     &                                  fac1*                           &
     &                                  BOUNDARY(ng)%tl_vbar_obc(i,     &
     &                                    inorth,it1,Linp)+             &
     &                                  fac2*                           &
     &                                  BOUNDARY(ng)%tl_vbar_obc(i,     &
     &                                    inorth,it2,Linp)
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE rp_obc_adjust_tile

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE rp_obc2d_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
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
      CALL wclock_on (ng, iRPM, 7)
#  endif
      CALL rp_obc2d_adjust_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Linp,                                  &
#  ifdef MASKING
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % vmask,                      &
#  endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % Hz_bry,                     &
     &                           GRID(ng) % tl_Hz,                      &
     &                           GRID(ng) % tl_Hz_bry)
#  ifdef PROFILE
      CALL wclock_off (ng, iRPM, 7)
#  endif

      RETURN
      END SUBROUTINE rp_obc2d_adjust
!
!***********************************************************************
      SUBROUTINE rp_obc2d_adjust_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBij, UBij,  &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 Linp,                            &
#  ifdef MASKING
     &                                 umask, vmask,                    &
#  endif
     &                                 Hz, Hz_bry,                      &
     &                                 tl_Hz, tl_Hz_bry)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
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
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz_bry(LBij:,:,:)

#  else

#   ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz_bry(LBij:UBij,N(ng),4)
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz_bry(LBij:UBij,N(ng),4)
#  endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j, k

      real(r8) :: fac, fac1, fac2
      real(r8) :: cff1, cff2, cff3, tl_cff1, tl_cff2

      real(r8), dimension(0:N(ng)) :: CF
      real(r8), dimension(0:N(ng)) :: DC

      real(r8), dimension(0:N(ng)) :: tl_CF
      real(r8), dimension(0:N(ng)) :: tl_DC

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust tangent linear open boundary fields with 4DVar increments.
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
!  2D U-momentum open boundaries: integrate 3D U-momentum at the open
!  boundaries.
!
      IF (tl_LBC(iwest,isUbar,ng)%acquire.and.                          &
     &      Lobc(iwest,isUbar,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        i=BOUNDS(ng)%edge(iwest,r2dvar)
        DO j=Jstr,Jend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_u_obc(j,k,iwest,it1,Linp)+        &
     &           fac2*BOUNDARY(ng)%tl_u_obc(j,k,iwest,it2,Linp)
            DC(k)=0.5_r8*(Hz_bry(j,k,iwest)+                            &
     &                    Hz(i+1,j,k))
            tl_DC(k)=0.5_r8*(tl_Hz_bry(j,k,iwest)+                      &
     &                       tl_Hz(i+1,j,k))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%u_west(j,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%u_west(j,k)+                 &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*umask(i,j)
#  endif
          BOUNDARY(ng)%tl_ubar_west(j)=BOUNDARY(ng)%tl_ubar_west(j)+    &
     &                                 tl_cff2
        END DO
      END IF

      IF (tl_LBC(ieast,isUbar,ng)%acquire.and.                          &
     &      Lobc(ieast,isUbar,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        i=BOUNDS(ng)%edge(ieast,r2dvar)
        DO j=Jstr,Jend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_u_obc(j,k,ieast,it1,Linp)+        &
     &           fac2*BOUNDARY(ng)%tl_u_obc(j,k,ieast,it2,Linp)
            DC(k)=0.5_r8*(Hz(i-1,j,k)+                                  &
     &                    Hz_bry(j,k,ieast))
            tl_DC(k)=0.5_r8*(tl_Hz(i-1,j,k)+                            &
     &                       tl_Hz_bry(j,k,ieast))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%u_east(j,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%u_east(j,k)+                 &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*umask(i,j)
#  endif
          BOUNDARY(ng)%tl_ubar_east(j)=BOUNDARY(ng)%tl_ubar_east(j)+    &
     &                                 tl_cff2
        END DO
      END IF

      IF (tl_LBC(isouth,isUbar,ng)%acquire.and.                         &
     &      Lobc(isouth,isUbar,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        j=BOUNDS(ng)%edge(isouth,r2dvar)
        DO i=Istr,Iend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_u_obc(i,k,isouth,it1,Linp)+       &
     &           fac2*BOUNDARY(ng)%tl_u_obc(i,k,isouth,it2,Linp)
            DC(k)=0.5_r8*(Hz_bry(i-1,k,isouth)+                         &
     &                    Hz_bry(i  ,k,isouth))
            tl_DC(k)=0.5_r8*(tl_Hz_bry(i-1,k,isouth)+                   &
     &                       tl_Hz_bry(i  ,k,isouth))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%u_south(i,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%u_south(i,k)+                &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*umask(i,j)
#  endif
          BOUNDARY(ng)%tl_ubar_south(i)=BOUNDARY(ng)%tl_ubar_south(i)+  &
     &                                  tl_cff2
        END DO
      END IF

      IF (tl_LBC(inorth,isUbar,ng)%acquire.and.                         &
     &      Lobc(inorth,isUbar,ng).and.                                 &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        j=BOUNDS(ng)%edge(inorth,r2dvar)
        DO i=Istr,Iend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_u_obc(i,k,inorth,it1,Linp)+       &
     &           fac2*BOUNDARY(ng)%tl_u_obc(i,k,inorth,it2,Linp)
            DC(k)=0.5_r8*(Hz_bry(i-1,k,inorth)+                         &
     &                    Hz_bry(i  ,k,inorth))
            tl_DC(k)=0.5_r8*(tl_Hz_bry(i-1,k,inorth)+                   &
     &                       tl_Hz_bry(i  ,k,inorth))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%u_north(i,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%u_north(i,k)+                &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*umask(i,j)
#  endif
          BOUNDARY(ng)%tl_ubar_north(i)=BOUNDARY(ng)%tl_ubar_north(i)+  &
     &                                  tl_cff2
        END DO
      END IF
!
!  2D V-momentum open boundaries: integrate 3D V-momentum at the open
!  boundaries.
!
      IF (tl_LBC(iwest,isVbar,ng)%acquire.and.                          &
     &      Lobc(iwest,isVbar,ng).and.                                  &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        i=BOUNDS(ng)%edge(iwest,r2dvar)
        DO j=JstrV,Jend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_v_obc(j,k,iwest,it1,Linp)+        &
     &           fac2*BOUNDARY(ng)%tl_v_obc(j,k,iwest,it2,Linp)
            DC(k)=0.5_r8*(Hz_bry(j-1,k,iwest)+                          &
     &                    Hz_bry(j  ,k,iwest))
            tl_DC(k)=0.5_r8*(tl_Hz_bry(j-1,k,iwest)+                    &
     &                       tl_Hz_bry(j  ,k,iwest))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%v_west(j,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%v_west(j,k)+                 &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*vmask(i,j)
#  endif
          BOUNDARY(ng)%tl_vbar_west(j)=BOUNDARY(ng)%tl_vbar_west(j)+    &
     &                                 tl_cff2
        END DO
      END IF

      IF (tl_LBC(ieast,isVbar,ng)%acquire.and.                          &
     &      Lobc(ieast,isVbar,ng).and.                                  &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        i=BOUNDS(ng)%edge(ieast,r2dvar)
        DO j=JstrV,Jend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_v_obc(j,k,ieast,it1,Linp)+        &
     &           fac2*BOUNDARY(ng)%tl_v_obc(j,k,ieast,it2,Linp)
            DC(k)=0.5_r8*(Hz_bry(j-1,k,ieast)+                          &
     &                    Hz_bry(j  ,k,ieast))
            tl_DC(k)=0.5_r8*(tl_Hz_bry(j-1,k,ieast)+                    &
     &                       tl_Hz_bry(j  ,k,ieast))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%v_east(j,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%v_east(j,k)+                 &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*vmask(i,j)
#  endif
          BOUNDARY(ng)%tl_vbar_east(j)=BOUNDARY(ng)%tl_vbar_east(j)+    &
     &                                 tl_cff2
        END DO
      END IF

      IF (tl_LBC(isouth,isVbar,ng)%acquire.and.                         &
     &      Lobc(isouth,isVbar,ng).and.                                 &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        j=BOUNDS(ng)%edge(isouth,r2dvar)
        DO i=Istr,Iend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_v_obc(i,k,isouth,it1,Linp)+       &
     &           fac2*BOUNDARY(ng)%tl_v_obc(i,k,isouth,it2,Linp)
            DC(k)=0.5_r8*(Hz_bry(i,k,isouth)+                           &
     &                    Hz(i+1,j,k))
            tl_DC(k)=0.5_r8*(tl_Hz_bry(i,k,isouth)+                     &
     &                       tl_Hz(i+1,j,k))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%v_south(i,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%v_south(i,k)+                &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*vmask(i,j)
#  endif
          BOUNDARY(ng)%tl_vbar_south(i)=BOUNDARY(ng)%tl_vbar_south(i)+  &
     &                                  tl_cff2
        END DO
      END IF

      IF (tl_LBC(inorth,isVbar,ng)%acquire.and.                         &
     &      Lobc(inorth,isVbar,ng).and.                                 &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        j=BOUNDS(ng)%edge(inorth,r2dvar)
        DO i=Istr,Iend
          DC(0)=0.0_r8
          tl_DC(0)=0.0_r8
          CF(0)=0.0_r8
          tl_CF(0)=0.0_r8
          DO k=1,N(ng)
            cff3=fac1*BOUNDARY(ng)%tl_v_obc(i,k,inorth,it1,Linp)+       &
     &           fac2*BOUNDARY(ng)%tl_v_obc(i,k,inorth,it2,Linp)
            DC(k)=0.5_r8*(Hz(i,j-1,k)+                                  &
     &                    Hz_bry(i,k,inorth))
            tl_DC(k)=0.5_r8*(tl_Hz(i,j-1,k)+                            &
     &                       tl_Hz_bry(i,k,inorth))
            DC(0)=DC(0)+DC(k)
            tl_DC(0)=tl_DC(0)+tl_DC(k)
            CF(0)=CF(0)+DC(k)*BOUNDARY(ng)%v_north(i,k)
            tl_CF(0)=tl_CF(0)+                                          &
     &               tl_DC(k)*BOUNDARY(ng)%v_north(i,k)+                &
     &               DC(k)*cff3
          END DO
          cff1=1.0_r8/DC(0)
          tl_cff1=-cff1*cff1*tl_DC(0)
          tl_cff2=tl_CF(0)*cff1+CF(0)*tl_cff1
#  ifdef MASKING
          tl_cff2=tl_cff2*vmask(i,j)
#  endif
          BOUNDARY(ng)%tl_vbar_north(i)=BOUNDARY(ng)%tl_vbar_north(i)+  &
     &                                  tl_cff2
        END DO
      END IF

      RETURN
      END SUBROUTINE rp_obc2d_adjust_tile
# endif
#endif
      END MODULE rp_obc_adjust_mod
