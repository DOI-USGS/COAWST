#include "cppdefs.h"
      MODULE obc_adjust_mod

#ifdef ADJUST_BOUNDARY
!
!svn $Id: obc_adjust.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!
      PRIVATE
      PUBLIC :: obc_adjust
      PUBLIC :: load_obc
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE obc_adjust (ng, tile, Linp)
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
      CALL wclock_on (ng, iNLM, 7)
# endif
      CALL obc_adjust_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Linp)
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
     &                            Linp)
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
      integer :: k, it
# endif
      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust nonlinear open boundaries with 4DVar increments.
!
!  HGA: need to set fac2=0, when the boundary condition in the NLM
!       are radiations until we know what to do.
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
      IF (LBC (iwest,isFsur,ng)%acquire.and.                            &
     &    Lobc(iwest,isFsur,ng).and.                                    &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%zeta_west(j)=BOUNDARY(ng)%zeta_west(j)+          &
     &                              fac1*                               &
     &                              BOUNDARY(ng)%tl_zeta_obc(j,         &
     &                                iwest,it1,Linp)+                  &
     &                              fac2*                               &
     &                              BOUNDARY(ng)%tl_zeta_obc(j,         &
     &                                iwest,it2,Linp)
        END DO
      END IF

      IF (LBC (ieast,isFsur,ng)%acquire.and.                            &
     &    Lobc(ieast,isFsur,ng).and.                                    &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%zeta_east(j)=BOUNDARY(ng)%zeta_east(j)+          &
     &                              fac1*                               &
     &                              BOUNDARY(ng)%tl_zeta_obc(j,         &
     &                                ieast,it1,Linp)+                  &
     &                              fac2*                               &
     &                              BOUNDARY(ng)%tl_zeta_obc(j,         &
     &                                ieast,it2,Linp)
        END DO
      END IF

      IF (LBC (isouth,isFsur,ng)%acquire.and.                           &
     &    Lobc(isouth,isFsur,ng).and.                                   &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%zeta_south(i)=BOUNDARY(ng)%zeta_south(i)+        &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_zeta_obc(i,        &
     &                                 isouth,it1,Linp)+                &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_zeta_obc(i,        &
     &                                 isouth,it2,Linp)
        END DO
      END IF

      IF (LBC (inorth,isFsur,ng)%acquire.and.                           &
     &    Lobc(inorth,isFsur,ng).and.                                   &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%zeta_north(i)=BOUNDARY(ng)%zeta_north(i)+        &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_zeta_obc(i,        &
     &                                 inorth,it1,Linp)+                &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_zeta_obc(i,        &
     &                                 inorth,it2,Linp)
        END DO
      END IF
!
!  2D U-momentum open boundaries.
!
      IF (LBC (iwest,isUbar,ng)%acquire.and.                            &
     &    Lobc(iwest,isUbar,ng).and.                                    &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%ubar_west(j)=BOUNDARY(ng)%ubar_west(j)+          &
     &                              fac1*                               &
     &                              BOUNDARY(ng)%tl_ubar_obc(j,         &
     &                                iwest,it1,Linp)+                  &
     &                              fac2*                               &
     &                              BOUNDARY(ng)%tl_ubar_obc(j,         &
     &                                iwest,it2,Linp)
        END DO
      END IF

      IF (LBC (ieast,isUbar,ng)%acquire.and.                            &
     &    Lobc(ieast,isUbar,ng).and.                                    &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          BOUNDARY(ng)%ubar_east(j)=BOUNDARY(ng)%ubar_east(j)+          &
     &                              fac1*                               &
     &                              BOUNDARY(ng)%tl_ubar_obc(j,         &
     &                                ieast,it1,Linp)+                  &
     &                              fac2*                               &
     &                              BOUNDARY(ng)%tl_ubar_obc(j,         &
     &                                ieast,it2,Linp)
        END DO
      END IF

      IF (LBC (isouth,isUbar,ng)%acquire.and.                           &
     &    Lobc(isouth,isUbar,ng).and.                                   &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrU,Iend
          BOUNDARY(ng)%ubar_south(i)=BOUNDARY(ng)%ubar_south(i)+        &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_ubar_obc(i,        &
     &                                 isouth,it1,Linp)+                &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_ubar_obc(i,        &
     &                                 isouth,it2,Linp)
        END DO
      END IF

      IF (LBC (inorth,isUbar,ng)%acquire.and.                           &
     &    Lobc(inorth,isUbar,ng).and.                                   &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrU,Iend
          BOUNDARY(ng)%ubar_north(i)=BOUNDARY(ng)%ubar_north(i)+        &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_ubar_obc(i,        &
     &                                 inorth,it1,Linp)+                &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_ubar_obc(i,        &
     &                                 inorth,it2,Linp)
        END DO
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (LBC (iwest,isVbar,ng)%acquire.and.                            &
     &    Lobc(iwest,isVbar,ng).and.                                    &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrV,Jend
          BOUNDARY(ng)%vbar_west(j)=BOUNDARY(ng)%vbar_west(j)+          &
     &                              fac1*                               &
     &                              BOUNDARY(ng)%tl_vbar_obc(j,         &
     &                                iwest,it1,Linp)+                  &
     &                              fac2*                               &
     &                              BOUNDARY(ng)%tl_vbar_obc(j,         &
     &                                iwest,it2,Linp)
        END DO
      END IF

      IF (LBC (ieast,isVbar,ng)%acquire.and.                            &
     &    Lobc(ieast,isVbar,ng).and.                                    &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrV,Jend
          BOUNDARY(ng)%vbar_east(j)=BOUNDARY(ng)%vbar_east(j)+          &
     &                              fac1*                               &
     &                              BOUNDARY(ng)%tl_vbar_obc(j,         &
     &                                ieast,it1,Linp)+                  &
     &                              fac2*                               &
     &                              BOUNDARY(ng)%tl_vbar_obc(j,         &
     &                                ieast,it2,Linp)
        END DO
      END IF

      IF (LBC (isouth,isVbar,ng)%acquire.and.                           &
     &    Lobc(isouth,isVbar,ng).and.                                   &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%vbar_south(i)=BOUNDARY(ng)%vbar_south(i)+        &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_vbar_obc(i,        &
     &                                 isouth,it1,Linp)+                &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_vbar_obc(i,        &
     &                                 isouth,it2,Linp)
        END DO
      END IF

      IF (LBC (inorth,isVbar,ng)%acquire.and.                           &
     &    Lobc(inorth,isVbar,ng).and.                                   &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
          BOUNDARY(ng)%vbar_north(i)=BOUNDARY(ng)%vbar_north(i)+        &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_vbar_obc(i,        &
     &                                 inorth,it1,Linp)+                &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_vbar_obc(i,        &
     &                                 inorth,it2,Linp)
        END DO
      END IF

# ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
      IF (LBC (iwest,isUvel,ng)%acquire.and.                            &
     &    Lobc(iwest,isUvel,ng).and.                                    &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            BOUNDARY(ng)%u_west(j,k)=BOUNDARY(ng)%u_west(j,k)+          &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_u_obc(j,k,         &
     &                                 iwest,it1,Linp)+                 &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_u_obc(j,k,         &
     &                                 iwest,it2,Linp)
          END DO
        END DO
      END IF

      IF (LBC (ieast,isUvel,ng)%acquire.and.                            &
     &    Lobc(ieast,isUvel,ng).and.                                    &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            BOUNDARY(ng)%u_east(j,k)=BOUNDARY(ng)%u_east(j,k)+          &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_u_obc(j,k,         &
     &                                 ieast,it1,Linp)+                 &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_u_obc(j,k,         &
     &                                 ieast,it2,Linp)
          END DO
        END DO
      END IF

      IF (LBC (isouth,isUvel,ng)%acquire.and.                           &
     &    Lobc(isouth,isUvel,ng).and.                                   &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BOUNDARY(ng)%u_south(i,k)=BOUNDARY(ng)%u_south(i,k)+        &
     &                                fac1*                             &
     &                                BOUNDARY(ng)%tl_u_obc(i,k,        &
     &                                  isouth,it1,Linp)+               &
     &                                fac2*                             &
     &                                BOUNDARY(ng)%tl_u_obc(i,k,        &
     &                                  isouth,it2,Linp)
          END DO
        END DO
      END IF

      IF (LBC (inorth,isUvel,ng)%acquire.and.                           &
     &    Lobc(inorth,isUvel,ng).and.                                   &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BOUNDARY(ng)%u_north(i,k)=BOUNDARY(ng)%u_north(i,k)+        &
     &                                fac1*                             &
     &                                BOUNDARY(ng)%tl_u_obc(i,k,        &
     &                                  inorth,it1,Linp)+               &
     &                                fac2*                             &
     &                                BOUNDARY(ng)%tl_u_obc(i,k,        &
     &                                  inorth,it2,Linp)
          END DO
        END DO
      END IF
!
!  3D V-momentum open boundaries.
!
      IF (LBC (iwest,isVvel,ng)%acquire.and.                            &
     &    Lobc(iwest,isVvel,ng).and.                                    &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            BOUNDARY(ng)%v_west(j,k)=BOUNDARY(ng)%v_west(j,k)+          &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_v_obc(j,k,         &
     &                                 iwest,it1,Linp)+                 &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_v_obc(j,k,         &
     &                                 iwest,it2,Linp)
          END DO
        END DO
      END IF

      IF (LBC (ieast,isVvel,ng)%acquire.and.                            &
     &    Lobc(ieast,isVvel,ng).and.                                    &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            BOUNDARY(ng)%v_east(j,k)=BOUNDARY(ng)%v_east(j,k)+          &
     &                               fac1*                              &
     &                               BOUNDARY(ng)%tl_v_obc(j,k,         &
     &                                 ieast,it1,Linp)+                 &
     &                               fac2*                              &
     &                               BOUNDARY(ng)%tl_v_obc(j,k,         &
     &                                 ieast,it2,Linp)
          END DO
        END DO
      END IF

      IF (LBC (isouth,isVvel,ng)%acquire.and.                           &
     &    Lobc(isouth,isVvel,ng).and.                                   &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            BOUNDARY(ng)%v_south(i,k)=BOUNDARY(ng)%v_south(i,k)+        &
     &                                fac1*                             &
     &                                BOUNDARY(ng)%tl_v_obc(i,k,        &
     &                                  isouth,it1,Linp)+               &
     &                                fac2*                             &
     &                                BOUNDARY(ng)%tl_v_obc(i,k,        &
     &                                  isouth,it2,Linp)
          END DO
        END DO
      END IF

      IF (LBC (inorth,isVvel,ng)%acquire.and.                           &
     &    Lobc(inorth,isVvel,ng).and.                                   &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            BOUNDARY(ng)%v_north(i,k)=BOUNDARY(ng)%v_north(i,k)+        &
     &                                fac1*                             &
     &                                BOUNDARY(ng)%tl_v_obc(i,k,        &
     &                                  inorth,it1,Linp)+               &
     &                                fac2*                             &
     &                                BOUNDARY(ng)%tl_v_obc(i,k,        &
     &                                  inorth,it2,Linp)
          END DO
        END DO
      END IF
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (LBC (iwest,isTvar(it),ng)%acquire.and.                      &
     &      Lobc(iwest,isTvar(it),ng).and.                              &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              BOUNDARY(ng)%t_west(j,k,it)=BOUNDARY(ng)%t_west(j,k,it)+  &
     &                                    fac1*                         &
     &                                    BOUNDARY(ng)%tl_t_obc(j,k,    &
     &                                      iwest,it1,Linp,it)+         &
     &                                    fac2*                         &
     &                                    BOUNDARY(ng)%tl_t_obc(j,k,    &
     &                                      iwest,it2,Linp,it)
            END DO
          END DO
        END IF

        IF (LBC (ieast,isTvar(it),ng)%acquire.and.                      &
     &      Lobc(ieast,isTvar(it),ng).and.                              &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              BOUNDARY(ng)%t_east(j,k,it)=BOUNDARY(ng)%t_east(j,k,it)+  &
     &                                    fac1*                         &
     &                                    BOUNDARY(ng)%tl_t_obc(j,k,    &
     &                                      ieast,it1,Linp,it)+         &
     &                                    fac2*                         &
     &                                    BOUNDARY(ng)%tl_t_obc(j,k,    &
     &                                      ieast,it2,Linp,it)
            END DO
          END DO
        END IF

        IF (LBC (isouth,isTvar(it),ng)%acquire.and.                     &
     &      Lobc(isouth,isTvar(it),ng).and.                             &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              BOUNDARY(ng)%t_south(i,k,it)=BOUNDARY(ng)%t_south(i,k,it)+&
     &                                     fac1*                        &
     &                                     BOUNDARY(ng)%tl_t_obc(i,k,   &
     &                                       isouth,it1,Linp,it)+       &
     &                                     fac2*                        &
     &                                     BOUNDARY(ng)%tl_t_obc(i,k,   &
     &                                       isouth,it2,Linp,it)
            END DO
          END DO
        END IF

        IF (LBC (inorth,isTvar(it),ng)%acquire.and.                     &
     &      Lobc(inorth,isTvar(it),ng).and.                             &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              BOUNDARY(ng)%t_north(i,k,it)=BOUNDARY(ng)%t_north(i,k,it)+&
     &                                     fac1*                        &
     &                                     BOUNDARY(ng)%tl_t_obc(i,k,   &
     &                                       inorth,it1,Linp,it)+       &
     &                                     fac2*                        &
     &                                     BOUNDARY(ng)%tl_t_obc(i,k,   &
     &                                       inorth,it2,Linp,it)
            END DO
          END DO
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
     &                    Lout)
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
     &                          Lout)
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
      integer, intent(in) :: Lout
!
!  Local variable declarations.
!
      integer :: i, ic, j
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
        ic=OBCcount(ng)
!
!  Free-surface open boundaries.
!
        IF (LBC (iwest,isFsur,ng)%acquire.and.                          &
     &      Lobc(iwest,isFsur,ng).and.                                  &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            BOUNDARY(ng)%zeta_obc(j,iwest,ic,Lout)=                     &
     &                              BOUNDARY(ng)%zeta_west(j)
          END DO
        END IF

        IF (LBC (ieast,isFsur,ng)%acquire.and.                          &
     &      Lobc(ieast,isFsur,ng).and.                                  &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            BOUNDARY(ng)%zeta_obc(j,ieast,ic,Lout)=                     &
     &                              BOUNDARY(ng)%zeta_east(j)
          END DO
        END IF

        IF (LBC (isouth,isFsur,ng)%acquire.and.                         &
     &      Lobc(isouth,isFsur,ng).and.                                 &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            BOUNDARY(ng)%zeta_obc(i,isouth,ic,Lout)=                    &
     &                              BOUNDARY(ng)%zeta_south(i)
          END DO
        END IF

        IF (LBC (inorth,isFsur,ng)%acquire.and.                         &
     &      Lobc(inorth,isFsur,ng).and.                                 &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            BOUNDARY(ng)%zeta_obc(i,inorth,ic,Lout)=                    &
     &                              BOUNDARY(ng)%zeta_north(i)
          END DO
        END IF
!
!  2D U-momentum open boundaries.
!
        IF (LBC (iwest,isUbar,ng)%acquire.and.                          &
     &      Lobc(iwest,isUbar,ng).and.                                  &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            BOUNDARY(ng)%ubar_obc(j,iwest,ic,Lout)=                     &
     &                              BOUNDARY(ng)%ubar_west(j)
          END DO
        END IF

        IF (LBC (ieast,isUbar,ng)%acquire.and.                          &
     &      Lobc(ieast,isUbar,ng).and.                                  &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            BOUNDARY(ng)%ubar_obc(j,ieast,ic,Lout)=                     &
     &                              BOUNDARY(ng)%ubar_east(j)
          END DO
        END IF

        IF (LBC (isouth,isUbar,ng)%acquire.and.                         &
     &      Lobc(isouth,isUbar,ng).and.                                 &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            BOUNDARY(ng)%ubar_obc(i,isouth,ic,Lout)=                    &
     &                              BOUNDARY(ng)%ubar_south(i)
          END DO
        END IF

        IF (LBC (inorth,isUbar,ng)%acquire.and.                         &
     &      Lobc(inorth,isUbar,ng).and.                                 &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            BOUNDARY(ng)%ubar_obc(i,inorth,ic,Lout)=                    &
     &                              BOUNDARY(ng)%ubar_north(i)
          END DO
        END IF
!
!  2D V-momentum open boundaries.
!
        IF (LBC (iwest,isVbar,ng)%acquire.and.                          &
     &      Lobc(iwest,isVbar,ng).and.                                  &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            BOUNDARY(ng)%vbar_obc(j,iwest,ic,Lout)=                     &
     &                              BOUNDARY(ng)%vbar_west(j)
          END DO
        END IF

        IF (LBC (ieast,isVbar,ng)%acquire.and.                          &
     &      Lobc(ieast,isVbar,ng).and.                                  &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            BOUNDARY(ng)%vbar_obc(j,ieast,ic,Lout)=                     &
     &                              BOUNDARY(ng)%vbar_east(j)
          END DO
        END IF

        IF (LBC (isouth,isVbar,ng)%acquire.and.                         &
     &      Lobc(isouth,isVbar,ng).and.                                 &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            BOUNDARY(ng)%vbar_obc(i,isouth,ic,Lout)=                    &
     &                              BOUNDARY(ng)%vbar_south(i)
          END DO
        END IF

        IF (LBC (inorth,isVbar,ng)%acquire.and.                         &
     &      Lobc(inorth,isVbar,ng).and.                                 &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            BOUNDARY(ng)%vbar_obc(i,inorth,ic,Lout)=                    &
     &                              BOUNDARY(ng)%vbar_north(i)
          END DO
        END IF

# ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
        IF (LBC (iwest,isUvel,ng)%acquire.and.                          &
     &      Lobc(iwest,isUvel,ng).and.                                  &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              BOUNDARY(ng)%u_obc(j,k,iwest,ic,Lout)=                    &
     &                               BOUNDARY(ng)%u_west(j,k)
            END DO
          END DO
        END IF

        IF (LBC (ieast,isUvel,ng)%acquire.and.                          &
     &      Lobc(ieast,isUvel,ng).and.                                  &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              BOUNDARY(ng)%u_obc(j,k,ieast,ic,Lout)=                    &
     &                               BOUNDARY(ng)%u_east(j,k)
            END DO
          END DO
        END IF

        IF (LBC (isouth,isUvel,ng)%acquire.and.                         &
     &      Lobc(isouth,isUvel,ng).and.                                 &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              BOUNDARY(ng)%u_obc(i,k,isouth,ic,Lout)=                   &
     &                               BOUNDARY(ng)%u_south(i,k)
            END DO
          END DO
        END IF

        IF (LBC (inorth,isUvel,ng)%acquire.and.                         &
     &      Lobc(inorth,isUvel,ng).and.                                 &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              BOUNDARY(ng)%u_obc(i,k,inorth,ic,Lout)=                   &
     &                               BOUNDARY(ng)%u_north(i,k)
            END DO
          END DO
        END IF
!
!  3D V-momentum open boundaries.
!
        IF (LBC (iwest,isVvel,ng)%acquire.and.                          &
     &      Lobc(iwest,isVvel,ng).and.                                  &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              BOUNDARY(ng)%v_obc(j,k,iwest,ic,Lout)=                    &
     &                               BOUNDARY(ng)%v_west(j,k)
            END DO
          END DO
        END IF

        IF (LBC (ieast,isVvel,ng)%acquire.and.                          &
     &      Lobc(ieast,isVvel,ng).and.                                  &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              BOUNDARY(ng)%v_obc(j,k,ieast,ic,Lout)=                    &
     &                               BOUNDARY(ng)%v_east(j,k)
            END DO
          END DO
        END IF

        IF (LBC (isouth,isVvel,ng)%acquire.and.                         &
     &      Lobc(isouth,isVvel,ng).and.                                 &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              BOUNDARY(ng)%v_obc(i,k,isouth,ic,Lout)=                   &
     &                               BOUNDARY(ng)%v_south(i,k)
            END DO
          END DO
        END IF

        IF (LBC (inorth,isVvel,ng)%acquire.and.                         &
     &      Lobc(inorth,isVvel,ng).and.                                 &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              BOUNDARY(ng)%v_obc(i,k,inorth,ic,Lout)=                   &
     &                               BOUNDARY(ng)%v_north(i,k)
            END DO
          END DO
        END IF
!
!  Tracers open boundaries.
!
        DO it=1,NT(ng)
          IF (LBC (iwest,isTvar(it),ng)%acquire.and.                    &
     &        Lobc(iwest,isTvar(it),ng).and.                            &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                BOUNDARY(ng)%t_obc(j,k,iwest,ic,Lout,it)=               &
     &                                 BOUNDARY(ng)%t_west(j,k,it)
              END DO
            END DO
          END IF

          IF (LBC (ieast,isTvar(it),ng)%acquire.and.                    &
     &        Lobc(ieast,isTvar(it),ng).and.                            &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                BOUNDARY(ng)%t_obc(j,k,ieast,ic,Lout,it)=               &
     &                                 BOUNDARY(ng)%t_east(j,k,it)
              END DO
            END DO
          END IF

          IF (LBC (isouth,isTvar(it),ng)%acquire.and.                   &
     &        Lobc(isouth,isTvar(it),ng).and.                           &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                BOUNDARY(ng)%t_obc(i,k,isouth,ic,Lout,it)=              &
     &                                 BOUNDARY(ng)%t_south(i,k,it)
              END DO
            END DO
          END IF

          IF (LBC (inorth,isTvar(it),ng)%acquire.and.                   &
     &        Lobc(inorth,isTvar(it),ng).and.                           &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                BOUNDARY(ng)%t_obc(i,k,inorth,ic,Lout,it)=              &
     &                                 BOUNDARY(ng)%t_north(i,k,it)
              END DO
            END DO
          END IF
        END DO
# endif
      END IF LOAD

      RETURN
      END SUBROUTINE load_obc_tile
#endif
      END MODULE obc_adjust_mod
