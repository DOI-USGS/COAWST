#include "cppdefs.h"
      MODULE state_copy_mod
!
!svn $Id: state_copy.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine makes a copies of the model state as follows:          !
!                                                                      !
!      s1_var(...,Lout) = s2_var(...,Linp)                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC  :: state_copy

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_copy (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Linp, Lout,                                &
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                       s1_t_obc, s2_t_obc,                        &
     &                       s1_u_obc, s2_u_obc,                        &
     &                       s1_v_obc, s2_v_obc,                        &
# endif
     &                       s1_ubar_obc, s2_ubar_obc,                  &
     &                       s1_vbar_obc, s2_vbar_obc,                  &
     &                       s1_zeta_obc, s2_zeta_obc,                  &
#endif
#ifdef ADJUST_WSTRESS
     &                       s1_sustr, s2_sustr,                        &
     &                       s1_svstr, s2_svstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       s1_tflux, s2_tflux,                        &
# endif
     &                       s1_t, s2_t,                                &
     &                       s1_u, s2_u,                                &
     &                       s1_v, s2_v,                                &
# if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       s1_ubar, s2_ubar,                          &
     &                       s1_vbar, s2_vbar,                          &
# endif
#else
     &                       s1_ubar, s2_ubar,                          &
     &                       s1_vbar, s2_vbar,                          &
#endif
     &                       s1_zeta, s2_zeta)
!***********************************************************************
!
      USE mod_param
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
      USE mod_ncparam
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Linp, Lout
!
#ifdef ASSUMED_SHAPE
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(in) :: s2_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: s2_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: s2_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(in) :: s2_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: s2_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: s2_zeta_obc(LBij:,:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: s2_sustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_svstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(in) :: s2_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(in) :: s2_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: s2_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_v(LBi:,LBj:,:,:)
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: s2_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_vbar(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: s2_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:,:)

# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s1_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s1_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s1_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: s1_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s1_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s1_zeta_obc(LBij:,:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_svstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s1_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s1_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_v(LBi:,LBj:,:,:)
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: s1_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s1_vbar(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(inout) :: s1_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s1_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s1_zeta(LBi:,LBj:,:)

#else

# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(in) :: s2_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: s2_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: s2_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(in) :: s2_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: s2_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: s2_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: s2_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: s2_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(in) :: s2_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(in) :: s2_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: s2_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: s2_v(LBi:UBi,LBj:UBj,N(ng),2)
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: s2_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: s2_vbar(LBi:UBi,LBj:UBj,3)
#  endif
# else
      real(r8), intent(in) :: s2_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: s2_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(in) :: s2_zeta(LBi:UBi,LBj:UBj,3)

# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s1_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s1_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s1_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: s1_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s1_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s1_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s1_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s1_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s1_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s1_v(LBi:UBi,LBj:UBj,N(ng),2)
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: s1_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s1_vbar(LBi:UBi,LBj:UBj,3)
#  endif
# else
      real(r8), intent(inout) :: s1_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s1_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s1_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ib, ir, it

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute the following operation between S1 and S2 model state
!  trajectories:
!                 S1(Lout) = fac1 * S1(Linp1) + fac2 * S2(Linp2)
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          s1_zeta(i,j,Lout)=s2_zeta(i,j,Linp)
        END DO
      END DO

#ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              s1_zeta_obc(j,ib,ir,Lout)=s2_zeta_obc(j,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              s1_zeta_obc(j,ib,ir,Lout)=s2_zeta_obc(j,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              s1_zeta_obc(i,ib,ir,Lout)=s2_zeta_obc(i,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              s1_zeta_obc(i,ib,ir,Lout)=s2_zeta_obc(i,ib,ir,Linp)
            END DO
          END IF
        END DO
      END IF
#endif

#if !defined SOLVE3D || (defined WEAK_CONSTRAINT && defined TIME_CONV)
!
!  2D U-momentum component.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          s1_ubar(i,j,Lout)=s2_ubar(i,j,Linp)
        END DO
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              s1_ubar_obc(j,ib,ir,Lout)=s2_ubar_obc(j,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              s1_ubar_obc(j,ib,ir,Lout)=s2_ubar_obc(j,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              s1_ubar_obc(i,ib,ir,Lout)=s2_ubar_obc(i,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              s1_ubar_obc(i,ib,ir,Lout)=s2_ubar_obc(i,ib,ir,Linp)
            END DO
          END IF
        END DO
      END IF
#endif

#if !defined SOLVE3D || (defined WEAK_CONSTRAINT && defined TIME_CONV)
!
!  2D V-momentum component.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          s1_vbar(i,j,Lout)=s2_vbar(i,j,Linp)
        END DO
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=JstrV,Jend
              s1_vbar_obc(j,ib,ir,Lout)=s2_vbar_obc(j,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              s1_vbar_obc(j,ib,ir,Lout)=s2_vbar_obc(j,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              s1_vbar_obc(i,ib,ir,Lout)=s2_vbar_obc(i,ib,ir,Linp)
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              s1_vbar_obc(i,ib,ir,Lout)=s2_vbar_obc(i,ib,ir,Linp)
            END DO
          END IF
        END DO
      END IF
#endif

#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO ir=1,Nfrec(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            s1_sustr(i,j,ir,Lout)=s2_sustr(i,j,ir,Linp)
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            s1_svstr(i,j,ir,Lout)=s2_svstr(i,j,ir,Linp)
          END DO
        END DO
      END DO
#endif

#ifdef SOLVE3D
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            s1_u(i,j,k,Lout)=s2_u(i,j,k,Linp)
          END DO
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
                s1_u_obc(j,k,ib,ir,Lout)=s2_u_obc(j,k,ib,ir,Linp)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                s1_u_obc(j,k,ib,ir,Lout)=s2_u_obc(j,k,ib,ir,Linp)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s1_u_obc(i,k,ib,ir,Lout)=s2_u_obc(i,k,ib,ir,Linp)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s1_u_obc(i,k,ib,ir,Lout)=s2_u_obc(i,k,ib,ir,Linp)
              END DO
            END DO
          END IF
        END DO
      END IF
# endif
!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            s1_v(i,j,k,Lout)=s2_v(i,j,k,Linp)
          END DO
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=JstrV,Jend
                s1_v_obc(j,k,ib,ir,Lout)=s2_v_obc(j,k,ib,ir,Linp)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                s1_v_obc(j,k,ib,ir,Lout)=s2_v_obc(j,k,ib,ir,Linp)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                s1_v_obc(i,k,ib,ir,Lout)=s2_v_obc(i,k,ib,ir,Linp)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                s1_v_obc(i,k,ib,ir,Lout)=s2_v_obc(i,k,ib,ir,Linp)
              END DO
            END DO
          END IF
        END DO
      END IF
# endif
!
!  Tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              s1_t(i,j,k,Lout,it)=s2_t(i,j,k,Linp,it)
            END DO
          END DO
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  s1_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                                       s2_t_obc(j,k,ib,ir,Linp,it)
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  s1_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                                       s2_t_obc(j,k,ib,ir,Linp,it)
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s1_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                                       s2_t_obc(i,k,ib,ir,Linp,it)
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s1_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                                       s2_t_obc(i,k,ib,ir,Linp,it)
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
# endif

# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          DO ir=1,Nfrec(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                s1_tflux(i,j,ir,Lout,it)=s2_tflux(i,j,ir,Linp,it)
              END DO
            END DO
          END DO
        END IF
      END DO
# endif

#endif

      RETURN
      END SUBROUTINE state_copy

      END MODULE state_copy_mod
