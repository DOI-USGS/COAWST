#include "cppdefs.h"
      MODULE state_scale_mod
!
!svn $Id: state_scale.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine scales a state variable by a constant:                 !
!                                                                      !
!     s_var(...,Lout) = fac * s_var(...,Linp)                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC  :: state_scale

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_scale (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        Linp, Lout, fac,                          &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                        s_t_obc,                                  &
     &                        s_u_obc, s_v_obc,                         &
# endif
     &                        s_ubar_obc, s_vbar_obc,                   &
     &                        s_zeta_obc,                               &
#endif
#ifdef ADJUST_WSTRESS
     &                        s_sustr, s_svstr,                         &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        s_tflux,                                  &
# endif
     &                        s_t, s_u, s_v,                            &
#else
     &                        s_ubar, s_vbar,                           &
#endif
     &                        s_zeta)
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
      real(r8), intent(in) :: fac
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_sustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_svstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ib, ir, it

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Scale model state variable by a constant.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          s_zeta(i,j,Lout)=fac*s_zeta(i,j,Linp)
#ifdef MASKING
          s_zeta(i,j,Lout)=s_zeta(i,j,Lout)*rmask(i,j)
#endif
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
              s_zeta_obc(j,ib,ir,Lout)=fac*s_zeta_obc(j,ib,ir,Linp)
# ifdef MASKING
              s_zeta_obc(j,ib,ir,Lout)=s_zeta_obc(j,ib,ir,Lout)*        &
     &                                 rmask(Istr-1,j)
# endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              s_zeta_obc(j,ib,ir,Lout)=fac*s_zeta_obc(j,ib,ir,Linp)
# ifdef MASKING
              s_zeta_obc(j,ib,ir,Lout)=s_zeta_obc(j,ib,ir,Lout)*        &
     &                                 rmask(Iend+1,j)
# endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              s_zeta_obc(i,ib,ir,Lout)=fac*s_zeta_obc(i,ib,ir,Linp)
# ifdef MASKING
              s_zeta_obc(i,ib,ir,Lout)=s_zeta_obc(i,ib,ir,Lout)*        &
     &                                 rmask(i,Jstr-1)
# endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              s_zeta_obc(i,ib,ir,Lout)=fac*s_zeta_obc(i,ib,ir,Linp)
# ifdef MASKING
              s_zeta_obc(i,ib,ir,Lout)=s_zeta_obc(i,ib,ir,Lout)*        &
     &                                 rmask(i,Jend+1)
# endif
            END DO
          END IF
        END DO
      END IF
#endif

#ifndef SOLVE3D
!
!  2D U-momentum component.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          s_ubar(i,j,Lout)=fac*s_ubar(i,j,Linp)
# ifdef MASKING
          s_ubar(i,j,Lout)=s_ubar(i,j,Lout)*umask(i,j)
# endif
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
              s_ubar_obc(j,ib,ir,Lout)=fac*s_ubar_obc(j,ib,ir,Linp)
# ifdef MASKING
              s_ubar_obc(j,ib,ir,Lout)=s_ubar_obc(j,ib,ir,Lout)*        &
     &                                 umask(Istr,j)
# endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              s_ubar_obc(j,ib,ir,Lout)=fac*s_ubar_obc(j,ib,ir,Linp)
# ifdef MASKING
              s_ubar_obc(j,ib,ir,Lout)=s_ubar_obc(j,ib,ir,Lout)*        &
     &                                 umask(Iend+1,j)
# endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              s_ubar_obc(i,ib,ir,Lout)=fac*s_ubar_obc(i,ib,ir,Linp)
# ifdef MASKING
              s_ubar_obc(i,ib,ir,Lout)=s_ubar_obc(i,ib,ir,Lout)*        &
     &                                 umask(i,Jstr-1)
# endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              s_ubar_obc(i,ib,ir,Lout)=fac*s_ubar_obc(i,ib,ir,Linp)
# ifdef MASKING
              s_ubar_obc(i,ib,ir,Lout)=s_ubar_obc(i,ib,ir,Lout)*        &
     &                                 umask(i,Jend+1)
# endif
            END DO
          END IF
        END DO
      END IF
#endif

#ifndef SOLVE3D
!
!  2D V-momentum component.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          s_vbar(i,j,Lout)=fac*s_vbar(i,j,Linp)
# ifdef MASKING
          s_vbar(i,j,Lout)=s_vbar(i,j,Lout)*vmask(i,j)
# endif
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
              s_vbar_obc(j,ib,ir,Lout)=fac*s_vbar_obc(j,ib,ir,Linp)
# ifdef MASKING
              s_vbar_obc(j,ib,ir,Lout)=s_vbar_obc(j,ib,ir,Lout)*        &
     &                                 vmask(Istr-1,j)
# endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              s_vbar_obc(j,ib,ir,Lout)=fac*s_vbar_obc(j,ib,ir,Linp)
# ifdef MASKING
              s_vbar_obc(j,ib,ir,Lout)=s_vbar_obc(j,ib,ir,Lout)*        &
     &                                 vmask(Iend+1,j)
# endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              s_vbar_obc(i,ib,ir,Lout)=fac*s_vbar_obc(i,ib,ir,Linp)
# ifdef MASKING
              s_vbar_obc(i,ib,ir,Lout)=s_vbar_obc(i,ib,ir,Lout)*        &
     &                                 vmask(i,Jstr)
# endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              s_vbar_obc(i,ib,ir,Lout)=fac*s_vbar_obc(i,ib,ir,Linp)
# ifdef MASKING
              s_vbar_obc(i,ib,ir,Lout)=s_vbar_obc(i,ib,ir,Lout)*        &
     &                                 vmask(i,Jend+1)
# endif
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
            s_sustr(i,j,ir,Lout)=fac*s_sustr(i,j,ir,Linp)
# ifdef MASKING
            s_sustr(i,j,ir,Lout)=s_sustr(i,j,ir,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            s_svstr(i,j,ir,Lout)=fac*s_svstr(i,j,ir,Linp)
# ifdef MASKING
            s_svstr(i,j,ir,Lout)=s_svstr(i,j,ir,Lout)*vmask(i,j)
# endif
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
            s_u(i,j,k,Lout)=fac*s_u(i,j,k,Linp)
# ifdef MASKING
            s_u(i,j,k,Lout)=s_u(i,j,k,Lout)*umask(i,j)
# endif
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
                s_u_obc(j,k,ib,ir,Lout)=fac*s_u_obc(j,k,ib,ir,Linp)
#  ifdef MASKING
                s_u_obc(j,k,ib,ir,Lout)=s_u_obc(j,k,ib,ir,Lout)*        &
     &                                  umask(Istr,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                s_u_obc(j,k,ib,ir,Lout)=fac*s_u_obc(j,k,ib,ir,Linp)
#  ifdef MASKING
                s_u_obc(j,k,ib,ir,Lout)=s_u_obc(j,k,ib,ir,Lout)*        &
     &                                  umask(Iend+1,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s_u_obc(i,k,ib,ir,Lout)=fac*s_u_obc(i,k,ib,ir,Linp)
#  ifdef MASKING
                s_u_obc(i,k,ib,ir,Lout)=s_u_obc(i,k,ib,ir,Lout)*        &
     &                                  umask(i,Jstr-1)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s_u_obc(i,k,ib,ir,Lout)=fac*s_u_obc(i,k,ib,ir,Linp)
#  ifdef MASKING
                s_u_obc(i,k,ib,ir,Lout)=s_u_obc(i,k,ib,ir,Lout)*        &
     &                                  umask(i,Jend+1)
#  endif
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
            s_v(i,j,k,Lout)=fac*s_v(i,j,k,Linp)
# ifdef MASKING
            s_v(i,j,k,Lout)=s_v(i,j,k,Lout)*vmask(i,j)
# endif
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
                s_v_obc(j,k,ib,ir,Lout)=fac*s_v_obc(j,k,ib,ir,Linp)
#  ifdef MASKING
                s_v_obc(j,k,ib,ir,Lout)=s_v_obc(j,k,ib,ir,Lout)*        &
     &                                  vmask(Istr-1,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                s_v_obc(j,k,ib,ir,Lout)=fac*s_v_obc(j,k,ib,ir,Linp)
#  ifdef MASKING
                s_v_obc(j,k,ib,ir,Lout)=s_v_obc(j,k,ib,ir,Lout)*        &
     &                                  vmask(Iend+1,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                s_v_obc(i,k,ib,ir,Lout)=fac*s_v_obc(i,k,ib,ir,Linp)
#  ifdef MASKING
                s_v_obc(i,k,ib,ir,Lout)=s_v_obc(i,k,ib,ir,Lout)*        &
     &                                  vmask(i,Jstr)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                s_v_obc(i,k,ib,ir,Lout)=fac*s_v_obc(i,k,ib,ir,Linp)
#  ifdef MASKING
                s_v_obc(i,k,ib,ir,Lout)=s_v_obc(i,k,ib,ir,Lout)*        &
     &                                  vmask(i,Jend+1)
#  endif
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
              s_t(i,j,k,Lout,it)=fac*s_t(i,j,k,Linp,it)
# ifdef MASKING
              s_t(i,j,k,Lout,it)=s_t(i,j,k,Lout,it)*rmask(i,j)
# endif
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
                  s_t_obc(j,k,ib,ir,Lout,it)=fac*                       &
     &                                       s_t_obc(j,k,ib,ir,Linp,it)
#  ifdef MASKING
                  s_t_obc(j,k,ib,ir,Lout,it)=s_t_obc(j,k,ib,ir,Lout,it)*&
     &                                       rmask(Istr-1,j)
#  endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  s_t_obc(j,k,ib,ir,Lout,it)=fac*                       &
     &                                       s_t_obc(j,k,ib,ir,Linp,it)
#  ifdef MASKING
                  s_t_obc(j,k,ib,ir,Lout,it)=s_t_obc(j,k,ib,ir,Lout,it)*&
     &                                       rmask(Iend+1,j)
#  endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s_t_obc(i,k,ib,ir,Lout,it)=fac*                       &
     &                                       s_t_obc(i,k,ib,ir,Linp,it)
#  ifdef MASKING
                  s_t_obc(i,k,ib,ir,Lout,it)=s_t_obc(i,k,ib,ir,Lout,it)*&
     &                                       rmask(i,Jstr-1)
#  endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s_t_obc(i,k,ib,ir,Lout,it)=fac*                       &
     &                                       s_t_obc(i,k,ib,ir,Linp,it)
#  ifdef MASKING
                  s_t_obc(i,k,ib,ir,Lout,it)=s_t_obc(i,k,ib,ir,Lout,it)*&
     &                                       rmask(i,Jend+1)
#  endif
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
                s_tflux(i,j,ir,Lout,it)=fac*s_tflux(i,j,ir,Linp,it)
#  ifdef MASKING
                s_tflux(i,j,ir,Lout,it)=s_tflux(i,j,ir,Lout,it)*        &
     &                                  rmask(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE state_scale

      END MODULE state_scale_mod
