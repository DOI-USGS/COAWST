#include "cppdefs.h"
      MODULE state_initialize_mod
!
!svn $Id: state_initialize.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initialize the state variable as follows:              !
!                                                                      !
!     s_var(...,Lout) = fac                                            !
!                                                                      !
!  where fac is scalar usually zero.                                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC  :: state_initialize

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_initialize (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             Lout, fac,                           &
#ifdef MASKING
     &                             rmask, umask, vmask,                 &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                             s_t_obc,                             &
     &                             s_u_obc, s_v_obc,                    &
# endif
     &                             s_ubar_obc, s_vbar_obc,              &
     &                             s_zeta_obc,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                             s_sustr, s_svstr,                    &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                             s_tflux,                             &
# endif
     &                             s_t, s_u, s_v,                       &
#else
     &                             s_ubar, s_vbar,                      &
#endif
     &                             s_zeta)
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
      integer, intent(in) :: Lout
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
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          s_zeta(i,j,Lout)=fac
        END DO
      END DO

#ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.WESTERN_EDGE) THEN
            DO j=Jstr,Jend
              s_zeta_obc(j,iwest,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
            DO j=Jstr,Jend
              s_zeta_obc(j,ieast,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO i=Istr,Iend
              s_zeta_obc(i,isouth,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
            DO i=Istr,Iend
              s_zeta_obc(i,inorth,ir,Lout)=fac
            END DO
          END IF
        END DO
      END IF
#endif

#ifndef SOLVE3D
!
!  2D U-momentum component.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          s_ubar(i,j,Lout)=fac
        END DO
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
            DO j=Jstr,Jend
              s_ubar_obc(j,iwest,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
            DO j=Jstr,Jend
              s_ubar_obc(j,ieast,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
            DO i=IstrU,Iend
              s_ubar_obc(i,isouth,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
            DO i=IstrU,Iend
              s_ubar_obc(i,inorth,ir,Lout)=fac
            END DO
          END IF
        END DO
      END IF
#endif

#ifndef SOLVE3D
!
!  2D V-momentum component.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          s_vbar(i,j,Lout)=fac
        END DO
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
            DO j=JstrV,Jend
              s_vbar_obc(j,iwest,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
            DO j=JstrV,Jend
              s_vbar_obc(j,ieast,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
            DO i=Istr,Iend
              s_vbar_obc(i,isouth,ir,Lout)=fac
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
            DO i=Istr,Iend
              s_vbar_obc(i,inorth,ir,Lout)=fac
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
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s_sustr(i,j,ir,Lout)=fac
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s_svstr(i,j,ir,Lout)=fac
          END DO
        END DO
      END DO
#endif

#ifdef SOLVE3D
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s_u(i,j,k,Lout)=fac
          END DO
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUvel,ng)).and.WESTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                s_u_obc(j,k,iwest,ir,Lout)=fac
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                s_u_obc(j,k,ieast,ir,Lout)=fac
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s_u_obc(i,k,isouth,ir,Lout)=fac
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s_u_obc(i,k,inorth,ir,Lout)=fac
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
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s_v(i,j,k,Lout)=fac
          END DO
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVvel,ng)).and.WESTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                s_v_obc(j,k,iwest,ir,Lout)=fac
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                s_v_obc(j,k,ieast,ir,Lout)=fac
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                s_v_obc(i,k,isouth,ir,Lout)=fac
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                s_v_obc(i,k,inorth,ir,Lout)=fac
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
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              s_t(i,j,k,Lout,it)=fac
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
            IF ((Lobc(iwest,isTvar(it),ng)).and.WESTERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  s_t_obc(j,k,iwest,ir,Lout,it)=fac
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  s_t_obc(j,k,ieast,ir,Lout,it)=fac
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s_t_obc(i,k,isouth,ir,Lout,it)=fac
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s_t_obc(i,k,inorth,ir,Lout,it)=fac
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
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                s_tflux(i,j,ir,Lout,it)=fac
              END DO
            END DO
          END DO
        END IF
      END DO
# endif

#endif

      RETURN
      END SUBROUTINE state_initialize

      END MODULE state_initialize_mod
