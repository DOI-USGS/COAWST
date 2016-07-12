#include "cppdefs.h"
      MODULE state_addition_mod
!
!svn $Id: state_addition.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the following model state addition:           !
!                                                                      !
!      s1_var(...,Lout) = fac1 * s1_var(...,Lin1) +                    !
!                         fac2 * s2_var(...,Lin2)                      !
!                                                                      !
!  where fac1 and fac2 are scalars.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC  :: state_addition

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_addition (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           Lin1, Lin2, Lout,                      &
     &                           fac1, fac2,                            &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                           s1_t_obc, s2_t_obc,                    &
     &                           s1_u_obc, s2_u_obc,                    &
     &                           s1_v_obc, s2_v_obc,                    &
# endif
     &                           s1_ubar_obc, s2_ubar_obc,              &
     &                           s1_vbar_obc, s2_vbar_obc,              &
     &                           s1_zeta_obc, s2_zeta_obc,              &
#endif
#ifdef ADJUST_WSTRESS
     &                           s1_sustr, s2_sustr,                    &
     &                           s1_svstr, s2_svstr,                    &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           s1_tflux, s2_tflux,                    &
# endif
     &                           s1_t, s2_t,                            &
     &                           s1_u, s2_u,                            &
     &                           s1_v, s2_v,                            &
# if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           s1_ubar, s2_ubar,                      &
     &                           s1_vbar, s2_vbar,                      &
# endif
#else
     &                           s1_ubar, s2_ubar,                      &
     &                           s1_vbar, s2_vbar,                      &
#endif
     &                           s1_zeta, s2_zeta)
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
      integer, intent(in) :: Lin1, Lin2, Lout
!
      real(r8), intent(in) :: fac1, fac2
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
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

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
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
!                 S1(Lout) = fac1 * S1(Lin1) + fac2 * S2(Lin2)
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          s1_zeta(i,j,Lout)=fac1*s1_zeta(i,j,Lin1)+                     &
     &                      fac2*s2_zeta(i,j,Lin2)
#ifdef MASKING
          s1_zeta(i,j,Lout)=s1_zeta(i,j,Lout)*rmask(i,j)
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
              s1_zeta_obc(j,ib,ir,Lout)=fac1*s1_zeta_obc(j,ib,ir,Lin1)+ &
     &                                  fac2*s2_zeta_obc(j,ib,ir,Lin2)
# ifdef MASKING
              s1_zeta_obc(j,ib,ir,Lout)=s1_zeta_obc(j,ib,ir,Lout)*      &
     &                                  rmask(Istr-1,j)
# endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              s1_zeta_obc(j,ib,ir,Lout)=fac1*s1_zeta_obc(j,ib,ir,Lin1)+ &
     &                                  fac2*s2_zeta_obc(j,ib,ir,Lin2)
# ifdef MASKING
              s1_zeta_obc(j,ib,ir,Lout)=s1_zeta_obc(j,ib,ir,Lout)*      &
     &                                  rmask(Iend+1,j)
# endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              s1_zeta_obc(i,ib,ir,Lout)=fac1*s1_zeta_obc(i,ib,ir,Lin1)+ &
     &                                  fac2*s2_zeta_obc(i,ib,ir,Lin2)
# ifdef MASKING
              s1_zeta_obc(i,ib,ir,Lout)=s1_zeta_obc(i,ib,ir,Lout)*      &
     &                                  rmask(i,Jstr-1)
# endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              s1_zeta_obc(i,ib,ir,Lout)=fac1*s1_zeta_obc(i,ib,ir,Lin1)+ &
     &                                  fac2*s2_zeta_obc(i,ib,ir,Lin2)
# ifdef MASKING
              s1_zeta_obc(i,ib,ir,Lout)=s1_zeta_obc(i,ib,ir,Lout)*      &
     &                                  rmask(i,Jend+1)
# endif
            END DO
          END IF
        END DO
      END IF
#endif

#if !defined SOLVE3D || (defined WEAK_CONSTRAINT && defined TIME_CONV )
!
!  2D U-momentum component.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          s1_ubar(i,j,Lout)=fac1*s1_ubar(i,j,Lin1)+                     &
     &                      fac2*s2_ubar(i,j,Lin2)
# ifdef MASKING
          s1_ubar(i,j,Lout)=s1_ubar(i,j,Lout)*umask(i,j)
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
              s1_ubar_obc(j,ib,ir,Lout)=fac1*s1_ubar_obc(j,ib,ir,Lin1)+ &
     &                                  fac2*s2_ubar_obc(j,ib,ir,Lin2)
# ifdef MASKING
              s1_ubar_obc(j,ib,ir,Lout)=s1_ubar_obc(j,ib,ir,Lout)*      &
     &                                  umask(Istr,j)
# endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              s1_ubar_obc(j,ib,ir,Lout)=fac1*s1_ubar_obc(j,ib,ir,Lin1)+ &
     &                                  fac2*s2_ubar_obc(j,ib,ir,Lin2)
# ifdef MASKING
              s1_ubar_obc(j,ib,ir,Lout)=s1_ubar_obc(j,ib,ir,Lout)*      &
     &                                  umask(Iend+1,j)
# endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              s1_ubar_obc(i,ib,ir,Lout)=fac1*s1_ubar_obc(i,ib,ir,Lin1)+ &
     &                                  fac2*s2_ubar_obc(i,ib,ir,Lin2)
# ifdef MASKING
              s1_ubar_obc(i,ib,ir,Lout)=s1_ubar_obc(i,ib,ir,Lout)*      &
     &                                  umask(i,Jstr-1)
# endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              s1_ubar_obc(i,ib,ir,Lout)=fac1*s1_ubar_obc(i,ib,ir,Lin1)+ &
     &                                  fac2*s2_ubar_obc(i,ib,ir,Lin2)
# ifdef MASKING
              s1_ubar_obc(i,ib,ir,Lout)=s1_ubar_obc(i,ib,ir,Lout)*      &
     &                                  umask(i,Jend+1)
# endif
            END DO
          END IF
        END DO
      END IF
#endif

#if !defined SOLVE3D || (defined WEAK_CONSTRAINT && defined TIME_CONV )
!
!  2D V-momentum component.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          s1_vbar(i,j,Lout)=fac1*s1_vbar(i,j,Lin1)+                     &
     &                      fac2*s2_vbar(i,j,Lin2)
# ifdef MASKING
          s1_vbar(i,j,Lout)=s1_vbar(i,j,Lout)*vmask(i,j)
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
              s1_vbar_obc(j,ib,ir,Lout)=fac1*s1_vbar_obc(j,ib,ir,Lin1)+ &
     &                                  fac2*s2_vbar_obc(j,ib,ir,Lin2)
# ifdef MASKING
              s1_vbar_obc(j,ib,ir,Lout)=s1_vbar_obc(j,ib,ir,Lout)*      &
     &                                  vmask(Istr-1,j)
# endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              s1_vbar_obc(j,ib,ir,Lout)=fac1*s1_vbar_obc(j,ib,ir,Lin1)+ &
     &                                  fac2*s2_vbar_obc(j,ib,ir,Lin2)
# ifdef MASKING
              s1_vbar_obc(j,ib,ir,Lout)=s1_vbar_obc(j,ib,ir,Lout)*      &
     &                                  vmask(Iend+1,j)
# endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              s1_vbar_obc(i,ib,ir,Lout)=fac1*s1_vbar_obc(i,ib,ir,Lin1)+ &
     &                                  fac2*s2_vbar_obc(i,ib,ir,Lin2)
# ifdef MASKING
              s1_vbar_obc(i,ib,ir,Lout)=s1_vbar_obc(i,ib,ir,Lout)*      &
     &                                  vmask(i,Jstr)
# endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              s1_vbar_obc(i,ib,ir,Lout)=fac1*s1_vbar_obc(i,ib,ir,Lin1)+ &
     &                                  fac2*s2_vbar_obc(i,ib,ir,Lin2)
# ifdef MASKING
              s1_vbar_obc(i,ib,ir,Lout)=s1_vbar_obc(i,ib,ir,Lout)*      &
     &                                  vmask(i,Jend+1)
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
            s1_sustr(i,j,ir,Lout)=fac1*s1_sustr(i,j,ir,Lin1)+           &
     &                            fac2*s2_sustr(i,j,ir,Lin2)
# ifdef MASKING
            s1_sustr(i,j,ir,Lout)=s1_sustr(i,j,ir,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            s1_svstr(i,j,ir,Lout)=fac1*s1_svstr(i,j,ir,Lin1)+           &
     &                            fac2*s2_svstr(i,j,ir,Lin2)
# ifdef MASKING
            s1_svstr(i,j,ir,Lout)=s1_svstr(i,j,ir,Lout)*vmask(i,j)
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
            s1_u(i,j,k,Lout)=fac1*s1_u(i,j,k,Lin1)+                     &
     &                       fac2*s2_u(i,j,k,Lin2)
# ifdef MASKING
            s1_u(i,j,k,Lout)=s1_u(i,j,k,Lout)*umask(i,j)
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
                s1_u_obc(j,k,ib,ir,Lout)=fac1*s1_u_obc(j,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_u_obc(j,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_u_obc(j,k,ib,ir,Lout)=s1_u_obc(j,k,ib,ir,Lout)*      &
     &                                   umask(Istr,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                s1_u_obc(j,k,ib,ir,Lout)=fac1*s1_u_obc(j,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_u_obc(j,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_u_obc(j,k,ib,ir,Lout)=s1_u_obc(j,k,ib,ir,Lout)*      &
     &                                   umask(Iend+1,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s1_u_obc(i,k,ib,ir,Lout)=fac1*s1_u_obc(i,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_u_obc(i,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_u_obc(i,k,ib,ir,Lout)=s1_u_obc(i,k,ib,ir,Lout)*      &
     &                                   umask(i,Jstr-1)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                s1_u_obc(i,k,ib,ir,Lout)=fac1*s1_u_obc(i,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_u_obc(i,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_u_obc(i,k,ib,ir,Lout)=s1_u_obc(i,k,ib,ir,Lout)*      &
     &                                   umask(i,Jend+1)
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
            s1_v(i,j,k,Lout)=fac1*s1_v(i,j,k,Lin1)+                     &
     &                       fac2*s2_v(i,j,k,Lin2)
# ifdef MASKING
            s1_v(i,j,k,Lout)=s1_v(i,j,k,Lout)*vmask(i,j)
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
                s1_v_obc(j,k,ib,ir,Lout)=fac1*s1_v_obc(j,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_v_obc(j,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_v_obc(j,k,ib,ir,Lout)=s1_v_obc(j,k,ib,ir,Lout)*      &
     &                                   vmask(Istr-1,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                s1_v_obc(j,k,ib,ir,Lout)=fac1*s1_v_obc(j,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_v_obc(j,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_v_obc(j,k,ib,ir,Lout)=s1_v_obc(j,k,ib,ir,Lout)*      &
     &                                   vmask(Iend+1,j)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                s1_v_obc(i,k,ib,ir,Lout)=fac1*s1_v_obc(i,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_v_obc(i,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_v_obc(i,k,ib,ir,Lout)=s1_v_obc(i,k,ib,ir,Lout)*      &
     &                                   vmask(i,Jstr)
#  endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                s1_v_obc(i,k,ib,ir,Lout)=fac1*s1_v_obc(i,k,ib,ir,Lin1)+ &
     &                                   fac2*s2_v_obc(i,k,ib,ir,Lin2)
#  ifdef MASKING
                s1_v_obc(i,k,ib,ir,Lout)=s1_v_obc(i,k,ib,ir,Lout)*      &
     &                                   vmask(i,Jend+1)
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
              s1_t(i,j,k,Lout,it)=fac1*s1_t(i,j,k,Lin1,it)+             &
     &                            fac2*s2_t(i,j,k,Lin2,it)
# ifdef MASKING
              s1_t(i,j,k,Lout,it)=s1_t(i,j,k,Lout,it)*rmask(i,j)
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
                  s1_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                               fac1*s1_t_obc(j,k,ib,ir,Lin1,it)+  &
     &                               fac2*s2_t_obc(j,k,ib,ir,Lin2,it)
#  ifdef MASKING
                  s1_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                    s1_t_obc(j,k,ib,ir,Lout,it)*rmask(Istr-1,j)
#  endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  s1_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                               fac1*s1_t_obc(j,k,ib,ir,Lin1,it)+  &
     &                               fac2*s2_t_obc(j,k,ib,ir,Lin2,it)
#  ifdef MASKING
                  s1_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                    s1_t_obc(j,k,ib,ir,Lout,it)*rmask(Iend+1,j)
#  endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s1_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                               fac1*s1_t_obc(i,k,ib,ir,Lin1,it)+  &
     &                               fac2*s2_t_obc(i,k,ib,ir,Lin2,it)
#  ifdef MASKING
                  s1_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                    s1_t_obc(i,k,ib,ir,Lout,it)*rmask(i,Jstr-1)
#  endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  s1_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                            fac1*s1_t_obc(i,k,ib,ir,Lin1,it)+     &
     &                            fac2*s2_t_obc(i,k,ib,ir,Lin2,it)
#  ifdef MASKING
                  s1_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                    s1_t_obc(i,k,ib,ir,Lout,it)*rmask(i,Jend+1)
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
                s1_tflux(i,j,ir,Lout,it)=fac1*s1_tflux(i,j,ir,Lin1,it)+ &
     &                                   fac2*s2_tflux(i,j,ir,Lin2,it)
#  ifdef MASKING
                s1_tflux(i,j,ir,Lout,it)=s1_tflux(i,j,ir,Lout,it)*      &
     &                                   rmask(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
      END DO
# endif

#endif

      RETURN
      END SUBROUTINE state_addition

      END MODULE state_addition_mod
