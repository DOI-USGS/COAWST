#include "cppdefs.h"
      MODULE state_product_mod
!
!svn $Id: state_product.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the product of each element of two model      !
!  states:                                                             !
!                                                                      !
!    s3(...) = s1(...) * s2(...)                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC  :: state_product

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_product (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                          s1_t_obc, s2_t_obc, s3_t_obc,           &
     &                          s1_u_obc, s2_u_obc, s3_u_obc,           &
     &                          s1_v_obc, s2_v_obc, s3_v_obc,           &
# endif
     &                          s1_ubar_obc, s2_ubar_obc, s3_ubar_obc,  &
     &                          s1_vbar_obc, s2_vbar_obc, s3_vbar_obc,  &
     &                          s1_zeta_obc, s2_zeta_obc, s3_zeta_obc,  &
#endif
#ifdef ADJUST_WSTRESS
     &                          s1_sustr, s2_sustr, s3_sustr,           &
     &                          s1_svstr, s2_svstr, s3_svstr,           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          s1_tflux, s2_tflux, s3_tflux,           &
# endif
     &                          s1_t, s2_t, s3_t,                       &
     &                          s1_u, s2_u, s3_u,                       &
     &                          s1_v, s2_v, s3_v,                       &
#else
     &                          s1_ubar, s2_ubar, s3_ubar,              &
     &                          s1_vbar, s2_vbar, s3_vbar,              &
#endif
     &                          s1_zeta, s2_zeta, s3_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
      USE mod_scalars
#endif
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(in) :: s1_t_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: s2_t_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: s1_u_obc(LBij:,:,:,:)
      real(r8), intent(in) :: s2_u_obc(LBij:,:,:,:)
      real(r8), intent(in) :: s1_v_obc(LBij:,:,:,:)
      real(r8), intent(in) :: s2_v_obc(LBij:,:,:,:)

      real(r8), intent(inout) :: s3_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s3_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s3_v_obc(LBij:,:,:,:)
#  endif
      real(r8), intent(in) :: s1_ubar_obc(LBij:,:,:)
      real(r8), intent(in) :: s2_ubar_obc(LBij:,:,:)
      real(r8), intent(in) :: s1_vbar_obc(LBij:,:,:)
      real(r8), intent(in) :: s2_vbar_obc(LBij:,:,:)
      real(r8), intent(in) :: s1_zeta_obc(LBij:,:,:)
      real(r8), intent(in) :: s2_zeta_obc(LBij:,:,:)

      real(r8), intent(inout) :: s3_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: s3_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: s3_zeta_obc(LBij:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: s1_sustr(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_sustr(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_svstr(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_svstr(LBi:,LBj:,:)

      real(r8), intent(inout) :: s3_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: s3_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(in) :: s1_tflux(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_tflux(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: s3_tflux(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: s1_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s1_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_v(LBi:,LBj:,:)

      real(r8), intent(inout) :: s3_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s3_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: s3_v(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: s1_ubar(LBi:,LBj:)
      real(r8), intent(in) :: s2_ubar(LBi:,LBj:)
      real(r8), intent(in) :: s1_vbar(LBi:,LBj:)
      real(r8), intent(in) :: s2_vbar(LBi:,LBj:)

      real(r8), intent(inout) :: s3_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: s3_vbar(LBi:,LBj:)
# endif
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)

      real(r8), intent(inout) :: s3_zeta(LBi:,LBj:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif

# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(in) :: s1_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),NT(ng))
      real(r8), intent(in) :: s2_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),NT(ng))
      real(r8), intent(in) :: s1_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(in) :: s2_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(in) :: s1_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(in) :: s2_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))

      real(r8), intent(inout) :: s3_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),NT(ng))
      real(r8), intent(inout) :: s3_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(inout) :: s3_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#  endif
      real(r8), intent(in) :: s1_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: s2_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: s1_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: s2_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: s1_zeta_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: s2_zeta_obc(LBij:UBij,4,Nbrec(ng))

      real(r8), intent(inout) :: s3_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: s3_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: s3_zeta_obc(LBij:UBij,4,Nbrec(ng))
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: s1_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(in) :: s2_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(in) :: s1_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(in) :: s2_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))

      real(r8), intent(inout) :: s3_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: s3_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(in) :: s1_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),NT(ng))
      real(r8), intent(in) :: s2_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),NT(ng))

      real(r8), intent(inout) :: s3_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),NT(ng))
#  endif
      real(r8), intent(in) :: s1_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: s2_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: s1_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: s2_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: s1_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: s2_v(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: s3_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: s3_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: s3_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(in) :: s1_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s1_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_vbar(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: s3_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: s3_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_zeta(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: s3_zeta(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      integer :: ir, it

      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute product between S1 and S2 model state trajectories and save
!  in S3.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff=s1_zeta(i,j)*s2_zeta(i,j)
#ifdef MASKING
          cff=cff*rmask(i,j)
#endif
          s3_zeta(i,j)=cff
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
            DO j=Jstr,Jend
              cff=s1_zeta_obc(j,iwest,ir)*s2_zeta_obc(j,iwest,ir)
# ifdef MASKING
              cff=cff*rmask(Istr-1,j)
# endif
              s3_zeta_obc(j,iwest,ir)=cff
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              cff=s1_zeta_obc(j,ieast,ir)*                              &
     &            s2_zeta_obc(j,ieast,ir)
# ifdef MASKING
              cff=cff*rmask(Iend+1,j)
# endif
              s3_zeta_obc(j,ieast,ir)=cff
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=s1_zeta_obc(i,isouth,ir)*                             &
     &            s2_zeta_obc(i,isouth,ir)
# ifdef MASKING
              cff=cff*rmask(i,Jstr-1)
# endif
              s3_zeta_obc(i,isouth,ir)=cff
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=s1_zeta_obc(i,inorth,ir)*                             &
     &            s2_zeta_obc(i,inorth,ir)
# ifdef MASKING
              cff=cff*rmask(i,Jend+1)
# endif
              s3_zeta_obc(i,inorth,ir)=cff
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
          cff=s1_ubar(i,j)*s2_ubar(i,j)
# ifdef MASKING
          cff=cff*umask(i,j)
# endif
          s3_ubar(i,j)=cff
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
            DO j=Jstr,Jend
              cff=s1_ubar_obc(j,iwest,ir)*                              &
     &            s2_ubar_obc(j,iwest,ir)
# ifdef MASKING
              cff=cff*umask(Istr,j)
# endif
              s3_ubar_obc(j,iwest,ir)=cff
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              cff=s1_ubar_obc(j,ieast,ir)*                              &
     &            s2_ubar_obc(j,ieast,ir)
# ifdef MASKING
              cff=cff*umask(Iend+1,j)
# endif
              s3_ubar_obc(j,ieast,ir)=cff
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
              cff=s1_ubar_obc(i,isouth,ir)*                             &
     &            s2_ubar_obc(i,isouth,ir)
# ifdef MASKING
              cff=cff*umask(i,Jstr-1)
# endif
              s3_ubar_obc(i,isouth,ir)=cff
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
              cff=s1_ubar_obc(i,inorth,ir)*                             &
     &            s2_ubar_obc(i,inorth,ir)
# ifdef MASKING
              cff=cff*umask(i,Jend+1)
# endif
              s3_ubar_obc(i,inorth,ir)=cff
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
          cff=s1_vbar(i,j)*s2_vbar(i,j)
# ifdef MASKING
          cff=cff*vmask(i,j)
# endif
          s3_vbar(i,j)=cff
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
            DO j=JstrV,Jend
              cff=s1_vbar_obc(j,iwest,ir)*                              &
     &            s2_vbar_obc(j,iwest,ir)
# ifdef MASKING
              cff=cff*vmask(Istr-1,j)
# endif
              s3_vbar_obc(j,iwest,ir)=cff
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
              cff=s1_vbar_obc(j,ieast,ir)*                              &
     &            s2_vbar_obc(j,ieast,ir)
# ifdef MASKING
              cff=cff*vmask(Iend+1,j)
# endif
              s3_vbar_obc(j,ieast,ir)=cff
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=s1_vbar_obc(i,isouth,ir)*                             &
     &            s2_vbar_obc(i,isouth,ir)
# ifdef MASKING
              cff=cff*vmask(i,Jstr)
# endif
              s3_vbar_obc(i,isouth,ir)=cff
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=s1_vbar_obc(i,inorth,ir)*                             &
     &            s2_vbar_obc(i,inorth,ir)
# ifdef MASKING
              cff=cff*vmask(i,Jend+1)
# endif
              s3_vbar_obc(i,inorth,ir)=cff
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
            cff=s1_sustr(i,j,ir)*s2_sustr(i,j,ir)
# ifdef MASKING
            cff=cff*umask(i,j)
# endif
            s3_sustr(i,j,ir)=cff
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            cff=s1_svstr(i,j,ir)*s2_svstr(i,j,ir)
# ifdef MASKING
            cff=cff*vmask(i,j)
# endif
            s3_svstr(i,j,ir)=cff
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
            cff=s1_u(i,j,k)*s2_u(i,j,k)
# ifdef MASKING
            cff=cff*umask(i,j)
# endif
            s3_u(i,j,k)=cff
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
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=s1_u_obc(j,k,iwest,ir)*                             &
     &              s2_u_obc(j,k,iwest,ir)
#  ifdef MASKING
                cff=cff*umask(Istr,j)
#  endif
                s3_u_obc(j,k,iwest,ir)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=s1_u_obc(j,k,ieast,ir)*                             &
     &              s2_u_obc(j,k,ieast,ir)
#  ifdef MASKING
                cff=cff*umask(Iend+1,j)
#  endif
                s3_u_obc(j,k,ieast,ir)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=s1_u_obc(i,k,isouth,ir)*                            &
     &              s2_u_obc(i,k,isouth,ir)
#  ifdef MASKING
                cff=cff*umask(i,Jstr-1)
#  endif
                s3_u_obc(i,k,isouth,ir)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=s1_u_obc(i,k,inorth,ir)*                            &
     &              s2_u_obc(i,k,inorth,ir)
#  ifdef MASKING
                cff=cff*umask(i,Jend+1)
#  endif
                s3_u_obc(i,k,inorth,ir)=cff
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
            cff=s1_v(i,j,k)*s2_v(i,j,k)
# ifdef MASKING
            cff=cff*vmask(i,j)
# endif
            s3_v(i,j,k)=cff
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
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=s1_v_obc(j,k,iwest,ir)*                             &
     &              s2_v_obc(j,k,iwest,ir)
#  ifdef MASKING
                cff=cff*vmask(Istr-1,j)
#  endif
                s3_v_obc(j,k,iwest,ir)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=s1_v_obc(j,k,ieast,ir)*                             &
     &              s2_v_obc(j,k,ieast,ir)
#  ifdef MASKING
                cff=cff*vmask(Iend+1,j)
#  endif
                s3_v_obc(j,k,ieast,ir)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=s1_v_obc(i,k,isouth,ir)*                            &
     &              s2_v_obc(i,k,isouth,ir)
#  ifdef MASKING
                cff=cff*vmask(i,Jstr)
#  endif
                s3_v_obc(i,k,isouth,ir)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=s1_v_obc(i,k,inorth,ir)*                            &
     &              s2_v_obc(i,k,inorth,ir)
#  ifdef MASKING
                cff=cff*vmask(i,Jend+1)
#  endif
                s3_v_obc(i,k,inorth,ir)=cff
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
              cff=s1_t(i,j,k,it)*s2_t(i,j,k,it)
# ifdef MASKING
              cff=cff*rmask(i,j)
# endif
              s3_t(i,j,k,it)=cff
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
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  cff=s1_t_obc(j,k,iwest,ir,it)*                        &
     &                s2_t_obc(j,k,iwest,ir,it)
#  ifdef MASKING
                  cff=cff*rmask(Istr-1,j)
#  endif
                  s3_t_obc(j,k,iwest,ir,it)=cff
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  cff=s1_t_obc(j,k,ieast,ir,it)*                        &
     &                s2_t_obc(j,k,ieast,ir,it)
#  ifdef MASKING
                  cff=cff*rmask(Iend+1,j)
#  endif
                  s3_t_obc(j,k,ieast,ir,it)=cff
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff=s1_t_obc(i,k,isouth,ir,it)*                       &
     &                s2_t_obc(i,k,isouth,ir,it)
#  ifdef MASKING
                  cff=cff*rmask(i,Jstr-1)
#  endif
                  s3_t_obc(i,k,isouth,ir,it)=cff
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff=s1_t_obc(i,k,inorth,ir,it)*                       &
     &                s2_t_obc(i,k,inorth,ir,it)
#  ifdef MASKING
                  cff=cff*rmask(i,Jend+1)
#  endif
                  s3_t_obc(i,k,inorth,ir,it)=cff
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
                cff=s1_tflux(i,j,ir,it)*s2_tflux(i,j,ir,it)
#  ifdef MASKING
                cff=cff*rmask(i,j)
#  endif
                s3_tflux(i,j,ir,it)=cff
              END DO
            END DO
          END DO
        END IF
      END DO
# endif

#endif

      RETURN
      END SUBROUTINE state_product

      END MODULE state_product_mod
