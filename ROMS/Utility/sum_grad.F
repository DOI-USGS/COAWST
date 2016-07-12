#include "cppdefs.h"
      MODULE sum_grad_mod
#if defined IS4DVAR || defined RPCG
!
!svn $Id: sum_grad.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the sum of the background cost function       !
!  gradients in v-space                                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC sum_grad

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sum_grad (ng, tile, Linp, Lout)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL sum_grad_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Linp, Lout,                                   &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    BOUNDARY(ng) % tl_t_obc,                      &
     &                    BOUNDARY(ng) % tl_u_obc,                      &
     &                    BOUNDARY(ng) % tl_v_obc,                      &
#  endif
     &                    BOUNDARY(ng) % tl_ubar_obc,                   &
     &                    BOUNDARY(ng) % tl_vbar_obc,                   &
     &                    BOUNDARY(ng) % tl_zeta_obc,                   &
# endif
# ifdef ADJUST_WSTRESS
     &                    FORCES(ng) % tl_ustr,                         &
     &                    FORCES(ng) % tl_vstr,                         &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    FORCES(ng) % tl_tflux,                        &
#  endif
     &                    OCEAN(ng) % tl_t,                             &
     &                    OCEAN(ng) % tl_u,                             &
     &                    OCEAN(ng) % tl_v,                             &
# else
     &                    OCEAN(ng) % tl_ubar,                          &
     &                    OCEAN(ng) % tl_vbar,                          &
# endif
     &                    OCEAN(ng) % tl_zeta)
      RETURN
      END SUBROUTINE sum_grad
!
!***********************************************************************
      SUBROUTINE sum_grad_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Linp, Lout,                             &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          tl_t_obc, tl_u_obc, tl_v_obc,           &
#  endif
     &                          tl_ubar_obc, tl_vbar_obc,               &
     &                          tl_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          tl_ustr, tl_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          tl_tflux,                               &
#  endif
     &                          tl_t, tl_u, tl_v,                       &
# else
     &                          tl_ubar, tl_vbar,                       &
# endif
     &                          tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
      USE mod_scalars
# endif
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ib, ir, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Sum of the background cost function gradients (v-space).
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Linp)+                          &
     &                      tl_zeta(i,j,Lout)
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Linp)+      &
     &                                  tl_zeta_obc(j,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Linp)+      &
     &                                  tl_zeta_obc(j,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Linp)+      &
     &                                  tl_zeta_obc(i,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Linp)+      &
     &                                  tl_zeta_obc(i,ib,ir,Lout)
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum component.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Linp)+                          &
     &                      tl_ubar(i,j,Lout)
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Linp)+      &
     &                                  tl_ubar_obc(j,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Linp)+      &
     &                                  tl_ubar_obc(j,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Linp)+      &
     &                                  tl_ubar_obc(i,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Linp)+      &
     &                                  tl_ubar_obc(i,ib,ir,Lout)
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Linp)+                          &
     &                      tl_vbar(i,j,Lout)
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=JstrV,Jend
              tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Linp)+      &
     &                                  tl_vbar_obc(j,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Linp)+      &
     &                                  tl_vbar_obc(j,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Linp)+      &
     &                                  tl_vbar_obc(i,ib,ir,Lout)
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Linp)+      &
     &                                  tl_vbar_obc(i,ib,ir,Lout)
            END DO
          END IF
        END DO
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Linp)+                    &
     &                          tl_ustr(i,j,k,Lout)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Linp)+                    &
     &                          tl_vstr(i,j,k,Lout)
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Linp)+                          &
     &                       tl_u(i,j,k,Lout)
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
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
                tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Linp)+      &
     &                                   tl_u_obc(j,k,ib,ir,Lout)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Linp)+      &
     &                                   tl_u_obc(j,k,ib,ir,Lout)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Linp)+      &
     &                                   tl_u_obc(i,k,ib,ir,Lout)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Linp)+      &
     &                                   tl_u_obc(i,k,ib,ir,Lout)
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Linp)+                          &
     &                       tl_v(i,j,k,Lout)
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
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
                tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Linp)+      &
     &                                   tl_v_obc(j,k,ib,ir,Lout)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Linp)+      &
     &                                   tl_v_obc(j,k,ib,ir,Lout)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Linp)+      &
     &                                   tl_v_obc(i,k,ib,ir,Lout)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Linp)+      &
     &                                   tl_v_obc(i,k,ib,ir,Lout)
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Linp,itrc)+              &
     &                              tl_t(i,j,k,Lout,itrc)
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isTvar(itrc),ng)).and.                      &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_t_obc(j,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(j,k,ib,ir,Linp,itrc)+     &
     &                               tl_t_obc(j,k,ib,ir,Lout,itrc)
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(itrc),ng)).and.                      &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_t_obc(j,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(j,k,ib,ir,Linp,itrc)+     &
     &                               tl_t_obc(j,k,ib,ir,Lout,itrc)
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(itrc),ng)).and.                     &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_t_obc(i,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(i,k,ib,ir,Linp,itrc)+     &
     &                               tl_t_obc(i,k,ib,ir,Lout,itrc)
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(itrc),ng)).and.                     &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_t_obc(i,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(i,k,ib,ir,Linp,itrc)+     &
     &                               tl_t_obc(i,k,ib,ir,Lout,itrc)
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          DO k=1,Nfrec(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Linp,itrc)+    &
     &                                    tl_tflux(i,j,k,Lout,itrc)
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE sum_grad_tile
#endif
      END MODULE sum_grad_mod
