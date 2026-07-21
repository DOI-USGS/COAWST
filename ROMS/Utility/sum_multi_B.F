#include "cppdefs.h"
      MODULE sum_multi_B_mod
#if defined FOUR_DVAR && defined MULTI_SCALE_B
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group           Andrew M. Moore    !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Multi-scale, Background-Error Covariace modeling:                   !
!                                                                      !
!  This routine computes the sum of the multiscale background cost     !
!  function gradients in v-space.                                      !
!                                                                      !
!  In the multiscale approach, expressed as B = SUM(Wi*Bi), represents !
!  a weighted sum of B values corresponding to various scales.         !
!  Typically, two to few different scales are combined in practical    !
!  applications. The weight coefficients Wi are constrained to sum     !
!  to unity.                                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
      USE roms_multiscale_mod, ONLY : multiscale         ! CLASS object
!
      implicit none
!
      PUBLIC  :: sum_multi_B
      PRIVATE :: sum_multi_B_tile
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sum_multi_B (ng, tile, ns, Linp, Lout)
!***********************************************************************
!
      USE roms_multiscale_mod, ONLY : MSB
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, ns, Linp, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL sum_multi_B_tile (MSB(ng), ng, tile,                         &
     &                       ns, Linp, Lout,                            &
     &                       LBi, UBi, LBj, UBj,                        &
# ifdef ADJUST_BOUNDARY
     &                       LBij, UBij,                                &
# endif
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef ADJUST_BOUNDARY
     &                       BOUNDARY(ng) % tl_ubar_obc,                &
     &                       BOUNDARY(ng) % tl_vbar_obc,                &
     &                       BOUNDARY(ng) % tl_zeta_obc,                &
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % tl_t_obc,                   &
     &                       BOUNDARY(ng) % tl_u_obc,                   &
     &                       BOUNDARY(ng) % tl_v_obc,                   &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % tl_ustr,                      &
     &                       FORCES(ng) % tl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       FORCES(ng) % tl_tflux,                     &
#  endif
     &                       OCEAN(ng) % tl_t,                          &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
# else
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
# endif
     &                       OCEAN(ng) % tl_zeta,                       &
# ifdef ADJUST_BOUNDARY
     &                       BOUNDARY(ng) % ad_ubar_obc,                &
     &                       BOUNDARY(ng) % ad_vbar_obc,                &
     &                       BOUNDARY(ng) % ad_zeta_obc,                &
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % ad_t_obc,                   &
     &                       BOUNDARY(ng) % ad_u_obc,                   &
     &                       BOUNDARY(ng) % ad_v_obc,                   &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % ad_ustr,                      &
     &                       FORCES(ng) % ad_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       FORCES(ng) % ad_tflux,                     &
#  endif
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
# else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
# endif
     &                       OCEAN(ng) % ad_zeta)
!
      RETURN
      END SUBROUTINE sum_multi_B
!
!***********************************************************************
      SUBROUTINE sum_multi_B_tile (self, ng, tile,                      &
     &                             ns, Linp, Lout,                      &
     &                             LBi, UBi, LBj, UBj,                  &
# ifdef ADJUST_BOUNDARY
     &                             LBij, UBij,                          &
# endif
     &                             IminS, ImaxS, JminS, JmaxS,          &
# ifdef ADJUST_BOUNDARY
     &                             tl_ubar_obc, tl_vbar_obc,            &
     &                             tl_zeta_obc,                         &
#  ifdef SOLVE3D
     &                             tl_t_obc, tl_u_obc, tl_v_obc,        &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                             tl_ustr, tl_vstr,                    &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                             tl_tflux,                            &
#  endif
     &                             tl_t, tl_u, tl_v,                    &
# else
     &                             tl_ubar, tl_vbar,                    &
# endif
     &                             tl_zeta,                             &
# ifdef ADJUST_BOUNDARY
     &                             ad_ubar_obc, ad_vbar_obc,            &
     &                             ad_zeta_obc,                         &
#  ifdef SOLVE3D
     &                             ad_t_obc, ad_u_obc, ad_v_obc,        &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                             ad_ustr, ad_vstr,                    &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                             ad_tflux,                            &
#  endif
     &                             ad_t, ad_u, ad_v,                    &
# else
     &                             ad_ubar, ad_vbar,                    &
# endif
     &                             ad_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: ns, Linp, Lout
      integer, intent(in) :: LBi, UBi, LBj, UBj
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,:)
# endif
!
!  Local variable declarations.
!
      logical :: Lboundary(4)
!
      integer :: i, ib, ir, j, k
# ifdef SOLVE3D
      integer :: ifield, itrc
# endif

# include "set_bounds.h"
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Combine weighted multiscle spatial scales.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ifdef ADJUST_BOUNDARY
!
!  Set switch to convolve boundary segments by the appropriate
!  tiles.
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
# endif
!
!  Free-surface.
!
      IF (ns.eq.1) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Lout)=0.0_r8
          END DO
        END DO
      END IF
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lout)=ad_zeta(i,j,Lout)+                          &
     &                      self%Bwgt(isFsur,ns)*tl_zeta(i,j,Linp)
        END DO
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isFsur,ng)) THEN
            IF (Lboundary(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  IF (ns.eq.1) THEN
                    DO j=JstrR,JendR
                      ad_zeta_obc(j,ib,ir,Lout)=0.0_r8
                    END DO
                  END IF
!
                  DO j=JstrR,JendR
                    ad_zeta_obc(j,ib,ir,Lout)=                          &
     &                               ad_zeta_obc(j,ib,ir,Lout)+         &
     &                               self%Bwgt(isFsur,ns)*              &
     &                               tl_zeta_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  IF (ns.eq.1) THEN
                    DO i=IstrR,IendR
                      ad_zeta_obc(i,ib,ir,Lout)=0.0_r8
                    END DO
                  END IF
!
                  DO i=IstrR,IendR
                    ad_zeta_obc(i,ib,ir,Lout)=                          &
                                     ad_zeta_obc(i,ib,ir,Lout)+         &
     &                               self%Bwgt(isFsur,ns)*              &
     &                               tl_zeta_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
# endif

# ifndef SOLVE3D
!
!  2D U-momentum component.
!
      IF (ns.eq.1) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,Lout)=0.0_r8
          END DO
        END DO
      END IF
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lout)=ad_ubar(i,j,Lout)+                          &
     &                      self%Bwgt(isUbar,ns)*tl_ubar(i,j,Linp)
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isUbar,ng)) THEN
            IF (Lboundary(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  IF (ns.eq.1) THEN
                    DO j=JstrR,JendR
                      ad_ubar_obc(j,ib,ir,Lout)=0.0_r8
                    END DO
                  END IF
!
                  DO j=JstrR,JendR
                    ad_ubar_obc(j,ib,ir,Lout)=                          &
     &                               ad_ubar_obc(j,ib,ir,Lout)+         &
     &                               self%Bwgt(isUbar,ns)*              &
     &                               tl_ubar_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  IF (ns.eq.1) THEN
                    DO i=Istr,IendR
                      ad_ubar_obc(i,ib,ir,Lout)=0.0_r8
                    END DO
                  END IF
!
                  DO i=Istr,IendR
                    ad_ubar_obc(i,ib,ir,Lout)=                          &
     &                               ad_ubar_obc(i,ib,ir,Lout)+         &
     &                               self%Bwgt(isUbar,ns)*              &
     &                               tl_ubar_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
# endif

# ifndef SOLVE3D
!
!  2D V-momentum component.
!
      IF (ns.eq.1) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Lout)=0.0_r8
          END DO
        END DO
      END IF
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lout)=ad_vbar(i,j,Lout)+                          &
     &                      self%Bwgt(isVbar,ns)*                       &
     &                      tl_vbar(i,j,Linp)
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!  2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isVbar,ng)) THEN
            IF (Lboundary(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  IF (ns.eq.1) THEN
                    DO j=Jstr,JendR
                      ad_vbar_obc(j,ib,ir,Lout)=0.0_r8
                    END DO
                  END IF
!
                  DO j=Jstr,JendR
                    ad_vbar_obc(j,ib,ir,Lout)=                          &
     &                               ad_vbar_obc(j,ib,ir,Lout)+         &
     &                               self%Bwgt(isVbar,ns)*              &
     &                               tl_vbar_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  IF (ns.eq.1) THEN
                    DO i=IstrR,IendR
                      ad_vbar_obc(i,ib,ir,Lout)=0.0_r8
                    END DO
                  END IF
!
                  DO i=IstrR,IendR
                    ad_vbar_obc(i,ib,ir,Lout)=                          &
     &                               ad_vbar_obc(i,ib,ir,Lout)+         &
     &                               self%Bwgt(isVbar,ns)*              &
     &                               tl_vbar_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        IF (ns.eq.1) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_ustr(i,j,k,Lout)=0.0_r8
            END DO
          END DO
!
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_vstr(i,j,k,Lout)=0.0_r8
            END DO
          END DO
        END IF
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ustr(i,j,k,Lout)=ad_ustr(i,j,k,Lout)+                    &
     &                          self%Bwgt(isUstr,ns)*                   &
     &                          tl_ustr(i,j,k,Linp)
          END DO
        END DO
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,k,Lout)=ad_vstr(i,j,k,Lout)+                    &
     &                          self%Bwgt(isVstr,ns)*                   &
     &                          tl_vstr(i,j,k,Linp)
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum component.
!
      IF (ns.eq.1) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_u(i,j,k,Lout)=0.0_r8
            END DO
          END DO
        END DO
      END IF
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lout)=ad_u(i,j,k,Lout)+                          &
     &                       self%Bwgt(isUvel,ns)*                      &
     &                       tl_u(i,j,k,Linp)
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isUvel,ng)) THEN
            IF (Lboundary(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  IF (ns.eq.1) THEN
                    DO k=1,N(ng)
                      DO j=JstrR,JendR
                        ad_u_obc(j,k,ib,ir,Lout)=0.0_r8
                      END DO
                    END DO
                  END IF
!
                  DO k=1,N(ng)
                    DO j=JstrR,JendR
                      ad_u_obc(j,k,ib,ir,Lout)=                         &
     &                                ad_u_obc(j,k,ib,ir,Lout)+         &
     &                                self%Bwgt(isUvel,ns)*             &
     &                                tl_u_obc(j,k,ib,ir,Linp)
                    END DO
                  END DO
                CASE (isouth, inorth)
                  IF (ns.eq.1) THEN
                    DO k=1,N(ng)
                      DO i=Istr,IendR
                        ad_u_obc(i,k,ib,ir,Lout)=0.0_r8
                      END DO
                    END DO
                  END IF
!
                  DO k=1,N(ng)
                    DO i=Istr,IendR
                      ad_u_obc(i,k,ib,ir,Lout)=                         &
     &                                ad_u_obc(i,k,ib,ir,Lout)+         &
     &                                self%Bwgt(isUvel,ns)*             &
     &                                tl_u_obc(i,k,ib,ir,Linp)
                    END DO
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
#  endif
!
!  3D V-momentum component.
!
      IF (ns.eq.1) THEN
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,Lout)=0.0_r8
            END DO
          END DO
        END DO
      END IF
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lout)=ad_v(i,j,k,Lout)+                          &
     &                       self%Bwgt(isVvel,ns)*                      &
     &                       tl_v(i,j,k,Linp)
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (Lobc(ib,isVvel,ng)) THEN
            IF (Lboundary(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  IF (ns.eq.1) THEN
                    DO k=1,N(ng)
                      DO j=Jstr,JendR
                        ad_v_obc(j,k,ib,ir,Lout)=0.0_r8
                      END DO
                    END DO
                  END IF
!
                  DO k=1,N(ng)
                    DO j=Jstr,JendR
                      ad_v_obc(j,k,ib,ir,Lout)=                         &
     &                                ad_v_obc(j,k,ib,ir,Lout)+         &
     &                                self%Bwgt(isVvel,ns)*             &
     &                                tl_v_obc(j,k,ib,ir,Linp)
                    END DO
                  END DO
                CASE (isouth, inorth)
                  IF (ns.eq.1) THEN
                    DO k=1,N(ng)
                      DO i=IstrR,IendR
                        ad_v_obc(i,k,ib,ir,Lout)=0.0_r8
                      END DO
                    END DO
                  END IF
!
                  DO k=1,N(ng)
                    DO i=IstrR,IendR
                      ad_v_obc(i,k,ib,ir,Lout)=                         &
     &                                ad_v_obc(i,k,ib,ir,Lout)+         &
     &                                self%Bwgt(isVvel,ns)*             &
     &                                tl_v_obc(i,k,ib,ir,Linp)
                    END DO
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
#  endif
!
!  Tracer variables.
!
      DO itrc=1,NT(ng)
        ifield=isTvar(itrc)
        IF (ns.eq.1) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_t(i,j,k,Lout,itrc)=0.0_r8
              END DO
            END DO
          END DO
        END IF
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lout,itrc)=ad_t(i,j,k,Lout,itrc)+              &
     &                              self%Bwgt(ifield,ns)*               &
     &                              tl_t(i,j,k,Linp,itrc)
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Tracers open boundaries.
!
      DO itrc=1,NT(ng)
        ifield=isTvar(itrc)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (ANY(Lobc(:,ifield,ng))) THEN
              IF (Lboundary(ib)) THEN
                SELECT CASE (ib)
                  CASE (iwest, ieast)
                    IF (ns.eq.1) THEN
                      DO k=1,N(ng)
                        DO j=JstrR,JendR
                          ad_t_obc(j,k,ib,ir,Lout,itrc)=0.0_r8
                        END DO
                      END DO
                    END IF
!
                    DO k=1,N(ng)
                      DO j=JstrR,JendR
                        ad_t_obc(j,k,ib,ir,Lout,itrc)=                  &
     &                                  ad_t_obc(j,k,ib,ir,Lout,itrc)+  &
     &                                  self%Bwgt(ifield,ns)*           &
     &                                  tl_t_obc(j,k,ib,ir,Linp,itrc)
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    IF (ns.eq.1) THEN
                      DO k=1,N(ng)
                        DO i=IstrR,IendR
                          ad_t_obc(i,k,ib,ir,Lout,itrc)=0.0_r8
                        END DO
                      END DO
                    END IF
!
                    DO k=1,N(ng)
                      DO i=IstrR,IendR
                        ad_t_obc(i,k,ib,ir,Lout,itrc)=                  &
     &                                  ad_t_obc(i,k,ib,ir,Lout,itrc)+  &
     &                                  self%Bwgt(ifield,ns)*           &
     &                                  tl_t_obc(i,k,ib,ir,Linp,itrc)
                      END DO
                    END DO
                END SELECT
              END IF
            END IF
          END DO
        END DO
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        ifield=isTsur(itrc)
        IF (Lstflux(itrc,ng)) THEN
          DO k=1,Nfrec(ng)
            IF (ns.eq.1) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ad_tflux(i,j,k,Lout,itrc)=0.0_r8
                END DO
              END DO
            END IF
!
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_tflux(i,j,k,Lout,itrc)=ad_tflux(i,j,k,Lout,itrc)+    &
     &                                    self%Bwgt(ifield,ns)*         &
     &                                    tl_tflux(i,j,k,Linp,itrc)
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif
!
      RETURN
      END SUBROUTINE sum_multi_B_tile
#endif
      END MODULE sum_multi_B_mod
