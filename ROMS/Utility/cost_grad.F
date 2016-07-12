#include "cppdefs.h"
      MODULE cost_grad_mod
#if defined IS4DVAR
!
!svn $Id: cost_grad.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the total cost function gradient by adding    !
!  background and observations contributions in v-space:               !
!                                                                      !
!      GRADv(J) = GRADv(Jb) + GRADv(Jo)                         (1)    !
!                                                                      !
!      GRADv(J) = SUM deltaV + B^(T/2) GRADx(Jo)                (2)    !
!                                                                      !
!      GRADv(J) = SUM deltaV + W^(-1/2) L^(T/2) G S             (3)    !
!                                                                      !
!  where                                                               !
!                                                                      !
!      deltaV     Increment vectors (TLM solution) in v-space          !
!      B          Background-error covariance matrix                   !
!      S          Background-error standard deviations                 !
!      G          Normalization coefficients matrix                    !
!      L          Self-adjoint filtering (diffusion) operators         !
!      W          Grid cell area or volume metric                      !
!                                                                      !
!  The last term in (3) is the result of the adjoint convolution.      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC cost_grad

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cost_grad (ng, tile, Linp1, Linp2, Lout)
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
      integer, intent(in) :: ng, tile, Linp1, Linp2, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL cost_grad_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Linp1, Linp2, Lout,                          &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % tl_t_obc,                     &
     &                     BOUNDARY(ng) % tl_u_obc,                     &
     &                     BOUNDARY(ng) % tl_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % tl_ubar_obc,                  &
     &                     BOUNDARY(ng) % tl_vbar_obc,                  &
     &                     BOUNDARY(ng) % tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
#  endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta,                         &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % ad_t_obc,                     &
     &                     BOUNDARY(ng) % ad_u_obc,                     &
     &                     BOUNDARY(ng) % ad_v_obc,                     &
#  endif
     &                     BOUNDARY(ng) % ad_ubar_obc,                  &
     &                     BOUNDARY(ng) % ad_vbar_obc,                  &
     &                     BOUNDARY(ng) % ad_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ad_ustr,                        &
     &                     FORCES(ng) % ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % ad_tflux,                       &
#  endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
# else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
# endif
     &                     OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE cost_grad
!
!***********************************************************************
      SUBROUTINE cost_grad_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Linp1, Linp2, Lout,                    &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           tl_t_obc, tl_u_obc, tl_v_obc,          &
#  endif
     &                           tl_ubar_obc, tl_vbar_obc,              &
     &                           tl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
#  endif
     &                           tl_t, tl_u, tl_v,                      &
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta,                               &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           ad_t_obc, ad_u_obc, ad_v_obc,          &
#  endif
     &                           ad_ubar_obc, ad_vbar_obc,              &
     &                           ad_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           ad_tflux,                              &
#  endif
     &                           ad_t, ad_u, ad_v,                      &
# else
     &                           ad_ubar, ad_vbar,                      &
# endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp1, Linp2, Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
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
!  Compute total cost function gradient: add background (Jb) and
!  observation contributions (Jo).  Over-write observation values
!  with total.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lout)=tl_zeta(i,j,Linp1)+                         &
     &                      tl_zeta(i,j,Linp2)+                         &
     &                      ad_zeta(i,j,Lout )
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
              ad_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Linp1)+     &
     &                                  tl_zeta_obc(j,ib,ir,Linp2)+     &
     &                                  ad_zeta_obc(j,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Linp1)+     &
     &                                  tl_zeta_obc(j,ib,ir,Linp2)+     &
     &                                  ad_zeta_obc(j,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Linp1)+     &
     &                                  tl_zeta_obc(i,ib,ir,Linp2)+     &
     &                                  ad_zeta_obc(i,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Linp1)+     &
     &                                  tl_zeta_obc(i,ib,ir,Linp2)+     &
     &                                  ad_zeta_obc(i,ib,ir,Lout )
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
          ad_ubar(i,j,Lout)=tl_ubar(i,j,Linp1)+                         &
     &                      tl_ubar(i,j,Linp2)+                         &
     &                      ad_ubar(i,j,Lout )
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
              ad_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Linp1)+     &
     &                                  tl_ubar_obc(j,ib,ir,Linp2)+     &
     &                                  ad_ubar_obc(j,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Linp1)+     &
     &                                  tl_ubar_obc(j,ib,ir,Linp2)+     &
     &                                  ad_ubar_obc(j,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Linp1)+     &
     &                                  tl_ubar_obc(i,ib,ir,Linp2)+     &
     &                                  ad_ubar_obc(i,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Linp1)+     &
     &                                  tl_ubar_obc(i,ib,ir,Linp2)+     &
     &                                  ad_ubar_obc(i,ib,ir,Lout )
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum component.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lout)=tl_vbar(i,j,Linp1)+                         &
     &                      tl_vbar(i,j,Linp2)+                         &
     &                      ad_vbar(i,j,Lout )
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
              ad_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Linp1)+     &
     &                                  tl_vbar_obc(j,ib,ir,Linp2)+     &
     &                                  ad_vbar_obc(j,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              ad_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Linp1)+     &
     &                                  tl_vbar_obc(j,ib,ir,Linp2)+     &
     &                                  ad_vbar_obc(j,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Linp1)+     &
     &                                  tl_vbar_obc(i,ib,ir,Linp2)+     &
     &                                  ad_vbar_obc(i,ib,ir,Lout )
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Linp1)+     &
     &                                  tl_vbar_obc(i,ib,ir,Linp2)+     &
     &                                  ad_vbar_obc(i,ib,ir,Lout )
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
            ad_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Linp1)+                   &
     &                          tl_ustr(i,j,k,Linp2)+                   &
     &                          ad_ustr(i,j,k,Lout )
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Linp1)+                   &
     &                          tl_vstr(i,j,k,Linp2)+                   &
     &                          ad_vstr(i,j,k,Lout )
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
            ad_u(i,j,k,Lout)=tl_u(i,j,k,Linp1)+                         &
     &                       tl_u(i,j,k,Linp2)+                         &
     &                       ad_u(i,j,k,Lout )
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
                ad_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Linp1)+     &
     &                                   tl_u_obc(j,k,ib,ir,Linp2)+     &
     &                                   ad_u_obc(j,k,ib,ir,Lout )
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                ad_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Linp1)+     &
     &                                   tl_u_obc(j,k,ib,ir,Linp2)+     &
     &                                   ad_u_obc(j,k,ib,ir,Lout )
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Linp1)+     &
     &                                   tl_u_obc(i,k,ib,ir,Linp2)+     &
     &                                   ad_u_obc(i,k,ib,ir,Lout )
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Linp1)+     &
     &                                   tl_u_obc(i,k,ib,ir,Linp2)+     &
     &                                   ad_u_obc(i,k,ib,ir,Lout )
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
            ad_v(i,j,k,Lout)=tl_v(i,j,k,Linp1)+                         &
     &                       tl_v(i,j,k,Linp2)+                         &
     &                       ad_v(i,j,k,Lout )
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
                ad_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Linp1)+     &
     &                                   tl_v_obc(j,k,ib,ir,Linp2)+     &
     &                                   ad_v_obc(j,k,ib,ir,Lout )
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                ad_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Linp1)+     &
     &                                   tl_v_obc(j,k,ib,ir,Linp2)+     &
     &                                   ad_v_obc(j,k,ib,ir,Lout )
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Linp1)+     &
     &                                   tl_v_obc(i,k,ib,ir,Linp2)+     &
     &                                   ad_v_obc(i,k,ib,ir,Lout )
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Linp1)+     &
     &                                   tl_v_obc(i,k,ib,ir,Linp2)+     &
     &                                   ad_v_obc(i,k,ib,ir,Lout )
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
              ad_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Linp1,itrc)+             &
     &                              tl_t(i,j,k,Linp2,itrc)+             &
     &                              ad_t(i,j,k,Lout ,itrc)
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
                  ad_t_obc(j,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(j,k,ib,ir,Linp1,itrc)+    &
     &                               tl_t_obc(j,k,ib,ir,Linp2,itrc)+    &
     &                               ad_t_obc(j,k,ib,ir,Lout ,itrc)
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(itrc),ng)).and.                      &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(j,k,ib,ir,Linp1,itrc)+    &
     &                               tl_t_obc(j,k,ib,ir,Linp2,itrc)+    &
     &                               ad_t_obc(j,k,ib,ir,Lout ,itrc)
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(itrc),ng)).and.                     &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(i,k,ib,ir,Linp1,itrc)+    &
     &                               tl_t_obc(i,k,ib,ir,Linp2,itrc)+    &
     &                               ad_t_obc(i,k,ib,ir,Lout ,itrc)
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(itrc),ng)).and.                     &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Lout,itrc)=                        &
     &                               tl_t_obc(i,k,ib,ir,Linp1,itrc)+    &
     &                               tl_t_obc(i,k,ib,ir,Linp2,itrc)+    &
     &                               ad_t_obc(i,k,ib,ir,Lout ,itrc)
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
                ad_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Linp1,itrc)+   &
     &                                    tl_tflux(i,j,k,Linp2,itrc)+   &
     &                                    ad_tflux(i,j,k,Lout ,itrc)
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE cost_grad_tile
#endif
      END MODULE cost_grad_mod
