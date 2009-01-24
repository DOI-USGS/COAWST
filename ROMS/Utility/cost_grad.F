#include "cppdefs.h"
      MODULE cost_grad_mod
#if defined IS4DVAR || defined IS4DVAR_OLD
!
!svn $Id: cost_grad.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the total cost function gradient by adding    !
!  background and observations contributions in v-space:               !
!                                                                      !
!      GRADv(J) = GRADv(Jb) + GRADv(Jo)                         (1)    !
!                                                                      !
!      GRADv(J) = deltaV + B^(T/2) GRADx(Jo)                    (2)    !
!                                                                      !
!      GRADv(J) = deltaV + W^(-1/2) L^(T/2) G S                 (3)    !
!                                                                      !
!  where                                                               !
!                                                                      !
!      deltaV     Increment vector (TLM solution) in v-space           !
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
      SUBROUTINE cost_grad (ng, tile, Linp, Lout)
!***********************************************************************
!
      USE mod_param
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
      CALL cost_grad_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Linp, Lout,                                  &
# ifdef ADJUST_WSTRESS
    &                      FORCES(ng) % tl_ustr,                        &
    &                      FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
    &                      FORCES(ng) % tl_tflux,                       &
#  endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta,                         &
# ifdef ADJUST_WSTRESS
    &                      FORCES(ng) % ad_ustr,                        &
    &                      FORCES(ng) % ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
    &                      FORCES(ng) % ad_tflux,                       &
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
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Linp, Lout,                            &
# ifdef ADJUST_WSTRESS
    &                            tl_ustr, tl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
    &                            tl_tflux,                              &
#  endif
     &                           tl_t, tl_u, tl_v,                      &
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta,                               &
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
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout
!
# ifdef ASSUMED_SHAPE
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
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
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
      integer :: i, j, k
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
!  Free-surface gradient norm.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lout)=tl_zeta(i,j,Linp)+                          &
     &                      ad_zeta(i,j,Lout)
        END DO
      END DO

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress gradient norm.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Linp)+                    &
     &                          ad_ustr(i,j,k,Lout)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Linp)+                    &
     &                          ad_vstr(i,j,k,Lout)
          END DO
        END DO
      END DO
# endif
# if !defined SOLVE3D
!
!  2D momentum gradient norm.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lout)=tl_ubar(i,j,Linp)+                          &
     &                      ad_ubar(i,j,Lout)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lout)=tl_vbar(i,j,Linp)+                          &
     &                      ad_vbar(i,j,Lout)
        END DO
      END DO
# else
!
!  3D momentum gradient norm.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lout)=tl_u(i,j,k,Linp)+                          &
     &                       ad_u(i,j,k,Lout)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lout)=tl_v(i,j,k,Linp)+                          &
     &                       ad_v(i,j,k,Lout)
          END DO
        END DO
      END DO
!
!  Tracers gradient norm.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Linp,itrc)+              &
     &                              ad_t(i,j,k,Lout,itrc)
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux gradient norm.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Linp,itrc)+      &
     &                                  ad_tflux(i,j,k,Lout,itrc)
            END DO
          END DO
        END DO
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE cost_grad_tile
#endif
      END MODULE cost_grad_mod
