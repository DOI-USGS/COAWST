#include "cppdefs.h"
      MODULE tl_ini_adjust_mod

#if defined IS4DVAR || defined IS4DVAR_OLD
!
!svn $Id: tl_ini_adjust.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine computes new tangent linear model initial conditions   !
!  after the minimization as:                                          !
!                                                                      !
!       deltaX(inner) = B^(1/2) deltaV - (Xi(outer) - Xb)              !
!                                                                      !
!  where                                                               !
!                                                                      !
!            deltaX : new increment vector (TLM initial conditions)    !
!    B^(1/2) deltaV : convolved increment after minimization           !
!                Xi : NLM initial conditions for current outer loop    !
!                Xb : background state                                 ! 
!                                                                      !
!  Here, the TLM state at index Lnew is adjusted by substracting the   !
!  current NLM initial condition departures from the background. The   !
!  convolved vector B^(1/2) deltaV, also at index Lnew, is distroyed.  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: tl_ini_adjust
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      PUBLIC :: tl_frc_adjust
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_ini_adjust (ng, tile, Lbck, Lini, Lnew)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck, Lini, Lnew
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
# endif
      CALL tl_ini_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Lbck, Lini, Lnew,                        &
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                         FORCES(ng) % tflux,                      &
     &                         FORCES(ng) % tl_tflux,                   &
# endif
# ifdef ADJUST_WSTRESS
     &                         FORCES(ng) % ustr,                       &
     &                         FORCES(ng) % tl_ustr,                    &
     &                         FORCES(ng) % vstr,                       &
     &                         FORCES(ng) % tl_vstr,                    &
# endif
# ifdef SOLVE3D
     &                         OCEAN(ng) % t,                           &
     &                         OCEAN(ng) % tl_t,                        &
     &                         OCEAN(ng) % u,                           &
     &                         OCEAN(ng) % tl_u,                        &
     &                         OCEAN(ng) % v,                           &
     &                         OCEAN(ng) % tl_v,                        &
# else
     &                         OCEAN(ng) % ubar,                        &
     &                         OCEAN(ng) % tl_ubar,                     &
     &                         OCEAN(ng) % vbar,                        &
     &                         OCEAN(ng) % tl_vbar,                     &
# endif
     &                         OCEAN(ng) % zeta,                        &
     &                         OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
# endif
      RETURN
      END SUBROUTINE tl_ini_adjust

!
!***********************************************************************
      SUBROUTINE tl_ini_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Lbck, Lini, Lnew,                  &
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                               tflux, tl_tflux,                   &
# endif
# ifdef ADJUST_WSTRESS
     &                               ustr, tl_ustr,                     &
     &                               vstr, tl_vstr,                     &
# endif
# ifdef SOLVE3D
     &                               t, tl_t,                           &
     &                               u, tl_u,                           &
     &                               v, tl_v,                           &
# else
     &                               ubar, tl_ubar,                     &
     &                               vbar, tl_vbar,                     &
# endif
     &                               zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck, Lini, Lnew
!
# ifdef ASSUMED_SHAPE
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tflux(LBi:UBi,LBj:UBj,                    &
     &                              Nfrec(ng),2,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
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
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

# include "set_bounds.h"
!      
!-----------------------------------------------------------------------
!  Adjust tangent linear initial conditions: substract nonlinear model
!  departures from background.
!-----------------------------------------------------------------------
!
!  Free-sruface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lnew)=tl_zeta(i,j,Lnew)-                          &
     &                      (zeta(i,j,Lini)-zeta(i,j,Lbck))
        END DO
      END DO

# ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Lnew)=tl_ubar(i,j,Lnew)-                          &
     &                      (ubar(i,j,Lini)-ubar(i,j,Lbck))
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lnew)=tl_vbar(i,j,Lnew)-                          &
     &                      (vbar(i,j,Lini)-vbar(i,j,Lbck))
        END DO
      END DO
# else
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_u(i,j,k,Lnew)=tl_u(i,j,k,Lnew)-                          &
     &                       (u(i,j,k,Lini)-u(i,j,k,Lbck))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lnew)=tl_v(i,j,k,Lnew)-                          &
     &                       (v(i,j,k,Lini)-v(i,j,k,Lbck))
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_t(i,j,k,Lnew,itrc)=tl_t(i,j,k,Lnew,itrc)-              &
     &                              (t(i,j,k,Lini,itrc)-                &
     &                               t(i,j,k,Lbck,itrc))
            END DO
          END DO          
        END DO
      END DO
# endif
# ifdef ADJUST_WSTRESS
!
!   Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ustr(i,j,k,Lnew)=tl_ustr(i,j,k,Lnew)-                    &
     &                          (ustr(i,j,k,Lini)-ustr(i,j,k,Lbck))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,k,Lnew)=tl_vstr(i,j,k,Lnew)-                    &
     &                          (vstr(i,j,k,Lini)-vstr(i,j,k,Lbck))
          END DO
        END DO
      END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Surface Tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_tflux(i,j,k,Lnew,itrc)=tl_tflux(i,j,k,Lnew,itrc)-      &
     &                                  (tflux(i,j,k,Lini,itrc)-        &
     &                                   tflux(i,j,k,Lbck,itrc))
            END DO
          END DO          
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE tl_ini_adjust_tile

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!***********************************************************************
      SUBROUTINE tl_frc_adjust (ng, tile, Lbck, Lini, Lnew)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck, Lini, Lnew
!
!  Local variable declarations.
!
#  include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL tl_frc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
#  if defined ADJUST_STFLUX && defined SOLVE3D
     &                         FORCES(ng) % tflux,                      &
     &                         FORCES(ng) % tl_tflux,                   &
#  endif
#  ifdef ADJUST_WSTRESS
     &                         FORCES(ng) % ustr,                       &
     &                         FORCES(ng) % tl_ustr,                    &
     &                         FORCES(ng) % vstr,                       &
     &                         FORCES(ng) % tl_vstr,                    &
#  endif
     &                         Lbck, Lini, Lnew)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif
      RETURN
      END SUBROUTINE tl_frc_adjust
!
!***********************************************************************
      SUBROUTINE tl_frc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
#  if defined ADJUST_STFLUX && defined SOLVE3D
     &                               tflux, tl_tflux,                   &
#  endif
#  ifdef ADJUST_WSTRESS
     &                               ustr, tl_ustr,                     &
     &                               vstr, tl_vstr,                     &
#  endif
     &                               Lbck, Lini, Lnew)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck, Lini, Lnew
!
#  ifdef ASSUMED_SHAPE
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tflux(LBi:,LBj:,:,:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: vstr(LBi:,LBj:,:,:)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#   endif
#  else
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tflux(LBi:UBi,LBj:UBj,                    &
     &                              Nfrec(ng),2,NT(ng))
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif

#  include "set_bounds.h"
!      
!-----------------------------------------------------------------------
!  Adjust tangent linear surface forcing fields: add nonlinear model
!  departures from background.
!-----------------------------------------------------------------------
!
#  ifdef ADJUST_WSTRESS
!
!   Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ustr(i,j,k,Lnew)=tl_ustr(i,j,k,Lnew)+                    &
     &                          (ustr(i,j,k,Lini)-ustr(i,j,k,Lbck))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,k,Lnew)=tl_vstr(i,j,k,Lnew)+                    &
     &                          (vstr(i,j,k,Lini)-vstr(i,j,k,Lbck))
          END DO
        END DO
      END DO
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Surface Tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_tflux(i,j,k,Lnew,itrc)=tl_tflux(i,j,k,Lnew,itrc)+      &
     &                                  (tflux(i,j,k,Lini,itrc)-        &
     &                                   tflux(i,j,k,Lbck,itrc))
            END DO
          END DO          
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_frc_adjust_tile
# endif
#endif       
      END MODULE tl_ini_adjust_mod
