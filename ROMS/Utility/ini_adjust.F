#include "cppdefs.h"
      MODULE ini_adjust_mod

#if defined TLM_CHECK  || defined OPT_PERTURBATION || \
    defined FORCING_SV || defined STOCHASTIC_OPT || \
    defined HESSIAN_SV || defined HESSIAN_SO || defined HESSIAN_FSV || \
    (defined ADJOINT   && defined TANGENT)
!
!svn $Id: ini_adjust.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines adjust and perturb state initial conditions.         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# if defined ADJOINT && defined TANGENT
#  if defined IS4DVAR || defined TL_W4DPSAS          || \
      defined W4DPSAS || defined W4DPSAS_SENSITIVITY
      PUBLIC :: ini_adjust
#  endif
# endif
# if defined TANGENT && defined ADJOINT
      PUBLIC :: load_ADtoTL
      PUBLIC :: load_TLtoAD
# endif
# ifdef TLM_CHECK
      PUBLIC :: ini_perturb
# endif
# if defined TL_W4DVAR          || defined W4DVAR      || \
     defined W4DVAR_SENSITIVITY || defined ARRAY_MODES || \
     defined CLIPPING
      PUBLIC :: rp_ini_adjust
# endif
# ifdef TLM_CHECK
      PUBLIC :: tl_ini_perturb
# endif
# if defined OPT_PERTURBATION || defined FORCING_SV || \
     defined STOCHASTIC_OPT   || defined HESSIAN_SV || \
     defined HESSIAN_SO       || defined HESSIAN_FSV
      PUBLIC :: ad_ini_perturb
# endif

      CONTAINS

# if defined ADJOINT && defined TANGENT
#  if defined IS4DVAR
      SUBROUTINE ini_adjust (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  This routine adds 4DVAR inner loop tangent linear increments to     !
!  nonlinear state initial conditions.  The boundary condition and     !
!  barotropic/baroclinic coupling, if any, are processed latter in     !
!  routine "ini_fields" before time-stepping.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Nonlinear state time index to update.                  !
!     tile      Domain partition.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
#   include "tile.h"
!
#   ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#   endif
      CALL ini_adjust_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Linp, Lout,                                 &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#   endif
#   ifdef SOLVE3D
     &                      OCEAN(ng) % tl_t,                           &
     &                      OCEAN(ng) % tl_u,                           &
     &                      OCEAN(ng) % tl_v,                           &
#   endif
     &                      OCEAN(ng) % tl_ubar,                        &
     &                      OCEAN(ng) % tl_vbar,                        &
     &                      OCEAN(ng) % tl_zeta,                        &
#   ifdef SOLVE3D
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
#   endif
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)
#   ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#   endif

      RETURN
      END SUBROUTINE ini_adjust
!
!***********************************************************************
      SUBROUTINE ini_adjust_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Linp, Lout,                           &
#   ifdef MASKING
     &                            rmask, umask, vmask,                  &
#   endif
#   ifdef SOLVE3D
     &                            tl_t, tl_u, tl_v,                     &
#   endif
     &                            tl_ubar, tl_vbar, tl_zeta,            &
#   ifdef SOLVE3D
     &                            t, u, v,                              &
#   endif
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
#   ifdef ASSUMED_SHAPE
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#    endif
#    ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#    ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#   else
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#    endif
#    ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#    endif
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#    ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#    endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#   endif
!
!  Local variable declarations.
!
      integer :: i, j
#   ifdef SOLVE3D
      integer :: itrc, k
#   endif

#   include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 2D state variables by adding tangent
!  linear increments from data assimilation.
!-----------------------------------------------------------------------
!
#   ifndef SOLVE3D
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ubar(i,j,Lout)=ubar(i,j,Lout)+tl_ubar(i,j,Linp)
#    ifdef MASKING
          ubar(i,j,Lout)=ubar(i,j,Lout)*umask(i,j)
#    endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vbar(i,j,Lout)=vbar(i,j,Lout)+tl_vbar(i,j,Linp)
#    ifdef MASKING
          vbar(i,j,Lout)=vbar(i,j,Lout)*vmask(i,j)
#    endif
        END DO
      END DO
#   endif

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          zeta(i,j,Lout)=zeta(i,j,Lout)+tl_zeta(i,j,Linp)
#   ifdef MASKING
          zeta(i,j,Lout)=zeta(i,j,Lout)*rmask(i,j)
#   endif
        END DO
      END DO

#   ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 3D state variables by adding tangent
!  linear increments from data assimilation.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            u(i,j,k,Lout)=u(i,j,k,Lout)+tl_u(i,j,k,Linp)
#    ifdef MASKING
            u(i,j,k,Lout)=u(i,j,k,Lout)*umask(i,j)
#    endif
          END DO
        END DO

        DO j=JstrP,JendT
          DO i=IstrT,IendT
            v(i,j,k,Lout)=v(i,j,k,Lout)+tl_v(i,j,k,Linp)
#    ifdef MASKING
            v(i,j,k,Lout)=v(i,j,k,Lout)*vmask(i,j)
#    endif
          END DO
        END DO
      END DO

      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)+                    &
     &                           tl_t(i,j,k,Linp,itrc)
#    ifdef MASKING
              t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)*rmask(i,j)
#    endif
            END DO
          END DO
        END DO
      END DO
#   endif

      RETURN
      END SUBROUTINE ini_adjust_tile
#  endif

#  if defined TL_W4DPSAS          || defined W4DPSAS || \
      defined W4DPSAS_SENSITIVITY
      SUBROUTINE ini_adjust (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  This routine adds convolved to adjoint solutions to nonlinear       !
!  conditions.  The boundary condition and barotropic/baroclinic       !
!  coupling, if any, are processed latter in routine "ini_fields"      !
!  before time-stepping.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Nonlinear state time index to update.                  !
!     tile      Domain partition.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
#   include "tile.h"
!
#   ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#   endif
      CALL ini_adjust_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Linp, Lout,                                 &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#   endif
#   ifdef SOLVE3D
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
#   endif
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
     &                      OCEAN(ng) % ad_zeta,                        &
#   ifdef SOLVE3D
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
#   endif
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)
#   ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#   endif

      RETURN
      END SUBROUTINE ini_adjust
!
!***********************************************************************
      SUBROUTINE ini_adjust_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Linp, Lout,                           &
#   ifdef MASKING
     &                            rmask, umask, vmask,                  &
#   endif
#   ifdef SOLVE3D
     &                            ad_t, ad_u, ad_v,                     &
#   endif
     &                            ad_ubar, ad_vbar, ad_zeta,            &
#   ifdef SOLVE3D
     &                            t, u, v,                              &
#   endif
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
#   ifdef ASSUMED_SHAPE
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#    endif
#    ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#    ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#   else
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#    endif
#    ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#    endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#    ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#    endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#   endif
!
!  Local variable declarations.
!
      integer :: i, j
#   ifdef SOLVE3D
      integer :: itrc, k
#   endif

#   include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 2D state variables by adding tangent
!  linear increments from data assimilation.
!-----------------------------------------------------------------------
!
#   ifndef SOLVE3D
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ubar(i,j,Lout)=ubar(i,j,Lout)+ad_ubar(i,j,Linp)
#    ifdef MASKING
          ubar(i,j,Lout)=ubar(i,j,Lout)*umask(i,j)
#    endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vbar(i,j,Lout)=vbar(i,j,Lout)+ad_vbar(i,j,Linp)
#    ifdef MASKING
          vbar(i,j,Lout)=vbar(i,j,Lout)*vmask(i,j)
#    endif
        END DO
      END DO
#   endif

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          zeta(i,j,Lout)=zeta(i,j,Lout)+ad_zeta(i,j,Linp)
#   ifdef MASKING
          zeta(i,j,Lout)=zeta(i,j,Lout)*rmask(i,j)
#   endif
        END DO
      END DO

#   ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 3D state variables by adding tangent
!  linear increments from data assimilation.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            u(i,j,k,Lout)=u(i,j,k,Lout)+ad_u(i,j,k,Linp)
#    ifdef MASKING
            u(i,j,k,Lout)=u(i,j,k,Lout)*umask(i,j)
#    endif
          END DO
        END DO

        DO j=JstrP,JendT
          DO i=IstrT,IendT
            v(i,j,k,Lout)=v(i,j,k,Lout)+ad_v(i,j,k,Linp)
#    ifdef MASKING
            v(i,j,k,Lout)=v(i,j,k,Lout)*vmask(i,j)
#    endif
          END DO
        END DO
      END DO

      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)+                    &
     &                           ad_t(i,j,k,Linp,itrc)
#    ifdef MASKING
              t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)*rmask(i,j)
#    endif
            END DO
          END DO
        END DO
      END DO
#   endif

      RETURN
      END SUBROUTINE ini_adjust_tile
#  endif

#  if defined TL_W4DVAR          || defined W4DVAR      || \
      defined W4DVAR_SENSITIVITY || defined ARRAY_MODES || \
      defined CLIPPING
      SUBROUTINE rp_ini_adjust (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  This routine adds weak constraint adjoint increments to nonlinear   !
!  state initial conditions  (background state)  at the end of inner   !
!  loop.  The boundary condition and barotropic/baroclinic coupling,   !
!  if any, are processed latter in routine "ini_fields".               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Nonlinear state time index to update.                  !
!     tile      Domain partition.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
#   include "tile.h"
!
#   ifdef PROFILE
      CALL wclock_on (ng, iRPM, 7)
#   endif
      CALL rp_ini_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp, Lout,                              &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
#   endif
#   ifdef SOLVE3D
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % ad_v,                        &
#   else
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ad_vbar,                     &
#   endif
     &                         OCEAN(ng) % ad_zeta,                     &
#   ifdef SOLVE3D
     &                         OCEAN(ng) % t,                           &
     &                         OCEAN(ng) % u,                           &
     &                         OCEAN(ng) % v,                           &
#   else
     &                         OCEAN(ng) % ubar,                        &
     &                         OCEAN(ng) % vbar,                        &
#   endif
     &                         OCEAN(ng) % zeta,                        &
#   ifdef SOLVE3D
     &                         OCEAN(ng) % tl_t,                        &
     &                         OCEAN(ng) % tl_u,                        &
     &                         OCEAN(ng) % tl_v,                        &
#   else
     &                         OCEAN(ng) % tl_ubar,                     &
     &                         OCEAN(ng) % tl_vbar,                     &
#   endif
     &                         OCEAN(ng) % tl_zeta)
#   ifdef PROFILE
      CALL wclock_off (ng, iRPM, 7)
#   endif

      RETURN
      END SUBROUTINE rp_ini_adjust
!
!***********************************************************************
      SUBROUTINE rp_ini_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp, Lout,                        &
#   ifdef MASKING
     &                               rmask, umask, vmask,               &
#   endif
#   ifdef SOLVE3D
     &                               ad_t, ad_u, ad_v,                  &
#   else
     &                               ad_ubar, ad_vbar,                  &
#   endif
     &                               ad_zeta,                           &
#   ifdef SOLVE3D
     &                               t, u, v,                           &
#   else
     &                               ubar, vbar,                        &
#   endif
     &                               zeta,                              &
#   ifdef SOLVE3D
     &                               tl_t, tl_u, tl_v,                  &
#   else
     &                               tl_ubar, tl_vbar,                  &
#   endif
     &                               tl_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
#   ifdef ASSUMED_SHAPE
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#    endif
#    ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#    endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#    ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
#    endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#    ifdef SOLVE3D
      real(r8), intent(out) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(out) :: tl_v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(out) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_vbar(LBi:,LBj:,:)
#    endif
      real(r8), intent(out) :: tl_zeta(LBi:,LBj:,:)
#   else
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#    endif
#    ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#    else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#    endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#    ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#    else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
#    endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#    ifdef SOLVE3D
      real(r8), intent(out) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(out) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#    else
      real(r8), intent(out) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#    endif
      real(r8), intent(out) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   endif
!
!  Local variable declarations.
!
      integer :: i, j
#   ifdef SOLVE3D
      integer :: itrc, k
#   endif

#   include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 2D state variables by adding adjoint
!  increments from weak constraint inner loop.
!-----------------------------------------------------------------------
!
#   ifndef SOLVE3D
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          tl_ubar(i,j,Lout)=ubar(i,j,Linp)+ad_ubar(i,j,Lout)
#    ifdef MASKING
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#    endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
          tl_vbar(i,j,Lout)=vbar(i,j,Linp)+ad_vbar(i,j,Lout)
#    ifdef MASKING
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#    endif
        END DO
      END DO
#   endif

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          tl_zeta(i,j,Lout)=zeta(i,j,Linp)+ad_zeta(i,j,Lout)
#   ifdef MASKING
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
#   endif
        END DO
      END DO

#   ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 3D state variables by adding adjoint
!  increments from weak constraint inner loop.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            tl_u(i,j,k,Lout)=u(i,j,k,Linp)+ad_u(i,j,k,Lout)
#    ifdef MASKING
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#    endif
          END DO
        END DO

        DO j=JstrP,JendT
          DO i=IstrT,IendT
            tl_v(i,j,k,Lout)=v(i,j,k,Linp)+ad_v(i,j,k,Lout)
#    ifdef MASKING
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#    endif
          END DO
        END DO
      END DO

      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              tl_t(i,j,k,Lout,itrc)=t(i,j,k,Linp,itrc)+                 &
     &                              ad_t(i,j,k,Lout,itrc)
#    ifdef MASKING
              tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
#    endif
            END DO
          END DO
        END DO
      END DO
#   endif

      RETURN
      END SUBROUTINE rp_ini_adjust_tile
#  endif
# endif

# if defined TANGENT && defined ADJOINT

      SUBROUTINE load_ADtoTL (ng, tile, Linp, Lout, add)
!
!=======================================================================
!                                                                      !
!  This routine loads or adds Linp adjoint state variables into Lout   !
!  Lout tangent linear state variables.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Nonlinear state time index to update.                  !
!     add       Logical switch to add to imported values.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef ADJUST_BOUNDARY
      USE mod_boundary
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
#  endif
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      logical, intent(in) :: add
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL load_ADtoTL_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Linp, Lout, add,                           &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       BOUNDARY(ng) % ad_t_obc,                   &
     &                       BOUNDARY(ng) % ad_u_obc,                   &
     &                       BOUNDARY(ng) % ad_v_obc,                   &
#   endif
     &                       BOUNDARY(ng) % ad_ubar_obc,                &
     &                       BOUNDARY(ng) % ad_vbar_obc,                &
     &                       BOUNDARY(ng) % ad_zeta_obc,                &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % ad_ustr,                      &
     &                       FORCES(ng) % ad_vstr,                      &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       FORCES(ng) % ad_tflux,                     &
#   endif
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
#   endif
#  else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
#  endif
     &                       OCEAN(ng) % ad_zeta,                       &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       BOUNDARY(ng) % tl_t_obc,                   &
     &                       BOUNDARY(ng) % tl_u_obc,                   &
     &                       BOUNDARY(ng) % tl_v_obc,                   &
#   endif
     &                       BOUNDARY(ng) % tl_ubar_obc,                &
     &                       BOUNDARY(ng) % tl_vbar_obc,                &
     &                       BOUNDARY(ng) % tl_zeta_obc,                &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % tl_ustr,                      &
     &                       FORCES(ng) % tl_vstr,                      &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       FORCES(ng) % tl_tflux,                     &
#   endif
     &                       OCEAN(ng) % tl_t,                          &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
#   endif
#  else
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
#  endif
     &                       OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif

      RETURN
      END SUBROUTINE load_ADtoTL
!
!***********************************************************************
      SUBROUTINE load_ADtoTL_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Linp, Lout, add,                     &
#  ifdef MASKING
     &                             rmask, umask, vmask,                 &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                             ad_t_obc, ad_u_obc, ad_v_obc,        &
#   endif
     &                             ad_ubar_obc, ad_vbar_obc,            &
     &                             ad_zeta_obc,                         &
#  endif
#  ifdef ADJUST_WSTRESS
     &                             ad_ustr, ad_vstr,                    &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                             ad_tflux,                            &
#   endif
     &                             ad_t, ad_u, ad_v,                    &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                             ad_ubar, ad_vbar,                    &
#   endif
#  else
     &                             ad_ubar, ad_vbar,                    &
#  endif
     &                             ad_zeta,                             &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                             tl_t_obc, tl_u_obc, tl_v_obc,        &
#   endif
     &                             tl_ubar_obc, tl_vbar_obc,            &
     &                             tl_zeta_obc,                         &
#  endif
#  ifdef ADJUST_WSTRESS
     &                             tl_ustr, tl_vstr,                    &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                             tl_tflux,                            &
#   endif
     &                             tl_t, tl_u, tl_v,                    &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                             tl_ubar, tl_vbar,                    &
#   endif
#  else
     &                             tl_ubar, tl_vbar,                    &
#  endif
     &                             tl_zeta)
!***********************************************************************
!
      USE mod_param
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#  endif
!
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
!
!  Imported variable declarations.
!
      logical, intent(in) :: add
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(in) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: ad_v_obc(LBij:,:,:,:,:)
#    endif
      real(r8), intent(in) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: ad_zeta_obc(LBij:,:,:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_vstr(LBi:,LBj:,:,:)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:,LBj:,:,:,:)
#    endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#    endif
#   else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#    endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#    endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#    endif
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(in) :: ad_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#    endif
      real(r8), intent(in) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#    endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#    endif
#   else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#    endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#    endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#    endif
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
      real(r8) :: fac1, fac2
!
!-----------------------------------------------------------------------
!  Load or add adjoint state variables into tangent linear state
!  variables.
!-----------------------------------------------------------------------
!
!  Add adjoint state to tangent linear state.
!
!    tl_var(Lout) = fac1 * tl_var(Lout) + fac2 * ad_var(Linp)
!
      IF (add) THEN

        fac1=1.0_r8
        fac2=1.0_r8

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lout, Linp, Lout, fac1, fac2,              &
#  ifdef MASKING
     &                       rmask, umask, vmask,                       &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       tl_t_obc, ad_t_obc,                        &
     &                       tl_u_obc, ad_u_obc,                        &
     &                       tl_v_obc, ad_v_obc,                        &
#   endif
     &                       tl_ubar_obc, ad_ubar_obc,                  &
     &                       tl_vbar_obc, ad_vbar_obc,                  &
     &                       tl_zeta_obc, ad_zeta_obc,                  &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       tl_ustr, ad_ustr,                          &
     &                       tl_vstr, ad_vstr,                          &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       tl_tflux, ad_tflux,                        &
#   endif
     &                       tl_t, ad_t,                                &
     &                       tl_u, ad_u,                                &
     &                       tl_v, ad_v,                                &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       tl_ubar, ad_ubar,                          &
     &                       tl_vbar, ad_vbar,                          &
#   endif
#  else
     &                       tl_ubar, ad_ubar,                          &
     &                       tl_vbar, ad_vbar,                          &
#  endif
     &                       tl_zeta, ad_zeta)
!
!  Otherwise, copy adjoint state into tangent linear state.
!
!    tl_var(Lout) = ad_var(Linp)
!
      ELSE

        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Linp, Lout,                                    &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                   tl_t_obc, ad_t_obc,                            &
     &                   tl_u_obc, ad_u_obc,                            &
     &                   tl_v_obc, ad_v_obc,                            &
#   endif
     &                   tl_ubar_obc, ad_ubar_obc,                      &
     &                   tl_vbar_obc, ad_vbar_obc,                      &
     &                   tl_zeta_obc, ad_zeta_obc,                      &
#  endif
#  ifdef ADJUST_WSTRESS
     &                   tl_ustr, ad_ustr,                              &
     &                   tl_vstr, ad_vstr,                              &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                   tl_tflux, ad_tflux,                            &
#   endif
     &                   tl_t, ad_t,                                    &
     &                   tl_u, ad_u,                                    &
     &                   tl_v, ad_v,                                    &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                   tl_ubar, ad_ubar,                              &
     &                   tl_vbar, ad_vbar,                              &
#   endif
#  else
     &                   tl_ubar, ad_ubar,                              &
     &                   tl_vbar, ad_vbar,                              &
#  endif
     &                   tl_zeta, ad_zeta)
      END IF

      RETURN
      END SUBROUTINE load_ADtoTL_tile

      SUBROUTINE load_TLtoAD (ng, tile, Linp, Lout, add)
!
!=======================================================================
!                                                                      !
!  This routine loads or adds Linp tangent linear state variables into !
!  Lout adjoint state variables.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Nonlinear state time index to update.                  !
!     add       Logical switch to add to imported values.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef ADJUST_BOUNDARY
      USE mod_boundary
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
#  endif
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      logical, intent(in) :: add
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      CALL load_TLtoAD_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Linp, Lout, add,                           &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       BOUNDARY(ng) % tl_t_obc,                   &
     &                       BOUNDARY(ng) % tl_u_obc,                   &
     &                       BOUNDARY(ng) % tl_v_obc,                   &
#   endif
     &                       BOUNDARY(ng) % tl_ubar_obc,                &
     &                       BOUNDARY(ng) % tl_vbar_obc,                &
     &                       BOUNDARY(ng) % tl_zeta_obc,                &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % tl_ustr,                      &
     &                       FORCES(ng) % tl_vstr,                      &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       FORCES(ng) % tl_tflux,                     &
#   endif
     &                       OCEAN(ng) % tl_t,                          &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
#   endif
#  else
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
#  endif
     &                       OCEAN(ng) % tl_zeta,                       &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       BOUNDARY(ng) % ad_t_obc,                   &
     &                       BOUNDARY(ng) % ad_u_obc,                   &
     &                       BOUNDARY(ng) % ad_v_obc,                   &
#   endif
     &                       BOUNDARY(ng) % ad_ubar_obc,                &
     &                       BOUNDARY(ng) % ad_vbar_obc,                &
     &                       BOUNDARY(ng) % ad_zeta_obc,                &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % ad_ustr,                      &
     &                       FORCES(ng) % ad_vstr,                      &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       FORCES(ng) % ad_tflux,                     &
#   endif
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
#   endif
#  else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
#  endif
     &                       OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif

      RETURN
      END SUBROUTINE load_TLtoAD
!
!***********************************************************************
      SUBROUTINE load_TLtoAD_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Linp, Lout, add,                     &
#  ifdef MASKING
     &                             rmask, umask, vmask,                 &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                             tl_t_obc, tl_u_obc, tl_v_obc,        &
#   endif
     &                             tl_ubar_obc, tl_vbar_obc,            &
     &                             tl_zeta_obc,                         &
#  endif
#  ifdef ADJUST_WSTRESS
     &                             tl_ustr, tl_vstr,                    &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                             tl_tflux,                            &
#   endif
     &                             tl_t, tl_u, tl_v,                    &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                             tl_ubar, tl_vbar,                    &
#   endif
#  else
     &                             tl_ubar, tl_vbar,                    &
#  endif
     &                             tl_zeta,                             &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                             ad_t_obc, ad_u_obc, ad_v_obc,        &
#   endif
     &                             ad_ubar_obc, ad_vbar_obc,            &
     &                             ad_zeta_obc,                         &
#  endif
#  ifdef ADJUST_WSTRESS
     &                             ad_ustr, ad_vstr,                    &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                             ad_tflux,                            &
#   endif
     &                             ad_t, ad_u, ad_v,                    &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                             ad_ubar, ad_vbar,                    &
#   endif
#  else
     &                             ad_ubar, ad_vbar,                    &
#  endif
     &                             ad_zeta)
!***********************************************************************
!
      USE mod_param
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#  endif
!
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
!
!  Imported variable declarations.
!
      logical, intent(in) :: add
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: tl_v_obc(LBij:,:,:,:,:)
#    endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_zeta_obc(LBij:,:,:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_vstr(LBi:,LBj:,:,:)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#    endif
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#    endif
#   else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#    endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#    endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#    endif
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#    endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(in) :: tl_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#    endif
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#    endif
#   else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#    endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#   ifdef SOLVE3D
#    ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#    endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#    endif
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
      real(r8) :: fac1, fac2
!
!-----------------------------------------------------------------------
!  Load or add tangent linear state variables into adjoint state
!  variables.
!-----------------------------------------------------------------------
!
!  Add tangent linear state to adjoint state.
!
!    ad_var(Lout) = fac1 * ad_var(Lout) + fac2 * tl_var(Linp)
!
      IF (add) THEN

        fac1=1.0_r8
        fac2=1.0_r8

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lout, Linp, Lout, fac1, fac2,              &
#  ifdef MASKING
     &                       rmask, umask, vmask,                       &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#   endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#   endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#   endif
#  else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#  endif
     &                       ad_zeta, tl_zeta)
!
!  Otherwise, copy tangent linear state into adjoint state.
!
!    ad_var(Lout) = tl_var(Linp)
!
      ELSE

        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Linp, Lout,                                    &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                   ad_t_obc, tl_t_obc,                            &
     &                   ad_u_obc, tl_u_obc,                            &
     &                   ad_v_obc, tl_v_obc,                            &
#   endif
     &                   ad_ubar_obc, tl_ubar_obc,                      &
     &                   ad_vbar_obc, tl_vbar_obc,                      &
     &                   ad_zeta_obc, tl_zeta_obc,                      &
#  endif
#  ifdef ADJUST_WSTRESS
     &                   ad_ustr, tl_ustr,                              &
     &                   ad_vstr, tl_vstr,                              &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                   ad_tflux, tl_tflux,                            &
#   endif
     &                   ad_t, tl_t,                                    &
     &                   ad_u, tl_u,                                    &
     &                   ad_v, tl_v,                                    &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                   ad_ubar, tl_ubar,                              &
     &                   ad_vbar, tl_vbar,                              &
#   endif
#  else
     &                   ad_ubar, tl_ubar,                              &
     &                   ad_vbar, tl_vbar,                              &
#  endif
     &                   ad_zeta, tl_zeta)
      END IF

      RETURN
      END SUBROUTINE load_TLtoAD_tile
# endif

# if defined TLM_CHECK
      SUBROUTINE ini_perturb (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  Add a perturbation to nonlinear state variables according to the    !
!  outer and inner loop iterations. The added term is a function of    !
!  the steepest descent direction (grad(J)) times the  perturbation    !
!  amplitude "p" (controlled by inner loop).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
#  if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#  endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Linp, Lout, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#  endif
      CALL ini_perturb_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng), Linp, Lout,            &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                       SEDBED(ng) % bed_thick,                    &
#   endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_u,                           &
#   ifdef ICESHELF
     &                       GRID(ng) % zice,                           &
#   endif
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
     &                       COUPLING(ng) % Zt_avg1,                    &
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
     &                       OCEAN(ng) % t,                             &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
#  endif
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
     &                       OCEAN(ng) % ad_zeta,                       &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#  endif

      RETURN
      END SUBROUTINE ini_perturb
!
!***********************************************************************
      SUBROUTINE ini_perturb_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew, Linp, Lout,              &
#  ifdef MASKING
     &                             rmask, umask, vmask,                 &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                             bed_thick,                           &
#   endif
     &                             Hz, h, om_v, on_u,                   &
#   ifdef ICESHELF
     &                             zice,                                &
#   endif
     &                             z_r, z_w, Zt_avg1,                   &
     &                             ad_t, ad_u, ad_v,                    &
     &                             t, u, v,                             &
#  endif
     &                             ad_ubar, ad_vbar, ad_zeta,           &
     &                             ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars
!
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#   ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#   endif
#  endif
#  ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
#  endif
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
      USE zetabc_mod, ONLY : zetabc_tile
#  ifdef SOLVE3D
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, Linp, Lout
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed_thick(LBi:,LBj:,:)
#    endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed_thick(LBi:UBi,LBj:UBj,2)
#    endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: i, ic, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      integer, dimension(NstateVar(ng)+1) :: StateVar

      real(r8) :: p, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variables: steepest descent
!  direction times the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
!  Set state variable to perturb according to outer loop index.
!
#  ifdef SOLVE3D
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
      StateVar(5)=isUvel
      StateVar(6)=isVvel
      DO i=1,NT(ng)
        StateVar(6+i)=isTvar(i)
      END DO
#  else
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
#  endif
!
!  Set perturbation amplitude as a function of the inner loop.
!
      p=10.0_r8**REAL(-inner,r8)
      scale=1.0_r8/SQRT(adDotProduct)
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        IF (Master) WRITE (stdout,10)
      END IF
!
!  Free-surface.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isFsur)) THEN
        DO j=JstrB,JendB
          DO i=IstrB,IendB
            zeta(i,j,Lout)=zeta(i,j,Lout)+p*ad_zeta(i,j,Linp)*scale
#  ifdef MASKING
            zeta(i,j,Lout)=zeta(i,j,Lout)*rmask(i,j)
#  endif
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isFsur)))
          END IF
        END IF
      END IF
!
!  2D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUbar)) THEN
        DO j=JstrB,JendB
          DO i=IstrM,IendB
            ubar(i,j,Lout)=ubar(i,j,Lout)+p*ad_ubar(i,j,Linp)*scale
#  ifdef MASKING
            ubar(i,j,Lout)=ubar(i,j,Lout)*umask(i,j)
#  endif
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUbar)))
          END IF
        END IF
      END IF
!
!  2D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVbar)) THEN
        DO j=JstrM,JendB
          DO i=IstrB,IendB
            vbar(i,j,Lout)=vbar(i,j,Lout)+p*ad_vbar(i,j,Linp)*scale
#  ifdef MASKING
            vbar(i,j,Lout)=vbar(i,j,Lout)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVbar)))
          END IF
        END IF
      END IF
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUvel)) THEN
        DO k=1,N(ng)
          DO j=JstrB,JendB
            DO i=IstrM,IendB
              u(i,j,k,Lout)=u(i,j,k,Lout)+p*ad_u(i,j,k,Linp)*scale
#   ifdef MASKING
              u(i,j,k,Lout)=u(i,j,k,Lout)*umask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUvel)))
          END IF
        END IF
      END IF
!
!  3D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVvel)) THEN
        DO k=1,N(ng)
          DO j=JstrM,JendB
            DO i=IstrB,IendB
              v(i,j,k,Lout)=v(i,j,k,Lout)+p*ad_v(i,j,k,Linp)*scale
#   ifdef MASKING
              v(i,j,k,Lout)=v(i,j,k,Lout)*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVvel)))
          END IF
        END IF
      END IF
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        IF ((StateVar(outer).eq.0).or.                                  &
     &      (StateVar(outer).eq.isTvar(itrc))) THEN
          DO k=1,N(ng)
            DO j=JstrB,JendB
              DO i=IstrB,IendB
                t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)+                  &
     &                             p*ad_t(i,j,k,Linp,itrc)*scale
#   ifdef MASKING
                t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) outer, inner,                           &
     &                          TRIM(Vname(1,idSvar(isTvar(itrc))))
            END IF
          END IF
        END IF
      END DO
#  endif
      IF (Master) WRITE (stdout,'(/)')
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
      CALL zetabc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  Lout, Lout, Lout,                               &
     &                  zeta)
#  ifndef SOLVE3D
      CALL u2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
#  endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,Lout))
#  ifndef SOLVE3D
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,Lout))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,Lout))
      END IF
#  endif

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta(:,:,Lout))
#   ifndef SOLVE3D
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar(:,:,Lout),                               &
     &                    vbar(:,:,Lout))
#   endif
#  endif
#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
#   ifdef ICESHELF
     &                     zice,                                        &
#   endif
#   if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
#   endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to perturbed 3D fields.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 Lout, Lout, u)
      CALL v3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 Lout, Lout, v)
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        CALL t3dbc_tile (ng, tile, itrc, ic,                            &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   Lout, Lout, t)
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,Lout))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,Lout))
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,Lout,itrc))
        END DO
      END IF

#   ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u(:,:,:,Lout),                                &
     &                    v(:,:,:,Lout))
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,Lout,:))
#   endif
#  endif
!
 10   FORMAT (/,' Perturbing Nonlinear model Initial Conditions:',/)
 20   FORMAT (' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,             &
     &        'Perturbing State Variable: ',a)

      RETURN
      END SUBROUTINE ini_perturb_tile
# endif

# if defined OPT_PERTURBATION || defined FORCING_SV || \
     defined STOCHASTIC_OPT   || defined HESSIAN_SV || \
     defined HESSIAN_SO       || defined HESSIAN_FSV

      SUBROUTINE ad_ini_perturb (ng, tile, Kinp, Kout, Ninp, Nout)
!
!=======================================================================
!                                                                      !
!  Initialize adjoint state variables with tangent linear state scaled !
!  by the energy norm, as required by the Generalized stability Theory !
!  propagators.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kinp, Kout, Ninp, Nout
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      CALL ad_ini_perturb_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Kinp, Kout, Ninp, Nout,                 &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
#  endif
#  ifdef BNORM
#    ifdef SOLVE3D
     &                          OCEAN(ng) % e_t,                        &
     &                          OCEAN(ng) % e_u,                        &
     &                          OCEAN(ng) % e_v,                        &
#    else
     &                          OCEAN(ng) % e_ubar,                     &
     &                          OCEAN(ng) % e_vbar,                     &
#    endif
     &                          OCEAN(ng) % e_zeta,                     &
#  endif
#  ifdef SOLVE3D
     &                          GRID(ng) % Hz,                          &
#  else
     &                          GRID(ng) % h,                           &
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
#  endif
     &                          OCEAN(ng) % tl_zeta,                    &
#  ifdef SOLVE3D
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#  else
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#  endif
     &                          OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif

      RETURN
      END SUBROUTINE ad_ini_perturb
!
!***********************************************************************
      SUBROUTINE ad_ini_perturb_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Kinp, Kout, Ninp, Nout,           &
#  ifdef MASKING
     &                                rmask, umask, vmask,              &
#  endif
#  ifdef BNORM
#    ifdef SOLVE3D
     &                                e_t, e_u, e_v,                    &
#    else
     &                                e_ubar, e_vbar,                   &
#    endif
     &                                e_zeta,                           &
#  endif
#  ifdef SOLVE3D
     &                                Hz,                               &
#  else
     &                                h,                                &
     &                                tl_ubar, tl_vbar,                 &
#  endif
     &                                tl_zeta,                          &
#  ifdef SOLVE3D
     &                                tl_t, tl_u, tl_v,                 &
     &                                ad_t, ad_u, ad_v,                 &
#  else
     &                                ad_ubar, ad_vbar,                 &
#  endif
     &                                ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Kout, Ninp, Nout
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   else
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   endif
#   ifdef BNORM
      real(r8), intent(in) :: e_zeta(LBi:,LBj:,:)
#    ifdef SOLVE3D
      real(r8), intent(in) :: e_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: e_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: e_v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(in) :: e_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: e_vbar(LBi:,LBj:,:)
#    endif
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#   ifdef BNORM
      real(r8), intent(in) :: e_zeta(LBi:UBi,LBj:UBj,NSA)
#    ifdef SOLVE3D
      real(r8), intent(in) :: e_t(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: e_u(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(in) :: e_v(LBi:UBi,LBj:UBj,N(ng),NSA)
#    else
      real(r8), intent(in) :: e_ubar(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: e_vbar(LBi:UBi,LBj:UBj,NSA)
#    endif
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: cff, scale
!
#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint state variables with tangent linear state scaled
!  by the energy norm.
!-----------------------------------------------------------------------

#  ifdef FULL_GRID
#   define I_R_RANGE IstrT,IendT
#   define I_U_RANGE IstrP,IendT
#   define J_R_RANGE JstrT,JendT
#   define J_V_RANGE JstrP,JendT
#  else
#   define I_R_RANGE Istr,Iend
#   define I_U_RANGE IstrU,Iend
#   define J_R_RANGE Jstr,Jend
#   define J_V_RANGE JstrV,Jend
#  endif
!
!  Free-surface.
!
#  ifndef BNORM
      scale=0.5_r8*g*rho0
#  endif
      DO j=J_R_RANGE
        DO i=I_R_RANGE
#  ifdef BNORM
          IF (e_zeta(i,j,1).gt.0.0_r8) THEN
            scale=1.0_r8/(e_zeta(i,j,1)*e_zeta(i,j,1))
          ELSE
            scale=1.0_r8
          END IF
#  endif
          ad_zeta(i,j,Kout)=scale*tl_zeta(i,j,Kinp)
#  ifdef MASKING
          ad_zeta(i,j,Kout)=ad_zeta(i,j,Kout)*rmask(i,j)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      cff=0.25_r8*rho0
      DO j=J_R_RANGE
        DO i=I_U_RANGE
#  ifdef BNORM
          IF (e_ubar(i,j,1).gt.0.0_r8) THEN
            scale=1.0_r8/(e_ubar(i,j,1)*e_ubar(i,j,1))
          ELSE
            scale=1.0_r8
          ENDIF
#  else
          scale=cff*(h(i-1,j)+h(i,j))
#  endif
          ad_ubar(i,j,Kout)=scale*tl_ubar(i,j,Kinp)
#   ifdef MASKING
          ad_ubar(i,j,Kout)=ad_ubar(i,j,Kout)*umask(i,j)
#   endif
        END DO
      END DO
!
!  2D v-momentum component.
!
      cff=0.25_r8*rho0
      DO j=J_V_RANGE
        DO i=I_R_RANGE
#  ifdef BNORM
          IF (e_vbar(i,j,1).gt.0.0_r8) THEN
            scale=1.0_r8/(e_vbar(i,j,1)*e_vbar(i,j,1))
          ELSE
            scale=1.0_r8
          ENDIF
#  else
          scale=cff*(h(i,j-1)+h(i,j))
#  endif
          ad_vbar(i,j,Kout)=scale*tl_vbar(i,j,Kinp)
#   ifdef MASKING
          ad_vbar(i,j,Kout)=ad_vbar(i,j,Kout)*vmask(i,j)
#   endif
        END DO
      END DO
#  else
!
!  3D u-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=J_R_RANGE
          DO i=I_U_RANGE
#  ifdef BNORM
            IF (e_u(i,j,k,1).gt.0.0_r8) THEN
              scale=1.0_r8/(e_u(i,j,k,1)*e_u(i,j,k,1))
            ELSE
              scale=1.0_r8
            END IF
#  else
            scale=cff*(Hz(i-1,j,k)+Hz(i,j,k))
#  endif
            ad_u(i,j,k,Nout)=scale*tl_u(i,j,k,Ninp)
#   ifdef MASKING
            ad_u(i,j,k,Nout)=ad_u(i,j,k,Nout)*umask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  3D v-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=J_V_RANGE
          DO i=I_R_RANGE
#  ifdef BNORM
            IF (e_v(i,j,k,1).gt.0.0_r8) THEN
              scale=1.0_r8/(e_v(i,j,k,1)*e_v(i,j,k,1))
            ELSE
              scale=1.0_r8
            END IF
#  else
            scale=cff*(Hz(i,j-1,k)+Hz(i,j,k))
#  endif
            ad_v(i,j,k,Nout)=scale*tl_v(i,j,k,Ninp)
#   ifdef MASKING
            ad_v(i,j,k,Nout)=ad_v(i,j,k,Nout)*vmask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  Tracers. For now, use salinity scale for passive tracers.
!
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
        DO k=1,N(ng)
          DO j=J_R_RANGE
            DO i=I_R_RANGE
#  ifdef BNORM
              IF (e_t(i,j,k,1,itrc).gt.0.0_r8) THEN
                scale=1.0_r8/(e_t(i,j,k,1,itrc)*e_t(i,j,k,1,itrc))
              ELSE
                scale=1.0_r8
              END IF
#  else
              scale=cff*Hz(i,j,k)
#  endif
              ad_t(i,j,k,Nout,itrc)=scale*tl_t(i,j,k,Ninp,itrc)
#   ifdef MASKING
              ad_t(i,j,k,Nout,itrc)=ad_t(i,j,k,Nout,itrc)*rmask(i,j)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

#  undef I_R_RANGE
#  undef I_U_RANGE
#  undef J_R_RANGE
#  undef J_V_RANGE

      RETURN
      END SUBROUTINE ad_ini_perturb_tile
# endif

# if defined TLM_CHECK

      SUBROUTINE tl_ini_perturb (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  Initialize tangent linear state variable according to the outer     !
!  and inner loop iterations.  The initial field is  a function of     !
!  the steepest descent direction (grad(J)) times the perturbation     !
!  amplitude "p" (controlled by inner loop).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
#  if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#  endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL tl_ini_perturb_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp(ng), nnew(ng), Linp, Lout,         &
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                          SEDBED(ng) % tl_bed_thick,              &
#   endif
     &                          GRID(ng) % tl_Hz,                       &
     &                          GRID(ng) % h,                           &
     &                          GRID(ng) % tl_h,                        &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % on_u,                        &
#   ifdef ICESHELF
     &                          GRID(ng) % zice,                        &
#   endif
     &                          GRID(ng) % tl_z_r,                      &
     &                          GRID(ng) % tl_z_w,                      &
     &                          COUPLING(ng) % Zt_avg1,                 &
     &                          COUPLING(ng) % tl_Zt_avg1,              &
#  endif
     &                          OCEAN(ng) % ubar,                       &
     &                          OCEAN(ng) % vbar,                       &
     &                          OCEAN(ng) % zeta,                       &
#  ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#  endif
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
     &                          OCEAN(ng) % ad_zeta,                    &
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
     &                          OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif

      RETURN
      END SUBROUTINE tl_ini_perturb
!
!***********************************************************************
      SUBROUTINE tl_ini_perturb_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                nstp, nnew, Linp, Lout,           &
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                                tl_bed_thick,                     &
#   endif
     &                                tl_Hz, h, tl_h,                   &
     &                                om_v, on_u,                       &
#   ifdef ICESHELF
     &                                zice,                             &
#   endif
     &                                tl_z_r, tl_z_w,                   &
     &                                Zt_avg1, tl_Zt_avg1,              &
#  endif
     &                                ubar, vbar, zeta,                 &
#  ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
     &                                tl_t, tl_u, tl_v,                 &
#  endif
     &                                ad_ubar, ad_vbar, ad_zeta,        &
     &                                tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars
!
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#   ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#   endif
#  endif
#  ifdef SOLVE3D
      USE tl_set_depth_mod, ONLY : tl_set_depth_tile
#  endif
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
#  ifdef SOLVE3D
      USE tl_t3dbc_mod, ONLY : tl_t3dbc_tile
      USE tl_u3dbc_mod, ONLY : tl_u3dbc_tile
      USE tl_v3dbc_mod, ONLY : tl_v3dbc_tile
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, Linp, Lout
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed_thick(LBi:,LBj:,:)
#    endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(out) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_w(LBi:,LBj:,0:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed_thick(LBi:UBi,LBj:UBj,2)
#    endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
       real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(out) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: i, ic, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      integer, dimension(NstateVar(ng)+1) :: StateVar

      real(r8) :: p, scale
!
#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize tangent linear with the steepest descent direction times
!  the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
!  Set state variable to perturb according to outer loop index.
!
#  ifdef SOLVE3D
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
      StateVar(5)=isUvel
      StateVar(6)=isVvel
      DO i=1,NT(ng)
        StateVar(6+i)=isTvar(i)
      END DO
#  else
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
#  endif
!
!  Set perturbation amplitude as a function of the inner loop.
!
      p=10.0_r8**REAL(-inner,r8)
      scale=1.0_r8/SQRT(adDotProduct)
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (Master) WRITE (stdout,10)
      END IF
!
!  Free-surface.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isFsur)) THEN
        DO j=JstrB,JendB
          DO i=IstrB,IendB
            tl_zeta(i,j,Lout)=p*ad_zeta(i,j,Linp)*scale
#  ifdef MASKING
            tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
#  endif
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isFsur)))
          END IF
        END IF
      END IF
!
!  2D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUbar)) THEN
        DO j=JstrB,JendB
          DO i=IstrM,IendB
            tl_ubar(i,j,Lout)=p*ad_ubar(i,j,Linp)*scale
#  ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#  endif
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUbar)))
          END IF
        END IF
      END IF
!
!  2D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVbar)) THEN
        DO j=JstrM,JendB
          DO i=IstrB,IendB
            tl_vbar(i,j,Lout)=p*ad_vbar(i,j,Linp)*scale
#  ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVbar)))
          END IF
        END IF
      END IF
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUvel)) THEN
        DO k=1,N(ng)
          DO j=JstrB,JendB
            DO i=IstrM,IendB
              tl_u(i,j,k,Lout)=p*ad_u(i,j,k,Linp)*scale
#   ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUvel)))
          END IF
        END IF
      END IF
!
!  3D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVvel)) THEN
        DO k=1,N(ng)
          DO j=JstrM,JendB
            DO i=IstrB,IendB
              tl_v(i,j,k,Lout)=p*ad_v(i,j,k,Linp)*scale
#   ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVvel)))
          END IF
        END IF
      END IF
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        IF ((StateVar(outer).eq.0).or.                                  &
     &      (StateVar(outer).eq.isTvar(itrc))) THEN
          DO k=1,N(ng)
            DO j=JstrB,JendB
              DO i=IstrB,IendB
                tl_t(i,j,k,Lout,itrc)=p*ad_t(i,j,k,Linp,itrc)*scale
#   ifdef MASKING
                tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) outer, inner,                           &
     &                          TRIM(Vname(1,idSvar(isTvar(itrc))))
            END IF
          END IF
        END IF
      END DO
#  endif
      IF (Master) WRITE (stdout,'(/)')
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
      CALL tl_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lout, Lout, Lout,                            &
     &                     zeta, tl_zeta)
#  ifndef SOLVE3D
      CALL tl_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lout, Lout, Lout,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lout, Lout, Lout,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
#  endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,Lout))
#  ifndef SOLVE3D
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_ubar(:,:,Lout))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_vbar(:,:,Lout))
#  endif
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,Lout))
#   ifndef SOLVE3D
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ubar(:,:,Lout),                            &
     &                    tl_vbar(:,:,Lout))
#   endif
#  endif
#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL tl_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp, nnew,                               &
     &                        h, tl_h,                                  &
#   ifdef ICESHELF
     &                        zice,                                     &
#   endif
#   if defined SEDIMENT && defined SED_MORPH
     &                        tl_bed_thick,                             &
#   endif
     &                        Zt_avg1, tl_Zt_avg1,                      &
     &                        tl_Hz, tl_z_r, tl_z_w)
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to perturbed 3D fields.
!-----------------------------------------------------------------------
!
      CALL tl_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lout, Lout, tl_u)
      CALL tl_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lout, Lout, tl_v)
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        CALL tl_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Lout, Lout, tl_t)
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_u(:,:,:,Lout))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_v(:,:,:,Lout))
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,Lout,itrc))
        END DO
      END IF

#   ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,Lout),                             &
     &                    tl_v(:,:,:,Lout))
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,Lout,:))
#   endif
#  endif
!
 10   FORMAT (/,' Perturbing Tangent Linear Initial Conditions:',/)
 20   FORMAT (' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,             &
     &        'Perturbing State Variable: ',a)

      RETURN
      END SUBROUTINE tl_ini_perturb_tile
# endif
#endif
      END MODULE ini_adjust_mod
