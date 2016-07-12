#include "cppdefs.h"

      MODULE posterior_var_mod

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F)

!
!svn $Id: posterior_var.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the full (diagonal) posterior analysis error  !
!  covariance matrix.                                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: posterior_var

      CONTAINS
!
!***********************************************************************
      SUBROUTINE posterior_var (ng, tile, model, outLoop)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, outLoop
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
!
      CALL posterior_var_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Lold(ng), Lnew(ng), outLoop,             &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         BOUNDARY(ng) % e_t_obc,                  &
     &                         BOUNDARY(ng) % e_u_obc,                  &
     &                         BOUNDARY(ng) % e_v_obc,                  &
#  endif
     &                         BOUNDARY(ng) % e_ubar_obc,               &
     &                         BOUNDARY(ng) % e_vbar_obc,               &
     &                         BOUNDARY(ng) % e_zeta_obc,               &
# endif
# ifdef ADJUST_WSTRESS
     &                         FORCES(ng) % e_sustr,                    &
     &                         FORCES(ng) % e_svstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                        FORCES(ng) % e_stflx,                     &
# endif
# ifdef SOLVE3D
     &                        OCEAN(ng) % e_t,                          &
     &                        OCEAN(ng) % e_u,                          &
     &                        OCEAN(ng) % e_v,                          &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                        OCEAN(ng) % e_ubar,                       &
     &                        OCEAN(ng) % e_vbar,                       &
#  endif
# else
     &                        OCEAN(ng) % e_ubar,                       &
     &                        OCEAN(ng) % e_vbar,                       &
# endif
     &                        OCEAN(ng) % e_zeta,                       &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        BOUNDARY(ng) % t_obc,                     &
     &                        BOUNDARY(ng) % u_obc,                     &
     &                        BOUNDARY(ng) % v_obc,                     &
#  endif
     &                        BOUNDARY(ng) % ubar_obc,                  &
     &                        BOUNDARY(ng) % vbar_obc,                  &
     &                        BOUNDARY(ng) % zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                        FORCES(ng) % ustr,                        &
     &                        FORCES(ng) % vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        FORCES(ng) % tflux,                       &
#  endif
     &                        OCEAN(ng) % t,                            &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                        OCEAN(ng) % ubar,                         &
     &                        OCEAN(ng) % vbar,                         &
#  endif
# else
     &                        OCEAN(ng) % ubar,                         &
     &                        OCEAN(ng) % vbar,                         &
# endif
     &                        OCEAN(ng) % zeta,                         &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        BOUNDARY(ng) % tl_t_obc,                  &
     &                        BOUNDARY(ng) % tl_u_obc,                  &
     &                        BOUNDARY(ng) % tl_v_obc,                  &
#  endif
     &                        BOUNDARY(ng) % tl_ubar_obc,               &
     &                        BOUNDARY(ng) % tl_vbar_obc,               &
     &                        BOUNDARY(ng) % tl_zeta_obc,               &
# endif
# ifdef ADJUST_WSTRESS
     &                        FORCES(ng) % tl_ustr,                     &
     &                        FORCES(ng) % tl_vstr,                     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        FORCES(ng) % tl_tflux,                    &
#  endif
     &                        OCEAN(ng) % tl_t,                         &
     &                        OCEAN(ng) % tl_u,                         &
     &                        OCEAN(ng) % tl_v,                         &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                        OCEAN(ng) % tl_ubar,                      &
     &                        OCEAN(ng) % tl_vbar,                      &
#  endif
# else
     &                        OCEAN(ng) % tl_ubar,                      &
     &                        OCEAN(ng) % tl_vbar,                      &
# endif
     &                        OCEAN(ng) % tl_zeta,                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        BOUNDARY(ng) % ad_t_obc,                  &
     &                        BOUNDARY(ng) % ad_u_obc,                  &
     &                        BOUNDARY(ng) % ad_v_obc,                  &
#  endif
     &                        BOUNDARY(ng) % ad_ubar_obc,               &
     &                        BOUNDARY(ng) % ad_vbar_obc,               &
     &                        BOUNDARY(ng) % ad_zeta_obc,               &
# endif
# ifdef ADJUST_WSTRESS
     &                        FORCES(ng) % ad_ustr,                     &
     &                        FORCES(ng) % ad_vstr,                     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        FORCES(ng) % ad_tflux,                    &
#  endif
     &                        OCEAN(ng) % ad_t,                         &
     &                        OCEAN(ng) % ad_u,                         &
     &                        OCEAN(ng) % ad_v,                         &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                        OCEAN(ng) % ad_ubar,                      &
     &                        OCEAN(ng) % ad_vbar,                      &
#  endif
# else
     &                        OCEAN(ng) % ad_ubar,                      &
     &                        OCEAN(ng) % ad_vbar,                      &
# endif
     &                        OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE posterior_var
!
!***********************************************************************
      SUBROUTINE posterior_var_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Lold, Lnew, outLoop,               &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                               t_obc_std, u_obc_std, v_obc_std,   &
#  endif
     &                               ubar_obc_std, vbar_obc_std,        &
     &                               zeta_obc_std,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                               sustr_std, svstr_std,              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                               stflx_std,                         &
# endif
# ifdef SOLVE3D
     &                               t_std, u_std, v_std,               &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                               ubar_std, vbar_std,                &
#  endif
# else
     &                               ubar_std, vbar_std,                &
# endif
     &                               zeta_std,                          &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                               nl_t_obc, nl_u_obc, nl_v_obc,      &
#  endif
     &                               nl_ubar_obc, nl_vbar_obc,          &
     &                               nl_zeta_obc,                       &
# endif
# ifdef ADJUST_WSTRESS
     &                               nl_ustr, nl_vstr,                  &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                               nl_tflux,                          &
#  endif
     &                               nl_t, nl_u, nl_v,                  &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                               nl_ubar, nl_vbar,                  &
#  endif
# else
     &                               nl_ubar, nl_vbar,                  &
# endif
     &                               nl_zeta,                           &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                               tl_t_obc, tl_u_obc, tl_v_obc,      &
#  endif
     &                               tl_ubar_obc, tl_vbar_obc,          &
     &                               tl_zeta_obc,                       &
# endif
# ifdef ADJUST_WSTRESS
     &                               tl_ustr, tl_vstr,                  &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                               tl_tflux,                          &
#  endif
     &                               tl_t, tl_u, tl_v,                  &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                               tl_ubar, tl_vbar,                  &
#  endif
# else
     &                               tl_ubar, tl_vbar,                  &
# endif
     &                               tl_zeta,                           &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                               ad_t_obc, ad_u_obc, ad_v_obc,      &
#  endif
     &                               ad_ubar_obc, ad_vbar_obc,          &
     &                               ad_zeta_obc,                       &
# endif
# ifdef ADJUST_WSTRESS
     &                               ad_ustr, ad_vstr,                  &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                               ad_tflux,                          &
#  endif
     &                               ad_t, ad_u, ad_v,                  &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                               ad_ubar, ad_vbar,                  &
#  endif
# else
     &                               ad_ubar, ad_vbar,                  &
# endif
     &                               ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
# endif
      USE state_copy_mod, ONLY : state_copy
      USE state_product_mod, ONLY : state_product
      USE state_addition_mod, ONLY : state_addition
      USE state_initialize_mod, ONLY : state_initialize
      USE posterior_mod, ONLY : read_state
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#   ifdef ADJUST_BOUNDARY
#    ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:,:,:,:)
      real(r8), intent(in) :: u_obc_std(LBij:,:,:)
      real(r8), intent(in) :: v_obc_std(LBij:,:,:)
#    endif
      real(r8), intent(in) :: ubar_obc_std(LBij:,:)
      real(r8), intent(in) :: vbar_obc_std(LBij:,:)
      real(r8), intent(in) :: zeta_obc_std(LBij:,:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:,LBj:,:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u_std(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v_std(LBi:,LBj:,:,:)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: ubar_std(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:,:)
#    endif
#   else
      real(r8), intent(in) :: ubar_std(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: zeta_std(LBi:,LBj:,:)
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
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: nl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: nl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
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
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(in) :: u_obc_std(LBij:UBij,N(ng),4)
      real(r8), intent(in) :: v_obc_std(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: vbar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: zeta_obc_std(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: u_std(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(in) :: v_std(LBi:UBi,LBj:UBj,N(ng),NSA)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj,NSA)
#   endif
#  else
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj,NSA)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:UBi,LBj:UBj,NSA)
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
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: nl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: nl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: nl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
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
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: L1 = 1
      integer :: L2 = 2

      integer :: Linp, Lout, Lscale, Lwrk, Lwrk1, i, j, ic
      integer :: info, ivec, jvec, rec

      real(r8) :: fac, fac1, fac2

      logical :: Lweak

      character (len=256) :: ncname

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Computes the full (diagonal) posterior analysis error covariance
# if defined POSTERIOR_ERROR_I
!  matrix at the initial time of the assimilation window.
# elif defined POSTERIOR_ERROR_F
!  matrix at the final time of the assimilation window.
# endif
!-----------------------------------------------------------------------
!
!  NOTE: In the case of weak constraint ("WEAK_CONSTRAINT") and
!        and time convolutions ("TIME_CONV"), the state arrays
!        tl_ubar, tl_vbar, ad_ubar, and ad_vbar are only passed
!        as required by the "state" operators routines but they
!        are not used in subsequent calculations.
!
      IF (Master) WRITE (stdout,10)
 10   FORMAT (/,' <<<< Full Posterior Error Covariance Matrix >>>>',/)
!
      Lweak=.FALSE.
!
!  Invert the tridiagonal matrix of inner-loop Lanczos vector
!  coefficients.
!
      zLanczos_coef=0.0_r8
      DO i=1,Ninner
        zLanczos_coef(i,i)=cg_delta(i,outLoop)
      END DO
      DO i=1,Ninner-1
        zLanczos_coef(i,i+1)=cg_beta(i+1,outLoop)
      END DO
      DO i=2,Ninner
        zLanczos_coef(i,i-1)=cg_beta(i,outLoop)
      END DO
      zLanczos_inv=zLanczos_coef
!
!  Invert "zLanczos_coef" using LAPACK routines DPOTRF and DPOTRI.
!  On exit "zLanczos_inv" contains the inverse of "zLanczos_coef".
!  The "zLanczos_coef" matrix is saved to check the inversion below.
!
      IF (Master) THEN
        CALL dpotrf ('U', Ninner, zLanczos_inv, Ninner, info)
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, info)
# endif
      IF (info.ne.0) THEN
        IF (Master) PRINT *,' Error in DPOTRF: info = ', info
        exit_flag=8
        RETURN
      END IF
!
      IF (Master) THEN
        CALL dpotri ('U', Ninner, zLanczos_inv, Ninner, info)
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, info)
      CALL mp_bcastf (ng, model, zLanczos_inv)
# endif
      IF (info.ne.0) THEN
        IF (Master) PRINT *,' Error in DPOTRI: info =', info
        exit_flag=8
        RETURN
      END IF
!
!  Copy the upper triangle of the inverse into the lower triangle
!  as required by the solution from DPOTRI.
!
      DO j=1,Ninner
        DO i=j+1,Ninner
          zLanczos_inv(i,j)=zLanczos_inv(j,i)
        END DO
      END DO
!
!  Test the inverse, we need to get the identity matrix within roundoff.
!
      DO j=1,Ninner
       DO i=1,Ninner
         zLanczos_err(i,j)=0.0_r8
         DO ic=1,Ninner
            zLanczos_err(i,j)=zLanczos_err(i,j)+                        &
     &                        zLanczos_inv(i,ic)*zLanczos_coef(ic,j)
         END DO
       END DO
      END DO
!
!  Now compute the diagonal of the posterior analysis error covariance
!  matrix. This will be stored in record L1 of the tl_var arrays.
!
      fac=0.0_r8
      CALL state_initialize (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       L1, fac,                                   &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       tl_t_obc, tl_u_obc, tl_v_obc,              &
#  endif
     &                       tl_ubar_obc, tl_vbar_obc,                  &
     &                       tl_zeta_obc,                               &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       tl_tflux,                                  &
#  endif
     &                       tl_t, tl_u, tl_v,                          &
# else
     &                       tl_ubar, tl_vbar,                          &
# endif
     &                       tl_zeta)
!
!  Read in turn each evolved and convolved Lanczos vector of the
!  stabilized representer matrix and compute the product with all
!  the other vectors. Use the ad_var arrays as temporary storage.
!
      ncname=HSS(ng)%name
      DO ivec=1,Ninner
        rec=ivec
!
!  Read vectors stored in hessian netcdf file using ad_var(L1) as
!  temporary storage.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L1, rec,                                       &
     &                   0, HSS(ng)%ncid, ncname,                       &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
#  endif
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
!
!  Compute state product of ad_var(L1) with itself using nl_var(L1)
!  as temporary storage.
!
        CALL state_product (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,L1,:),                     &
     &                      ad_t_obc(:,:,:,:,L1,:),                     &
     &                      nl_t_obc(:,:,:,:,L1,:),                     &
     &                      ad_u_obc(:,:,:,:,L1),                       &
     &                      ad_u_obc(:,:,:,:,L1),                       &
     &                      nl_u_obc(:,:,:,:,L1),                       &
     &                      ad_v_obc(:,:,:,:,L1),                       &
     &                      ad_v_obc(:,:,:,:,L1),                       &
     &                      nl_v_obc(:,:,:,:,L1),                       &
#  endif
     &                      ad_ubar_obc(:,:,:,L1),                      &
     &                      ad_ubar_obc(:,:,:,L1),                      &
     &                      nl_ubar_obc(:,:,:,L1),                      &
     &                      ad_vbar_obc(:,:,:,L1),                      &
     &                      ad_vbar_obc(:,:,:,L1),                      &
     &                      nl_vbar_obc(:,:,:,L1),                      &
     &                      ad_zeta_obc(:,:,:,L1),                      &
     &                      ad_zeta_obc(:,:,:,L1),                      &
     &                      nl_zeta_obc(:,:,:,L1),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,L1),                          &
     &                      ad_ustr(:,:,:,L1),                          &
     &                      nl_ustr(:,:,:,L1),                          &
     &                      ad_vstr(:,:,:,L1),                          &
     &                      ad_vstr(:,:,:,L1),                          &
     &                      nl_vstr(:,:,:,L1),                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,L1,:),                       &
     &                      ad_tflux(:,:,:,L1,:),                       &
     &                      nl_tflux(:,:,:,L1,:),                       &
#  endif
     &                      ad_t(:,:,:,L1,:),                           &
     &                      ad_t(:,:,:,L1,:),                           &
     &                      nl_t(:,:,:,L1,:),                           &
     &                      ad_u(:,:,:,L1),                             &
     &                      ad_u(:,:,:,L1),                             &
     &                      nl_u(:,:,:,L1),                             &
     &                      ad_v(:,:,:,L1),                             &
     &                      ad_v(:,:,:,L1),                             &
     &                      nl_v(:,:,:,L1),                             &
# else
     &                      ad_ubar(:,:,L1),                            &
     &                      ad_ubar(:,:,L1),                            &
     &                      nl_ubar(:,:,L1),                            &
     &                      ad_vbar(:,:,L1),                            &
     &                      ad_vbar(:,:,L1),                            &
     &                      nl_vbar(:,:,L1),                            &
# endif
     &                      ad_zeta(:,:,L1),                            &
     &                      ad_zeta(:,:,L1),                            &
     &                      nl_zeta(:,:,L1))
!
!  tl_var(L1) = fac1 * tl_var(L1) + fac2 * nl_var(L1).  See NOTE above.
!
        fac1=1.0_r8
        fac2=zLanczos_inv(ivec,ivec)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       L1, L1, L1, fac1, fac2,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       tl_t_obc, nl_t_obc,                        &
     &                       tl_u_obc, nl_u_obc,                        &
     &                       tl_v_obc, nl_v_obc,                        &
#  endif
     &                       tl_ubar_obc, nl_ubar_obc,                  &
     &                       tl_vbar_obc, nl_vbar_obc,                  &
     &                       tl_zeta_obc, nl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, nl_ustr,                          &
     &                       tl_vstr, nl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       tl_tflux, nl_tflux,                        &
#  endif
     &                       tl_t, nl_t,                                &
     &                       tl_u, nl_u,                                &
     &                       tl_v, nl_v,                                &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       tl_ubar, nl_ubar,                          &
     &                       tl_vbar, nl_vbar,                          &
#  endif
# else
     &                       tl_ubar, nl_ubar,                          &
     &                       tl_vbar, nl_vbar,                          &
# endif
     &                       tl_zeta, nl_zeta)
!
        DO jvec=ivec+1,Ninner
          rec=jvec
!
!  Read vectors stored in hessian netcdf file using ad_var(L2) as
!  temporary storage.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     L2, rec,                                     &
     &                     0, HSS(ng)%ncid, ncname,                     &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     ad_t_obc, ad_u_obc, ad_v_obc,                &
#  endif
     &                     ad_ubar_obc, ad_vbar_obc,                    &
     &                     ad_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     ad_ustr, ad_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     ad_tflux,                                    &
#  endif
     &                     ad_t, ad_u, ad_v,                            &
# else
     &                     ad_ubar, ad_vbar,                            &
# endif
     &                     ad_zeta)
!
!  Compute state product of ad_var(L2) with ad_var(L1) using nl_var(L1)
!  as temporary storage.
!
          CALL state_product (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        ad_t_obc(:,:,:,:,L1,:),                   &
     &                        ad_t_obc(:,:,:,:,L2,:),                   &
     &                        nl_t_obc(:,:,:,:,L1,:),                   &
     &                        ad_u_obc(:,:,:,:,L1),                     &
     &                        ad_u_obc(:,:,:,:,L2),                     &
     &                        nl_u_obc(:,:,:,:,L1),                     &
     &                        ad_v_obc(:,:,:,:,L1),                     &
     &                        ad_v_obc(:,:,:,:,L2),                     &
     &                        nl_v_obc(:,:,:,:,L1),                     &
#  endif
     &                        ad_ubar_obc(:,:,:,L1),                    &
     &                        ad_ubar_obc(:,:,:,L2),                    &
     &                        nl_ubar_obc(:,:,:,L1),                    &
     &                        ad_vbar_obc(:,:,:,L1),                    &
     &                        ad_vbar_obc(:,:,:,L2),                    &
     &                        nl_vbar_obc(:,:,:,L1),                    &
     &                        ad_zeta_obc(:,:,:,L1),                    &
     &                        ad_zeta_obc(:,:,:,L2),                    &
     &                        nl_zeta_obc(:,:,:,L1),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                        ad_ustr(:,:,:,L1),                        &
     &                        ad_ustr(:,:,:,L2),                        &
     &                        nl_ustr(:,:,:,L1),                        &
     &                        ad_vstr(:,:,:,L1),                        &
     &                        ad_vstr(:,:,:,L2),                        &
     &                        nl_vstr(:,:,:,L1),                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        ad_tflux(:,:,:,L1,:),                     &
     &                        ad_tflux(:,:,:,L2,:),                     &
     &                        nl_tflux(:,:,:,L1,:),                     &
#  endif
     &                        ad_t(:,:,:,L1,:),                         &
     &                        ad_t(:,:,:,L2,:),                         &
     &                        nl_t(:,:,:,L1,:),                         &
     &                        ad_u(:,:,:,L1),                           &
     &                        ad_u(:,:,:,L2),                           &
     &                        nl_u(:,:,:,L1),                           &
     &                        ad_v(:,:,:,L1),                           &
     &                        ad_v(:,:,:,L2),                           &
     &                        nl_v(:,:,:,L1),                           &
# else
     &                        ad_ubar(:,:,L1),                          &
     &                        ad_ubar(:,:,L2),                          &
     &                        nl_ubar(:,:,L1),                          &
     &                        ad_vbar(:,:,L1),                          &
     &                        ad_vbar(:,:,L2),                          &
     &                        nl_vbar(:,:,L1),                          &
# endif
     &                        ad_zeta(:,:,L1),                          &
     &                        ad_zeta(:,:,L2),                          &
     &                        nl_zeta(:,:,L1))
!
!  tl_var(L1) = fac1 * tl_var(L1) + fac2 * nl_var(L1). See NOTE above.
!
          fac1=1.0_r8
          fac2=2.0_r8*zLanczos_inv(ivec,jvec)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         L1, L1, L1, fac1, fac2,                  &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         tl_t_obc, nl_t_obc,                      &
     &                         tl_u_obc, nl_u_obc,                      &
     &                         tl_v_obc, nl_v_obc,                      &
#  endif
     &                         tl_ubar_obc, nl_ubar_obc,                &
     &                         tl_vbar_obc, nl_vbar_obc,                &
     &                         tl_zeta_obc, nl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         tl_ustr, nl_ustr,                        &
     &                         tl_vstr, nl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         tl_tflux, nl_tflux,                      &
#  endif
     &                         tl_t, nl_t,                              &
     &                         tl_u, nl_u,                              &
     &                         tl_v, nl_v,                              &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                         tl_ubar, nl_ubar,                        &
     &                         tl_vbar, nl_vbar,                        &
#  endif
# else
     &                         tl_ubar, nl_ubar,                        &
     &                         tl_vbar, nl_vbar,                        &
# endif
     &                         tl_zeta, nl_zeta)
        END DO
      END DO
!
!  Finally subtract the result in tl_var(L1) from the background
!  error variances to get the posterior/analysis error variance.
!
      CALL analysis_var (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L1, Lweak,                                     &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   t_obc_std, u_obc_std, v_obc_std,               &
#  endif
     &                   ubar_obc_std, vbar_obc_std, zeta_obc_std,      &
# endif
# ifdef ADJUST_WSTRESS
     &                   sustr_std, svstr_std,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   stflx_std,                                     &
#  endif
     &                   t_std, u_std, v_std,                           &
# else
     &                   ubar_std, vbar_std,                            &
# endif
     &                   zeta_std,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc, tl_zeta_obc,         &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
#  endif
     &                   tl_t, tl_u, tl_v,                              &
# else
     &                   tl_ubar, tl_vbar,                              &
# endif
     &                   tl_zeta)

      RETURN
      END SUBROUTINE posterior_var_tile
!
!***********************************************************************
      SUBROUTINE analysis_var (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lout, Lweak,                             &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         t_obc_std, u_obc_std, v_obc_std,         &
#  endif
     &                         ubar_obc_std, vbar_obc_std,              &
     &                         zeta_obc_std,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                         sustr_std, svstr_std,                    &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         stflx_std,                               &
#  endif
     &                         t_std, u_std, v_std,                     &
# else
     &                         ubar_std, vbar_std,                      &
# endif
     &                         zeta_std,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         s1_t_obc, s1_u_obc, s1_v_obc,            &
#  endif
     &                         s1_ubar_obc, s1_vbar_obc,                &
     &                         s1_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         s1_sustr, s1_svstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         s1_tflux,                                &
#  endif
     &                         s1_t, s1_u, s1_v,                        &
# else
     &                         s1_ubar, s1_vbar,                        &
# endif
     &                         s1_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
# if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
     defined ADJUST_WSTRESS
      USE mod_scalars
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Lout
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      logical, intent(in) :: Lweak
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s1_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s1_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s1_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: s1_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s1_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s1_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_svstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: s1_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s1_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: s1_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s1_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: s1_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:,:,:,:)
      real(r8), intent(in) :: u_obc_std(LBij:,:,:)
      real(r8), intent(in) :: v_obc_std(LBij:,:,:)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:,:)
      real(r8), intent(in) :: vbar_obc_std(LBij:,:)
      real(r8), intent(in) :: zeta_obc_std(LBij:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u_std(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v_std(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ubar_std(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif

#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s1_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s1_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s1_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: s1_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s1_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s1_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s1_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: s1_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s1_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s1_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: s1_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s1_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: s1_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(in) :: u_obc_std(LBij:UBij,N(ng),4)
      real(r8), intent(in) :: v_obc_std(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: vbar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: zeta_obc_std(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: u_std(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(in) :: v_std(LBi:UBi,LBj:UBj,N(ng),NSA)
#  else
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj,NSA)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:UBi,LBj:UBj,NSA)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      integer :: ir, it, rec

      real(r8) :: cff

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine error covariance standard deviation factors to use.
!-----------------------------------------------------------------------
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint
      ELSE
        rec=1                        ! strong constraint
      END IF
!
!-----------------------------------------------------------------------
!  Compute S1=STD*STD-S1.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff=zeta_std(i,j,rec)*zeta_std(i,j,rec)-                      &
     &        s1_zeta(i,j,Lout)
# ifdef MASKING
          cff=cff*rmask(i,j)
# endif
          s1_zeta(i,j,Lout)=cff
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
            DO j=Jstr,Jend
              cff=zeta_obc_std(j,iwest)*zeta_obc_std(j,iwest)-          &
     &            s1_zeta_obc(j,iwest,ir,Lout)
#  ifdef MASKING
              cff=cff*rmask(Istr-1,j)
#  endif
              s1_zeta_obc(j,iwest,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              cff=zeta_obc_std(j,ieast)*zeta_obc_std(j,ieast)-          &
     &            s1_zeta_obc(j,ieast,ir,Lout)
#  ifdef MASKING
              cff=cff*rmask(Iend+1,j)
#  endif
              s1_zeta_obc(j,ieast,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=zeta_obc_std(i,isouth)*zeta_obc_std(i,isouth)-        &
     &            s1_zeta_obc(i,isouth,ir,Lout)
#  ifdef MASKING
              cff=cff*rmask(i,Jstr-1)
#  endif
              s1_zeta_obc(i,isouth,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=zeta_obc_std(i,inorth)*zeta_obc_std(i,inorth)-        &
     &            s1_zeta_obc(i,inorth,ir,Lout)
#  ifdef MASKING
              cff=cff*rmask(i,Jend+1)
#  endif
              s1_zeta_obc(i,inorth,ir,Lout)=cff
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum component.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          cff=ubar_std(i,j,rec)*ubar_std(i,j,rec)-                      &
     &        s1_ubar(i,j,Lout)
#  ifdef MASKING
          cff=cff*umask(i,j)
#  endif
          s1_ubar(i,j,Lout)=cff
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
            DO j=Jstr,Jend
              cff=ubar_obc_std(j,iwest)*ubar_obc_std(j,iwest)-          &
     &            s1_ubar_obc(j,iwest,ir,Lout)
#  ifdef MASKING
              cff=cff*umask(Istr,j)
#  endif
              s1_ubar_obc(j,iwest,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              cff=ubar_obc_std(j,ieast)*ubar_obc_std(j,ieast)-          &
     &            s1_ubar_obc(j,ieast,ir,Lout)
#  ifdef MASKING
              cff=cff*umask(Iend+1,j)
#  endif
              s1_ubar_obc(j,ieast,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
              cff=ubar_obc_std(i,isouth)*ubar_obc_std(i,isouth)-        &
     &            s1_ubar_obc(i,isouth,ir,Lout)
#  ifdef MASKING
              cff=cff*umask(i,Jstr-1)
#  endif
              s1_ubar_obc(i,isouth,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
              cff=ubar_obc_std(i,inorth)*ubar_obc_std(i,inorth)-        &
     &            s1_ubar_obc(i,inorth,ir,Lout)
#  ifdef MASKING
              cff=cff*umask(i,Jend+1)
#  endif
              s1_ubar_obc(i,inorth,ir,Lout)=cff
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum component.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          cff=vbar_std(i,j,rec)*vbar_std(i,j,rec)-                      &
     &        s1_vbar(i,j,Lout)
#  ifdef MASKING
          cff=cff*vmask(i,j)
#  endif
          s1_vbar(i,j,Lout)=cff
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
            DO j=JstrV,Jend
              cff=vbar_obc_std(j,iwest)*vbar_obc_std(j,iwest)-          &
     &            s1_vbar_obc(j,iwest,ir,Lout)
#  ifdef MASKING
              cff=cff*vmask(Istr-1,j)
#  endif
              s1_vbar_obc(j,iwest,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
              cff=vbar_obc_std(j,ieast)*vbar_obc_std(j,ieast)-          &
     &            s1_vbar_obc(j,ieast,ir,Lout)
#  ifdef MASKING
              cff=cff*vmask(Iend+1,j)
#  endif
              s1_vbar_obc(j,ieast,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=vbar_obc_std(i,isouth)*vbar_obc_std(i,isouth)-        &
     &            s1_vbar_obc(i,isouth,ir,Lout)
#  ifdef MASKING
              cff=cff*vmask(i,Jstr)
#  endif
              s1_vbar_obc(i,isouth,ir,Lout)=cff
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              cff=vbar_obc_std(i,inorth)*vbar_obc_std(i,inorth)-        &
     &            s1_vbar_obc(i,inorth,ir,Lout)
#  ifdef MASKING
              cff=cff*vmask(i,Jend+1)
#  endif
              s1_vbar_obc(i,inorth,ir,Lout)=cff
            END DO
          END IF
        END DO
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO ir=1,Nfrec(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            cff=sustr_std(i,j)*sustr_std(i,j)-                          &
     &          s1_sustr(i,j,ir,Lout)
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
            s1_sustr(i,j,ir,Lout)=cff
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            cff=svstr_std(i,j)*svstr_std(i,j)-                          &
     &          s1_svstr(i,j,ir,Lout)
#  ifdef MASKING
            cff=cff*vmask(i,j)
#  endif
            s1_svstr(i,j,ir,Lout)=cff
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            cff=u_std(i,j,k,rec)*u_std(i,j,k,rec)-                      &
     &          s1_u(i,j,k,Lout)
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
            s1_u(i,j,k,Lout)=cff
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
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=u_obc_std(j,k,iwest)*u_obc_std(j,k,iwest)-          &
     &              s1_u_obc(j,k,iwest,ir,Lout)
#   ifdef MASKING
                cff=cff*umask(Istr,j)
#   endif
                s1_u_obc(j,k,iwest,ir,Lout)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                cff=u_obc_std(j,k,ieast)*u_obc_std(j,k,ieast)-          &
     &              s1_u_obc(j,k,ieast,ir,Lout)
#   ifdef MASKING
                cff=cff*umask(Iend+1,j)
#   endif
                s1_u_obc(j,k,ieast,ir,Lout)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=u_obc_std(i,k,isouth)*u_obc_std(i,k,isouth)-        &
     &              s1_u_obc(i,k,isouth,ir,Lout)
#   ifdef MASKING
                cff=cff*umask(i,Jstr-1)
#   endif
                s1_u_obc(i,k,isouth,ir,Lout)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=IstrU,Iend
                cff=u_obc_std(i,k,inorth)*u_obc_std(i,k,inorth)-        &
     &              s1_u_obc(i,k,inorth,ir,Lout)
#   ifdef MASKING
                cff=cff*umask(i,Jend+1)
#   endif
                s1_u_obc(i,k,inorth,ir,Lout)=cff
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
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            cff=v_std(i,j,k,rec)*v_std(i,j,k,rec)-                      &
     &          s1_v(i,j,k,Lout)
#  ifdef MASKING
            cff=cff*vmask(i,j)
#  endif
            s1_v(i,j,k,Lout)=cff
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
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=v_obc_std(j,k,iwest)*v_obc_std(j,k,iwest)-          &
     &              s1_v_obc(j,k,iwest,ir,Lout)
#   ifdef MASKING
                cff=cff*vmask(Istr-1,j)
#   endif
                s1_v_obc(j,k,iwest,ir,Lout)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO j=JstrV,Jend
                cff=v_obc_std(j,k,ieast)*v_obc_std(j,k,ieast)-          &
     &              s1_v_obc(j,k,ieast,ir,Lout)
#   ifdef MASKING
                cff=cff*vmask(Iend+1,j)
#   endif
                s1_v_obc(j,k,ieast,ir,Lout)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=v_obc_std(i,k,isouth)*v_obc_std(i,k,isouth)-        &
     &              s1_v_obc(i,k,isouth,ir,Lout)
#   ifdef MASKING
                cff=cff*vmask(i,Jstr)
#   endif
                s1_v_obc(i,k,isouth,ir,Lout)=cff
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=v_obc_std(i,k,inorth)*v_obc_std(i,k,inorth)-        &
     &              s1_v_obc(i,k,inorth,ir,Lout)
#   ifdef MASKING
                cff=cff*vmask(i,Jend+1)
#   endif
                s1_v_obc(i,k,inorth,ir,Lout)=cff
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  Tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              cff=t_std(i,j,k,rec,it)*t_std(i,j,k,rec,it)-              &
     &            s1_t(i,j,k,Lout,it)
#  ifdef MASKING
              cff=cff*rmask(i,j)
#  endif
              s1_t(i,j,k,Lout,it)=cff
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_BOUNDARY
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
                  cff=t_obc_std(j,k,iwest,it)*t_obc_std(j,k,iwest,it)-  &
     &                s1_t_obc(j,k,iwest,ir,Lout,it)
#   ifdef MASKING
                  cff=cff*rmask(Istr-1,j)
#   endif
                  s1_t_obc(j,k,iwest,ir,Lout,it)=cff
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  cff=t_obc_std(j,k,ieast,it)*t_obc_std(j,k,ieast,it)-  &
     &                s1_t_obc(j,k,ieast,ir,Lout,it)
#   ifdef MASKING
                  cff=cff*rmask(Iend+1,j)
#   endif
                  s1_t_obc(j,k,ieast,ir,Lout,it)=cff
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff=t_obc_std(i,k,isouth,it)*t_obc_std(i,k,isouth,it)-&
     &                s1_t_obc(i,k,isouth,ir,Lout,it)
#   ifdef MASKING
                  cff=cff*rmask(i,Jstr-1)
#   endif
                  s1_t_obc(i,k,isouth,ir,Lout,it)=cff
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff=t_obc_std(i,k,inorth,it)*t_obc_std(i,k,inorth,it)-&
     &                s1_t_obc(i,k,inorth,ir,Lout,it)
#   ifdef MASKING
                  cff=cff*rmask(i,Jend+1)
#   endif
                  s1_t_obc(i,k,inorth,ir,Lout,it)=cff
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
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          DO ir=1,Nfrec(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                cff=stflx_std(i,j,it)*stflx_std(i,j,it)-                &
     &              s1_tflux(i,j,ir,Lout,it)
#   ifdef MASKING
                cff=cff*rmask(i,j)
#   endif
                s1_tflux(i,j,ir,Lout,it)=cff
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif

# endif

      RETURN
      END SUBROUTINE analysis_var
#endif
      END MODULE posterior_var_mod
