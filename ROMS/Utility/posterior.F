#include "cppdefs.h"
      MODULE posterior_mod

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
    defined POSTERIOR_ERROR_F)
!
!svn $Id: posterior.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module computes the Lanczos vectors and eigenvectors of the    !
!  posterior analysis error covariance matrix.                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: posterior
      PUBLIC :: read_state

      CONTAINS
!
!***********************************************************************
      SUBROUTINE posterior (ng, tile, model, innLoop, outLoop, Ltrace)
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
      integer, intent(in) :: ng, tile, model, innLoop, outLoop
      logical, intent(in) :: Ltrace
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
!
      CALL posterior_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lold(ng), Lnew(ng),                          &
     &                     innLoop, outLoop, Ltrace,                    &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % t_obc,                        &
     &                     BOUNDARY(ng) % u_obc,                        &
     &                     BOUNDARY(ng) % v_obc,                        &
#  endif
     &                     BOUNDARY(ng) % ubar_obc,                     &
     &                     BOUNDARY(ng) % vbar_obc,                     &
     &                     BOUNDARY(ng) % zeta_obc,                     &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ustr,                           &
     &                     FORCES(ng) % vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tflux,                          &
#  endif
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
#  endif
# else
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
# endif
     &                     OCEAN(ng) % zeta,                            &
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#  endif
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta,                         &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     BOUNDARY(ng) % d_t_obc,                      &
     &                     BOUNDARY(ng) % d_u_obc,                      &
     &                     BOUNDARY(ng) % d_v_obc,                      &
#  endif
     &                     BOUNDARY(ng) % d_ubar_obc,                   &
     &                     BOUNDARY(ng) % d_vbar_obc,                   &
     &                     BOUNDARY(ng) % d_zeta_obc,                   &
# endif
# ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % d_sustr,                        &
     &                     FORCES(ng) % d_svstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     FORCES(ng) % d_stflx,                        &
#  endif
     &                     OCEAN(ng) % d_t,                             &
     &                     OCEAN(ng) % d_u,                             &
     &                     OCEAN(ng) % d_v,                             &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                     OCEAN(ng) % d_ubar,                          &
     &                     OCEAN(ng) % d_vbar,                          &
#  endif
# else
     &                     OCEAN(ng) % d_ubar,                          &
     &                     OCEAN(ng) % d_vbar,                          &
# endif
     &                     OCEAN(ng) % d_zeta,                          &
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
#  endif
# else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
# endif
     &                     OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE posterior
!
!***********************************************************************
      SUBROUTINE posterior_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lold, Lnew,                            &
     &                           innLoop, outLoop, Ltrace,              &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           nl_t_obc, nl_u_obc, nl_v_obc,          &
#  endif
     &                           nl_ubar_obc, nl_vbar_obc,              &
     &                           nl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           nl_ustr, nl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           nl_tflux,                              &
#  endif
     &                           nl_t, nl_u, nl_v,                      &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           nl_ubar, nl_vbar,                      &
#  endif
# else
     &                           nl_ubar, nl_vbar,                      &
# endif
     &                           nl_zeta,                               &
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           tl_ubar, tl_vbar,                      &
#  endif
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta,                               &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           d_t_obc, d_u_obc, d_v_obc,             &
#  endif
     &                           d_ubar_obc, d_vbar_obc,                &
     &                           d_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                           d_sustr, d_svstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           d_stflx,                               &
#  endif
     &                           d_t, d_u, d_v,                         &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           d_ubar, d_vbar,                        &
#  endif
# else
     &                           d_ubar, d_vbar,                        &
# endif
     &                           d_zeta,                                &
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           ad_ubar, ad_vbar,                      &
#  endif
# else
     &                           ad_ubar, ad_vbar,                      &
# endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
# endif
      USE state_copy_mod, ONLY : state_copy
# ifdef BEOFS_ONLY
      USE state_dotprod_mod, ONLY : state_dotprod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew
      integer, intent(in) :: innLoop, outLoop
      logical, intent(in) :: Ltrace
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
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
      real(r8), intent(inout) :: d_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: d_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: d_v_obc(LBij:,:,:,:)
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_zeta_obc(LBij:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
#   endif
#  else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
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
      real(r8), intent(inout) :: d_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),NT(ng))
      real(r8), intent(inout) :: d_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(inout) :: d_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_zeta_obc(LBij:UBij,4,Nbrec(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#   endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
#   endif
#  else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
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
      logical :: Ltrans

      integer :: L1 = 1
      integer :: L2 = 2

      integer :: Linp, Lout, Lscale, Lwrk, Lwrk1, i, j, ic
      integer :: info, itheta1

      real(r8) :: norm, zbeta

      real(r8), dimension(2*NpostI-2) :: work
# ifdef BEOFS_ONLY
      integer :: LAMM
      real(r8), dimension(0:NstateVar(ng)) :: dot
# endif

      character (len=13) :: string

# include "set_bounds.h"
# ifdef BEOFS_ONLY
!
!-----------------------------------------------------------------------
!  Compute the action of the background error covariance matrix.
!-----------------------------------------------------------------------
!
!  NOTE: In the case of WEAK_CONSTRAINT and TIME_CONV, tl_ubar, tl_vbar
!        ad_ubar and ad_vbar are only passed as required but are not
!        used in subsequent calculations.
!
!  Copy tl_var(Lnew) into ad_var(Lnew).
!
      IF (innLoop.eq.0) THEN
        LAMM=Lold
      ELSE
        LAMM=Lnew
      END IF

      CALL state_copy (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 LAMM, Lnew,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                 ad_t_obc, tl_t_obc,                              &
     &                 ad_u_obc, tl_u_obc,                              &
     &                 ad_v_obc, tl_v_obc,                              &
#  endif
     &                 ad_ubar_obc, tl_ubar_obc,                        &
     &                 ad_vbar_obc, tl_vbar_obc,                        &
     &                 ad_zeta_obc, tl_zeta_obc,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                 ad_ustr, tl_ustr,                                &
     &                 ad_vstr, tl_vstr,                                &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                 ad_tflux, tl_tflux,                              &
#  endif
     &                 ad_t, tl_t,                                      &
     &                 ad_u, tl_u,                                      &
     &                 ad_v, tl_v,                                      &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                 ad_ubar, tl_ubar,                                &
     &                 ad_vbar, tl_vbar,                                &
#  endif
# else
     &                 ad_ubar, tl_ubar,                                &
     &                 ad_vbar, tl_vbar,                                &
# endif
     &                 ad_zeta, tl_zeta)
!
      Lwrk=1
      IF ((innLoop.gt.0).or.Ltrace) THEN
        Linp=1
        Lout=2
!
!-----------------------------------------------------------------------
!  Compute norm Delta(k) as the dot-product between the new vector
!  and previous Lanczos vector.
!-----------------------------------------------------------------------
!
!  Compute current iteration norm Delta(k) used to compute tri-diagonal
!  matrix T(k) in the Lanczos recurrence.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc(:,:,:,:,Lold,:),                     &
     &                    tl_t_obc(:,:,:,:,Lnew,:),                     &
     &                    tl_u_obc(:,:,:,:,Lold),                       &
     &                    tl_u_obc(:,:,:,:,Lnew),                       &
     &                    tl_v_obc(:,:,:,:,Lold),                       &
     &                    tl_v_obc(:,:,:,:,Lnew),                       &
#  endif
     &                    tl_ubar_obc(:,:,:,Lold),                      &
     &                    tl_ubar_obc(:,:,:,Lnew),                      &
     &                    tl_vbar_obc(:,:,:,Lold),                      &
     &                    tl_vbar_obc(:,:,:,Lnew),                      &
     &                    tl_zeta_obc(:,:,:,Lold),                      &
     &                    tl_zeta_obc(:,:,:,Lnew),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr(:,:,:,Lold), tl_ustr(:,:,:,Lnew),     &
     &                    tl_vstr(:,:,:,Lold), tl_vstr(:,:,:,Lnew),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    tl_tflux(:,:,:,Lold,:),                       &
     &                    tl_tflux(:,:,:,Lnew,:),                       &
#  endif
     &                    tl_t(:,:,:,Lold,:), tl_t(:,:,:,Lnew,:),       &
     &                    tl_u(:,:,:,Lold), tl_u(:,:,:,Lnew),           &
     &                    tl_v(:,:,:,Lold), tl_v(:,:,:,Lnew),           &
# else
     &                    tl_ubar(:,:,Lold), tl_ubar(:,:,Lnew),         &
     &                    tl_vbar(:,:,Lold), tl_vbar(:,:,Lnew),         &
# endif
     &                    tl_zeta(:,:,Lold), tl_zeta(:,:,Lnew))

      ae_delta(innLoop,outLoop)=dot(0)
# else
!
!-----------------------------------------------------------------------
!  Compute the action of the act analysis error covariance matrix.
!-----------------------------------------------------------------------
!
!  NOTE: In the case of weak constraint ("WEAK_CONSTRAINT") and
!        and time convolutions ("TIME_CONV"), the state arrays
!        tl_ubar, tl_vbar, ad_ubar, and ad_vbar are only passed
!        as required by the "state" operators routines but they
!        are not used in subsequent calculations.
!
!  Copy tl_var(Lold) into ad_var(Lnew). See NOTE above.
!
      CALL state_copy (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 Lold, Lnew,                                      &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                 ad_t_obc, tl_t_obc,                              &
     &                 ad_u_obc, tl_u_obc,                              &
     &                 ad_v_obc, tl_v_obc,                              &
#   endif
     &                 ad_ubar_obc, tl_ubar_obc,                        &
     &                 ad_vbar_obc, tl_vbar_obc,                        &
     &                 ad_zeta_obc, tl_zeta_obc,                        &
#  endif
#  ifdef ADJUST_WSTRESS
     &                 ad_ustr, tl_ustr,                                &
     &                 ad_vstr, tl_vstr,                                &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                 ad_tflux, tl_tflux,                              &
#   endif
     &                 ad_t, tl_t,                                      &
     &                 ad_u, tl_u,                                      &
     &                 ad_v, tl_v,                                      &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                 ad_ubar, tl_ubar,                                &
     &                 ad_vbar, tl_vbar,                                &
#   endif
#  else
     &                 ad_ubar, tl_ubar,                                &
     &                 ad_vbar, tl_vbar,                                &
#  endif
     &                 ad_zeta, tl_zeta)
!
      Lwrk=1
      IF ((innLoop.gt.0).or.Ltrace) THEN
        Linp=1
        Lout=2
        CALL analysis_error (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Linp, Lout, Lwrk,                          &
     &                       innLoop, outLoop,                          &
#  ifdef MASKING
     &                       rmask, umask, vmask,                       &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       ad_t_obc, ad_u_obc, ad_v_obc,              &
#   endif
     &                       ad_ubar_obc, ad_vbar_obc,                  &
     &                       ad_zeta_obc,                               &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       ad_ustr, ad_vstr,                          &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       ad_tflux,                                  &
#   endif
     &                       ad_t, ad_u, ad_v,                          &
#  else
     &                       ad_ubar, ad_vbar,                          &
#  endif
     &                       ad_zeta,                                   &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       tl_t_obc, tl_u_obc, tl_v_obc,              &
#   endif
     &                       tl_ubar_obc, tl_vbar_obc,                  &
     &                       tl_zeta_obc,                               &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       tl_ustr, tl_vstr,                          &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       tl_tflux,                                  &
#   endif
     &                       tl_t, tl_u, tl_v,                          &
#  else
     &                       tl_ubar, tl_vbar,                          &
#  endif
     &                       tl_zeta,                                   &
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                       nl_t_obc, nl_u_obc, nl_v_obc,              &
#   endif
     &                       nl_ubar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc,                               &
#  endif
#  ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_vstr,                          &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                       nl_tflux,                                  &
#   endif
     &                       nl_t, nl_u, nl_v,                          &
#  else
     &                       nl_ubar, nl_vbar,                          &
#  endif
     &                       nl_zeta)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Check that the analysis error covariance matrix is positive
!  definite, or report the trace estimate.
!
        IF (Ltrace) THEN
          ae_trace(innLoop+1)=ae_delta(innLoop,outLoop)

          CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'ae_trace',    &
     &                          ae_trace(innLoop+1:),                   &
     &                          (/innLoop+1/), (/1/),                   &
     &                          ncid = HSS(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN

          IF (Master) THEN
            WRITE (stdout,10) outLoop, innLoop,                         &
     &                        ae_trace(innLoop+1)
 10         FORMAT (1x,'(',i3.3,',',i3.3,'): ',                         &
     &              'Analysis Error Trace Estimate, ae_trace  = ',      &
     &              1p,e14.7)
          END IF
          RETURN
        END IF
        IF (ae_delta(innLoop,outLoop).le.0.0_r8) THEN
          PRINT *, ' AE_DELTA not positive.'
          PRINT *, ' AE_DELTA = ', ae_delta(innLoop,outLoop),           &
     &             ', outer = ', outLoop, ', inner = ', innLoop
          exit_flag=8
          RETURN
        END IF
      END IF
!
!  Apply the Lanczos recurrence and orthonormalize.
!
      Linp=1
      Lout=2
      CALL lanczos (ng, tile, model,                                    &
     &              LBi, UBi, LBj, UBj, LBij, UBij,                     &
     &              IminS, ImaxS, JminS, JmaxS,                         &
     &              Linp, Lout, Lwrk,                                   &
     &              innLoop, outLoop,                                   &
# ifdef MASKING
     &              rmask, umask, vmask,                                &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &              tl_t_obc, tl_u_obc, tl_v_obc,                       &
#  endif
     &              tl_ubar_obc, tl_vbar_obc,                           &
     &              tl_zeta_obc,                                        &
# endif
# ifdef ADJUST_WSTRESS
     &              tl_ustr, tl_vstr,                                   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &              tl_tflux,                                           &
#  endif
     &              tl_t, tl_u, tl_v,                                   &
# else
     &              tl_ubar, tl_vbar,                                   &
# endif
     &              tl_zeta,                                            &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &              ad_t_obc, ad_u_obc, ad_v_obc,                       &
#  endif
     &              ad_ubar_obc, ad_vbar_obc,                           &
     &              ad_zeta_obc,                                        &
# endif
# ifdef ADJUST_WSTRESS
     &              ad_ustr, ad_vstr,                                   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &              ad_tflux,                                           &
#  endif
     &              ad_t, ad_u, ad_v,                                   &
# else
     &              ad_ubar, ad_vbar,                                   &
# endif
     &              ad_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!  Compute new direction, d(k+1).
!
      Linp=1
      Lout=2
      CALL new_direction (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Linp, Lout,                                   &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc, ad_u_obc, ad_v_obc,                 &
#  endif
     &                    ad_ubar_obc, ad_vbar_obc,                     &
     &                    ad_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
#  endif
     &                    ad_t, ad_u, ad_v,                             &
# else
     &                    ad_ubar, ad_vbar,                             &
# endif
     &                    ad_zeta,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    d_t_obc, d_u_obc, d_v_obc,                    &
#  endif
     &                    d_ubar_obc, d_vbar_obc,                       &
     &                    d_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                    d_sustr, d_svstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    d_stflx,                                      &
#  endif
     &                    d_t, d_u, d_v,                                &
# else
     &                    d_ubar, d_vbar,                               &
# endif
     &                    d_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Determine the eigenvalues and eigenvectors of the tridiagonal matrix.
!  These will be used on the last inner-loop to compute the eigenvectors
!  of the Hessian.
!-----------------------------------------------------------------------
!
      IF (innLoop.gt.0) THEN
        DO i=1,innLoop
          ae_Ritz(i,outLoop)=ae_delta(i,outLoop)
        END DO
        DO i=1,innLoop-1
          ae_Tmatrix(i,1)=ae_beta(i+1,outLoop)
        END DO
!
!  Use the LAPACK routine DSTEQR to compute the eigenvectors and
!  eigenvalues of the tridiagonal matrix. If applicable, the
!  eigenpairs is computed by master thread only. Notice that on
!  exit, the matrix "ae_Tmatrix" is destroyed.
!
        IF (Master) THEN
          CALL DSTEQR ('I', innLoop, ae_Ritz(1,outLoop), ae_Tmatrix,    &
     &                 ae_zv, NpostI, work, info)
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, info)
# endif
        IF (info.ne.0) THEN
          PRINT *, ' Error in DSTEQR: info=', info
          exit_flag=8
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, ae_Ritz(:,outLoop))
        CALL mp_bcastf (ng, model, ae_zv)
# endif
!
!  Estimate the Ritz value error bounds.
!
        DO i=1,innLoop
          ae_RitzErr(i,outLoop)=ABS(ae_beta(innLoop+1,outLoop)*         &
     &                              ae_zv(innLoop,i))
        END DO
!
!  Check for exploding or negative Ritz values.
!
        DO i=1,innLoop
          IF (ae_Ritz(i,outLoop).lt.0.0_r8) THEN
            PRINT *, ' Negative Ritz value found.'
            exit_flag=8
            RETURN
          END IF
        END DO
!
!  Calculate the converged eigenvectors of the Hessian.
!
        IF (innLoop.eq.NpostI) THEN
          RitzMaxErr=HevecErr
          DO i=1,innLoop
            ae_RitzErr(i,outLoop)=ae_RitzErr(i,outLoop)/                &
     &                            ae_Ritz(NpostI,outLoop)
          END DO
          Lwrk=2
          Linp=1
          Lout=2
          CALL posterior_eofs (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp, Lout, Lwrk,                        &
     &                         innLoop, outLoop,                        &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         nl_t_obc, nl_u_obc, nl_v_obc,            &
#  endif
     &                         nl_ubar_obc, nl_vbar_obc,                &
     &                         nl_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         nl_ustr, nl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         nl_tflux,                                &
#  endif
     &                         nl_t, nl_u, nl_v,                        &
# else
     &                         nl_ubar, nl_vbar,                        &
# endif
     &                         nl_zeta,                                 &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         tl_t_obc, tl_u_obc, tl_v_obc,            &
#  endif
     &                         tl_ubar_obc, tl_vbar_obc,                &
     &                         tl_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         tl_ustr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         tl_tflux,                                &
#  endif
     &                         tl_t, tl_u, tl_v,                        &
# else
     &                         tl_ubar, tl_vbar,                        &
# endif
     &                         tl_zeta,                                 &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, ad_u_obc, ad_v_obc,            &
#  endif
     &                         ad_ubar_obc, ad_vbar_obc,                &
     &                         ad_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
#  endif
     &                         ad_t, ad_u, ad_v,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
     &                         ad_zeta)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master.and.(nConvRitz.eq.0)) THEN
            PRINT *, ' No converged Hesssian eigenvectors found.'
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set TLM initial conditions for next inner loop, X(k+1).
!-----------------------------------------------------------------------
!
!   X(k+1) = tau(k+1) * d(k+1)
!
      Linp=1
      Lout=2
      CALL tl_new_vector (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Linp, Lout,                                   &
     &                    innLoop, outLoop,                             &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    d_t_obc, d_u_obc, d_v_obc,                    &
#  endif
     &                    d_ubar_obc, d_vbar_obc,                       &
     &                    d_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                    d_sustr, d_svstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    d_stflx,                                      &
#  endif
     &                    d_t, d_u, d_v,                                &
# else
     &                    d_ubar, d_vbar,                               &
# endif
     &                    d_zeta,                                       &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc, tl_u_obc, tl_v_obc,                 &
#  endif
     &                    tl_ubar_obc, tl_vbar_obc,                     &
     &                    tl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
#  endif
     &                    tl_t, tl_u, tl_v,                             &
# else
     &                    tl_ubar, tl_vbar,                             &
# endif
     &                    tl_zeta,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc, ad_u_obc, ad_v_obc,                 &
#  endif
     &                    ad_ubar_obc, ad_vbar_obc,                     &
     &                    ad_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
#  endif
     &                    ad_t, ad_u, ad_v,                             &
# else
     &                    ad_ubar, ad_vbar,                             &
# endif
     &                    ad_zeta)
!
!-----------------------------------------------------------------------
!  Report posterior analysis error covariance estimation parameters.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        IF (inner.eq.0) THEN
          WRITE (stdout,20) outLoop, innLoop,                           &
     &                      ae_Gnorm(outLoop)
 20       FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                         &
     &            'Analysis Error gradient norm, ae_Gnorm  = ',         &
     &            1p,e14.7,/)
        END IF
        IF (innLoop.gt.0) THEN
          WRITE (stdout,30) RitzMaxErr
 30       FORMAT (/,' Ritz Eigenvalues and relative accuracy: ',        &
     &             'RitzMaxErr = ',1p,e14.7,/)
          ic=0
          DO i=1,innLoop
            IF (ae_RitzErr(i,outLoop).le.RitzMaxErr) THEN
              string='converged'
              ic=ic+1
              WRITE (stdout,40) i, ae_Ritz(i,outLoop),                  &
     &                          ae_RitzErr(i,outLoop),                  &
     &                          TRIM(ADJUSTL(string)), ic
 40           FORMAT(5x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a,2x,           &
     &               '(Good='i3.3,')')
            ELSE
              string='not converged'
              WRITE (stdout,50) i, ae_Ritz(i,outLoop),                  &
     &                          ae_RitzErr(i,outLoop),                  &
     &                          TRIM(ADJUSTL(string))
 50           FORMAT(5x,i3.3,2x,1p,e14.7,2x,1p,e14.7,2x,a)
            END IF
          END DO
          WRITE (stdout,'(/)')
        END IF
      END IF

      RETURN
      END SUBROUTINE posterior_tile

!
!***********************************************************************
      SUBROUTINE tl_new_vector (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Linp, Lout,                             &
     &                          innLoop, outLoop,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          d_t_obc, d_u_obc, d_v_obc,              &
#  endif
     &                          d_ubar_obc, d_vbar_obc,                 &
     &                          d_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                          d_sustr, d_svstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          d_stflx,                                &
#  endif
     &                          d_t, d_u, d_v,                          &
# else
     &                          d_ubar, d_vbar,                         &
# endif
     &                          d_zeta,                                 &
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
     &                          tl_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          ad_t_obc, ad_u_obc, ad_v_obc,           &
#  endif
     &                          ad_ubar_obc, ad_vbar_obc,               &
     &                          ad_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
#  endif
     &                          ad_t, ad_u, ad_v,                       &
# else
     &                          ad_ubar, ad_vbar,                       &
# endif
     &                          ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars
      USE mod_iounits
!
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
      USE state_initialize_mod, ONLY : state_initialize
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: d_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: d_v_obc(LBij:,:,:,:)
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_zeta_obc(LBij:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(in) :: d_svstr(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: d_stflx(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: d_u(LBi:,LBj:,:)
      real(r8), intent(in) :: d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: d_ubar(LBi:,LBj:)
      real(r8), intent(in) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: d_zeta(LBi:,LBj:)
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

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: d_t_obc(LBij:UBij,N(ng),4,                &
     &                                Nbrec(ng),NT(ng))
      real(r8), intent(in) :: d_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(in) :: d_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#   endif
      real(r8), intent(in) :: d_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: d_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(in) :: d_zeta_obc(LBij:UBij,4,Nbrec(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(in) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: d_stflx(LBi:UBi,LBj:UBj,                  &
     &                                Nfrec(ng),NT(ng))
#   endif
      real(r8), intent(in) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: d_zeta(LBi:UBi,LBj:UBj)
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
      integer :: i, j, k, lstr, rec
      integer :: ib, ir, it

      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          tl_zeta(i,j,Lout)=d_zeta(i,j)
# ifdef MASKING
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
# endif
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
              tl_zeta_obc(j,ib,ir,Lout)=d_zeta_obc(j,ib,ir)
#  ifdef MASKING
              tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Lout)*      &
     &                                  rmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              tl_zeta_obc(j,ib,ir,Lout)=d_zeta_obc(j,ib,ir)
#  ifdef MASKING
              tl_zeta_obc(j,ib,ir,Lout)=tl_zeta_obc(j,ib,ir,Lout)*      &
     &                                  rmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              tl_zeta_obc(i,ib,ir,Lout)=d_zeta_obc(i,ib,ir)
#  ifdef MASKING
              tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Lout)*      &
     &                                  rmask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              tl_zeta_obc(i,ib,ir,Lout)=d_zeta_obc(i,ib,ir)
#  ifdef MASKING
              tl_zeta_obc(i,ib,ir,Lout)=tl_zeta_obc(i,ib,ir,Lout)*      &
     &                                  rmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          tl_ubar(i,j,Lout)=d_ubar(i,j)
#  ifdef MASKING
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#  endif
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
              tl_ubar_obc(j,ib,ir,Lout)=d_ubar_obc(j,ib,ir)
#  ifdef MASKING
              tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Lout)*      &
     &                                  umask(Istr,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              tl_ubar_obc(j,ib,ir,Lout)=d_ubar_obc(j,ib,ir)
#  ifdef MASKING
              tl_ubar_obc(j,ib,ir,Lout)=tl_ubar_obc(j,ib,ir,Lout)*      &
     &                                  umask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              tl_ubar_obc(i,ib,ir,Lout)=d_ubar_obc(i,ib,ir)
#  ifdef MASKING
              tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Lout)*      &
     &                                  umask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              tl_ubar_obc(i,ib,ir,Lout)=d_ubar_obc(i,ib,ir)
#  ifdef MASKING
              tl_ubar_obc(i,ib,ir,Lout)=tl_ubar_obc(i,ib,ir,Lout)*      &
     &                                  umask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          tl_vbar(i,j,Lout)=d_vbar(i,j)
#  ifdef MASKING
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#  endif
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
              tl_vbar_obc(j,ib,ir,Lout)=d_vbar_obc(j,ib,ir)
#  ifdef MASKING
              tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Lout)*      &
     &                                  vmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              tl_vbar_obc(j,ib,ir,Lout)=d_vbar_obc(j,ib,ir)
#  ifdef MASKING
              tl_vbar_obc(j,ib,ir,Lout)=tl_vbar_obc(j,ib,ir,Lout)*      &
     &                                  vmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              tl_vbar_obc(i,ib,ir,Lout)=d_vbar_obc(i,ib,ir)
#  ifdef MASKING
              tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Lout)*      &
     &                                  vmask(i,Jstr)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              tl_vbar_obc(i,ib,ir,Lout)=d_vbar_obc(i,ib,ir)
#  ifdef MASKING
              tl_vbar_obc(i,ib,ir,Lout)=tl_vbar_obc(i,ib,ir,Lout)*      &
     &                                  vmask(i,Jend+1)
#  endif
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
            tl_ustr(i,j,ir,Lout)=d_sustr(i,j,ir)
#  ifdef MASKING
            tl_ustr(i,j,ir,Lout)=tl_ustr(i,j,ir,Lout)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            tl_vstr(i,j,ir,Lout)=d_svstr(i,j,ir)
#  ifdef MASKING
            tl_vstr(i,j,ir,Lout)=tl_vstr(i,j,ir,Lout)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            tl_u(i,j,k,Lout)=d_u(i,j,k)
#  ifdef MASKING
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#  endif
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
                tl_u_obc(j,k,ib,ir,Lout)=d_u_obc(j,k,ib,ir)
#   ifdef MASKING
                tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Lout)*      &
     &                                   umask(Istr,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                tl_u_obc(j,k,ib,ir,Lout)=d_u_obc(j,k,ib,ir)
#   ifdef MASKING
                tl_u_obc(j,k,ib,ir,Lout)=tl_u_obc(j,k,ib,ir,Lout)*      &
     &                                   umask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                tl_u_obc(i,k,ib,ir,Lout)=d_u_obc(i,k,ib,ir)
#   ifdef MASKING
                tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Lout)*      &
     &                                   umask(i,Jstr-1)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                tl_u_obc(i,k,ib,ir,Lout)=d_u_obc(i,k,ib,ir)
#   ifdef MASKING
                tl_u_obc(i,k,ib,ir,Lout)=tl_u_obc(i,k,ib,ir,Lout)*      &
     &                                   umask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum.
!
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            tl_v(i,j,k,Lout)=d_v(i,j,k)
#  ifdef MASKING
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#  endif
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
                tl_v_obc(j,k,ib,ir,Lout)=d_v_obc(j,k,ib,ir)
#   ifdef MASKING
                tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Lout)*      &
     &                                   vmask(Istr-1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                tl_v_obc(j,k,ib,ir,Lout)=d_v_obc(j,k,ib,ir)
#   ifdef MASKING
                tl_v_obc(j,k,ib,ir,Lout)=tl_v_obc(j,k,ib,ir,Lout)*      &
     &                                   vmask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_v_obc(i,k,ib,ir,Lout)=d_v_obc(i,k,ib,ir)
#   ifdef MASKING
                tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Lout)*      &
     &                                   vmask(i,Jstr)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_v_obc(i,k,ib,ir,Lout)=d_v_obc(i,k,ib,ir)
#   ifdef MASKING
                tl_v_obc(i,k,ib,ir,Lout)=tl_v_obc(i,k,ib,ir,Lout)*      &
     &                                   vmask(i,Jend+1)
#   endif
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
              tl_t(i,j,k,Lout,it)=d_t(i,j,k,it)
#  ifdef MASKING
              tl_t(i,j,k,Lout,it)=tl_t(i,j,k,Lout,it)*rmask(i,j)
#  endif
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
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_t_obc(j,k,ib,ir,Lout,it)=d_t_obc(j,k,ib,ir,it)
#   ifdef MASKING
                  tl_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                     tl_t_obc(j,k,ib,ir,Lout,it)*rmask(Istr-1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_t_obc(j,k,ib,ir,Lout,it)=d_t_obc(j,k,ib,ir,it)
#   ifdef MASKING
                  tl_t_obc(j,k,ib,ir,Lout,it)=                          &
     &                     tl_t_obc(j,k,ib,ir,Lout,it)*rmask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_t_obc(i,k,ib,ir,Lout,it)=d_t_obc(i,k,ib,ir,it)
#   ifdef MASKING
                  tl_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                     tl_t_obc(i,k,ib,ir,Lout,it)*rmask(i,Jstr-1)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_t_obc(i,k,ib,ir,Lout,it)=d_t_obc(i,k,ib,ir,it)
#   ifdef MASKING
                  tl_t_obc(i,k,ib,ir,Lout,it)=                          &
     &                     tl_t_obc(i,k,ib,ir,Lout,it)*rmask(i,Jend+1)
#   endif
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
                tl_tflux(i,j,ir,Lout,it)=d_stflx(i,j,ir,it)
#   ifdef MASKING
                tl_tflux(i,j,ir,Lout,it)=tl_tflux(i,j,ir,Lout,it)*      &
     &                                   rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif

# endif

      RETURN
      END SUBROUTINE tl_new_vector
!
!***********************************************************************
      SUBROUTINE read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       s_t_obc, s_u_obc, s_v_obc,                 &
#  endif
     &                       s_ubar_obc, s_vbar_obc,                    &
     &                       s_zeta_obc,                                &
# endif
# ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
#  endif
     &                       s_t, s_u, s_v,                             &
# else
     &                       s_ubar, s_vbar,                            &
# endif
     &                       s_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
#  ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
#  endif
# endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod, ONLY : nf_fread3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lwrk, rec, ndef

      integer, intent(inout) :: ncfileid

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ifield, it
      integer :: gtype, ncid, status, varid

      integer, dimension(4) :: Vsize

      real(r8) :: Fmin, Fmax, scale

# include "set_bounds.h"
!
      SourceFile='posterior_lanczos.h, read_state'
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_open (ng, model, ncname, 0, ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
        ncfileid=ncid
      ELSE
        ncid=ncfileid
      END IF

      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idFsur), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  rmask,                                          &
# endif
     &                  s_zeta(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        ifield=idSbry(isFsur)
        gtype=r2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_zeta_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifndef SOLVE3D
!
!  Read in 2D U-momentum component.
!
      gtype=u2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ubar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in 2D V-momentum component.
!
      gtype=v2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vbar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        ifield=idSbry(isUbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_ubar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Read in 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        ifield=idSbry(isVbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_vbar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ustr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vstr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif

# ifdef SOLVE3D
!
!  Read in 3D U-momentum component.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_u(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        ifield=idSbry(isUvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_u_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Read in 3D U-momentum component.
!
      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_v(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        ifield=idSbry(isVvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_v_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Read in tracers.
!
      gtype=r3dvar
      scale=1.0_r8
      DO it=1,NT(ng)
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idTvar(it)),  &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d(ng, model, ncname, ncid,                      &
     &                    Vname(1,idTvar(it)), varid,                   &
     &                    rec, gtype, Vsize,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask,                                        &
#  endif
     &                    s_t(:,:,:,Lwrk,it))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(it))), rec,           &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Read in tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
          ifield=idSbry(isTvar(it))
          gtype=r3dvar
          scale=1.0_r8
          CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),    &
     &                           ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d_bry (ng, model, ncname, ncid,               &
     &                           Vname(1,ifield), varid,                &
     &                           rec, gtype,                            &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale, Fmin, Fmax,                     &
     &                           s_t_obc(:,:,:,:,Lwrk,it))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      scale=1.0_r8
      DO it=1,NT(ng)
        IF (Lstflux(it,ng)) THEN
          CALL netcdf_inq_varid (ng, model, ncname,                     &
     &                           Vname(1,idTsur(it)), ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d(ng, model, ncname, ncid,                    &
     &                      Vname(1,idTsur(it)), varid,                 &
     &                      rec, gtype, Vsize,                          &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      rmask,                                      &
#   endif
     &                      s_tflux(:,:,:,Lwrk,it))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(it))), rec,         &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
!
!  If multiple files, close current file.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state

!
!***********************************************************************
      SUBROUTINE new_direction (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lold, Lnew,                             &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          ad_t_obc, ad_u_obc, ad_v_obc,           &
#  endif
     &                          ad_ubar_obc, ad_vbar_obc,               &
     &                          ad_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
#  endif
     &                          ad_t, ad_u, ad_v,                       &
# else
     &                          ad_ubar, ad_vbar,                       &
# endif
     &                          ad_zeta,                                &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          d_t_obc, d_u_obc, d_v_obc,              &
#  endif
     &                          d_ubar_obc, d_vbar_obc,                 &
     &                          d_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                          d_sustr, d_svstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                          d_stflx,                                &
#  endif
     &                          d_t, d_u, d_v,                          &
# else
     &                          d_ubar, d_vbar,                         &
# endif
     &                          d_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_parallel
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
    defined ADJUST_BOUNDARY
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: d_u_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: d_v_obc(LBij:,:,:,:)
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_vbar_obc(LBij:,:,:)
      real(r8), intent(inout) :: d_zeta_obc(LBij:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: ad_vstr(LBi:UBI,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(in) :: ad_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: d_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),NT(ng))
      real(r8), intent(inout) :: d_u_obc(LBij:UBij,N(ng),4,Nbrec(ng))
      real(r8), intent(inout) :: d_v_obc(LBij:UBij,N(ng),4,Nbrec(ng))
#   endif
      real(r8), intent(inout) :: d_ubar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_vbar_obc(LBij:UBij,4,Nbrec(ng))
      real(r8), intent(inout) :: d_zeta_obc(LBij:UBij,4,Nbrec(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBI,LBj:UBj,Nfrec(ng))
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
#   endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ib, ir, it

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new conjugate descent direction, d(k+1). Notice that the old
!  descent direction is overwritten.
!-----------------------------------------------------------------------
!
!  Free-sruface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          d_zeta(i,j)=ad_zeta(i,j,Lnew)
# ifdef MASKING
          d_zeta(i,j)=d_zeta(i,j)*rmask(i,j)
# endif
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
              d_zeta_obc(j,ib,ir)=ad_zeta_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(j,ib,ir)=d_zeta_obc(j,ib,ir)*                  &
     &                            rmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              d_zeta_obc(j,ib,ir)=ad_zeta_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(j,ib,ir)=d_zeta_obc(j,ib,ir)*                  &
     &                            rmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              d_zeta_obc(i,ib,ir)=ad_zeta_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(i,ib,ir)=d_zeta_obc(i,ib,ir)*                  &
     &                            rmask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              d_zeta_obc(i,ib,ir)=ad_zeta_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_zeta_obc(i,ib,ir)=d_zeta_obc(i,ib,ir)*                  &
     &                            rmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          d_ubar(i,j)=ad_ubar(i,j,Lnew)
#  ifdef MASKING
          d_ubar(i,j)=d_ubar(i,j)*umask(i,j)
#  endif
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
              d_ubar_obc(j,ib,ir)=ad_ubar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(j,ib,ir)=d_ubar_obc(j,ib,ir)*                  &
     &                            umask(Istr,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              d_ubar_obc(j,ib,ir)=ad_ubar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(j,ib,ir)=d_ubar_obc(j,ib,ir)*                  &
     &                            umask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              d_ubar_obc(i,ib,ir)=ad_ubar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(i,ib,ir)=d_ubar_obc(i,ib,ir)*                  &
     &                            umask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              d_ubar_obc(i,ib,ir)=ad_ubar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_ubar_obc(i,ib,ir)=d_ubar_obc(i,ib,ir)*                  &
     &                            umask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          d_vbar(i,j)=ad_vbar(i,j,Lnew)
#  ifdef MASKING
          d_vbar(i,j)=d_vbar(i,j)*vmask(i,j)
#  endif
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
              d_vbar_obc(j,ib,ir)=ad_vbar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(j,ib,ir)=d_vbar_obc(j,ib,ir)*                  &
     &                            vmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              d_vbar_obc(j,ib,ir)=ad_vbar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(j,ib,ir)=d_vbar_obc(j,ib,ir)*                  &
     &                            vmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              d_vbar_obc(i,ib,ir)=ad_vbar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(i,ib,ir)=d_vbar_obc(i,ib,ir)*                  &
     &                            vmask(i,Jstr)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              d_vbar_obc(i,ib,ir)=ad_vbar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              d_vbar_obc(i,ib,ir)=d_vbar_obc(i,ib,ir)*                  &
     &                            vmask(i,Jend+1)
#  endif
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
            d_sustr(i,j,ir)=ad_ustr(i,j,ir,Lnew)
#  ifdef MASKING
            d_sustr(i,j,ir)=d_sustr(i,j,ir)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            d_svstr(i,j,ir)=ad_vstr(i,j,ir,Lnew)
#  ifdef MASKING
            d_svstr(i,j,ir)=d_svstr(i,j,ir)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            d_u(i,j,k)=ad_u(i,j,k,Lnew)
#  ifdef MASKING
            d_u(i,j,k)=d_u(i,j,k)*umask(i,j)
#  endif
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
                d_u_obc(j,k,ib,ir)=ad_u_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(j,k,ib,ir)=d_u_obc(j,k,ib,ir)*                  &
     &                             umask(Istr,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                d_u_obc(j,k,ib,ir)=ad_u_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(j,k,ib,ir)=d_u_obc(j,k,ib,ir)*                  &
     &                             umask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                d_u_obc(i,k,ib,ir)=ad_u_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(i,k,ib,ir)=d_u_obc(i,k,ib,ir)*                  &
     &                             umask(i,Jstr-1)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                d_u_obc(i,k,ib,ir)=ad_u_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_u_obc(i,k,ib,ir)=d_u_obc(i,k,ib,ir)*                  &
     &                             umask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum.
!
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            d_v(i,j,k)=ad_v(i,j,k,Lnew)
#  ifdef MASKING
            d_v(i,j,k)=d_v(i,j,k)*vmask(i,j)
#  endif
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
                d_v_obc(j,k,ib,ir)=ad_v_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(j,k,ib,ir)=d_v_obc(j,k,ib,ir)*                  &
     &                             vmask(Istr-1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                d_v_obc(j,k,ib,ir)=ad_v_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(j,k,ib,ir)=d_v_obc(j,k,ib,ir)*                  &
     &                             vmask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                d_v_obc(i,k,ib,ir)=ad_v_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(i,k,ib,ir)=d_v_obc(i,k,ib,ir)*                  &
     &                             vmask(i,Jstr)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                d_v_obc(i,k,ib,ir)=ad_v_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                d_v_obc(i,k,ib,ir)=d_v_obc(i,k,ib,ir)*                  &
     &                             vmask(i,Jend+1)
#   endif
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
              d_t(i,j,k,it)=ad_t(i,j,k,Lnew,it)
#  ifdef MASKING
              d_t(i,j,k,it)=d_t(i,j,k,it)*rmask(i,j)
#  endif
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
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  d_t_obc(j,k,ib,ir,it)=ad_t_obc(j,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(j,k,ib,ir,it)=d_t_obc(j,k,ib,ir,it)*          &
     &                                  rmask(Istr-1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  d_t_obc(j,k,ib,ir,it)=ad_t_obc(j,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(j,k,ib,ir,it)=d_t_obc(j,k,ib,ir,it)*          &
     &                                  rmask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  d_t_obc(i,k,ib,ir,it)=ad_t_obc(i,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(i,k,ib,ir,it)=d_t_obc(i,k,ib,ir,it)*          &
     &                                  rmask(i,Jstr-1)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  d_t_obc(i,k,ib,ir,it)=ad_t_obc(i,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  d_t_obc(i,k,ib,ir,it)=d_t_obc(i,k,ib,ir,it)*          &
     &                                  rmask(i,Jend+1)
#   endif
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
                d_stflx(i,j,ir,it)=ad_tflux(i,j,ir,Lnew,it)
#   ifdef MASKING
                d_stflx(i,j,ir,it)=d_stflx(i,j,ir,it)*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE new_direction
!
!***********************************************************************
      SUBROUTINE analysis_error (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lold, Lnew, Lwrk,                      &
     &                           innLoop, outLoop,                      &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           ad_ubar, ad_vbar,                      &
#  endif
# else
     &                           ad_ubar, ad_vbar,                      &
# endif
     &                           ad_zeta,                               &
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           tl_ubar, tl_vbar,                      &
#  endif
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta,                               &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           nl_t_obc, nl_u_obc, nl_v_obc,          &
#  endif
     &                           nl_ubar_obc, nl_vbar_obc,              &
     &                           nl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           nl_ustr, nl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           nl_tflux,                              &
#  endif
     &                           nl_t, nl_u, nl_v,                      &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           nl_ubar, nl_vbar,                      &
#  endif
# else
     &                           nl_ubar, nl_vbar,                      &
# endif
     &                           nl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_copy_mod, ONLY : state_copy
      USE state_addition_mod, ONLY : state_addition
      USE state_initialize_mod, ONLY : state_initialize
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
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

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
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
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, lstr
      integer :: ib, ir, it, ivec, Ltmp1, Ltmp2, rec

      real(r8) :: fac, fac1, fac2, zbet

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(1:Ninner) :: DotProd

      real(r8), dimension(Ninner) :: zu, zgam

      character (len=256) :: ncname

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute analysis error covariance matrix.
!-----------------------------------------------------------------------
!
!  NOTE: In the case of weak constraint ("WEAK_CONSTRAINT") and
!        and time convolutions ("TIME_CONV"), the state arrays
!        tl_ubar, tl_vbar, ad_ubar, and ad_vbar are only passed
!        as required by the "state" operators routines but they
!        are not used in subsequent calculations.
!
      Ltmp1=1
      Ltmp2=2
!
!  Compute the dot-product of ad_var(Lnew) with each evolved and
!  convolved Lanczos vector of the stabilized representer matrix
!  which are temporarily stored in the hessian netcdf file.
!
      ncname=HSS(ng)%name
      DO ivec=1,Ninner
        rec=ivec
!
!  Read vectors stored in hessian netcdf file using nl_var(Ltmp1) as
!  temporary storage.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Ltmp1, rec,                                    &
     &                   0, HSS(ng)%ncid, ncname,                       &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, nl_u_obc, nl_v_obc,                  &
#  endif
     &                   nl_ubar_obc, nl_vbar_obc,                      &
     &                   nl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux,                                      &
#  endif
     &                   nl_t, nl_u, nl_v,                              &
# else
     &                   nl_ubar, nl_vbar,                              &
# endif
     &                   nl_zeta)
!
!  Compute dot product.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      nl_t_obc(:,:,:,:,Ltmp1,:),                  &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      nl_u_obc(:,:,:,:,Ltmp1),                    &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      nl_v_obc(:,:,:,:,Ltmp1),                    &
#  endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      nl_ubar_obc(:,:,:,Ltmp1),                   &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      nl_vbar_obc(:,:,:,Ltmp1),                   &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      nl_zeta_obc(:,:,:,Ltmp1),                   &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), nl_ustr(:,:,:,Ltmp1),  &
     &                      ad_vstr(:,:,:,Lnew), nl_vstr(:,:,:,Ltmp1),  &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      nl_tflux(:,:,:,Ltmp1,:),                    &
#  endif
     &                      ad_t(:,:,:,Lnew,:), nl_t(:,:,:,Ltmp1,:),    &
     &                      ad_u(:,:,:,Lnew), nl_u(:,:,:,Ltmp1),        &
     &                      ad_v(:,:,:,Lnew), nl_v(:,:,:,Ltmp1),        &
# else
     &                      ad_ubar(:,:,Lnew), nl_ubar(:,:,Ltmp1),      &
     &                      ad_vbar(:,:,Lnew), nl_vbar(:,:,Ltmp1),      &
# endif
     &                      ad_zeta(:,:,Lnew), nl_zeta(:,:,Ltmp1))

        DotProd(ivec)=dot(0)
      END DO
!
!   Now multiply the result by the inverse tridiagonal matrix of
!   stabilized representer matrix Lanczos vector coefficients.
!
      zbet=cg_delta(1,outLoop)
      zu(1)=DotProd(1)/zbet
      DO ivec=2,Ninner
        zgam(ivec)=cg_beta(ivec,outLoop)/zbet
        zbet=cg_delta(ivec,outLoop)-cg_beta(ivec,outLoop)*zgam(ivec)
        zu(ivec)=(DotProd(ivec)-cg_beta(ivec,outLoop)*                  &
     &            zu(ivec-1))/zbet
      END DO

      DO ivec=Ninner-1,1,-1
        zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
      END DO
!
!   Now form the weighted sum of the evolved and convolved Lanczos vector
!   of the stabilized representer matrix. Use nl_var(2) as temporary
!   storage.
!
!  Initialize nonl-linear state arrays: nl_var(Ltmp2) = fac
!
      fac=0.0_r8

      CALL state_initialize (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Ltmp2, fac,                                &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_u_obc, nl_v_obc,              &
#  endif
     &                       nl_ubar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc,                               &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       nl_tflux,                                  &
#  endif
     &                       nl_t, nl_u, nl_v,                          &
# else
     &                       nl_ubar, nl_vbar,                          &
# endif
     &                       nl_zeta)
!
      ncname=HSS(ng)%name
      DO ivec=1,Ninner
        rec=ivec
!
!  Read vectors stored in hessian netcdf file using nl_var(Ltmp1) as
!  temporary storage.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Ltmp1, rec,                                    &
     &                   0, HSS(ng)%ncid, ncname,                       &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, nl_u_obc, nl_v_obc,                  &
#  endif
     &                   nl_ubar_obc, nl_vbar_obc,                      &
     &                   nl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux,                                      &
#  endif
     &                   nl_t, nl_u, nl_v,                              &
# else
     &                   nl_ubar, nl_vbar,                              &
# endif
     &                   nl_zeta)
!
!  Compute the sum: (See NOTE above)
!
!    nl_var(Ltmp2) = fac1 * nl_var(Ltmp2) + fac2 * nl_var(Ltmp1)
!
        fac1=1.0_r8
        fac2=zu(ivec)
!
        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Ltmp2, Ltmp1, Ltmp2, fac1, fac2,           &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       nl_t_obc, nl_t_obc,                        &
     &                       nl_u_obc, nl_u_obc,                        &
     &                       nl_v_obc, nl_v_obc,                        &
#  endif
     &                       nl_ubar_obc, nl_ubar_obc,                  &
     &                       nl_vbar_obc, nl_vbar_obc,                  &
     &                       nl_zeta_obc, nl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_ustr,                          &
     &                       nl_vstr, nl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       nl_tflux, nl_tflux,                        &
#  endif
     &                       nl_t, nl_t,                                &
     &                       nl_u, nl_u,                                &
     &                       nl_v, nl_v,                                &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       nl_ubar, nl_ubar,                          &
     &                       nl_vbar, nl_vbar,                          &
#  endif
# else
     &                       nl_ubar, nl_ubar,                          &
     &                       nl_vbar, nl_vbar,                          &
# endif
     &                       nl_zeta, nl_zeta)
      END DO
!
!  Copy nl_var(Ltmp2) into ad_var(Lnew).  See NOTE above.
!
      CALL state_copy (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, LBij, UBij,                  &
     &                 Ltmp2, Lnew,                                     &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                 ad_t_obc, nl_t_obc,                              &
     &                 ad_u_obc, nl_u_obc,                              &
     &                 ad_v_obc, nl_v_obc,                              &
#  endif
     &                 ad_ubar_obc, nl_ubar_obc,                        &
     &                 ad_vbar_obc, nl_vbar_obc,                        &
     &                 ad_zeta_obc, nl_zeta_obc,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                 ad_ustr, nl_ustr,                                &
     &                 ad_vstr, nl_vstr,                                &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                 ad_tflux, nl_tflux,                              &
#  endif
     &                 ad_t, nl_t,                                      &
     &                 ad_u, nl_u,                                      &
     &                 ad_v, nl_v,                                      &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                 ad_ubar, nl_ubar,                                &
     &                 ad_vbar, nl_vbar,                                &
#  endif
# else
     &                 ad_ubar, nl_ubar,                                &
     &                 ad_vbar, nl_vbar,                                &
# endif
     &                 ad_zeta, nl_zeta)
!
!  Now form the final analysis error covariance matrix-vector product:
!
!    y=(D-DG'VT^-1V'GD)x
!
!  where tl_var(Lnew)=Dx and ad_var(Lnew)=DG'VT^-1V'GDx
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta(i,j,Lnew)=tl_zeta(i,j,Lnew)-ad_zeta(i,j,Lnew)
# ifdef MASKING
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)*rmask(i,j)
# endif
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
              ad_zeta_obc(j,ib,ir,Lnew)=tl_zeta_obc(j,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              ad_zeta_obc(j,ib,ir,Lnew)=ad_zeta_obc(j,ib,ir,Lnew)*      &
     &                                  rmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_zeta_obc(j,ib,ir,Lnew)=tl_zeta_obc(j,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              ad_zeta_obc(j,ib,ir,Lnew)=ad_zeta_obc(j,ib,ir,Lnew)*      &
     &                                  rmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Lnew)=tl_zeta_obc(i,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              ad_zeta_obc(i,ib,ir,Lnew)=ad_zeta_obc(i,ib,ir,Lnew)*      &
     &                                  rmask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Lnew)=tl_zeta_obc(i,ib,ir,Lnew)-      &
     &                                  ad_zeta_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              ad_zeta_obc(i,ib,ir,Lnew)=ad_zeta_obc(i,ib,ir,Lnew)*      &
     &                                  rmask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D U-momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar(i,j,Lnew)=tl_ubar(i,j,Lnew)-ad_ubar(i,j,Lnew)
#  ifdef MASKING
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)*umask(i,j)
#  endif
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
              ad_ubar_obc(j,ib,ir,Lnew)=tl_ubar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              ad_ubar_obc(j,ib,ir,Lnew)=ad_ubar_obc(j,ib,ir,Lnew)*      &
     &                                  umask(Istr,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_ubar_obc(j,ib,ir,Lnew)=tl_ubar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              ad_ubar_obc(j,ib,ir,Lnew)=ad_ubar_obc(j,ib,ir,Lnew)*      &
     &                                  umask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Lnew)=tl_ubar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              ad_ubar_obc(i,ib,ir,Lnew)=ad_ubar_obc(i,ib,ir,Lnew)*      &
     &                                  umask(i,Jstr-1)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Lnew)=tl_ubar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_ubar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              ad_ubar_obc(i,ib,ir,Lnew)=ad_ubar_obc(i,ib,ir,Lnew)*      &
     &                                  umask(i,Jend+1)
#  endif
            END DO
          END IF
        END DO
      END IF
# endif

# ifndef SOLVE3D
!
!  2D V-momentum.
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          ad_vbar(i,j,Lnew)=tl_vbar(i,j,Lnew)-ad_vbar(i,j,Lnew)
#  ifdef MASKING
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)*vmask(i,j)
#  endif
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
              ad_vbar_obc(j,ib,ir,Lnew)=tl_vbar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              ad_vbar_obc(j,ib,ir,Lnew)=ad_vbar_obc(j,ib,ir,Lnew)*      &
     &                                  vmask(Istr-1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              ad_vbar_obc(j,ib,ir,Lnew)=tl_vbar_obc(j,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(j,ib,ir,Lnew)
#  ifdef MASKING
              ad_vbar_obc(j,ib,ir,Lnew)=ad_vbar_obc(j,ib,ir,Lnew)*      &
     &                                  vmask(Iend+1,j)
#  endif
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Lnew)=tl_vbar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              ad_vbar_obc(i,ib,ir,Lnew)=ad_vbar_obc(i,ib,ir,Lnew)*      &
     &                                  vmask(i,Jstr)
#  endif
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Lnew)=tl_vbar_obc(i,ib,ir,Lnew)-      &
     &                                  ad_vbar_obc(i,ib,ir,Lnew)
#  ifdef MASKING
              ad_vbar_obc(i,ib,ir,Lnew)=ad_vbar_obc(i,ib,ir,Lnew)*      &
     &                                  vmask(i,Jend+1)
#  endif
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
            ad_ustr(i,j,ir,Lnew)=tl_ustr(i,j,ir,Lnew)-                  &
     &                           ad_ustr(i,j,ir,Lnew)
#  ifdef MASKING
            ad_ustr(i,j,ir,Lnew)=ad_ustr(i,j,ir,Lnew)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            ad_vstr(i,j,ir,Lnew)=tl_vstr(i,j,ir,Lnew)-                  &
     &                           ad_vstr(i,j,ir,Lnew)
#  ifdef MASKING
            ad_vstr(i,j,ir,Lnew)=ad_vstr(i,j,ir,Lnew)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
!
!  3D U-momentum.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            ad_u(i,j,k,Lnew)=tl_u(i,j,k,Lnew)-ad_u(i,j,k,Lnew)
#  ifdef MASKING
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)*umask(i,j)
#  endif
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
                ad_u_obc(j,k,ib,ir,Lnew)=tl_u_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_u_obc(j,k,ib,ir,Lnew)=ad_u_obc(j,k,ib,ir,Lnew)*      &
     &                                   umask(Istr,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                ad_u_obc(j,k,ib,ir,Lnew)=tl_u_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_u_obc(j,k,ib,ir,Lnew)=ad_u_obc(j,k,ib,ir,Lnew)*      &
     &                                   umask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Lnew)=tl_u_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_u_obc(i,k,ib,ir,Lnew)=ad_u_obc(i,k,ib,ir,Lnew)*      &
     &                                   umask(i,Jstr-1)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Lnew)=tl_u_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_u_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_u_obc(i,k,ib,ir,Lnew)=ad_u_obc(i,k,ib,ir,Lnew)*      &
     &                                   umask(i,Jend+1)
#   endif
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
!
!  3D V-momentum.
!
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            ad_v(i,j,k,Lnew)=tl_v(i,j,k,Lnew)-ad_v(i,j,k,Lnew)
#  ifdef MASKING
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)*vmask(i,j)
#  endif
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
                ad_v_obc(j,k,ib,ir,Lnew)=tl_v_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_v_obc(j,k,ib,ir,Lnew)=ad_v_obc(j,k,ib,ir,Lnew)*      &
     &                                   vmask(Istr-1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                ad_v_obc(j,k,ib,ir,Lnew)=tl_v_obc(j,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(j,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_v_obc(j,k,ib,ir,Lnew)=ad_v_obc(j,k,ib,ir,Lnew)*      &
     &                                   vmask(Iend+1,j)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Lnew)=tl_v_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_v_obc(i,k,ib,ir,Lnew)=ad_v_obc(i,k,ib,ir,Lnew)*      &
     &                                   vmask(i,Jstr)
#   endif
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Lnew)=tl_v_obc(i,k,ib,ir,Lnew)-      &
     &                                   ad_v_obc(i,k,ib,ir,Lnew)
#   ifdef MASKING
                ad_v_obc(i,k,ib,ir,Lnew)=ad_v_obc(i,k,ib,ir,Lnew)*      &
     &                                   vmask(i,Jend+1)
#   endif
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
              ad_t(i,j,k,Lnew,it)=tl_t(i,j,k,Lnew,it)-                  &
     &                            ad_t(i,j,k,Lnew,it)
#  ifdef MASKING
              ad_t(i,j,k,Lnew,it)=ad_t(i,j,k,Lnew,it)*rmask(i,j)
#  endif
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
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                                    tl_t_obc(j,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(j,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(j,k,ib,ir,Lnew,it)*rmask(Istr-1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                                    tl_t_obc(j,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(j,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  ad_t_obc(j,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(j,k,ib,ir,Lnew,it)*rmask(Iend+1,j)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                                    tl_t_obc(i,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(i,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(i,k,ib,ir,Lnew,it)*rmask(i,Jstr-1)
#   endif
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                                    tl_t_obc(i,k,ib,ir,Lnew,it)-  &
     &                                    ad_t_obc(i,k,ib,ir,Lnew,it)
#   ifdef MASKING
                  ad_t_obc(i,k,ib,ir,Lnew,it)=                          &
     &                    ad_t_obc(i,k,ib,ir,Lnew,it)*rmask(i,Jend+1)
#   endif
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
                ad_tflux(i,j,ir,Lnew,it)=tl_tflux(i,j,ir,Lnew,it)-      &
     &                                   ad_tflux(i,j,ir,Lnew,it)
#   ifdef MASKING
                ad_tflux(i,j,ir,Lnew,it)=ad_tflux(i,j,ir,Lnew,it)*      &
     &                                   rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute norm Delta(k) as the dot-product between the new vector
!  and previous Lanczos vector.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Compute current iteration norm Delta(k) used to compute tri-diagonal
!  matrix T(k) in the Lanczos recurrence.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc(:,:,:,:,Lnew,:),                     &
     &                    tl_t_obc(:,:,:,:,Lwrk,:),                     &
     &                    ad_u_obc(:,:,:,:,Lnew),                       &
     &                    tl_u_obc(:,:,:,:,Lwrk),                       &
     &                    ad_v_obc(:,:,:,:,Lnew),                       &
     &                    tl_v_obc(:,:,:,:,Lwrk),                       &
#  endif
     &                    ad_ubar_obc(:,:,:,Lnew),                      &
     &                    tl_ubar_obc(:,:,:,Lwrk),                      &
     &                    ad_vbar_obc(:,:,:,Lnew),                      &
     &                    tl_vbar_obc(:,:,:,Lwrk),                      &
     &                    ad_zeta_obc(:,:,:,Lnew),                      &
     &                    tl_zeta_obc(:,:,:,Lwrk),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),     &
     &                    ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lnew,:),                       &
     &                    tl_tflux(:,:,:,Lwrk,:),                       &
#  endif
     &                    ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),       &
     &                    ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),           &
     &                    ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),           &
# else
     &                    ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),         &
     &                    ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),         &
# endif
     &                    ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))

      ae_delta(innLoop,outLoop)=dot(0)

      RETURN
      END SUBROUTINE analysis_error
!
!***********************************************************************
      SUBROUTINE lanczos (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lold, Lnew, Lwrk,                             &
     &                    innLoop, outLoop,                             &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc, tl_u_obc, tl_v_obc,                 &
#  endif
     &                    tl_ubar_obc, tl_vbar_obc,                     &
     &                    tl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
#  endif
     &                    tl_t, tl_u, tl_v,                             &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                    tl_ubar, tl_vbar,                             &
#  endif
# else
     &                    tl_ubar, tl_vbar,                             &
# endif
     &                    tl_zeta,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc, ad_u_obc, ad_v_obc,                 &
#  endif
     &                    ad_ubar_obc, ad_vbar_obc,                     &
     &                    ad_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
#  endif
     &                    ad_t, ad_u, ad_v,                             &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                    ad_ubar, ad_vbar,                             &
#  endif
# else
     &                    ad_ubar, ad_vbar,                             &
# endif
     &                    ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
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
      integer :: i, j, lstr, rec

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:NpostI) :: DotProd, dot_new, dot_old

      character (len=256) :: ncname

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Calculate the new Lanczos vector, q(k+1) using reccurence equation
!  for the gradient vectors:
!
!     H q(k+1) = Gamma(k+1) q(k+2) + Delta(k+1) q(k+1) + Gamma(k) q(k)
!
!  where  Gamma(k) = - SQRT ( Beta(k+1) ) / Alpha(k)
!-----------------------------------------------------------------------
!
!  NOTE: In the case of weak constraint ("WEAK_CONSTRAINT") and
!        and time convolutions ("TIME_CONV"), the state arrays
!        tl_ubar, tl_vbar, ad_ubar, and ad_vbar are only passed
!        as required by the "state" operators routines but they
!        are not used in subsequent calculations.
!
!  At this point, the previous orthonormal Lanczos vector is still in
!  tangent linear state arrays (index Lwrk) - it was read in the
!  routine anerror.
!
      IF (innLoop.gt.0) THEN
!
!  Compute new Lanczos vector:  (See NOTE above)
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
        fac1=1.0_r8
        fac2=-ae_delta(innLoop,outLoop)
!
        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#  endif
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
     &                       ad_zeta, tl_zeta)
      END IF
!
!  Substract previous orthonormal Lanczos vector.
!
      IF (innLoop.gt.1) THEN
!
!  Determine adjoint file to process.
!
        ncname=ADM(ng)%name
!
!  Read in the previous (innLoop-1) orthonormal Lanczos vector.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, innLoop-1,                               &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Substract previous orthonormal Lanczos vector:  (See NOTE above)
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
        fac1=1.0_r8
        fac2=-ae_beta(innLoop,outLoop)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#  endif
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
     &                       ad_zeta, tl_zeta)
      END IF
!
!-----------------------------------------------------------------------
!  Orthogonalize current gradient, q(k+1), against all previous
!  gradients (reverse order) using Gramm-Schmidt procedure.
!-----------------------------------------------------------------------
!
!  We can overwrite adjoint arrays at index Lnew each time around the
!  the following loop because the preceding gradient vectors that we
!  read are orthogonal to each other. The reversed order of the loop
!  is important for the Lanczos vector calculations.
!
      ncname=ADM(ng)%name
!
      DO rec=innLoop,1,-1
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product <q(k+1), q(rec)>.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#  endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
#  endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
# else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
# endif
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
!
!  Compute Gramm-Schmidt scaling coefficient.
!
        DotProd(rec)=dot(0)
!
!  Gramm-Schmidt orthonormalization, free-surface.  See NOTE above.
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
        fac1=1.0_r8
        fac2=-DotProd(rec)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       ad_t_obc, tl_t_obc,                        &
     &                       ad_u_obc, tl_u_obc,                        &
     &                       ad_v_obc, tl_v_obc,                        &
#  endif
     &                       ad_ubar_obc, tl_ubar_obc,                  &
     &                       ad_vbar_obc, tl_vbar_obc,                  &
     &                       ad_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
#  endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#  endif
# else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
# endif
     &                       ad_zeta, tl_zeta)
      END DO
!
!-----------------------------------------------------------------------
!  Normalize current orthogonal Lanczos vector.
!-----------------------------------------------------------------------
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    ad_t_obc(:,:,:,:,Lnew,:),                     &
     &                    ad_t_obc(:,:,:,:,Lnew,:),                     &
     &                    ad_u_obc(:,:,:,:,Lnew),                       &
     &                    ad_u_obc(:,:,:,:,Lnew),                       &
     &                    ad_v_obc(:,:,:,:,Lnew),                       &
     &                    ad_v_obc(:,:,:,:,Lnew),                       &
#  endif
     &                    ad_ubar_obc(:,:,:,Lnew),                      &
     &                    ad_ubar_obc(:,:,:,Lnew),                      &
     &                    ad_vbar_obc(:,:,:,Lnew),                      &
     &                    ad_vbar_obc(:,:,:,Lnew),                      &
     &                    ad_zeta_obc(:,:,:,Lnew),                      &
     &                    ad_zeta_obc(:,:,:,Lnew),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),     &
     &                    ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),     &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    ad_tflux(:,:,:,Lnew,:),                       &
     &                    ad_tflux(:,:,:,Lnew,:),                       &
#  endif
     &                    ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),       &
     &                    ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),           &
     &                    ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),           &
# else
     &                    ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),         &
     &                    ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),         &
# endif
     &                    ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
!
!  Compute normalization factor.
!
      IF (innLoop.eq.0) THEN
        ae_Gnorm(outLoop)=SQRT(dot(0))
      ELSE
        ae_beta(innLoop+1,outLoop)=SQRT(dot(0))
      END IF
!
!  Normalize the vector: ad_var(Lnew) = fac * ad_var(Lnew)
!
      fac=1.0_r8/SQRT(dot(0))

      CALL state_scale (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, LBij, UBij,                 &
     &                  Lnew, Lnew, fac,                                &
# ifdef MASKING
     &                  rmask, umask, vmask,                            &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                  ad_t_obc, ad_u_obc, ad_v_obc,                   &
#  endif
     &                  ad_ubar_obc, ad_vbar_obc,                       &
     &                  ad_zeta_obc,                                    &
# endif
# ifdef ADJUST_WSTRESS
     &                  ad_ustr, ad_vstr,                               &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                  ad_tflux,                                       &
#  endif
     &                  ad_t, ad_u, ad_v,                               &
# else
     &                  ad_ubar, ad_vbar,                               &
# endif
     &                  ad_zeta)
!
# ifdef TEST_ORTHOGONALIZATION
!
!-----------------------------------------------------------------------
!  Test orthogonal properties of the new gradient.
!-----------------------------------------------------------------------
!
!  Determine adjoint file to process.
!
      ncname=ADM(ng)%name
!
      DO rec=innLoop,1,-1
!
!  Read in each previous gradient state solutions, q(0) to q(k), and
!  compute its associated dot angaint orthogonalized q(k+1). Again,
!  each gradient solution is loaded into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ADM(ng)%ncid, ncname,             &
#  ifdef MASKING
     &                   rmask, umask, vmask,                           &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#   endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
#  endif
#  ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
#   endif
     &                   tl_t, tl_u, tl_v,                              &
#  else
     &                   tl_ubar, tl_vbar,                              &
#  endif
     &                   tl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
#  ifdef MASKING
     &                      rmask, umask, vmask,                        &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#   endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
#  endif
#  ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
#   endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#  else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
#  endif
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
        dot_new(rec)=dot(0)
      END DO
!
!  Report dot products. If everything is working correctly, at the
!  end of the orthogonalization dot_new(rec) << dot_old(rec).
!
      IF (Master) THEN
        WRITE (stdout,20) outLoop, innLoop
        DO rec=innLoop,1,-1
          WRITE (stdout,30) DotProd(rec), rec-1
        END DO
        WRITE (stdout,*) ' '
        DO rec=innLoop,1,-1
          WRITE (stdout,40) innLoop, rec-1, dot_new(rec),               &
     &                      rec-1, rec-1, dot_old(rec)
        END DO
 20     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'Gramm-Schmidt Orthogonalization:',/)
 30     FORMAT (12x,'Orthogonalization Factor = ',1p,e19.12,3x,         &
     &          '(Iter=',i3.3,')')
 40     FORMAT (2x,'Ortho Test: ',                                      &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8,1x,             &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8)
      END IF
# endif

      RETURN
      END SUBROUTINE lanczos
!
!***********************************************************************
      SUBROUTINE posterior_eofs (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lold, Lnew, Lwrk,                      &
     &                           innLoop, outLoop,                      &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           nl_t_obc, nl_u_obc, nl_v_obc,          &
#  endif
     &                           nl_ubar_obc, nl_vbar_obc,              &
     &                           nl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           nl_ustr, nl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           nl_tflux,                              &
#  endif
     &                           nl_t, nl_u, nl_v,                      &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           nl_ubar, nl_vbar,                      &
#  endif
# else
     &                           nl_ubar, nl_vbar,                      &
# endif
     &                           nl_zeta,                               &
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           tl_ubar, tl_vbar,                      &
#  endif
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
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                           ad_ubar, ad_vbar,                      &
#  endif
# else
     &                           ad_ubar, ad_vbar,                      &
# endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
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
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
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
#  else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ingood, lstr, rec, nvec, status, varid
      integer :: L1
      integer :: start(4), total(4)

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(NpostI) :: RitzErr

      real(r8), dimension(0:NstateVar(ng)) :: dot

      character (len=256) :: ncname

# include "set_bounds.h"
!
      SourceFile='posterior.F, posterior_eofs'
!
!-----------------------------------------------------------------------
!  Calculate converged eigenvectors of the Hessian.
!-----------------------------------------------------------------------
!
!  NOTE: In the case of weak constraint ("WEAK_CONSTRAINT") and
!        and time convolutions ("TIME_CONV"), the state arrays
!        tl_ubar, tl_vbar, ad_ubar, and ad_vbar are only passed
!        as required by the "state" operators routines but they
!        are not used in subsequent calculations.
!
!  Count and collect the converged eigenvalues.
!
      ingood=0
      DO i=innLoop,1,-1
        ingood=ingood+1
        Ritz(ingood)=ae_Ritz(i,outLoop)
        RitzErr(ingood)=ae_RitzErr(i,outLoop)
      END DO
      nConvRitz=ingood
!
!  Write out number of converged eigenvalues.
!
      CALL netcdf_put_ivar (ng, model, HSS(ng)%name, 'nConvRitz',       &
     &                      nConvRitz, (/0/), (/0/),                    &
     &                      ncid = HSS(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  First, premultiply the eigenvectors of the tridiagonal
!  matrix T(k) by the matrix of Lanczos vectors Q(k).  Use tangent
!  linear (index Lwrk) and adjoint (index Lold) state arrays as
!  temporary storage.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,10)
!
      HSS(ng)%Rindex=0
!
      COLUMNS : DO nvec=innLoop,1,-1
!
!  Initialize adjoint state arrays: ad_var(Lold) = fac
!
        fac=0.0_r8

        CALL state_initialize (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lold, fac,                               &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, ad_u_obc, ad_v_obc,            &
#  endif
     &                         ad_ubar_obc, ad_vbar_obc,                &
     &                         ad_zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
#  endif
     &                         ad_t, ad_u, ad_v,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
     &                         ad_zeta)
!
!  Compute Hessian eigenvectors.
!
        ncname=ADM(ng)%name
!
        ROWS : DO rec=1,innLoop
!
!  Read gradient solution and load it into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     Lwrk, rec,                                   &
     &                     ndefADJ(ng), ADM(ng)%ncid, ncname,           &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     tl_t_obc, tl_u_obc, tl_v_obc,                &
#  endif
     &                     tl_ubar_obc, tl_vbar_obc,                    &
     &                     tl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
#  endif
     &                     tl_t, tl_u, tl_v,                            &
# else
     &                     tl_ubar, tl_vbar,                            &
# endif
     &                     tl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute Hessian eigenvectors:  (See NOTE above)
!
!    ad_var(Lold) = fac1 * ad_var(Lold) + fac2 * tl_var(Lwrk)
!
          fac1=1.0_r8
          fac2=ae_zv(rec,nvec)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lold, Lwrk, Lold, fac1, fac2,            &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         ad_t_obc, tl_t_obc,                      &
     &                         ad_u_obc, tl_u_obc,                      &
     &                         ad_v_obc, tl_v_obc,                      &
#  endif
     &                         ad_ubar_obc, tl_ubar_obc,                &
     &                         ad_vbar_obc, tl_vbar_obc,                &
     &                         ad_zeta_obc, tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         ad_ustr, tl_ustr,                        &
     &                         ad_vstr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         ad_tflux, tl_tflux,                      &
#  endif
     &                         ad_t, tl_t,                              &
     &                         ad_u, tl_u,                              &
     &                         ad_v, tl_v,                              &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                         ad_ubar, tl_ubar,                        &
     &                         ad_vbar, tl_vbar,                        &
#  endif
# else
     &                         ad_ubar, tl_ubar,                        &
     &                         ad_vbar, tl_vbar,                        &
# endif
     &                         ad_zeta, tl_zeta)
        END DO ROWS
!
!  Write eigenvectors into Hessian NetCDF.
!
        LwrtState2d(ng)=.TRUE.
        CALL wrt_hessian (ng, Lold, Lold)
        LwrtState2d(ng)=.FALSE.
        IF (exit_flag.ne.NoERRor) RETURN

      END DO COLUMNS
!
!-----------------------------------------------------------------------
!  Second, orthonormalize the converged Hessian vectors against each
!  other. Use tangent linear state arrays (index Lwrk) as temporary
!  storage.
!-----------------------------------------------------------------------
!
!  Use nl_var(1) as temporary storage since we need to preserve
!  ad_var(Lnew).
!
      ncname=HSS(ng)%name
      IF (Master) WRITE (stdout,30)
!
      DO nvec=1,innLoop
!
!  Read in just computed Hessian eigenvectors into adjoint state array
!  index Lold.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lold, nvec,                                    &
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize nonlinear state arrays index L1 with just read Hessian
!  vector in index Lold (initialize the summation): (See NOTE above)
!
!  Copy ad_var(Lold) into nl_var(L1).
!
        L1=1

        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lold, L1,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   nl_t_obc, ad_t_obc,                            &
     &                   nl_u_obc, ad_u_obc,                            &
     &                   nl_v_obc, ad_v_obc,                            &
#  endif
     &                   nl_ubar_obc, ad_ubar_obc,                      &
     &                   nl_vbar_obc, ad_vbar_obc,                      &
     &                   nl_zeta_obc, ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   nl_ustr, ad_ustr,                              &
     &                   nl_vstr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   nl_tflux, ad_tflux,                            &
#  endif
     &                   nl_t, ad_t,                                    &
     &                   nl_u, ad_u,                                    &
     &                   nl_v, ad_v,                                    &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                   nl_ubar, ad_ubar,                              &
     &                   nl_vbar, ad_vbar,                              &
#  endif
# else
     &                   nl_ubar, ad_ubar,                              &
     &                   nl_vbar, ad_vbar,                              &
# endif
     &                   nl_zeta, ad_zeta)
!
!  Orthogonalize Hessian eigenvectors against each other.
!
        DO rec=1,nvec-1
!
!  Read in gradient just computed Hessian eigenvectors into tangent
!  linear state array index Lwrk.
!
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj, LBij, UBij,              &
     &                     Lwrk, rec,                                   &
     &                     0, HSS(ng)%ncid, ncname,                     &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                     tl_t_obc, tl_u_obc, tl_v_obc,                &
#  endif
     &                     tl_ubar_obc, tl_vbar_obc,                    &
     &                     tl_zeta_obc,                                 &
# endif
# ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
#  endif
     &                     tl_t, tl_u, tl_v,                            &
# else
     &                     tl_ubar, tl_vbar,                            &
# endif
     &                     tl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product.
!
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        NstateVar(ng), dot(0:),                   &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        ad_t_obc(:,:,:,:,Lold,:),                 &
     &                        tl_t_obc(:,:,:,:,Lwrk,:),                 &
     &                        ad_u_obc(:,:,:,:,Lold),                   &
     &                        tl_u_obc(:,:,:,:,Lwrk),                   &
     &                        ad_v_obc(:,:,:,:,Lold),                   &
     &                        tl_v_obc(:,:,:,:,Lwrk),                   &
#  endif
     &                        ad_ubar_obc(:,:,:,Lold),                  &
     &                        tl_ubar_obc(:,:,:,Lwrk),                  &
     &                        ad_vbar_obc(:,:,:,Lold),                  &
     &                        tl_vbar_obc(:,:,:,Lwrk),                  &
     &                        ad_zeta_obc(:,:,:,Lold),                  &
     &                        tl_zeta_obc(:,:,:,Lwrk),                  &
# endif
# ifdef ADJUST_WSTRESS
     &                        ad_ustr(:,:,:,Lold), tl_ustr(:,:,:,Lwrk), &
     &                        ad_vstr(:,:,:,Lold), tl_vstr(:,:,:,Lwrk), &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        ad_tflux(:,:,:,Lold,:),                   &
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
#  endif
     &                        ad_t(:,:,:,Lold,:), tl_t(:,:,:,Lwrk,:),   &
     &                        ad_u(:,:,:,Lold), tl_u(:,:,:,Lwrk),       &
     &                        ad_v(:,:,:,Lold), tl_v(:,:,:,Lwrk),       &
# else
     &                        ad_ubar(:,:,Lold), tl_ubar(:,:,Lwrk),     &
     &                        ad_vbar(:,:,Lold), tl_vbar(:,:,Lwrk),     &
# endif
     &                        ad_zeta(:,:,Lold), tl_zeta(:,:,Lwrk))
!
!  Orthogonalize Hessian eigenvectors: (See NOTE above)
!
!    nl_var(L1) = fac1 * nl_var(L1) + fac2 * tl_var(Lwrk)
!
          fac1=1.0_r8
          fac2=-dot(0)

          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         L1, Lwrk, L1, fac1, fac2,                &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         nl_t_obc, tl_t_obc,                      &
     &                         nl_u_obc, tl_u_obc,                      &
     &                         nl_v_obc, tl_v_obc,                      &
#  endif
     &                         nl_ubar_obc, tl_ubar_obc,                &
     &                         nl_vbar_obc, tl_vbar_obc,                &
     &                         nl_zeta_obc, tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         nl_ustr, tl_ustr,                        &
     &                         nl_vstr, tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         nl_tflux, tl_tflux,                      &
#  endif
     &                         nl_t, tl_t,                              &
     &                         nl_u, tl_u,                              &
     &                         nl_v, tl_v,                              &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                         nl_ubar, tl_ubar,                        &
     &                         nl_vbar, tl_vbar,                        &
#  endif
# else
     &                         nl_ubar, tl_ubar,                        &
     &                         nl_vbar, tl_vbar,                        &
# endif
     &                         nl_zeta, tl_zeta)
        END DO
!
!  Compute normalization factor.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      nl_t_obc(:,:,:,:,L1,:),                     &
     &                      nl_t_obc(:,:,:,:,L1,:),                     &
     &                      nl_u_obc(:,:,:,:,L1),                       &
     &                      nl_u_obc(:,:,:,:,L1),                       &
     &                      nl_v_obc(:,:,:,:,L1),                       &
     &                      nl_v_obc(:,:,:,:,L1),                       &
#  endif
     &                      nl_ubar_obc(:,:,:,L1),                      &
     &                      nl_ubar_obc(:,:,:,L1),                      &
     &                      nl_vbar_obc(:,:,:,L1),                      &
     &                      nl_vbar_obc(:,:,:,L1),                      &
     &                      nl_zeta_obc(:,:,:,L1),                      &
     &                      nl_zeta_obc(:,:,:,L1),                      &
# endif
# ifdef ADJUST_WSTRESS
     &                      nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L1),       &
     &                      nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L1),       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      nl_tflux(:,:,:,L1,:),                       &
     &                      nl_tflux(:,:,:,L1,:),                       &
#  endif
     &                      nl_t(:,:,:,L1,:), nl_t(:,:,:,L1,:),         &
     &                      nl_u(:,:,:,L1), nl_u(:,:,:,L1),             &
     &                      nl_v(:,:,:,L1), nl_v(:,:,:,L1),             &
# else
     &                      nl_ubar(:,:,L1), nl_ubar(:,:,L1),           &
     &                      nl_vbar(:,:,L1), nl_vbar(:,:,L1),           &
# endif
     &                      nl_zeta(:,:,L1), nl_zeta(:,:,L1))
!
!  Normalize Hessian eigenvectors:
!
!    nl_var(L1) = fac * nl_var(L1)
!
        fac=1.0_r8/SQRT(dot(0))

        CALL state_scale (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    L1, L1, fac,                                  &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    nl_t_obc, nl_u_obc, nl_v_obc,                 &
#  endif
     &                    nl_ubar_obc, nl_vbar_obc,                     &
     &                    nl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    nl_ustr, nl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    nl_tflux,                                     &
#  endif
     &                    nl_t, nl_u, nl_v,                             &
# else
     &                    nl_ubar, nl_vbar,                             &
# endif
     &                    nl_zeta)
!
!  Copy nl_var(L1) into  ad_var(Lold).  See NOTE above.
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   L1, Lold,                                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, nl_t_obc,                            &
     &                   ad_u_obc, nl_u_obc,                            &
     &                   ad_v_obc, nl_v_obc,                            &
#  endif
     &                   ad_ubar_obc, nl_ubar_obc,                      &
     &                   ad_vbar_obc, nl_vbar_obc,                      &
     &                   ad_zeta_obc, nl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, nl_ustr,                              &
     &                   ad_vstr, nl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux, nl_tflux,                            &
#  endif
     &                   ad_t, nl_t,                                    &
     &                   ad_u, nl_u,                                    &
     &                   ad_v, nl_v,                                    &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                   ad_ubar, nl_ubar,                              &
     &                   ad_vbar, nl_vbar,                              &
#  endif
# else
     &                   ad_ubar, nl_ubar,                              &
     &                   ad_vbar, nl_vbar,                              &
# endif
     &                   ad_zeta, nl_zeta)
!
!  Write out converged Ritz eigenvalues and is associated accuracy.
!
        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'Ritz',          &
     &                        Ritz(nvec:), (/nvec/), (/1/),             &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, HSS(ng)%name, 'Ritz_error',    &
     &                        RitzErr(nvec:), (/nvec/), (/1/),          &
     &                        ncid = HSS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Replace record "nvec" of Hessian eigenvectors NetCDF with the
!  normalized value in adjoint state arrays at index Lold.
!
        HSS(ng)%Rindex=nvec-1
        LwrtState2d(ng)=.TRUE.
        CALL wrt_hessian (ng, Lold, Lold)
        LwrtState2d(ng)=.FALSE.
        IF (exit_flag.ne.NoERRor) RETURN

      END DO

  10  FORMAT (/,' Computing converged analysis error eofs...',/)
  20  FORMAT (a,'_',i3.3,'.nc')
  30  FORMAT (/,' Orthonormalizing converged analysis error eofs...',/)

      RETURN
      END SUBROUTINE posterior_eofs
#endif

      END MODULE posterior_mod
