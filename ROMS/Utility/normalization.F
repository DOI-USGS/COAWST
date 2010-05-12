#include "cppdefs.h"
      MODULE normalization_mod

#if defined FOUR_DVAR && defined CONVOLVE
# ifdef EW_PERIODIC
#  define IR_RANGE Istr-1,Iend+1
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JR_RANGE Jstr-1,Jend+1
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: normalization.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the background covariance, B, normalization   !
!  factors using the exact or approximated approach.  These  factors   !
!  ensure that the diagonal elements of B are equal to unity. Notice   !
!  that in applications with land/sea masking, it will produce large   !
!  changes in the covariance structures near the boundary.             !
!                                                                      !
!  The exact method is very expensive: the normalization factors are   !
!  computed by perturbing each model grid cell with a delta function   !
!  scaled by the area (2D factors) or volume (3D factors),  and then   !
!  convolving with the squared-root  adjoint  and  tangent diffusion   !
!  operators.                                                          !
!                                                                      !
!  The approximated method is cheaper: the normalization factors are   !
!  computed using the randomization approach of  Fisher and Courtier   !
!  (1995).  The factors are initialized with randon numbers having a   !
!  uniform distribution  (drawn from a normal distribution with zero   !
!  mean and unity variance).  Then, they scaled the inverse, squared   !
!  root cell area (2D factors) or volume (3D factors) and convoluted   !
!  with the squared-root adjoint and tangent  diffuse operators over   !
!  a specified number of iterations, Nrandom.                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fisher, M. and. P. Courtier, 1995:  Estimating the covariance     !
!      matrices of analysis and forecast error in variational data     !
!      assimilation, ECMWF Technical Memo N. 220, ECMWF,  Reading,     !
!      UK.                                                             !
!                                                                      !
!      www.ecmwf.int/publications/library/ecpublications/_pdf/tm/      !
!                                                001-300/tm220.pdf     !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteo.    !
!      Soc, 127, 1815-1846.                                            !
!                                                                      !
!======================================================================!
!
      USE mod_kinds

      implicit none

      PRIVATE
      PUBLIC :: normalization
      PUBLIC :: wrt_norm2d
      PUBLIC :: wrt_norm3d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE normalization (ng, tile, ifac)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
      USE mod_grid
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
# endif
      USE mod_fourdvar
      USE mod_mixing
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, ifac
!
!  Local variable declarations.
!
# include "tile.h"
!
!  Compute background error covariance normalization factors using
!  the very expensive exact method.
!
      IF (Nmethod(ng).eq.0) THEN
        CALL normalization_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp(ng), nnew(ng), ifac,              &
     &                           GRID(ng) % pm,                         &
     &                           GRID(ng) % om_r,                       &
     &                           GRID(ng) % om_u,                       &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % pn,                         &
     &                           GRID(ng) % on_r,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % on_v,                       &
     &                           GRID(ng) % pmon_p,                     &
     &                           GRID(ng) % pmon_r,                     &
     &                           GRID(ng) % pmon_u,                     &
     &                           GRID(ng) % pnom_p,                     &
     &                           GRID(ng) % pnom_r,                     &
     &                           GRID(ng) % pnom_v,                     &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % vmask,                      &
# endif
# ifdef SOLVE3D
     &                           GRID(ng) % h,                          &
#  ifdef ICESHELF
     &                           GRID(ng) % zice,                       &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                           SEDBED(ng) % bed_thick,                &
#  endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % z_r,                        &
     &                           GRID(ng) % z_w,                        &
# endif
     &                           MIXING(ng) % Kh,                       &
# ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           BOUNDARY(ng) % b_t_obc,                &
     &                           BOUNDARY(ng) % b_u_obc,                &
     &                           BOUNDARY(ng) % b_v_obc,                &
#  endif
     &                           BOUNDARY(ng) % b_ubar_obc,             &
     &                           BOUNDARY(ng) % b_vbar_obc,             &
     &                           BOUNDARY(ng) % b_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                           FORCES(ng) % b_stflx,                  &
# endif
# ifdef SOLVE3D
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
# endif
     &                           OCEAN(ng) % b_zeta,                    &
     &                           OCEAN(ng) % b_ubar,                    &
     &                           OCEAN(ng) % b_vbar)
!
!  Compute background error covariance normalization factors using
!  the approximated randomization method.
!
      ELSE IF (Nmethod(ng).eq.1) THEN
        CALL randomization_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp(ng), nnew(ng), ifac,              &
     &                           GRID(ng) % pm,                         &
     &                           GRID(ng) % om_r,                       &
     &                           GRID(ng) % om_u,                       &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % pn,                         &
     &                           GRID(ng) % on_r,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % on_v,                       &
     &                           GRID(ng) % pmon_p,                     &
     &                           GRID(ng) % pmon_r,                     &
     &                           GRID(ng) % pmon_u,                     &
     &                           GRID(ng) % pnom_p,                     &
     &                           GRID(ng) % pnom_r,                     &
     &                           GRID(ng) % pnom_v,                     &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % vmask,                      &
# endif
# ifdef SOLVE3D
     &                           GRID(ng) % h,                          &
#  ifdef ICESHELF
     &                           GRID(ng) % zice,                       &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                           SEDBED(ng) % bed_thick,                &
#  endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % z_r,                        &
     &                           GRID(ng) % z_w,                        &
# endif
     &                           MIXING(ng) % Kh,                       &
# ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           BOUNDARY(ng) % b_t_obc,                &
     &                           BOUNDARY(ng) % b_u_obc,                &
     &                           BOUNDARY(ng) % b_v_obc,                &
#  endif
     &                           BOUNDARY(ng) % b_ubar_obc,             &
     &                           BOUNDARY(ng) % b_vbar_obc,             &
     &                           BOUNDARY(ng) % b_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                           FORCES(ng) % b_stflx,                  &
# endif
# ifdef SOLVE3D
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
# endif
     &                           OCEAN(ng) % b_zeta,                    &
     &                           OCEAN(ng) % b_ubar,                    &
     &                           OCEAN(ng) % b_vbar)
      END IF

      RETURN
      END SUBROUTINE normalization
!
!***********************************************************************
      SUBROUTINE normalization_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nstp, nnew, ifac,                  &
     &                               pm, om_r, om_u, om_v,              &
     &                               pn, on_r, on_u, on_v,              &
     &                               pmon_p, pmon_r, pmon_u,            &
     &                               pnom_p, pnom_r, pnom_v,            &
# ifdef MASKING
     &                               rmask, pmask, umask, vmask,        &
# endif
# ifdef SOLVE3D
     &                               h,                                 &
#  ifdef ICESHELF
     &                               zice,                              &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                               bed_thick,                         &
#  endif
     &                               Hz, z_r, z_w,                      &
# endif
     &                               Kh,                                &
# ifdef SOLVE3D
     &                               Kv,                                &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                               VnormRobc, VnormUobc, VnormVobc,   &
#  endif
     &                               HnormRobc, HnormUobc, HnormVobc,   &
# endif
# ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                               HnormSTF,                          &
# endif
# ifdef SOLVE3D
     &                               VnormR, VnormU, VnormV,            &
# endif
     &                               HnormR, HnormU, HnormV)
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
      USE bc_2d_mod
      USE ad_conv_2d_mod
      USE tl_conv_2d_mod
# ifdef SOLVE3D
      USE bc_3d_mod
      USE ad_conv_3d_mod
      USE tl_conv_3d_mod
# endif
# ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod
      USE ad_conv_bry2d_mod
      USE tl_conv_bry2d_mod
#  ifdef SOLVE3D
      USE bc_bry3d_mod
      USE ad_conv_bry3d_mod
      USE tl_conv_bry3d_mod
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE set_depth_mod
# ifdef DISTRIBUTE
#  ifdef ADJUST_BOUNDARY
      USE distribute_mod, ONLY : mp_collect
#  endif
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, ifac
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Kh(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent(out) :: VnormUobc(LBij:,:,:)
      real(r8), intent(out) :: VnormVobc(LBij:,:,:)
#   endif
      real(r8), intent(out) :: HnormRobc(LBij:,:)
      real(r8), intent(out) :: HnormUobc(LBij:,:)
      real(r8), intent(out) :: HnormVobc(LBij:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
#  endif
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick0(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(out) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent(out) :: VnormVobc(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent(out) :: HnormRobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormUobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormVobc(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng),NSA)
#  endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj,NSA)
#  ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
# ifdef SOLVE3D
      logical :: Ldiffer, Lsame
# endif
# ifdef ADJUST_BOUNDARY
      logical :: bounded
      logical :: Lconvolve(4)
# endif
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ic, ifile, is, j, jc, rec
# ifdef SOLVE3D
      integer :: UBt, itrc, k, kc, ntrc
# endif
# ifdef ADJUST_BOUNDARY
      integer :: Bmin, Bmax, IJlen, IJKlen
      integer :: bry, ib, ifield, kb, tindex

      real(r8), parameter :: Aspv = 0.0_r8
# endif
      real(r8) :: cff, compute

      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
# endif
# ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
      real(r8), dimension(LBij:UBij) :: HscaleB
#  ifdef SOLVE3D
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
      real(r8), dimension(LBij:UBij,1:N(ng)) :: VscaleB
#  endif
# endif

      character (len=40) :: Text
      character (len=80) :: ncname

# include "set_bounds.h"
!
      SourceFile='normalization.F, normalization_tile'

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
#  endif
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
# endif
!
!-----------------------------------------------------------------------
!  Compute initial conditions and model error covariance, B,
!  normalization factors using the exact method. It involves
!  computing the filter variance (convolution) at each point
!  independenly.  That is, each point is perturbed with a delta
!  function, scaled by the inverse squared root of the area (2D)
!  or volume (3D), and then convoluted.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,10)

      FILE_LOOP : DO ifile=1,NSA

        IF (LwrtNRM(ifile,ng)) THEN
          IF (ifile.eq.1) THEN
            Text='initial conditions'
          ELSE IF (ifile.eq.2) THEN
            Text='model'
          END IF
!
!  Set time record index to write in normalization NetCDF file.
!
          ncname=NRMname(ifile,ng)
          tNRMindx(ifile,ng)=tNRMindx(ifile,ng)+1
          NrecNRM(ifile,ng)=NrecNRM(ifile,ng)+1
!
!  Write out model time (s).
!
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,idtime),      &
     &                          time(ng),                               &
     &                          start = (/tNRMindx(ifile,ng)/),         &
     &                          total = (/1/),                          &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,idtime,ng))
          IF (exit_flag.ne.NoError) RETURN
!
!  2D norm at RHO-points.
!
          IF (Cnorm(ifile,isFsur)) THEN
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at RHO-points'
              CALL my_flush (stdout)
            END IF
            DO j=JR_RANGE
              DO i=IR_RANGE
                Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
              END DO
            END DO
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
# ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (rmask(ic,jc).gt.0) compute=1.0_r8
                END IF
#  ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                compute=1.0_r8
# endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=Hscale(ic,jc)
                  END IF
                  CALL ad_conv_r2d_tile (ng, tile, iADM,                &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,isFsur)/ifac,    &
     &                                   DTsizeH(ifile,isFsur),         &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_u, pnom_v,        &
# ifdef MASKING
     &                                   rmask, umask, vmask,           &
# endif
     &                                   A2d)
                  DO j=JR_RANGE
                    DO i=IR_RANGE
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
                  CALL tl_conv_r2d_tile (ng, tile, iTLM,                &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,isFsur)/ifac,    &
     &                                   DTsizeH(ifile,isFsur),         &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_u, pnom_v,        &
# ifdef MASKING
     &                                   rmask, umask, vmask,           &
# endif
     &                                   A2d)
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    cff=1.0_r8/SQRT(A2d(ic,jc))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormR(ic,jc,ifile)=cff
                END IF
              END DO
            END DO
            CALL bc_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormR(:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          HnormR(:,:,ifile))
# endif
            CALL wrt_norm2d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, idFsur,                &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idFsur,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
# ifdef MASKING
     &                       rmask,                                     &
# endif
     &                       HnormR(:,:,ifile))
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  2D norm at U-points.
!
          IF (Cnorm(ifile,isUbar)) THEN
# ifdef EW_PERIODIC
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
# else
            Imin=2
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
# endif
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at   U-points'
              CALL my_flush (stdout)
            END IF
            DO j=JU_RANGE
              DO i=IU_RANGE
                Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
              END DO
            END DO
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
# ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (umask(ic,jc).gt.0) compute=1.0_r8
                END IF
#  ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                compute=1.0_r8
# endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=Hscale(ic,jc)
                  END IF
                  CALL ad_conv_u2d_tile (ng, tile, iADM,                &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,isUbar)/ifac,    &
     &                                   DTsizeH(ifile,isUbar),         &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_r, pnom_p,        &
# ifdef MASKING
     &                                   umask, pmask,                  &
# endif
     &                                   A2d)
                  DO j=JU_RANGE
                    DO i=IU_RANGE
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
                  CALL tl_conv_u2d_tile (ng, tile, iTLM,                &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,isUbar)/ifac,    &
     &                                   DTsizeH(ifile,isUbar),         &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_r, pnom_p,        &
# ifdef MASKING
     &                                   umask, pmask,                  &
# endif
     &                                   A2d)
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    cff=1.0_r8/SQRT(A2d(ic,jc))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormU(ic,jc,ifile)=cff
                END IF
              END DO
            END DO
            CALL bc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormU(:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          HnormU(:,:,ifile))
# endif
            CALL wrt_norm2d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, idUbar,                &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idUbar,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
# ifdef MASKING
     &                       umask,                                     &
# endif
     &                       HnormU(:,:,ifile))
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  2D norm at V-points.
!
          IF (Cnorm(ifile,isVbar)) THEN
# ifdef NS_PERIODIC
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
# else
            Imin=1
            Imax=Lm(ng)
            Jmin=2
            Jmax=Mm(ng)
# endif
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at   V-points'
              CALL my_flush (stdout)
            END IF
            DO j=JV_RANGE
              DO i=IV_RANGE
                Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
              END DO
            END DO
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
# ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (vmask(ic,jc).gt.0) compute=1.0_r8
                END IF
#  ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                compute=1.0_r8
# endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=Hscale(ic,jc)
                  END IF
                  CALL ad_conv_v2d_tile (ng, tile, iADM,                &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,isVbar)/ifac,    &
     &                                   DTsizeH(ifile,isVbar),         &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_p, pnom_r,        &
# ifdef MASKING
     &                                   vmask, pmask,                  &
# endif
     &                                   A2d)
                  DO j=JV_RANGE
                    DO i=IV_RANGE
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
                  CALL tl_conv_v2d_tile (ng, tile, iTLM,                &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,isVbar)/ifac,    &
     &                                   DTsizeH(ifile,isVbar),         &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_p, pnom_r,        &
# ifdef MASKING
     &                                   vmask, pmask,                  &
# endif
     &                                   A2d)
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    cff=1.0_r8/SQRT(A2d(ic,jc))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormV(ic,jc,ifile)=cff
                END IF
              END DO
            END DO
            CALL bc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormV(:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          HnormV(:,:,ifile))
# endif
            CALL wrt_norm2d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, idVbar,                &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idVbar,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
# ifdef MASKING
     &                       vmask,                                     &
# endif
     &                       HnormV(:,:,ifile))
            IF (exit_flag.ne.NoError) RETURN
          END IF

# ifdef SOLVE3D
!
!  3D norm U-points.
!
          IF (Cnorm(ifile,isUvel)) THEN
# ifdef EW_PERIODIC
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
# else
            Imin=2
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
# endif
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at   U-points'
              CALL my_flush (stdout)
            END IF
            DO j=JU_RANGE
              DO i=IU_RANGE
                cff=om_u(i,j)*on_u(i,j)*0.5_r8
                DO k=1,N(ng)
                  Vscale(i,j,k)=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
                END DO
              END DO
            END DO
            DO kc=1,N(ng)
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#  ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (umask(ic,jc).gt.0) compute=1.0_r8
                  END IF
#   ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                  compute=1.0_r8
#  endif
                  IF (compute.gt.0.0_r8) THEN
                    DO k=1,N(ng)
                      DO j=LBj,UBj
                        DO i=LBi,UBi
                          A3d(i,j,k)=0.0_r8
                        END DO
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A3d(ic,jc,kc)=Vscale(ic,jc,kc)
                    END IF
                    CALL ad_conv_u3d_tile (ng, tile, iADM,              &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     1, N(ng),                    &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHsteps(ifile,isUvel)/ifac,  &
     &                                     NVsteps(ifile,isUvel)/ifac,  &
     &                                     DTsizeH(ifile,isUvel),       &
     &                                     DTsizeV(ifile,isUvel),       &
     &                                     Kh, Kv,                      &
     &                                     pm, pn, pmon_r, pnom_p,      &
#  ifdef MASKING
     &                                     umask, pmask,                &
#  endif
     &                                     Hz, z_r,                     &
     &                                     A3d)
                    DO k=1,N(ng)
                      DO j=JU_RANGE
                        DO i=IU_RANGE
                          A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                        END DO
                      END DO
                    END DO
                    CALL tl_conv_u3d_tile (ng, tile, iTLM,              &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     1, N(ng),                    &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHsteps(ifile,isUvel)/ifac,  &
     &                                     NVsteps(ifile,isUvel)/ifac,  &
     &                                     DTsizeH(ifile,isUvel),       &
     &                                     DTsizeV(ifile,isUvel),       &
     &                                     Kh, Kv,                      &
     &                                     pm, pn, pmon_r, pnom_p,      &
#  ifdef MASKING
     &                                     umask, pmask,                &
#  endif
     &                                     Hz, z_r,                     &
     &                                     A3d)
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      cff=1.0_r8/SQRT(A3d(ic,jc,kc))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    VnormU(ic,jc,kc,ifile)=cff
                  END IF
                END DO
              END DO
            END DO
            CALL bc_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        VnormU(:,:,:,ifile))
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          VnormU(:,:,:,ifile))
#  endif
            CALL wrt_norm3d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), idUvel,      &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idUvel,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
#  ifdef MASKING
     &                       umask,                                     &
#  endif
     &                       VnormU(:,:,:,ifile))
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  3D norm at V-points.
!
          IF (Cnorm(ifile,isVvel)) THEN
# ifdef NS_PERIODIC
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
# else
            Imin=1
            Imax=Lm(ng)
            Jmin=2
            Jmax=Mm(ng)
# endif
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at   V-points'
              CALL my_flush (stdout)
            END IF
            DO j=JV_RANGE
              DO i=IV_RANGE
                cff=om_v(i,j)*on_v(i,j)*0.5_r8
                DO k=1,N(ng)
                  Vscale(i,j,k)=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
                END DO
              END DO
            END DO
            DO kc=1,N(ng)
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#  ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (vmask(ic,jc).gt.0) compute=1.0_r8
                  END IF
#   ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                  compute=1.0_r8
#  endif
                  IF (compute.gt.0.0_r8) THEN
                    DO k=1,N(ng)
                      DO j=LBj,UBj
                        DO i=LBi,UBi
                          A3d(i,j,k)=0.0_r8
                        END DO
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A3d(ic,jc,kc)=Vscale(ic,jc,kc)
                    END IF
                    CALL ad_conv_v3d_tile (ng, tile, iADM,              &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     1, N(ng),                    &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHsteps(ifile,isVvel)/ifac,  &
     &                                     NVsteps(ifile,isVvel)/ifac,  &
     &                                     DTsizeH(ifile,isVvel),       &
     &                                     DTsizeV(ifile,isVvel),       &
     &                                     Kh, Kv,                      &
     &                                     pm, pn, pmon_p, pnom_r,      &
#  ifdef MASKING
     &                                     vmask, pmask,                &
#  endif
     &                                     Hz, z_r,                     &
     &                                     A3d)
                    DO k=1,N(ng)
                      DO j=JV_RANGE
                        DO i=IV_RANGE
                          A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                        END DO
                      END DO
                    END DO
                    CALL tl_conv_v3d_tile (ng, tile, iTLM,              &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     1, N(ng),                    &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHsteps(ifile,isVvel)/ifac,  &
     &                                     NVsteps(ifile,isVvel)/ifac,  &
     &                                     DTsizeH(ifile,isVvel),       &
     &                                     DTsizeV(ifile,isVvel),       &
     &                                     Kh, Kv,                      &
     &                                     pm, pn, pmon_p, pnom_r,      &
#  ifdef MASKING
     &                                     vmask, pmask,                &
#  endif
     &                                     Hz, z_r,                     &
     &                                     A3d)
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      cff=1.0_r8/SQRT(A3d(ic,jc,kc))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    VnormV(ic,jc,kc,ifile)=cff
                  END IF
                END DO
              END DO
            END DO
            CALL bc_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        VnormV(:,:,:,ifile))
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          VnormV(:,:,:,ifile))
#  endif
            CALL wrt_norm3d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), idVvel,      &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idVvel,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
#  ifdef MASKING
     &                       vmask,                                     &
#  endif
     &                       VnormV(:,:,:,ifile))
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  3D norm at RHO-points.
!
          IF (Master) THEN
            Lsame=.FALSE.
            DO itrc=1,NT(ng)
              is=isTvar(itrc)
              IF (Cnorm(ifile,is)) Lsame=.TRUE.
            END DO
            IF (Lsame) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at RHO-points'
              CALL my_flush (stdout)
            END IF
          END IF
!
!  Check if the decorrelation scales for all the tracers are different.
!  If not, just compute the normalization factors for the first tracer
!  and assign the same value to the rest.  Recall that this computation
!  is very expensive.
!
          Ldiffer=.FALSE.
          Imin=1
          Imax=Lm(ng)
          Jmin=1
          Jmax=Mm(ng)
          DO itrc=2,NT(ng)
            IF ((Hdecay(ifile,isTvar(itrc  ),ng).ne.                    &
     &           Hdecay(ifile,isTvar(itrc-1),ng)).or.                   &
     &          (Vdecay(ifile,isTvar(itrc  ),ng).ne.                    &
     &           Vdecay(ifile,isTvar(itrc-1),ng))) THEN
              Ldiffer=.TRUE.
            END IF
          END DO
          IF (.not.Ldiffer) THEN
            Lsame=.TRUE.
            UBt=1
          ELSE
            Lsame=.FALSE.
            UBt=NT(ng)
          END IF
!
          DO j=JR_RANGE
            DO i=IR_RANGE
              cff=om_r(i,j)*on_r(i,j)
              DO k=1,N(ng)
                Vscale(i,j,k)=1.0_r8/SQRT(cff*Hz(i,j,k))
              END DO
            END DO
          END DO
          DO itrc=1,UBt
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              DO kc=1,N(ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
#  ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (rmask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#   ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                    compute=1.0_r8
#  endif
                    IF (compute.gt.0.0_r8) THEN
                      DO k=1,N(ng)
                        DO j=LBj,UBj
                          DO i=LBi,UBi
                            A3d(i,j,k)=0.0_r8
                          END DO
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A3d(ic,jc,kc)=Vscale(ic,jc,kc)
                      END IF
                      CALL ad_conv_r3d_tile (ng, tile, iADM,            &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHsteps(ifile,is)/ifac,    &
     &                                       NVsteps(ifile,is)/ifac,    &
     &                                       DTsizeH(ifile,is),         &
     &                                       DTsizeV(ifile,is),         &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_u, pnom_v,    &
#  ifdef MASKING
     &                                       rmask, umask, vmask,       &
#  endif
     &                                       Hz, z_r,                   &
     &                                       A3d)
                      DO k=1,N(ng)
                        DO j=JR_RANGE
                          DO i=IR_RANGE
                            A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                          END DO
                        END DO
                      END DO
                      CALL tl_conv_r3d_tile (ng, tile, iTLM,            &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHsteps(ifile,is)/ifac,    &
     &                                       NVsteps(ifile,is)/ifac,    &
     &                                       DTsizeH(ifile,is),         &
     &                                       DTsizeV(ifile,is),         &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_u, pnom_v,    &
#  ifdef MASKING
     &                                       rmask, umask, vmask,       &
#  endif
     &                                       Hz, z_r,                   &
     &                                       A3d)
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        cff=1.0_r8/SQRT(A3d(ic,jc,kc))
                      END IF
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (Lsame) THEN
                        DO ntrc=1,NT(ng)
                          VnormR(ic,jc,kc,ifile,ntrc)=cff
                        END DO
                      ELSE
                        VnormR(ic,jc,kc,ifile,itrc)=cff
                      END IF
                    END IF
                  END DO
                END DO
              END DO
            END IF
          END DO
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              CALL bc_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormR(:,:,:,ifile,itrc))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints, EWperiodic, NSperiodic, &
     &                            VnormR(:,:,:,ifile,itrc))
#  endif
              CALL wrt_norm3d (ng, tile, iTLM, ncname,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         idTvar(itrc), ncNRMid(ifile,ng),         &
     &                         nrmVid(ifile,idTvar(itrc),ng),           &
     &                         tNRMindx(ifile,ng),                      &
#  ifdef MASKING
     &                         rmask,                                   &
#  endif
     &                         VnormR(:,:,:,ifile,itrc))
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END DO
# endif
        END IF
      END DO FILE_LOOP

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Compute open boundaries error covariance, B, normalization factors
!  using the exact method.
!-----------------------------------------------------------------------
!
      ifile=3
      IF (LwrtNRM(ifile,ng)) THEN
        Text='boundary conditions'
        IJlen=UBij-LBij+1
#  ifdef SOLVE3D
        IJKlen=IJlen*N(ng)
#  endif
        Lconvolve(iwest )=WESTERN_EDGE
        Lconvolve(ieast )=EASTERN_EDGE
        Lconvolve(isouth)=SOUTHERN_EDGE
        Lconvolve(inorth)=NORTHERN_EDGE
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRMname(ifile,ng)
        tNRMindx(ifile,ng)=tNRMindx(ifile,ng)+1
        NrecNRM(ifile,ng)=NrecNRM(ifile,ng)+1
        tindex=tNRMindx(ifile,ng)
!
!  Write out model time (s).
!
        CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,idtime),        &
     &                        time(ng),                                 &
     &                        start = (/tindex/),                       &
     &                        total = (/1/),                            &
     &                        ncid = ncNRMid(ifile,ng),                 &
     &                        varid = nrmVid(ifile,idtime,ng))
        IF (exit_flag.ne.NoError) RETURN
!
!  2D boundary norm at RHO-points.
!
        HnormRobc=Aspv

        IF (Master.and.ANY(CnormB(isFsur,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at RHO-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          HscaleB=0.0_r8
          IF (CnormB(isFsur,bry)) THEN
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,r2dvar)
              Bmin=1
              Bmax=Mm(ng)
              IF (Lconvolve(bry)) THEN
                DO j=JR_RANGE
                  HscaleB(j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,r2dvar)
              Bmin=1
              Bmax=Lm(ng)
              IF (Lconvolve(bry)) THEN
                DO i=IR_RANGE
                  HscaleB(i)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
                END DO
              END IF
            END IF
            DO ib=Bmin,Bmax
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                bounded=Lconvolve(bry).and.                             &
     &                  ((Jstr.le.ib).and.(ib.le.Jend))
                j=ib
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                bounded=Lconvolve(bry).and.                             &
     &                  ((Istr.le.ib).and.(ib.le.Iend))
                i=ib
              END IF
#  ifdef MASKING
              IF (bounded) THEN
                compute=rmask(i,j)
              ELSE
                compute=0.0_r8
              END IF
#   ifdef DISTRIBUTE
              CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
              compute=1.0_r8
#  endif
              IF (compute.gt.0.0_r8) THEN
                B2d=0.0_r8
                IF (bounded) THEN
                  B2d(ib)=HscaleB(ib)
                END IF
                CALL ad_conv_r2d_bry_tile (ng, tile, iADM, bry,         &
     &                                     BOUNDS(ng)%edge(:,r2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,isFsur)/ifac,   &
     &                                     DTsizeHB(bry,isFsur),        &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_u, pnom_v,      &
#  ifdef MASKING
     &                                     rmask, umask, vmask,         &
#  endif
     &                                     B2d)
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO j=JR_RANGE
                    B2d(j)=B2d(j)*HscaleB(j)
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO i=IR_RANGE
                    B2d(i)=B2d(i)*HscaleB(i)
                  END DO
                END IF
                CALL tl_conv_r2d_bry_tile (ng, tile, iTLM, bry,         &
     &                                     BOUNDS(ng)%edge(:,r2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,isFsur)/ifac,   &
     &                                     DTsizeHB(bry,isFsur),        &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_u, pnom_v,      &
#  ifdef MASKING
     &                                     rmask, umask, vmask,         &
#  endif
     &                                     B2d)
                IF (bounded) THEN
                  cff=1.0_r8/SQRT(B2d(ib))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (bounded) THEN
                HnormRobc(ib,bry)=cff
              END IF
            END DO
            CALL bc_r2d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij,                           &
     &                            HnormRobc(:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormRobc(LBij:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isFsur,:))) THEN
          ifield=idSbry(isFsur)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          HnormRobc(LBij:,:),                     &
     &                          start = (/1,1,tindex/),                 &
     &                          total = (/IJlen,4,1/),                  &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  2D boundary norm at U-points.
!
        HnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   U-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          HscaleB=0.0_r8
          IF (CnormB(isUbar,bry)) THEN
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,u2dvar)
              Bmin=1
              Bmax=Mm(ng)
              IF (Lconvolve(bry)) THEN
                DO j=JU_RANGE
                  HscaleB(j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,u2dvar)
#  ifdef EW_PERIODIC
              Bmin=1
              Bmax=Lm(ng)
#  else
              Bmin=2
              Bmax=Lm(ng)
#  endif
              IF (Lconvolve(bry)) THEN
                DO i=IU_RANGE
                  HscaleB(i)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
                END DO
              END IF
            END IF
            DO ib=Bmin,Bmax
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                bounded=Lconvolve(bry).and.                             &
     &                  ((Jstr.le.ib).and.(ib.le.Jend))
                j=ib
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                bounded=Lconvolve(bry).and.                             &
     &                  ((Istr.le.ib).and.(ib.le.Iend))
                i=ib
              END IF
#  ifdef MASKING
              IF (bounded) THEN
                compute=umask(i,j)
              ELSE
                compute=0.0_r8
              END IF
#   ifdef DISTRIBUTE
              CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
              compute=1.0_r8
#  endif
              IF (compute.gt.0.0_r8) THEN
                B2d=0.0_r8
                IF (bounded) THEN
                  B2d(ib)=HscaleB(ib)
                END IF
                CALL ad_conv_u2d_bry_tile (ng, tile, iADM, bry,         &
     &                                     BOUNDS(ng)%edge(:,u2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,isUbar)/ifac,   &
     &                                     DTsizeHB(bry,isUbar),        &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_r, pnom_p,      &
#  ifdef MASKING
     &                                     umask, pmask,                &
#  endif
     &                                     B2d)
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO j=JU_RANGE
                    B2d(j)=B2d(j)*HscaleB(j)
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO i=IU_RANGE
                    B2d(i)=B2d(i)*HscaleB(i)
                  END DO
                END IF
                CALL tl_conv_u2d_bry_tile (ng, tile, iTLM, bry,         &
     &                                     BOUNDS(ng)%edge(:,u2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,isUbar)/ifac,   &
     &                                     DTsizeHB(bry,isUbar),        &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_r, pnom_p,      &
#  ifdef MASKING
     &                                     umask, pmask,                &
#  endif
     &                                     B2d)
                IF (bounded) THEN
                  cff=1.0_r8/SQRT(B2d(ib))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (bounded) THEN
                HnormUobc(ib,bry)=cff
              END IF
            END DO
            CALL bc_u2d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij,                           &
     &                            HnormUobc(:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormUobc(LBij:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isUbar,:))) THEN
          ifield=idSbry(isUbar)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          HnormUobc(LBij:,:),                     &
     &                          start = (/1,1,tindex/),                 &
     &                          total = (/IJlen,4,1/),                  &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  2D boundary norm at V-points.
!
        HnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   V-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          HscaleB=0.0_r8
          IF (CnormB(isVbar,bry)) THEN
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,v2dvar)
#  ifdef NS_PERIODIC
              Bmin=1
              Bmax=Mm(ng)
#  else
              Bmin=2
              Bmax=Mm(ng)
#  endif
              IF (Lconvolve(bry)) THEN
                DO j=JR_RANGE
                  HscaleB(j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,v2dvar)
              Bmin=1
              Bmax=Lm(ng)
              IF (Lconvolve(bry)) THEN
                DO i=IR_RANGE
                  HscaleB(i)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
                END DO
              END IF
            END IF
            DO ib=Bmin,Bmax
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                bounded=Lconvolve(bry).and.                             &
     &                  ((Jstr.le.ib).and.(ib.le.Jend))
                j=ib
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                bounded=Lconvolve(bry).and.                             &
     &                  ((Istr.le.ib).and.(ib.le.Iend))
                i=ib
              END IF
#  ifdef MASKING
              IF (bounded) THEN
                compute=vmask(i,j)
              ELSE
                compute=0.0_r8
              END IF
#   ifdef DISTRIBUTE
              CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
              compute=1.0_r8
#  endif
              IF (compute.gt.0.0_r8) THEN
                B2d=0.0_r8
                IF (bounded) THEN
                  B2d(ib)=HscaleB(ib)
                END IF
                CALL ad_conv_v2d_bry_tile (ng, tile, iADM, bry,         &
     &                                     BOUNDS(ng)%edge(:,v2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,isVbar)/ifac,   &
     &                                     DTsizeHB(bry,isVbar),        &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_p, pnom_r,      &
#  ifdef MASKING
     &                                     vmask, pmask,                &
#  endif
     &                                     B2d)
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO j=JV_RANGE
                    B2d(j)=B2d(j)*HscaleB(j)
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO i=IV_RANGE
                    B2d(i)=B2d(i)*HscaleB(i)
                  END DO
                END IF
                CALL tl_conv_v2d_bry_tile (ng, tile, iTLM, bry,         &
     &                                     BOUNDS(ng)%edge(:,v2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,isVbar)/ifac,   &
     &                                     DTsizeHB(bry,isVbar),        &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_p, pnom_r,      &
#  ifdef MASKING
     &                                     vmask, pmask,                &
#  endif
     &                                     B2d)
                IF (bounded) THEN
                  cff=1.0_r8/SQRT(B2d(ib))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (bounded) THEN
                HnormVobc(ib,bry)=cff
              END IF
            END DO
            CALL bc_v2d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij,                           &
     &                            HnormVobc(:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormVobc(LBij:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isVbar,:))) THEN
          ifield=idSbry(isVbar)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          HnormVobc(LBij:,:),                     &
     &                          start = (/1,1,tindex/),                 &
     &                          total = (/IJlen,4,1/),                  &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  ifdef SOLVE3D
!
!  3D boundary norm at U-points.
!
        VnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   U-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          VscaleB=0.0_r8
          IF (CnormB(isUvel,bry)) THEN
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,u2dvar)
              Bmin=1
              Bmax=Mm(ng)
              IF (Lconvolve(bry)) THEN
                DO j=JU_RANGE
                  cff=om_u(i,j)*on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(j,k)=1.0_r8/                                &
     &                           SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,u2dvar)
#   ifdef EW_PERIODIC
              Bmin=1
              Bmax=Lm(ng)
#   else
              Bmin=2
              Bmax=Lm(ng)
#   endif
              IF (Lconvolve(bry)) THEN
                DO i=IU_RANGE
                  cff=om_u(i,j)*on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(i,k)=1.0_r8/                                &
     &                           SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            END IF
            DO kb=1,N(ng)
              DO ib=Bmin,Bmax
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  bounded=Lconvolve(bry).and.                           &
     &                    ((Jstr.le.ib).and.(ib.le.Jend))
                  j=ib
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  bounded=Lconvolve(bry).and.                           &
     &                    ((Istr.le.ib).and.(ib.le.Iend))
                  i=ib
                END IF
#   ifdef MASKING
                IF (bounded) THEN
                  compute=umask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#    ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                compute=1.0_r8
#   endif
                IF (compute.gt.0.0_r8) THEN
                  B3d=0.0_r8
                  IF (bounded) THEN
                    B3d(ib,kb)=VscaleB(ib,kb)
                  END IF
                  CALL ad_conv_u3d_bry_tile (ng, tile, iADM, bry,       &
     &                                       BOUNDS(ng)%edge(:,u2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(bry,isUvel)/ifac, &
     &                                       NVstepsB(bry,isUvel)/ifac, &
     &                                       DTsizeHB(bry,isUvel),      &
     &                                       DTsizeVB(bry,isUvel),      &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_r, pnom_p,    &
#   ifdef MASKING
     &                                       umask, pmask,              &
#   endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
                  IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                    DO k=1,N(ng)
                      DO j=JU_RANGE
                        B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                      END DO
                    END DO
                  ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                    DO k=1,N(ng)
                      DO i=IU_RANGE
                        B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                      END DO
                    END DO
                  END IF
                  CALL tl_conv_u3d_bry_tile (ng, tile, iTLM, bry,       &
     &                                       BOUNDS(ng)%edge(:,u2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(bry,isUvel)/ifac, &
     &                                       NVstepsB(bry,isUvel)/ifac, &
     &                                       DTsizeHB(bry,isUvel),      &
     &                                       DTsizeVB(bry,isUvel),      &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_r, pnom_p,    &
#   ifdef MASKING
     &                                       umask, pmask,              &
#   endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B3d(ib,kb))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  VnormUobc(ib,kb,bry)=cff
                END IF
              END DO
            END DO
            CALL bc_u3d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormUobc(:,:,bry))
#   ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJKlen, Aspv,                    &
     &                       VnormUobc(LBij:,:,bry))
#   endif
          END IF
        END DO
        IF (ANY(CnormB(isUvel,:))) THEN
          ifield=idSbry(isUvel)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          VnormUobc(LBij:,:,:),                   &
     &                          start = (/1,1,1,tindex/),               &
     &                          total = (/IJlen,N(ng),4,1/),            &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  3D boundary norm at V-points.
!
        VnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   V-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          VscaleB=0.0_r8
          IF (CnormB(isVvel,bry)) THEN
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,v2dvar)
#   ifdef NS_PERIODIC
              Bmin=1
              Bmax=Mm(ng)
#   else
              Bmin=2
              Bmax=Mm(ng)
#   endif
              IF (Lconvolve(bry)) THEN
                DO j=JV_RANGE
                  cff=om_v(i,j)*on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(j,k)=1.0_r8/                                &
     &                           SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,v2dvar)
              Bmin=1
              Bmax=Lm(ng)
              IF (Lconvolve(bry)) THEN
                DO i=IV_RANGE
                  cff=om_v(i,j)*on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(i,k)=1.0_r8/                                &
     &                           SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            END IF
            DO kb=1,N(ng)
              DO ib=Bmin,Bmax
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  bounded=Lconvolve(bry).and.                           &
     &                    ((Jstr.le.ib).and.(ib.le.Jend))
                  j=ib
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  bounded=Lconvolve(bry).and.                           &
     &                    ((Istr.le.ib).and.(ib.le.Iend))
                  i=ib
                END IF
#   ifdef MASKING
                IF (bounded) THEN
                  compute=vmask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#    ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                compute=1.0_r8
#   endif
                IF (compute.gt.0.0_r8) THEN
                  B3d=0.0_r8
                  IF (bounded) THEN
                    B3d(ib,kb)=VscaleB(ib,kb)
                  END IF
                  CALL ad_conv_v3d_bry_tile (ng, tile, iADM, bry,       &
     &                                       BOUNDS(ng)%edge(:,v2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(bry,isVvel)/ifac, &
     &                                       NVstepsB(bry,isVvel)/ifac, &
     &                                       DTsizeHB(bry,isVvel),      &
     &                                       DTsizeVB(bry,isVvel),      &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_p, pnom_r,    &
#   ifdef MASKING
     &                                       vmask, pmask,              &
#   endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
                  IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                    DO k=1,N(ng)
                      DO j=JV_RANGE
                        B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                      END DO
                    END DO
                  ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                    DO k=1,N(ng)
                      DO i=IV_RANGE
                        B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                      END DO
                    END DO
                  END IF
                  CALL tl_conv_v3d_bry_tile (ng, tile, iTLM, bry,       &
     &                                       BOUNDS(ng)%edge(:,v2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(bry,isVvel)/ifac, &
     &                                       NVstepsB(bry,isVvel)/ifac, &
     &                                       DTsizeHB(bry,isVvel),      &
     &                                       DTsizeVB(bry,isVvel),      &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_p, pnom_r,    &
#   ifdef MASKING
     &                                       vmask, pmask,              &
#   endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B3d(ib,kb))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  VnormVobc(ib,kb,bry)=cff
                END IF
              END DO
            END DO
            CALL bc_v3d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormVobc(:,:,bry))
#   ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJKlen, Aspv,                    &
     &                       VnormVobc(LBij:,:,bry))
#   endif
          END IF
        END DO
        IF (ANY(CnormB(isVvel,:))) THEN
          ifield=idSbry(isVvel)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          VnormVobc(LBij:,:,:),                   &
     &                          start = (/1,1,1,tindex/),               &
     &                          total = (/IJlen,N(ng),4,1/),            &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  3D boundary norm at RHO-points.
!
        IF (Master) THEN
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (ANY(CnormB(is,:))) THEN
              Lsame=.TRUE.
              EXIT
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '3D normalization factors at RHO-points'
            CALL my_flush (stdout)
          END IF
        END IF

        DO itrc=1,NT(ng)
          VnormRobc=Aspv
          is=isTvar(itrc)
          DO bry=1,4
            VscaleB=0.0_r8
            IF (CnormB(is,bry)) THEN
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(bry,r2dvar)
                Bmin=2
                Bmax=Mm(ng)
                IF (Lconvolve(bry)) THEN
                  DO j=JR_RANGE
                    cff=om_r(i,j)*on_r(i,j)
                    DO k=1,N(ng)
                      VscaleB(j,k)=1.0_r8/SQRT(cff*Hz(i,j,k))
                    END DO
                  END DO
                END IF
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(bry,r2dvar)
                Bmin=1
                Bmax=Lm(ng)
                IF (Lconvolve(bry)) THEN
                  DO i=IR_RANGE
                    cff=om_r(i,j)*on_r(i,j)
                    DO k=1,N(ng)
                      VscaleB(i,k)=1.0_r8/SQRT(cff*Hz(i,j,k))
                    END DO
                  END DO
                END IF
              END IF
              DO kb=1,N(ng)
                DO ib=Bmin,Bmax
                  IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                    bounded=Lconvolve(bry).and.                         &
     &                      ((Jstr.le.ib).and.(ib.le.Jend))
                    j=ib
                  ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                    bounded=Lconvolve(bry).and.                         &
     &                      ((Istr.le.ib).and.(ib.le.Iend))
                    i=ib
                  END IF
#   ifdef MASKING
                  IF (bounded) THEN
                    compute=rmask(i,j)
                  ELSE
                    compute=0.0_r8
                  END IF
#    ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                  compute=1.0_r8
#   endif
                  IF (compute.gt.0.0_r8) THEN
                    B3d=0.0_r8
                    IF (bounded) THEN
                      B3d(ib,kb)=VscaleB(ib,kb)
                    END IF
                    CALL ad_conv_r3d_bry_tile (ng, tile, iADM, bry,     &
     &                                         BOUNDS(ng)%edge(:,       &
     &                                                         r2dvar), &
     &                                         LBij, UBij,              &
     &                                         LBi, UBi, LBj, UBj,      &
     &                                         1, N(ng),                &
     &                                         IminS, ImaxS,            &
     &                                         JminS, JmaxS,            &
     &                                         NghostPoints,            &
     &                                         NHstepsB(bry,is)/ifac,   &
     &                                         NVstepsB(bry,is)/ifac,   &
     &                                         DTsizeHB(bry,is),        &
     &                                         DTsizeVB(bry,is),        &
     &                                         Kh, Kv,                  &
     &                                         pm, pn, pmon_u, pnom_v,  &
#   ifdef MASKING
     &                                         rmask, umask, vmask,     &
#   endif
     &                                         Hz, z_r,                 &
     &                                         B3d)
                    IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                      DO k=1,N(ng)
                        DO j=JR_RANGE
                          B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                        END DO
                      END DO
                    ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                      DO k=1,N(ng)
                        DO i=IR_RANGE
                          B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                        END DO
                      END DO
                    END IF
                    CALL tl_conv_r3d_bry_tile (ng, tile, iTLM, bry,     &
     &                                         BOUNDS(ng)%edge(:,       &
     &                                                         r2dvar), &
     &                                         LBij, UBij,              &
     &                                         LBi, UBi, LBj, UBj,      &
     &                                         1, N(ng),                &
     &                                         IminS, ImaxS,            &
     &                                         JminS, JmaxS,            &
     &                                         NghostPoints,            &
     &                                         NHstepsB(bry,is)/ifac,   &
     &                                         NVstepsB(bry,is)/ifac,   &
     &                                         DTsizeHB(bry,is),        &
     &                                         DTsizeVB(bry,is),        &
     &                                         Kh, Kv,                  &
     &                                         pm, pn, pmon_u, pnom_v,  &
#   ifdef MASKING
     &                                         rmask, umask, vmask,     &
#   endif
     &                                         Hz, z_r,                 &
     &                                         B3d)
                    IF (bounded) THEN
                      cff=1.0_r8/SQRT(B3d(ib,kb))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (bounded) THEN
                    VnormRobc(ib,kb,bry,itrc)=cff
                  END IF
                END DO
              END DO
              CALL bc_r3d_bry_tile (ng, tile, bry,                      &
     &                              LBij, UBij, 1, N(ng),               &
     &                              VnormRobc(:,:,bry,itrc))
#   ifdef DISTRIBUTE
              CALL mp_collect (ng, iTLM, IJKlen, Aspv,                  &
     &                         VnormRobc(LBij:,:,bry,itrc))
#   endif
            END IF
          END DO
          IF (ANY(CnormB(is,:))) THEN
            ifield=idSbry(isTvar(itrc))
            CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),    &
     &                            VnormRobc(LBij:,:,:,itrc),            &
     &                            start = (/1,1,1,tindex/),             &
     &                            total = (/IJlen,N(ng),4,1/),          &
     &                            ncid = ncNRMid(ifile,ng),             &
     &                            varid = nrmVid(ifile,ifield,ng))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
!
!  Synchronize open boundaries normalization NetCDF file to disk to
!  allow other processes to access data immediately after it is
!  written.
!
        CALL netcdf_sync (ng, iTLM, ncname, ncNRMid(ifile,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Compute surface forcing error covariance, B, normalization factors
!  using the exact method.
!-----------------------------------------------------------------------
!
      ifile=4
      IF (LwrtNRM(ifile,ng)) THEN
        rec=1
        Text='surface forcing'
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRMname(ifile,ng)
        tNRMindx(ifile,ng)=tNRMindx(ifile,ng)+1
        NrecNRM(ifile,ng)=NrecNRM(ifile,ng)+1
!
!  Write out model time (s).
!
        CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,idtime),        &
     &                        time(ng),                                 &
     &                        start = (/tNRMindx(ifile,ng)/),           &
     &                        total = (/1/),                            &
     &                        ncid = ncNRMid(ifile,ng),                 &
     &                        varid = nrmVid(ifile,idtime,ng))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef ADJUST_WSTRESS
!
!  2D norm at U-stress points.
!
        IF (Cnorm(rec,isUstr)) THEN
#   ifdef EW_PERIODIC
          Imin=1
          Imax=Lm(ng)
          Jmin=1
          Jmax=Mm(ng)
#   else
          Imin=2
          Imax=Lm(ng)
          Jmin=1
          Jmax=Mm(ng)
#   endif
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                     '2D normalization factors at U-stress points'
            CALL my_flush (stdout)
          END IF
          DO j=JU_RANGE
            DO i=IU_RANGE
              Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
            END DO
          END DO
          DO jc=Jmin,Jmax
            DO ic=Imin,Imax
#   ifdef MASKING
              compute=0.0_r8
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                IF (umask(ic,jc).gt.0) compute=1.0_r8
              END IF
#    ifdef DISTRIBUTE
              CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
              compute=1.0_r8
#   endif
              IF (compute.gt.0.0_r8) THEN
                DO j=LBj,UBj
                  DO i=LBi,UBi
                    A2d(i,j)=0.0_r8
                  END DO
                END DO
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  A2d(ic,jc)=Hscale(ic,jc)
                END IF
                CALL ad_conv_u2d_tile (ng, tile, iADM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(rec,isUstr)/ifac,        &
     &                                 DTsizeH(rec,isUstr),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_r, pnom_p,          &
#   ifdef MASKING
     &                                 umask, pmask,                    &
#   endif
     &                                 A2d)
                DO j=JU_RANGE
                  DO i=IU_RANGE
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
                CALL tl_conv_u2d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(rec,isUstr)/ifac,        &
     &                                 DTsizeH(rec,isUstr),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_r, pnom_p,          &
#   ifdef MASKING
     &                                 umask, pmask,                    &
#   endif
     &                                 A2d)
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  cff=1.0_r8/SQRT(A2d(ic,jc))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                HnormSUS(ic,jc)=cff
              END IF
            END DO
          END DO
          CALL bc_u2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      HnormSUS)
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        HnormSUS)
#   endif
          CALL wrt_norm2d (ng, tile, iTLM, ncname,                      &
     &                     LBi, UBi, LBj, UBj, idUsms,                  &
     &                     ncNRMid(ifile,ng),                           &
     &                     nrmVid(ifile,idUsms,ng),                     &
     &                     tNRMindx(ifile,ng),                          &
#   ifdef MASKING
     &                     umask,                                       &
#   endif
     &                     HnormSUS)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  2D norm at V-stress points.
!
        IF (Cnorm(rec,isVstr)) THEN
#   ifdef NS_PERIODIC
          Imin=1
          Imax=Lm(ng)
          Jmin=1
          Jmax=Mm(ng)
#   else
          Imin=1
          Imax=Lm(ng)
          Jmin=2
          Jmax=Mm(ng)
#   endif
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                     '2D normalization factors at V-stress points'
            CALL my_flush (stdout)
          END IF
          DO j=JV_RANGE
            DO i=IV_RANGE
              Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
            END DO
          END DO
          DO jc=Jmin,Jmax
            DO ic=Imin,Imax
#   ifdef MASKING
              compute=0.0_r8
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                IF (vmask(ic,jc).gt.0) compute=1.0_r8
              END IF
#    ifdef DISTRIBUTE
              CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
              compute=1.0_r8
#   endif
              IF (compute.gt.0.0_r8) THEN
                DO j=LBj,UBj
                  DO i=LBi,UBi
                    A2d(i,j)=0.0_r8
                  END DO
                END DO
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  A2d(ic,jc)=Hscale(ic,jc)
                END IF
                CALL ad_conv_v2d_tile (ng, tile, iADM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(rec,isVstr)/ifac,        &
     &                                 DTsizeH(rec,isVstr),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_p, pnom_r,          &
#   ifdef MASKING
     &                                 vmask, pmask,                    &
#   endif
     &                                 A2d)
                DO j=JV_RANGE
                  DO i=IV_RANGE
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
                CALL tl_conv_v2d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(rec,isVstr)/ifac,        &
     &                                 DTsizeH(rec,isVstr),             &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_p, pnom_r,          &
#   ifdef MASKING
     &                                 vmask, pmask,                    &
#   endif
     &                                 A2d)
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  cff=1.0_r8/SQRT(A2d(ic,jc))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                HnormSVS(ic,jc)=cff
              END IF
            END DO
          END DO
          CALL bc_v2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      HnormSVS)
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        HnormSVS)
#   endif
          CALL wrt_norm2d (ng, tile, iTLM, ncname,                      &
     &                     LBi, UBi, LBj, UBj, idVsms,                  &
     &                     ncNRMid(ifile,ng),                           &
     &                     nrmVid(ifile,idVsms,ng),                     &
     &                     tNRMindx(ifile,ng),                          &
#   ifdef MASKING
     &                     vmask,                                       &
#   endif
     &                     HnormSVS)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  2D norm at surface tracer fluxes points.
!
        IF (Master) THEN
          Lsame=.FALSE.
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              is=isTsur(itrc)
              IF (Cnorm(rec,is)) Lsame=.TRUE.
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
                              '2D normalization factors at RHO-points'
            CALL my_flush (stdout)
          END IF
        END IF
!
!  Check if the decorrelation scales for all the surface tracer fluxes
!  are different. If not, just compute the normalization factors for the
!  first tracer and assign the same value to the rest.  Recall that this
!  computation is very expensive.
!
        Ldiffer=.FALSE.
        Imin=1
        Imax=Lm(ng)
        Jmin=1
        Jmax=Mm(ng)
        DO itrc=2,NT(ng)
          IF (Hdecay(rec,isTsur(itrc  ),ng).ne.                         &
     &        Hdecay(rec,isTsur(itrc-1),ng)) THEN
            Ldiffer=.TRUE.
          END IF
        END DO
        IF (.not.Ldiffer) THEN
          Lsame=.TRUE.
          UBt=1
        ELSE
          Lsame=.FALSE.
          UBt=NT(ng)
        END IF
!
        DO j=JR_RANGE
          DO i=IR_RANGE
            Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
          END DO
        END DO
        DO itrc=1,UBt
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#   ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (rmask(ic,jc).gt.0) compute=1.0_r8
                  END IF
#    ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                  compute=1.0_r8
#   endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=Hscale(ic,jc)
                    END IF
                    CALL ad_conv_r2d_tile (ng, tile, iADM,              &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHsteps(rec,is)/ifac,        &
     &                                     DTsizeH(rec,is),             &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_u, pnom_v,      &
#   ifdef MASKING
     &                                     rmask, umask, vmask,         &
#   endif
     &                                     A2d)
                    DO j=JR_RANGE
                      DO i=IR_RANGE
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
                    CALL tl_conv_r2d_tile (ng, tile, iTLM,              &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHsteps(rec,is)/ifac,        &
     &                                     DTsizeH(rec,is),             &
     &                                     Kh,                          &
     &                                     pm, pn, pmon_u, pnom_v,      &
#   ifdef MASKING
     &                                     rmask, umask, vmask,         &
#   endif
     &                                     A2d)
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      cff=1.0_r8/SQRT(A2d(ic,jc))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (Lsame) THEN
                      DO ntrc=1,NT(ng)
                        IF (Lstflux(ntrc,ng)) THEN
                          HnormSTF(ic,jc,ntrc)=cff
                        END IF
                      END DO
                    ELSE
                      HnormSTF(ic,jc,itrc)=cff
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
        END DO
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              CALL bc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormSTF(:,:,itrc))
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, EWperiodic, NSperiodic, &
     &                            HnormSTF(:,:,itrc))
#   endif
              CALL wrt_norm2d (ng, tile, iTLM, ncname,                  &
     &                         LBi, UBi, LBj, UBj, idTsur(itrc),        &
     &                         ncNRMid(ifile,ng),                       &
     &                         nrmVid(ifile,idTsur(itrc),ng),           &
     &                         tNRMindx(ifile,ng),                      &
#   ifdef MASKING
     &                         rmask,                                   &
#   endif
     &                         HnormSTF(:,:,itrc))
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END IF
        END DO
#  endif
      END IF
# endif
!
      IF (Master) THEN
        WRITE (stdout,30)
      END IF

 10   FORMAT (/,' Error Covariance Normalization Factors: ',            &
     &        'Exact Method',/)
 20   FORMAT (4x,'Computing',1x,a,1x,a)
 30   FORMAT (/)

      RETURN
      END SUBROUTINE normalization_tile

!
!***********************************************************************
      SUBROUTINE randomization_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nstp, nnew, ifac,                  &
     &                               pm, om_r, om_u, om_v,              &
     &                               pn, on_r, on_u, on_v,              &
     &                               pmon_p, pmon_r, pmon_u,            &
     &                               pnom_p, pnom_r, pnom_v,            &
# ifdef MASKING
     &                               rmask, pmask, umask, vmask,        &
# endif
# ifdef SOLVE3D
     &                               h,                                 &
#  ifdef ICESHELF
     &                               zice,                              &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                               bed_thick,                         &
#  endif
     &                               Hz, z_r, z_w,                      &
# endif
     &                               Kh,                                &
# ifdef SOLVE3D
     &                               Kv,                                &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                               VnormRobc, VnormUobc, VnormVobc,   &
#  endif
     &                               HnormRobc, HnormUobc, HnormVobc,   &
# endif
# ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                               HnormSTF,                          &
# endif
# ifdef SOLVE3D
     &                               VnormR, VnormU, VnormV,            &
# endif
     &                               HnormR, HnormU, HnormV)
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
      USE bc_2d_mod
      USE tl_conv_2d_mod
# ifdef SOLVE3D
      USE bc_3d_mod
      USE tl_conv_3d_mod
# endif
# ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod
      USE tl_conv_bry2d_mod
#  ifdef SOLVE3D
      USE bc_bry3d_mod
      USE tl_conv_bry3d_mod
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE set_depth_mod
      USE white_noise_mod
# ifdef DISTRIBUTE
#  ifdef ADJUST_BOUNDARY
      USE distribute_mod, ONLY : mp_collect
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, ifac
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Kh(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent(out) :: VnormUobc(LBij:,:,:)
      real(r8), intent(out) :: VnormVobc(LBij:,:,:)
#   endif
      real(r8), intent(out) :: HnormRobc(LBij:,:)
      real(r8), intent(out) :: HnormUobc(LBij:,:)
      real(r8), intent(out) :: HnormVobc(LBij:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
#  endif
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(out) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent(out) :: VnormVobc(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent(out) :: HnormRobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormUobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormVobc(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng),NSA)
#  endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj,NSA)
#  ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
# ifdef SOLVE3D
      logical :: Ldiffer, Lsame
# endif
# ifdef ADJUST_BOUNDARY
      logical :: Lconvolve(4)
# endif
      integer :: i, ifile, is, iter, j, rec
# ifdef SOLVE3D
      integer :: UBt, itrc, k
# endif
# ifdef ADJUST_BOUNDARY
      integer :: IJlen, IJKlen, bry, ifield, tindex
# endif
      integer :: start(4), total(4)

      real(r8) :: Aavg, Amax, Amin, FacAvg, FacSqr, cff, val
# ifdef ADJUST_BOUNDARY
      real(r8) :: Bavg, Bmin, Bmax

      real(r8), parameter :: Aspv = 0.0_r8
# endif

      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2davg
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
# ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
      real(r8), dimension(LBij:UBij) :: B2davg
      real(r8), dimension(LBij:UBij) :: B2dsqr
      real(r8), dimension(LBij:UBij) :: HscaleB
# endif
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3davg
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
#  ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3davg
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3dsqr
      real(r8), dimension(LBij:UBij,1:N(ng)) :: VscaleB
#  endif
# endif

      character (len=40) :: Text
      character (len=80) :: ncname

# include "set_bounds.h"
!
      SourceFile='normalization.F, randomization_tile'

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
#  endif
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
# endif
!
!-----------------------------------------------------------------------
!  Compute initial conditions and model erro covariance, B,
!  normalization factors using the randomization approach of Fisher
!  and Courtier (1995). These factors ensure that the diagonal
!  elements of B are equal to unity. Notice that in applications
!  with land/sea masking, the boundary conditions will produce
!  large changes in the covariance structures near the boundary.
!
!  Initialize factors with randon numbers ("white-noise") having an
!  uniform distribution (zero mean and unity variance). Then, scale
!  by the inverse squared root area (2D) or volume (3D) and "color"
!  with the diffusion operator. Iterate this step over a specified
!  number of ensamble members, Nrandom.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,10)

      FILE_LOOP : DO ifile=1,NSA

        IF (LwrtNRM(ifile,ng)) THEN
          IF (ifile.eq.1) THEN
            Text='initial conditions'
          ELSE IF (ifile.eq.2) THEN
            Text='model'
          END IF
!
!  Set randomization summation factors.
!
          FacAvg=1.0_r8/REAL(Nrandom,r8)
          FacSqr=SQRT(REAL(Nrandom-1,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
          ncname=NRMname(ifile,ng)
          tNRMindx(ifile,ng)=tNRMindx(ifile,ng)+1
          NrecNRM(ifile,ng)=NrecNRM(ifile,ng)+1
!
!  Write out model time (s).
!
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,idtime),      &
     &                          time(ng),                               &
     &                          start = (/tNRMindx(ifile,ng)/),         &
     &                          total = (/1/),                          &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,idtime,ng))
          IF (exit_flag.ne.NoError) RETURN
!
!  2D norm at RHO-points.
!
          IF (Cnorm(ifile,isFsur)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at RHO-points'
              CALL my_flush (stdout)
            END IF
            DO j=JR_RANGE
              DO i=IR_RANGE
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),        &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
              DO j=JR_RANGE
                DO i=IR_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_r2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(ifile,isFsur)/ifac,        &
     &                               DTsizeH(ifile,isFsur),             &
     &                               Kh,                                &
     &                               pm, pn, pmon_u, pnom_v,            &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
     &                               A2d)
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=Istr,Iend
                Aavg=FacAvg*A2davg(i,j)
# ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
!!                HnormR(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                  HnormR(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j))
                ELSE
                  HnormR(i,j,ifile)=0.0_r8
                END IF
# else
!!              HnormR(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                HnormR(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j))
# endif
              END DO
            END DO
            CALL bc_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormR(:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          HnormR(:,:,ifile))
# endif
            CALL wrt_norm2d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, idFsur,                &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idFsur,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
# ifdef MASKING
     &                       rmask,                                     &
# endif
     &                       HnormR(:,:,ifile))
          END IF
!
!  2D norm at U-points.
!
          IF (Cnorm(ifile,isUbar)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at   U-points'
              CALL my_flush (stdout)
            END IF
            DO j=JU_RANGE
              DO i=IU_RANGE
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),        &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
              DO j=JU_RANGE
                DO i=IU_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_u2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(ifile,isUbar)/ifac,        &
     &                               DTsizeH(ifile,isUbar),             &
     &                               Kh,                                &
     &                               pm, pn, pmon_r, pnom_p,            &
# ifdef MASKING
     &                               umask, pmask,                      &
# endif
     &                               A2d)
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                Aavg=FacAvg*A2davg(i,j)
# ifdef MASKING
                IF (umask(i,j).gt.0.0_r8) THEN
!!                HnormU(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                  HnormU(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j))
                ELSE
                  HnormU(i,j,ifile)=0.0_r8
                END IF
# else
!!              HnormU(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                HnormU(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j))
# endif
              END DO
            END DO
            CALL bc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormU(:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          HnormU(:,:,ifile))
# endif
            CALL wrt_norm2d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, idUbar,                &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idUbar,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
# ifdef MASKING
     &                       umask,                                     &
# endif
     &                       HnormU(:,:,ifile))
          END IF
!
!  2D norm at V-points.
!
          IF (Cnorm(ifile,isVbar)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at   V-points'
              CALL my_flush (stdout)
            END IF
            DO j=JV_RANGE
              DO i=IV_RANGE
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),        &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
              DO j=JV_RANGE
                DO i=IV_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_v2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(ifile,isVbar)/ifac,        &
     &                               DTsizeH(ifile,isVbar),             &
     &                               Kh,                                &
     &                               pm, pn, pmon_p, pnom_r,            &
# ifdef MASKING
     &                               vmask, pmask,                      &
# endif
     &                               A2d)
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO
            DO j=JstrV,Jend
              DO i=Istr,Iend
                Aavg=FacAvg*A2davg(i,j)
# ifdef MASKING
                IF (vmask(i,j).gt.0.0_r8) THEN
!!                HnormV(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                  HnormV(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j))
                ELSE
                  HnormV(i,j,ifile)=0.0_r8
                END IF
# else
!!              HnormV(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                HnormV(i,j,ifile)=FacSqr/SQRT(A2dsqr(i,j))
# endif
              END DO
            END DO
            CALL bc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormV(:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          HnormV(:,:,ifile))
# endif
            CALL wrt_norm2d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, idVbar,                &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idVbar,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
# ifdef MASKING
     &                       vmask,                                     &
# endif
     &                       HnormV(:,:,ifile))
          END IF

# ifdef SOLVE3D
!
!  3D norm U-points.
!
          IF (Cnorm(ifile,isUvel)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at   U-points'
              CALL my_flush (stdout)
            END IF
            DO j=JU_RANGE
              DO i=IU_RANGE
                val=om_u(i,j)*on_u(i,j)*0.5_r8
                DO k=1,N(ng)
                  A3davg(i,j,k)=0.0_r8
                  A3dsqr(i,j,k)=0.0_r8
                  Vscale(i,j,k)=1.0_r8/SQRT(val*(Hz(i-1,j,k)+Hz(i,j,k)))
                END DO
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise3d (ng, iTLM, u3dvar, Rscheme(ng),        &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            Amin, Amax, A3d)
              DO k=1,N(ng)
                DO j=JU_RANGE
                  DO i=IU_RANGE
                    A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                  END DO
                END DO
              END DO
              CALL tl_conv_u3d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(ifile,isUvel)/ifac,        &
     &                               NVsteps(ifile,isUvel)/ifac,        &
     &                               DTsizeH(ifile,isUvel),             &
     &                               DTsizeV(ifile,isUvel),             &
     &                               Kh, Kv,                            &
     &                               pm, pn, pmon_r, pnom_p,            &
#  ifdef MASKING
     &                               umask, pmask,                      &
#  endif
     &                               Hz, z_r,                           &
     &                               A3d)
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=IstrU,Iend
                    A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                    A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                  END DO
                END DO
              END DO
            END DO
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  Aavg=FacAvg*A3davg(i,j,k)
#  ifdef MASKING
                  IF (umask(i,j).gt.0.0_r8) THEN
!!                  VnormU(i,j,k,ifile)=FacSqr/                         &
!!   &                                  SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                    VnormU(i,j,k,ifile)=FacSqr/SQRT(A3dsqr(i,j,k))
                  ELSE
                    VnormU(i,j,k,ifile)=0.0_r8
                  END IF
#  else
!!                VnormU(i,j,k,ifile)=FacSqr/                           &
!!   &                                SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                  VnormU(i,j,k,ifile)=FacSqr/SQRT(A3dsqr(i,j,k))
#  endif
                END DO
              END DO
            END DO
            CALL bc_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        VnormU(:,:,:,ifile))
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          VnormU(:,:,:,ifile))
#  endif
            CALL wrt_norm3d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), idUvel,      &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idUvel,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
#  ifdef MASKING
     &                       umask,                                     &
#  endif
     &                       VnormU(:,:,:,ifile))
          END IF
!
!  3D norm at V-points.
!
          IF (Cnorm(ifile,isVvel)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at   V-points'
              CALL my_flush (stdout)
            END IF
            DO j=JV_RANGE
              DO i=IV_RANGE
                val=om_v(i,j)*on_v(i,j)*0.5_r8
                DO k=1,N(ng)
                  A3davg(i,j,k)=0.0_r8
                  A3dsqr(i,j,k)=0.0_r8
                  Vscale(i,j,k)=1.0_r8/SQRT(val*(Hz(i,j-1,k)+Hz(i,j,k)))
                END DO
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise3d (ng, iTLM, v3dvar, Rscheme(ng),        &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            Amin, Amax, A3d)
              DO k=1,N(ng)
                DO j=JV_RANGE
                  DO i=IV_RANGE
                    A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                  END DO
                END DO
              END DO
              CALL tl_conv_v3d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(ifile,isVvel)/ifac,        &
     &                               NVsteps(ifile,isVvel)/ifac,        &
     &                               DTsizeH(ifile,isVvel),             &
     &                               DTsizeV(ifile,isVvel),             &
     &                               Kh, Kv,                            &
     &                               pm, pn, pmon_p, pnom_r,            &
#  ifdef MASKING
     &                               vmask, pmask,                      &
#  endif
     &                               Hz, z_r,                           &
     &                               A3d)
              DO k=1,N(ng)
                DO j=JstrV,Jend
                  DO i=Istr,Iend
                    A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                    A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                  END DO
                END DO
              END DO
            END DO
            DO k=1,N(ng)
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  Aavg=FacAvg*A3davg(i,j,k)
#  ifdef MASKING
                  IF (vmask(i,j).gt.0.0_r8) THEN
!!                  VnormV(i,j,k,ifile)=FacSqr/                         &
!!   &                                  SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                    VnormV(i,j,k,ifile)=FacSqr/SQRT(A3dsqr(i,j,k))
                  ELSE
                    VnormV(i,j,k,ifile)=0.0_r8
                  END IF
#  else
!!                VnormV(i,j,k,ifile)=FacSqr/                           &
!!   &                                SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                  VnormV(i,j,k,ifile)=FacSqr/SQRT(A3dsqr(i,j,k))
#  endif
                END DO
              END DO
            END DO
            CALL bc_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        VnormV(:,:,:,ifile))
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          VnormV(:,:,:,ifile))
#  endif
            CALL wrt_norm3d (ng, tile, iTLM, ncname,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), idVvel,      &
     &                       ncNRMid(ifile,ng),                         &
     &                       nrmVid(ifile,idVvel,ng),                   &
     &                       tNRMindx(ifile,ng),                        &
#  ifdef MASKING
     &                       vmask,                                     &
#  endif
     &                       VnormV(:,:,:,ifile))
          END IF
!
!  3D norm at RHO-points.
!
          IF (Master) THEN
            Lsame=.FALSE.
            DO itrc=1,NT(ng)
              is=isTvar(itrc)
              IF (Cnorm(ifile,is)) Lsame=.TRUE.
            END DO
            IF (Lsame) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at RHO-points'
              CALL my_flush (stdout)
            END IF
          END IF
!
!  Check if the decorrelation scales for all the tracers are different.
!  If not, just compute the normalization factors for the first tracer
!  and assign the same value to the rest.  Recall that this computation
!  is very expensive.
!
          Ldiffer=.FALSE.
          DO itrc=2,NT(ng)
            IF ((Hdecay(ifile,isTvar(itrc  ),ng).ne.                    &
     &           Hdecay(ifile,isTvar(itrc-1),ng)).or.                   &
     &          (Vdecay(ifile,isTvar(itrc  ),ng).ne.                    &
     &           Vdecay(ifile,isTvar(itrc-1),ng))) THEN
              Ldiffer=.TRUE.
            END IF
          END DO
          IF (.not.Ldiffer) THEN
            Lsame=.TRUE.
            UBt=1
          ELSE
            Lsame=.FALSE.
            UBt=NT(ng)
          END IF
!
          DO j=JR_RANGE
            DO i=IR_RANGE
              val=om_r(i,j)*on_r(i,j)
              DO k=1,N(ng)
                Vscale(i,j,k)=1.0_r8/SQRT(val*Hz(i,j,k))
              END DO
            END DO
          END DO
          DO itrc=1,UBt
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              DO k=1,N(ng)
                DO j=JR_RANGE
                  DO i=IR_RANGE
                    A3davg(i,j,k)=0.0_r8
                    A3dsqr(i,j,k)=0.0_r8
                  END DO
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise3d (ng, iTLM, r3dvar, Rscheme(ng),      &
     &                              IstrR, IendR, JstrR, JendR,         &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              Amin, Amax, A3d)
                DO k=1,N(ng)
                  DO j=JR_RANGE
                    DO i=IR_RANGE
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
                CALL tl_conv_r3d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(ifile,is)/ifac,          &
     &                                 NVsteps(ifile,is)/ifac,          &
     &                                 DTsizeH(ifile,is),               &
     &                                 DTsizeV(ifile,is),               &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_u, pnom_v,          &
#  ifdef MASKING
     &                                 rmask, umask, vmask,             &
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    DO i=Istr,Iend
                      A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                      A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                    END DO
                  END DO
                END DO
              END DO
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Aavg=FacAvg*A3davg(i,j,k)
#  ifdef MASKING
                    IF (rmask(i,j).gt.0.0_r8) THEN
!!                    VnormR(i,j,k,ifile,itrc)=FacSqr/                  &
!!   &                                         SQRT(A3dsqr(i,j,k)-      &
!!   &                                              Aavg*Aavg)
                      VnormR(i,j,k,ifile,itrc)=FacSqr/                  &
     &                                         SQRT(A3dsqr(i,j,k))
                    ELSE
                      VnormR(i,j,k,ifile,itrc)=0.0_r8
                    END IF
#  else
!!                  VnormR(i,j,k,ifile,itrc)=FacSqr/                    &
!!   &                                       SQRT(A3dsqr(i,j,k)-        &
!!   &                                            Aavg*Aavg)
                    VnormR(i,j,k,ifile,itrc)=FacSqr/SQRT(A3dsqr(i,j,k))
#  endif
                  END DO
                END DO
              END DO
            END IF
          END DO
          IF (Lsame) THEN
            DO itrc=2,NT(ng)
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    VnormR(i,j,k,ifile,itrc)=VnormR(i,j,k,ifile,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              CALL bc_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormR(:,:,:,ifile,itrc))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints, EWperiodic, NSperiodic, &
     &                            VnormR(:,:,:,ifile,itrc))
#  endif
              CALL wrt_norm3d (ng, tile, iTLM, ncname,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         idTvar(itrc), ncNRMid(ifile,ng),         &
     &                         nrmVid(ifile,idTvar(itrc),ng),           &
     &                         tNRMindx(ifile,ng),                      &
#  ifdef MASKING
     &                         rmask,                                   &
#  endif
     &                         VnormR(:,:,:,ifile,itrc))
            END IF
          END DO
# endif
        END IF
      END DO FILE_LOOP

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Compute open boundaries error covariance, B, normalization factors
!  using the randomization approach of Fisher and Courtier (1995).
!-----------------------------------------------------------------------
!
      ifile=3
      IF (LwrtNRM(ifile,ng)) THEN
        Text='boundary conditions'
        IJlen=UBij-LBij+1
#  ifdef SOLVE3D
        IJKlen=IJlen*N(ng)
#  endif
        Lconvolve(iwest )=WESTERN_EDGE
        Lconvolve(ieast )=EASTERN_EDGE
        Lconvolve(isouth)=SOUTHERN_EDGE
        Lconvolve(inorth)=NORTHERN_EDGE
!
!  Set randomization summation factors.
!
        FacAvg=1.0_r8/REAL(Nrandom,r8)
        FacSqr=SQRT(REAL(Nrandom-1,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRMname(ifile,ng)
        tNRMindx(ifile,ng)=tNRMindx(ifile,ng)+1
        NrecNRM(ifile,ng)=NrecNRM(ifile,ng)+1
        tindex=tNRMindx(ifile,ng)
!
!  Write out model time (s).
!
        CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,idtime),        &
     &                        time(ng),                                 &
     &                        start = (/tindex/),                       &
     &                        total = (/1/),                            &
     &                        ncid = ncNRMid(ifile,ng),                 &
     &                        varid = nrmVid(ifile,idtime,ng))
        IF (exit_flag.ne.NoError) RETURN
!
!  2D boundary norm at RHO-points.
!
        HnormRobc=Aspv

        IF (Master.and.ANY(CnormB(isFsur,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at RHO-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          IF (CnormB(isFsur,bry)) THEN
            HscaleB=0.0_r8
            B2davg=0.0_r8
            B2dsqr=0.0_r8
            B2d=0.0_r8
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,r2dvar)
              IF (Lconvolve(bry)) THEN
                DO j=JR_RANGE
                  HscaleB(j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,r2dvar)
              IF (Lconvolve(bry)) THEN
                DO i=IR_RANGE
                  HscaleB(i)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
                END DO
              END IF
            END IF
            DO iter=1,Nrandom
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                CALL white_noise2d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  JstrR, JendR,                   &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                DO j=JR_RANGE
                  B2d(j)=B2d(j)*HscaleB(j)
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                CALL white_noise2d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  IstrR, IendR,                   &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                DO i=IR_RANGE
                  B2d(i)=B2d(i)*HscaleB(i)
                END DO
              END IF
              CALL tl_conv_r2d_bry_tile (ng, tile, iTLM, bry,           &
     &                                   BOUNDS(ng)%edge(:,r2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(bry,isFsur)/ifac,     &
     &                                   DTsizeHB(bry,isFsur),          &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_u, pnom_v,        &
#  ifdef MASKING
     &                                   rmask, umask, vmask,           &
#  endif
     &                                   B2d)
              IF (Lconvolve(bry)) THEN
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO j=Jstr,Jend
                    B2davg(j)=B2davg(j)+B2d(j)
                    B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO i=Istr,Iend
                    B2davg(i)=B2davg(i)+B2d(i)
                    B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                  END DO
                END IF
              END IF
            END DO
            IF (Lconvolve(bry)) THEN
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                DO j=Jstr,Jend
                  Bavg=FacAvg*B2davg(j)
#  ifdef MASKING
                  IF (rmask(i,j).gt.0.0_r8) THEN
!!                  HnormRobc(j,bry)=FacSqr/SQRT(B2dsqr(j)-Bavg*Bavg)
                    HnormRobc(j,bry)=FacSqr/SQRT(B2dsqr(j))
                  ELSE
                    HnormRobc(j,bry)=0.0_r8
                  END IF
#  else
!!                HnormRobc(j,bry)=FacSqr/SQRT(B2dsqr(j)-Bavg*Bavg)
                  HnormRobc(j,bry)=FacSqr/SQRT(B2dsqr(j))
#  endif
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                DO i=Istr,Iend
                  Bavg=FacAvg*B2davg(i)
#  ifdef MASKING
                  IF (rmask(i,j).gt.0.0_r8) THEN
!!                  HnormRobc(i,bry)=FacSqr/SQRT(B2dsqr(i)-Bavg*Bavg)
                    HnormRobc(i,bry)=FacSqr/SQRT(B2dsqr(i))
                  ELSE
                    HnormRobc(i,bry)=0.0_r8
                  END IF
#  else
!!                HnormRobc(i,bry)=FacSqr/SQRT(B2dsqr(i)-Bavg*Bavg)
                  HnormRobc(i,bry)=FacSqr/SQRT(B2dsqr(i))
#  endif
                END DO
              END IF
            END IF
            CALL bc_r2d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij,                           &
     &                            HnormRobc(:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormRobc(LBij:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isFsur,:))) THEN
          ifield=idSbry(isFsur)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          HnormRobc(LBij:,:),                     &
     &                          start = (/1,1,tindex/),                 &
     &                          total = (/IJlen,4,1/),                  &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  2D boundary norm at U-points.
!
        HnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   U-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          IF (CnormB(isUbar,bry)) THEN
            HscaleB=0.0_r8
            B2davg=0.0_r8
            B2dsqr=0.0_r8
            B2d=0.0_r8
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,u2dvar)
              IF (Lconvolve(bry)) THEN
                DO j=JU_RANGE
                  HscaleB(j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,u2dvar)
              IF (Lconvolve(bry)) THEN
                DO i=IU_RANGE
                  HscaleB(i)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
                END DO
              END IF
            END IF
            DO iter=1,Nrandom
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                CALL white_noise2d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  JstrR, JendR,                   &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                DO j=JU_RANGE
                  B2d(j)=B2d(j)*HscaleB(j)
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                CALL white_noise2d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  Istr, IendR,                    &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                DO i=IU_RANGE
                  B2d(i)=B2d(i)*HscaleB(i)
                END DO
              END IF
              CALL tl_conv_u2d_bry_tile (ng, tile, iTLM, bry,           &
     &                                   BOUNDS(ng)%edge(:,u2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(bry,isUbar)/ifac,     &
     &                                   DTsizeHB(bry,isUbar),          &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_r, pnom_p,        &
#  ifdef MASKING
     &                                   umask, pmask,                  &
#  endif
     &                                   B2d)
              IF (Lconvolve(bry)) THEN
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO j=Jstr,Jend
                    B2davg(j)=B2davg(j)+B2d(j)
                    B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO i=IstrU,Iend
                    B2davg(i)=B2davg(i)+B2d(i)
                    B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                  END DO
                END IF
              END IF
            END DO
            IF (Lconvolve(bry)) THEN
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                DO j=Jstr,Jend
                  Bavg=FacAvg*B2davg(j)
#  ifdef MASKING
                  IF (umask(i,j).gt.0.0_r8) THEN
!!                  HnormUobc(j,bry)=FacSqr/SQRT(B2dsqr(j)-Bavg*Bavg)
                    HnormUobc(j,bry)=FacSqr/SQRT(B2dsqr(j))
                  ELSE
                    HnormUobc(j,bry)=0.0_r8
                  END IF
#  else
!!                HnormUobc(j,bry)=FacSqr/SQRT(B2dsqr(j)-Bavg*Bavg)
                  HnormUobc(j,bry)=FacSqr/SQRT(B2dsqr(j))
#  endif
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                DO i=IstrU,Iend
                  Bavg=FacAvg*B2davg(i)
#  ifdef MASKING
                  IF (umask(i,j).gt.0.0_r8) THEN
!!                  HnormUobc(i,bry)=FacSqr/SQRT(B2dsqr(i)-Bavg*Bavg)
                    HnormUobc(i,bry)=FacSqr/SQRT(B2dsqr(i))
                  ELSE
                    HnormUobc(i,bry)=0.0_r8
                  END IF
#  else
!!                HnormUobc(i,bry)=FacSqr/SQRT(B2dsqr(i)-Bavg*Bavg)
                  HnormUobc(i,bry)=FacSqr/SQRT(B2dsqr(i))
#  endif
                END DO
              END IF
            END IF
            CALL bc_u2d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij,                           &
     &                            HnormUobc(:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormUobc(LBij:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isUbar,:))) THEN
          ifield=idSbry(isUbar)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          HnormUobc(LBij:,:),                     &
     &                          start = (/1,1,tindex/),                 &
     &                          total = (/IJlen,4,1/),                  &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  2D boundary norm at V-points.
!
        HnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   V-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          IF (CnormB(isVbar,bry)) THEN
            HscaleB=0.0_r8
            B2davg=0.0_r8
            B2dsqr=0.0_r8
            B2d=0.0_r8
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,v2dvar)
              IF (Lconvolve(bry)) THEN
                DO j=JV_RANGE
                  HscaleB(j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,v2dvar)
              IF (Lconvolve(bry)) THEN
                DO i=IV_RANGE
                  HscaleB(i)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
                END DO
              END IF
            END IF
            DO iter=1,Nrandom
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                CALL white_noise2d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  Jstr, JendR,                    &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                DO j=JV_RANGE
                  B2d(j)=B2d(j)*HscaleB(j)
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                CALL white_noise2d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  IstrR, IendR,                   &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                DO i=IV_RANGE
                  B2d(i)=B2d(i)*HscaleB(i)
                END DO
              END IF
              CALL tl_conv_v2d_bry_tile (ng, tile, iTLM, bry,           &
     &                                   BOUNDS(ng)%edge(:,v2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(bry,isFsur)/ifac,     &
     &                                   DTsizeHB(bry,isFsur),          &
     &                                   Kh,                            &
     &                                   pm, pn, pmon_p, pnom_r,        &
#  ifdef MASKING
     &                                   vmask, pmask,                  &
#  endif
     &                                   B2d)
              IF (Lconvolve(bry)) THEN
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO j=JstrV,Jend
                    B2davg(j)=B2davg(j)+B2d(j)
                    B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO i=Istr,Iend
                    B2davg(i)=B2davg(i)+B2d(i)
                    B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                  END DO
                END IF
              END IF
            END DO
            IF (Lconvolve(bry)) THEN
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                DO j=JstrV,Jend
                  Bavg=FacAvg*B2davg(j)
#  ifdef MASKING
                  IF (vmask(i,j).gt.0.0_r8) THEN
!!                  HnormVobc(j,bry)=FacSqr/SQRT(B2dsqr(j)-Bavg*Bavg)
                    HnormVobc(j,bry)=FacSqr/SQRT(B2dsqr(j))
                  ELSE
                    HnormVobc(j,bry)=0.0_r8
                  END IF
#  else
!!                HnormVobc(j,bry)=FacSqr/SQRT(B2dsqr(j)-Bavg*Bavg)
                  HnormVobc(j,bry)=FacSqr/SQRT(B2dsqr(j))
#  endif
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                DO i=Istr,Iend
                  Bavg=FacAvg*B2davg(i)
#  ifdef MASKING
                  IF (vmask(i,j).gt.0.0_r8) THEN
!!                  HnormVobc(i,bry)=FacSqr/SQRT(B2dsqr(i)-Bavg*Bavg)
                    HnormVobc(i,bry)=FacSqr/SQRT(B2dsqr(i))
                  ELSE
                    HnormVobc(i,bry)=0.0_r8
                  END IF
#  else
!!                HnormVobc(i,bry)=FacSqr/SQRT(B2dsqr(i)-Bavg*Bavg)
                  HnormVobc(i,bry)=FacSqr/SQRT(B2dsqr(i))
#  endif
                END DO
              END IF
            END IF
            CALL bc_v2d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij,                           &
     &                            HnormVobc(:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormVobc(LBij:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isVbar,:))) THEN
          ifield=idSbry(isVbar)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          HnormVobc(LBij:,:),                     &
     &                          start = (/1,1,tindex/),                 &
     &                          total = (/IJlen,4,1/),                  &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  ifdef SOLVE3D
!
!  3D boundary norm at U-points.
!
        VnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   U-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          IF (CnormB(isUvel,bry)) THEN
            VscaleB=0.0_r8
            B3davg=0.0_r8
            B3dsqr=0.0_r8
            B3d=0.0_r8
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,u2dvar)
              IF (Lconvolve(bry)) THEN
                DO j=JU_RANGE
                  val=om_u(i,j)*on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(j,k)=1.0_r8/                                &
     &                           SQRT(val*(Hz(i-1,j,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,u2dvar)
              IF (Lconvolve(bry)) THEN
                DO i=IU_RANGE
                  val=om_u(i,j)*on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(i,k)=1.0_r8/                                &
     &                           SQRT(val*(Hz(i-1,j,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            END IF
            DO iter=1,Nrandom
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                CALL white_noise3d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  JstrR, JendR,                   &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  Bmin, Bmax, B3d)
                DO k=1,N(ng)
                  DO j=JU_RANGE
                    B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                  END DO
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                CALL white_noise3d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  Istr, IendR,                    &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  Bmin, Bmax, B3d)
                DO k=1,N(ng)
                  DO i=IU_RANGE
                    B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                  END DO
                END DO
              END IF
              CALL tl_conv_u3d_bry_tile (ng, tile, iTLM, bry,           &
     &                                   BOUNDS(ng)%edge(:,u2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(bry,isUvel)/ifac,     &
     &                                   NVstepsB(bry,isUvel)/ifac,     &
     &                                   DTsizeHB(bry,isUvel),          &
     &                                   DTsizeVB(bry,isUvel),          &
     &                                   Kh, Kv,                        &
     &                                   pm, pn, pmon_r, pnom_p,        &
#  ifdef MASKING
     &                                   umask, pmask,                  &
#  endif
     &                                   Hz, z_r,                       &
     &                                   B3d)
              IF (Lconvolve(bry)) THEN
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO k=1,N(ng)
                    DO j=Jstr,Jend
                      B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                      B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                    END DO
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO k=1,N(ng)
                    DO i=IstrU,Iend
                      B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                      B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                    END DO
                  END DO
                END IF
              END IF
            END DO
            IF (Lconvolve(bry)) THEN
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    Bavg=FacAvg*B3davg(j,k)
#  ifdef MASKING
                    IF (umask(i,j).gt.0.0_r8) THEN
!!                    VnormUobc(j,k,bry)=FacSqr/                        &
!!   &                                   SQRT(B3dsqr(j,k)-Bavg*Bavg)
                      VnormUobc(j,k,bry)=FacSqr/SQRT(B3dsqr(j,k))
                    ELSE
                      VnormUobc(j,k,bry)=0.0_r8
                    END IF
#  else
!!                  VnormUobc(j,k,bry)=FacSqr/                          &
!!   &                                 SQRT(B3dsqr(j,k)-Bavg*Bavg)
                    VnormUobc(j,k,bry)=FacSqr/SQRT(B3dsqr(j,k))
#  endif
                  END DO
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                DO k=1,N(ng)
                  DO i=IstrU,Iend
                    Bavg=FacAvg*B3davg(i,k)
#  ifdef MASKING
                    IF (umask(i,j).gt.0.0_r8) THEN
!!                    VnormUobc(i,k,bry)=FacSqr/                        &
!!   &                                   SQRT(B3dsqr(i,k)-Bavg*Bavg)
                      VnormUobc(i,k,bry)=FacSqr/SQRT(B3dsqr(i,k))
                    ELSE
                      VnormUobc(i,k,bry)=0.0_r8
                    END IF
#  else
!!                  VnormUobc(i,k,bry)=FacSqr/                          &
!!   &                                 SQRT(B3dsqr(i,k)-Bavg*Bavg)
                    VnormUobc(i,k,bry)=FacSqr/SQRT(B3dsqr(i,k))
#  endif
                  END DO
                END DO
              END IF
            END IF
            CALL bc_u3d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormUobc(:,:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJKlen, Aspv,                    &
     &                       VnormUobc(LBij:,:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isUvel,:))) THEN
          ifield=idSbry(isUvel)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          VnormUobc(LBij:,:,:),                   &
     &                          start = (/1,1,1,tindex/),               &
     &                          total = (/IJlen,N(ng),4,1/),            &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  3D boundary norm at V-points.
!
        VnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   V-points'
          CALL my_flush (stdout)
        END IF

        DO bry=1,4
          IF (CnormB(isVvel,bry)) THEN
            VscaleB=0.0_r8
            B3davg=0.0_r8
            B3dsqr=0.0_r8
            B3d=0.0_r8
            IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(bry,v2dvar)
              IF (Lconvolve(bry)) THEN
                DO j=JV_RANGE
                  val=om_v(i,j)*on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(j,k)=1.0_r8/                                &
     &                           SQRT(val*(Hz(i,j-1,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(bry,v2dvar)
              IF (Lconvolve(bry)) THEN
                DO i=IV_RANGE
                  val=om_v(i,j)*on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    VscaleB(i,k)=1.0_r8/                                &
     &                           SQRT(val*(Hz(i,j-1,k)+Hz(i,j,k)))
                  END DO
                END DO
              END IF
            END IF
            DO iter=1,Nrandom
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                CALL white_noise3d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  Jstr, JendR,                    &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  Bmin, Bmax, B3d)
                DO k=1,N(ng)
                  DO j=JV_RANGE
                    B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                  END DO
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                CALL white_noise3d_bry (ng, tile, iTLM, bry,            &
     &                                  Rscheme(ng),                    &
     &                                  IstrR, IendR,                   &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  Bmin, Bmax, B3d)
                DO k=1,N(ng)
                  DO i=IV_RANGE
                    B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                  END DO
                END DO
              END IF
              CALL tl_conv_v3d_bry_tile (ng, tile, iTLM, bry,           &
     &                                   BOUNDS(ng)%edge(:,v2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(bry,isVvel)/ifac,     &
     &                                   NVstepsB(bry,isVvel)/ifac,     &
     &                                   DTsizeHB(bry,isVvel),          &
     &                                   DTsizeVB(bry,isVvel),          &
     &                                   Kh, Kv,                        &
     &                                   pm, pn, pmon_p, pnom_r,        &
#  ifdef MASKING
     &                                   vmask, pmask,                  &
#  endif
     &                                   Hz, z_r,                       &
     &                                   B3d)
              IF (Lconvolve(bry)) THEN
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO k=1,N(ng)
                    DO j=JstrV,Jend
                      B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                      B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                    END DO
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO k=1,N(ng)
                    DO i=Istr,Iend
                      B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                      B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                    END DO
                  END DO
                END IF
              END IF
            END DO
            IF (Lconvolve(bry)) THEN
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                DO k=1,N(ng)
                  DO j=JstrV,Jend
                    Bavg=FacAvg*B3davg(j,k)
#  ifdef MASKING
                    IF (vmask(i,j).gt.0.0_r8) THEN
!!                    VnormVobc(j,k,bry)=FacSqr/                        &
!!   &                                   SQRT(B3dsqr(j,k)-Bavg*Bavg)
                      VnormVobc(j,k,bry)=FacSqr/SQRT(B3dsqr(j,k))
                    ELSE
                      VnormVobc(j,k,bry)=0.0_r8
                    END IF
#  else
!!                  VnormVobc(j,k,bry)=FacSqr/                          &
!!   &                                 SQRT(B3dsqr(j,k)-Bavg*Bavg)
                    VnormVobc(j,k,bry)=FacSqr/SQRT(B3dsqr(j,k))
#  endif
                  END DO
                END DO
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bavg=FacAvg*B3davg(i,k)
#  ifdef MASKING
                    IF (vmask(i,j).gt.0.0_r8) THEN
!!                    VnormVobc(i,k,bry)=FacSqr/                        &
!!   &                                   SQRT(B3dsqr(i,k)-Bavg*Bavg)
                      VnormVobc(i,k,bry)=FacSqr/SQRT(B3dsqr(i,k))
                    ELSE
                      VnormVobc(i,k,bry)=0.0_r8
                    END IF
#  else
!!                  VnormVobc(i,k,bry)=FacSqr/                          &
!!   &                                 SQRT(B3dsqr(i,k)-Bavg*Bavg)
                    VnormVobc(i,k,bry)=FacSqr/SQRT(B3dsqr(i,k))
#  endif
                  END DO
                END DO
              END IF
            END IF
            CALL bc_v3d_bry_tile (ng, tile, bry,                        &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormVobc(:,:,bry))
#  ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJKlen, Aspv,                    &
     &                       VnormVobc(LBij:,:,bry))
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isVvel,:))) THEN
          ifield=idSbry(isVvel)
          CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),      &
     &                          VnormVobc(LBij:,:,:),                   &
     &                          start = (/1,1,1,tindex/),               &
     &                          total = (/IJlen,N(ng),4,1/),            &
     &                          ncid = ncNRMid(ifile,ng),               &
     &                          varid = nrmVid(ifile,ifield,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  3D boundary norm at RHO-points.
!
        IF (Master) THEN
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (ANY(CnormB(is,:))) THEN
              Lsame=.TRUE.
              EXIT
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '3D normalization factors at RHO-points'
            CALL my_flush (stdout)
          END IF
        END IF

        DO itrc=1,NT(ng)
          VnormRobc=Aspv
          is=isTvar(itrc)
          DO bry=1,4
            IF (CnormB(is,bry)) THEN
              VscaleB=0.0_r8
              B3davg=0.0_r8
              B3dsqr=0.0_r8
              B3d=0.0_r8
              IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(bry,r2dvar)
                IF (Lconvolve(bry)) THEN
                  DO j=JR_RANGE
                    val=om_r(i,j)*on_r(i,j)
                    DO k=1,N(ng)
                      VscaleB(j,k)=1.0_r8/SQRT(val*Hz(i,j,k))
                    END DO
                  END DO
                END IF
              ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(bry,r2dvar)
                IF (Lconvolve(bry)) THEN
                  DO i=IR_RANGE
                    val=om_r(i,j)*on_r(i,j)
                    DO k=1,N(ng)
                      VscaleB(i,k)=1.0_r8/SQRT(val*Hz(i,j,k))
                    END DO
                  END DO
                END IF
              END IF
              DO iter=1,Nrandom
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  CALL white_noise3d_bry (ng, tile, iTLM, bry,          &
     &                                    Rscheme(ng),                  &
     &                                    JstrR, JendR,                 &
     &                                    LBij, UBij, 1, N(ng),         &
     &                                    Bmin, Bmax, B3d)
                  DO k=1,N(ng)
                    DO j=JR_RANGE
                      B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                    END DO
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  CALL white_noise3d_bry (ng, tile, iTLM, bry,          &
     &                                    Rscheme(ng),                  &
     &                                    IstrR, IendR,                 &
     &                                    LBij, UBij, 1, N(ng),         &
     &                                    Bmin, Bmax, B3d)
                  DO k=1,N(ng)
                    DO i=IR_RANGE
                      B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                    END DO
                  END DO
                END IF
                CALL tl_conv_r3d_bry_tile (ng, tile, iTLM, bry,         &
     &                                     BOUNDS(ng)%edge(:,r2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     1, N(ng),                    &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(bry,is)/ifac,       &
     &                                     NVstepsB(bry,is)/ifac,       &
     &                                     DTsizeHB(bry,is),            &
     &                                     DTsizeVB(bry,is),            &
     &                                     Kh, Kv,                      &
     &                                     pm, pn, pmon_u, pnom_v,      &
#  ifdef MASKING
     &                                     rmask, umask, vmask,         &
#  endif
     &                                     Hz, z_r,                     &
     &                                     B3d)
                IF (Lconvolve(bry)) THEN
                  IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                    DO k=1,N(ng)
                      DO j=Jstr,Jend
                        B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                        B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                      END DO
                    END DO
                  ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                    DO k=1,N(ng)
                      DO i=Istr,Iend
                        B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                        B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                      END DO
                    END DO
                  END IF
                END IF
              END DO
              IF (Lconvolve(bry)) THEN
                IF ((bry.eq.iwest).or.(bry.eq.ieast)) THEN
                  DO k=1,N(ng)
                    DO j=Jstr,Jend
                      Bavg=FacAvg*B3davg(j,k)
#  ifdef MASKING
                      IF (rmask(i,j).gt.0.0_r8) THEN
!!                      VnormRobc(j,k,bry,itrc)=FacSqr/                 &
!!   &                                          SQRT(B3dsqr(j,k)-       &
!!   &                                               Bavg*Bavg)
                        VnormRobc(j,k,bry,itrc)=FacSqr/SQRT(B3dsqr(j,k))
                      ELSE
                        VnormRobc(j,k,bry,itrc)=0.0_r8
                      END IF
#  else
!!                    VnormRobc(j,k,bry,itrc)=FacSqr/                   &
!!   &                                        SQRT(B3dsqr(j,k)-         &
!!   &                                             Bavg*Bavg)
                      VnormRobc(j,k,bry,itrc)=FacSqr/SQRT(B3dsqr(j,k))
#  endif
                    END DO
                  END DO
                ELSE IF ((bry.eq.isouth).or.(bry.eq.inorth)) THEN
                  DO k=1,N(ng)
                    DO i=Istr,Iend
                      Bavg=FacAvg*B3davg(i,k)
#  ifdef MASKING
                      IF (rmask(i,j).gt.0.0_r8) THEN
!!                      VnormRobc(i,k,bry,itrc)=FacSqr/                 &
!!   &                                          SQRT(B3dsqr(i,k)-       &
!!   &                                               Bavg*Bavg)
                        VnormRobc(i,k,bry,itrc)=FacSqr/SQRT(B3dsqr(i,k))
                      ELSE
                        VnormRobc(i,k,bry,itrc)=0.0_r8
                      END IF
#  else
!!                    VnormRobc(i,k,bry,itrc)=FacSqr/                   &
!!   &                                        SQRT(B3dsqr(i,k)-
!!   &                                             Bavg*Bavg)
                      VnormRobc(i,k,bry,itrc)=FacSqr/SQRT(B3dsqr(i,k))
#  endif
                    END DO
                  END DO
                END IF
              END IF
              CALL bc_r3d_bry_tile (ng, tile, bry,                      &
     &                              LBij, UBij, 1, N(ng),               &
     &                              VnormRobc(:,:,bry,itrc))
#  ifdef DISTRIBUTE
              CALL mp_collect (ng, iTLM, IJKlen, Aspv,                  &
     &                         VnormRobc(LBij:,:,bry,itrc))
#  endif
            END IF
          END DO
          IF (ANY(CnormB(is,:))) THEN
            ifield=idSbry(isTvar(itrc))
            CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,ifield),    &
     &                            VnormRobc(LBij:,:,:,itrc),            &
     &                            start = (/1,1,1,tindex/),             &
     &                            total = (/IJlen,N(ng),4,1/),          &
     &                            ncid = ncNRMid(ifile,ng),             &
     &                            varid = nrmVid(ifile,ifield,ng))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
!
!  Synchronize open boundaries normalization NetCDF file to disk to
!  allow other processes to access data immediately after it is
!  written.
!
        CALL netcdf_sync (ng, iTLM, ncname, ncNRMid(ifile,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Compute surface forcing error covariance, B, normalization factors
!  using the randomization approach of Fisher and Courtier (1995).
!-----------------------------------------------------------------------
!
      ifile=4
      IF (LwrtNRM(ifile,ng)) THEN
        rec=1
        Text='surface forcing'
!
!  Set randomization summation factors.
!
        FacAvg=1.0_r8/REAL(Nrandom,r8)
        FacSqr=SQRT(REAL(Nrandom-1,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRMname(ifile,ng)
        tNRMindx(ifile,ng)=tNRMindx(ifile,ng)+1
        NrecNRM(ifile,ng)=NrecNRM(ifile,ng)+1
!
!  Write out model time (s).
!
        CALL netcdf_put_fvar (ng, iTLM, ncname, Vname(1,idtime),        &
     &                        time(ng),                                 &
     &                        start = (/tNRMindx(ifile,ng)/),           &
     &                        total = (/1/),                            &
     &                        ncid = ncNRMid(ifile,ng),                 &
     &                        varid = nrmVid(ifile,idtime,ng))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef ADJUST_WSTRESS
!
!  2D norm at U-stress points.
!
        IF (Cnorm(rec,isUstr)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '2D normalization factors at U-points'
            CALL my_flush (stdout)
          END IF
          DO j=JU_RANGE
            DO i=IU_RANGE
              A2davg(i,j)=0.0_r8
              A2dsqr(i,j)=0.0_r8
              Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
            END DO
          END DO
          DO iter=1,Nrandom
            CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),          &
     &                          Istr, IendR, JstrR, JendR,              &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Amin, Amax, A2d)
            DO j=JU_RANGE
              DO i=IU_RANGE
                A2d(i,j)=A2d(i,j)*Hscale(i,j)
              END DO
            END DO
            CALL tl_conv_u2d_tile (ng, tile, iTLM,                      &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             NghostPoints,                        &
     &                             NHsteps(rec,isUstr)/ifac,            &
     &                             DTsizeH(rec,isUstr),                 &
     &                             Kh,                                  &
     &                             pm, pn, pmon_r, pnom_p,              &
#   ifdef MASKING
     &                             umask, pmask,                        &
#   endif
     &                             A2d)
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
              END DO
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              Aavg=FacAvg*A2davg(i,j)
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
!!              HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
              ELSE
                HnormSUS(i,j)=0.0_r8
              END IF
#   else
!!            HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
#   endif
            END DO
          END DO
          CALL bc_u2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      HnormSUS)
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        HnormSUS)
#   endif
          CALL wrt_norm2d (ng, tile, iTLM, ncname,                      &
     &                     LBi, UBi, LBj, UBj, idUsms,                  &
     &                     ncNRMid(ifile,ng),                           &
     &                     nrmVid(ifile,idUsms,ng),                     &
     &                     tNRMindx(ifile,ng),                          &
#   ifdef MASKING
     &                     umask,                                       &
#   endif
     &                     HnormSUS)
        END IF
!
!  2D norm at V-stress points.
!
        IF (Cnorm(rec,isVstr)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '2D normalization factors at V-points'
            CALL my_flush (stdout)
          END IF
          DO j=JV_RANGE
            DO i=IV_RANGE
              A2davg(i,j)=0.0_r8
              A2dsqr(i,j)=0.0_r8
              Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
            END DO
          END DO
          DO iter=1,Nrandom
            CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),          &
     &                          IstrR, IendR, Jstr, JendR,              &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Amin, Amax, A2d)
            DO j=JV_RANGE
              DO i=IV_RANGE
                A2d(i,j)=A2d(i,j)*Hscale(i,j)
              END DO
            END DO
            CALL tl_conv_v2d_tile (ng, tile, iTLM,                      &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             NghostPoints,                        &
     &                             NHsteps(rec,isVstr)/ifac,            &
     &                             DTsizeH(rec,isVstr),                 &
     &                             Kh,                                  &
     &                             pm, pn, pmon_p, pnom_r,              &
#   ifdef MASKING
     &                             vmask, pmask,                        &
#   endif
     &                             A2d)
            DO j=JstrV,Jend
              DO i=Istr,Iend
                A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
              END DO
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              Aavg=FacAvg*A2davg(i,j)
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
!!              HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
              ELSE
                HnormSVS(i,j)=0.0_r8
              END IF
#   else
!!            HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
#   endif
            END DO
          END DO
          CALL bc_v2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      HnormSVS)
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        HnormSVS)
#   endif
          CALL wrt_norm2d (ng, tile, iTLM, ncname,                      &
     &                     LBi, UBi, LBj, UBj, idVsms,                  &
     &                     ncNRMid(ifile,ng),                           &
     &                     nrmVid(ifile,idVsms,ng),                     &
     &                     tNRMindx(ifile,ng),                          &
#   ifdef MASKING
     &                     vmask,                                       &
#   endif
     &                     HnormSVS)
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  2D norm at surface treace flux points.
!
        IF (Master) THEN
          Lsame=.FALSE.
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              is=isTsur(itrc)
              IF (Cnorm(rec,is)) Lsame=.TRUE.
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '2D normalization factors at RHO-points'
            CALL my_flush (stdout)
          END IF
        END IF
!
!  Check if the decorrelation scales for all the surface tracer fluxes
!  are different. If not, just compute the normalization factors for the
!  first tracer and assign the same value to the rest.  Recall that this
!  computation is very expensive.
!
        Ldiffer=.FALSE.
        DO itrc=2,NT(ng)
          IF (Hdecay(rec,isTvar(itrc  ),ng).ne.                         &
     &        Hdecay(rec,isTvar(itrc-1),ng)) THEN
            Ldiffer=.TRUE.
          END IF
        END DO
        IF (.not.Ldiffer) THEN
          Lsame=.TRUE.
          UBt=1
        ELSE
          Lsame=.FALSE.
          UBt=NT(ng)
        END IF
!
        DO j=JR_RANGE
          DO i=IR_RANGE
            Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
          END DO
        END DO
        DO itrc=1,UBt
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              DO j=JR_RANGE
                DO i=IR_RANGE
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),      &
     &                              IstrR, IendR, JstrR, JendR,         &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Amin, Amax, A2d)
                DO j=JR_RANGE
                  DO i=IR_RANGE
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
                CALL tl_conv_r2d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(rec,is)/ifac,            &
     &                                 DTsizeH(rec,is),                 &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_u, pnom_v,          &
#   ifdef MASKING
     &                                 rmask, umask, vmask,             &
#   endif
     &                                 A2d)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Aavg=FacAvg*A2davg(i,j)
#   ifdef MASKING
                  IF (rmask(i,j).gt.0.0_r8) THEN
!!                  HnormSTF(i,j,itrc)=FacSqr/                          &
!!   &                                 SQRT(A2dsqr(i,j)-Aavg*Aavg)
                    HnormSTF(i,j,itrc)=FacSqr/SQRT(A2dsqr(i,j))
                  ELSE
                    HnormSTF(i,j,itrc)=0.0_r8
                  END IF
#   else
!!                HnormSTF(i,j,itrc)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
                  HnormSTF(i,j,itrc)=FacSqr/SQRT(A2dsqr(i,j))
#   endif
                END DO
              END DO
            END IF
          END IF
        END DO
        IF (Lsame) THEN
          DO itrc=2,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  HnormSTF(i,j,itrc)=HnormSTF(i,j,1)
                END DO
              END DO
            END IF
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              CALL bc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormSTF(:,:,itrc))
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints, EWperiodic, NSperiodic, &
     &                            HnormSTF(:,:,itrc))
#   endif
              CALL wrt_norm2d (ng, tile, iTLM, ncname,                  &
     &                         LBi, UBi, LBj, UBj, idTsur(itrc),        &
     &                         ncNRMid(ifile,ng),                       &
     &                         nrmVid(ifile,idTsur(itrc),ng),           &
     &                         tNRMindx(ifile,ng),                      &
#   ifdef MASKING
     &                         rmask,                                   &
#   endif
     &                         HnormSTF(:,:,itrc))
            END IF
          END IF
        END DO
#  endif
      END IF
# endif
!
      IF (Master) THEN
        WRITE (stdout,30)
      END IF

 10   FORMAT (/,' Error Covariance Factors: Randomization Method',/)
 20   FORMAT (4x,'Computing',1x,a,1x,a)
 30   FORMAT (/)

      RETURN
      END SUBROUTINE randomization_tile

!
!***********************************************************************
      SUBROUTINE wrt_norm2d (ng, tile, model, ncname,                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ifield, ncid, ncvarid, tindex,             &
# ifdef MASKING
     &                       Amask,                                     &
# endif
     &                       A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, status

      real(r8) :: scale

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 2D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
        gfactor=-1
# else
        gfactor=1
# endif
!
!  Write out 2D normalization field.
!
      gtype=gfactor*Iinfo(1,ifield,ng)
      scale=1.0_r8
      status=nf_fwrite2d(ng, model, ncid, ncvarid, tindex, gtype,       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   Amask,                                         &
# endif
     &                   A)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex
      CALL my_flush (stdout)
!
  10  FORMAT (/,' WRT_NORM2D - error while writing variable: ',a,/,11x, &
     &        'into normalization NetCDF file for time record: ',i4)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i7.7)

      END SUBROUTINE wrt_norm2d

!
!***********************************************************************
      SUBROUTINE wrt_norm3d (ng, tile, model, ncname,                   &
     &                       LBi, UBi, LBj, UBj, LBk, UBk,              &
     &                       ifield, ncid, ncvarid, tindex,             &
# ifdef MASKING
     &                       Amask,                                     &
# endif
     &                       A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, status

      real(r8) :: scale

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 3D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
        gfactor=-1
# else
        gfactor=1
# endif
!
!  Write out 3D normalization field.
!
      gtype=gfactor*Iinfo(1,ifield,ng)
      scale=1.0_r8
      status=nf_fwrite3d(ng, model, ncid, ncvarid, tindex, gtype,       &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, scale,           &
# ifdef MASKING
     &                   Amask,                                         &
# endif
     &                   A)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex
      CALL my_flush (stdout)
!
  10  FORMAT (/,' WRT_NORM3D - error while writing variable: ',a,/,11x, &
     &        'into normalization NetCDF file for time record: ',i4)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i7.7)
      RETURN

      END SUBROUTINE wrt_norm3d

# undef IR_RANGE
# undef JR_RANGE
# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE

#endif
      END MODULE normalization_mod
