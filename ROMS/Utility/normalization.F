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
!svn $Id: normalization.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE mod_grid
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
# endif
      USE mod_fourdvar
      USE mod_mixing
      USE mod_ocean
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
     &                           LBi, UBi, LBj, UBj,                    &
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
     &                           GRID(ng) % bed_thick,                  &
#  endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % z_r,                        &
     &                           GRID(ng) % z_w,                        &
# endif
     &                           MIXING(ng) % Kh,                       &
# ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
#  ifdef ADJUST_STFLUX
     &                           FORCES(ng) % b_stflx,                  &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
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
     &                           LBi, UBi, LBj, UBj,                    &
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
     &                           GRID(ng) % bed_thick,                  &
#  endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % z_r,                        &
     &                           GRID(ng) % z_w,                        &
# endif
     &                           MIXING(ng) % Kh,                       &
# ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
#  ifdef ADJUST_STFLUX
     &                           FORCES(ng) % b_stflx,                   &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
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
     &                               LBi, UBi, LBj, UBj,                &
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
     &                               VnormR, VnormU, VnormV,            &
#  ifdef ADJUST_STFLUX
     &                               HnormSTF,                          &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
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
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE set_depth_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
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
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:)
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
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng))
#   ifdef ADJUST_STFLUX
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj)
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
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ic, is, j, jc, status
# ifdef SOLVE3D
      integer :: UBt, itrc, k, kc, ntrc
# endif
      integer :: start(4), total(4)

      real(r8) :: cff, compute

      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d 
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
# endif

# include "set_bounds.h"

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
!  Compute background/model error covariance, B, normalization factors
!  using the exact method. It involves computing the filter variance
!  (convolution) at each point independenly.  That is, each point is
!  perturbed with a delta function, scaled by the inverse squared
!  root of the area (2D) or volume (3D), and then convoluted.  
!-----------------------------------------------------------------------
!
!  Set time record index to write in normalization NetCDF file.
!
      tNRMindx(ng)=tNRMindx(ng)+1
      NrecNRM(ng)=NrecNRM(ng)+1
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        start(1)=tNRMindx(ng)
        total(1)=1
        status=nf90_put_var(ncNRMid(ng), nrmVid(idtime,ng),             &
     &                      time(ng:), start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tNRMindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Master) WRITE (stdout,20)
!
!  2D norm at RHO-points.
!
      IF (Cnorm(isFsur)) THEN
        Imin=1
        Imax=Lm(ng)
        Jmin=1
        Jmax=Mm(ng)
        IF (Master) THEN
          WRITE (stdout,30) '2D normalization factors at RHO-points'
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
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                A2d(ic,jc)=Hscale(ic,jc)
              END IF
              CALL ad_conv_r2d_tile (ng, tile, iADM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isFsur)/ifac,              &
     &                               DTsizeH(isFsur),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_u, pnom_v,            &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
     &                               A2d)
              DO j=JR_RANGE
                DO i=IR_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_r2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isFsur)/ifac,              &
     &                               DTsizeH(isFsur),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_u, pnom_v,            &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
     &                               A2d)
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                cff=1.0_r8/SQRT(A2d(ic,jc))
              END IF
            ELSE
              cff=0.0_r8
            END IF
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              HnormR(ic,jc)=cff
            END IF
          END DO
        END DO
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormR)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormR)
# endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idFsur,                    &
     &                   ncNRMid(ng), nrmVid(idFsur,ng), tNRMindx(ng),  &
# ifdef MASKING
     &                   rmask,                                         &
# endif
     &                   HnormR)
      END IF
!
!  2D norm at U-points.
!
      IF (Cnorm(isUbar)) THEN
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
          WRITE (stdout,30) '2D normalization factors at   U-points'
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
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                A2d(ic,jc)=Hscale(ic,jc)
              END IF
              CALL ad_conv_u2d_tile (ng, tile, iADM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isUbar)/ifac,              &
     &                               DTsizeH(isUbar),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_r, pnom_p,            &
# ifdef MASKING
     &                               umask, pmask,                      &
# endif
     &                               A2d)
              DO j=JU_RANGE
                DO i=IU_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_u2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isUbar)/ifac,              &
     &                               DTsizeH(isUbar),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_r, pnom_p,            &
# ifdef MASKING
     &                               umask, pmask,                      &
# endif
     &                               A2d)
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                cff=1.0_r8/SQRT(A2d(ic,jc))
              END IF
            ELSE
              cff=0.0_r8
            END IF
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              HnormU(ic,jc)=cff
            END IF
          END DO
        END DO
        CALL bc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormU)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormU)
# endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idUbar,                    &
     &                   ncNRMid(ng), nrmVid(idUbar,ng), tNRMindx(ng),  &
# ifdef MASKING
     &                   umask,                                         &
# endif
     &                   HnormU)
      END IF
!
!  2D norm at V-points.
!
      IF (Cnorm(isVbar)) THEN
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
          WRITE (stdout,30) '2D normalization factors at   V-points'
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
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                A2d(ic,jc)=Hscale(ic,jc)
              END IF
              CALL ad_conv_v2d_tile (ng, tile, iADM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isVbar)/ifac,              &
     &                               DTsizeH(isVbar),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_p, pnom_r,            &
# ifdef MASKING
     &                               vmask, pmask,                      &
# endif
     &                               A2d)
              DO j=JV_RANGE
                DO i=IV_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_v2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isVbar)/ifac,              &
     &                               DTsizeH(isVbar),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_p, pnom_r,            &
# ifdef MASKING
     &                               vmask, pmask,                      &
# endif
     &                               A2d)
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                cff=1.0_r8/SQRT(A2d(ic,jc))
              END IF
            ELSE
              cff=0.0_r8
            END IF
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              HnormV(ic,jc)=cff
            END IF
          END DO
        END DO
        CALL bc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormV)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormV)
# endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idVbar,                    &
     &                   ncNRMid(ng), nrmVid(idVbar,ng), tNRMindx(ng),  &
# ifdef MASKING
     &                   vmask,                                         &
# endif
     &                   HnormV)
      END IF

# ifdef ADJUST_WSTRESS
!
!  2D norm at U-stress points.
!
      IF (Cnorm(isUstr)) THEN
#  ifdef EW_PERIODIC
        Imin=1
        Imax=Lm(ng)
        Jmin=1
        Jmax=Mm(ng)
#  else
        Imin=2
        Imax=Lm(ng)
        Jmin=1
        Jmax=Mm(ng)
#  endif
        IF (Master) THEN
          WRITE (stdout,30)                                             &
     &          '2D normalization factors at SUS-stress points'
          CALL my_flush (stdout)
        END IF
        DO j=JU_RANGE
          DO i=IU_RANGE
            Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
          END DO
        END DO
        DO jc=Jmin,Jmax
          DO ic=Imin,Imax
#  ifdef MASKING
            compute=0.0_r8
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              IF (umask(ic,jc).gt.0) compute=1.0_r8
            END IF
#   ifdef DISTRIBUTE
            CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
            compute=1.0_r8
#  endif
            IF (compute.gt.0.0_r8) THEN
              DO j=LBj,UBj
                DO i=LBi,UBi
                  A2d(i,j)=0.0_r8
                END DO
              END DO
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                A2d(ic,jc)=Hscale(ic,jc)
              END IF
              CALL ad_conv_u2d_tile (ng, tile, iADM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isUstr)/ifac,              &
     &                               DTsizeH(isUstr),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_r, pnom_p,            &
#  ifdef MASKING
     &                               umask, pmask,                      &
#  endif
     &                               A2d)
              DO j=JU_RANGE
                DO i=IU_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_u2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isUstr)/ifac,              &
     &                               DTsizeH(isUstr),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_r, pnom_p,            &
#  ifdef MASKING
     &                               umask, pmask,                      &
#  endif
     &                               A2d)
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                cff=1.0_r8/SQRT(A2d(ic,jc))
              END IF
            ELSE
              cff=0.0_r8
            END IF
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              HnormSUS(ic,jc)=cff
            END IF
          END DO
        END DO
        CALL bc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormSUS)
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormSUS)
#  endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idUsms,                    &
     &                   ncNRMid(ng), nrmVid(idUsms,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   umask,                                         &
#  endif
     &                   HnormSUS)
      END IF
!
!  2D norm at V-stress points.
!
      IF (Cnorm(isVstr)) THEN
#  ifdef NS_PERIODIC
        Imin=1
        Imax=Lm(ng)
        Jmin=1
        Jmax=Mm(ng)
#  else
        Imin=1
        Imax=Lm(ng)
        Jmin=2
        Jmax=Mm(ng)
#  endif
        IF (Master) THEN
          WRITE (stdout,30)                                             &
     &          '2D normalization factors at SVS-stress points'
          CALL my_flush (stdout)
        END IF
        DO j=JV_RANGE
          DO i=IV_RANGE
            Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
          END DO
        END DO
        DO jc=Jmin,Jmax
          DO ic=Imin,Imax
#  ifdef MASKING
            compute=0.0_r8
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              IF (vmask(ic,jc).gt.0) compute=1.0_r8
            END IF
#   ifdef DISTRIBUTE
            CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
            compute=1.0_r8
#  endif
            IF (compute.gt.0.0_r8) THEN
              DO j=LBj,UBj
                DO i=LBi,UBi
                  A2d(i,j)=0.0_r8
                END DO
              END DO
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                A2d(ic,jc)=Hscale(ic,jc)
              END IF
              CALL ad_conv_v2d_tile (ng, tile, iADM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isVstr)/ifac,              &
     &                               DTsizeH(isVstr),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_p, pnom_r,            &
#  ifdef MASKING
     &                               vmask, pmask,                      &
#  endif
     &                               A2d)
              DO j=JV_RANGE
                DO i=IV_RANGE
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_conv_v2d_tile (ng, tile, iTLM,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               NHsteps(isVstr)/ifac,              &
     &                               DTsizeH(isVstr),                   &
     &                               Kh,                                &
     &                               pm, pn, pmon_p, pnom_r,            &
#  ifdef MASKING
     &                               vmask, pmask,                      &
#  endif
     &                               A2d)
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                cff=1.0_r8/SQRT(A2d(ic,jc))
              END IF
            ELSE
              cff=0.0_r8
            END IF
            IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                    &
     &          ((Istr.le.ic).and.(ic.le.Iend))) THEN
              HnormSVS(ic,jc)=cff
            END IF
          END DO
        END DO
        CALL bc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormSVS)
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormSVS)
#  endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idVsms,                    &
     &                   ncNRMid(ng), nrmVid(idVsms,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   vmask,                                         &
#  endif
     &                   HnormSVS)
      END IF
# endif

# ifdef SOLVE3D
!
!  3D norm U-points.
!
      IF (Cnorm(isUvel)) THEN
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
          WRITE (stdout,30) '3D normalization factors at   U-points'
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
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  A3d(ic,jc,kc)=Vscale(ic,jc,kc)
                END IF
                CALL ad_conv_u3d_tile (ng, tile, iADM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 NghostPoints,                    &
     &                                 NHsteps(isUvel)/ifac,            &
     &                                 NVsteps(isUvel)/ifac,            &
     &                                 DTsizeH(isUvel),                 &
     &                                 DTsizeV(isUvel),                 &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_r, pnom_p,          &
#  ifdef MASKING
     &                                 umask, pmask,                    &
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                DO k=1,N(ng)
                  DO j=JU_RANGE
                    DO i=IU_RANGE
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
                CALL tl_conv_u3d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 NghostPoints,                    &
     &                                 NHsteps(isUvel)/ifac,            &
     &                                 NVsteps(isUvel)/ifac,            &
     &                                 DTsizeH(isUvel),                 &
     &                                 DTsizeV(isUvel),                 &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_r, pnom_p,          &
#  ifdef MASKING
     &                                 umask, pmask,                    &
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  cff=1.0_r8/SQRT(A3d(ic,jc,kc))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                VnormU(ic,jc,kc)=cff
              END IF
            END DO
          END DO
        END DO
        CALL bc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    VnormU)
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      VnormU)
#  endif
        CALL wrt_norm3d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), idUvel,          &
     &                   ncNRMid(ng), nrmVid(idUvel,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   umask,                                         &
#  endif
     &                   VnormU)
      END IF
!
!  3D norm at V-points.
!
      IF (Cnorm(isVvel)) THEN
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
          WRITE (stdout,30) '3D normalization factors at   V-points'
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
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  A3d(ic,jc,kc)=Vscale(ic,jc,kc)
                END IF
                CALL ad_conv_v3d_tile (ng, tile, iADM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 NghostPoints,                    &
     &                                 NHsteps(isVvel)/ifac,            &
     &                                 NVsteps(isVvel)/ifac,            &
     &                                 DTsizeH(isVvel),                 &
     &                                 DTsizeV(isVvel),                 &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_p, pnom_r,          &
#  ifdef MASKING
     &                                 vmask, pmask,                    &
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                DO k=1,N(ng)
                  DO j=JV_RANGE
                    DO i=IV_RANGE
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
                CALL tl_conv_v3d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 NghostPoints,                    &
     &                                 NHsteps(isVvel)/ifac,            &
     &                                 NVsteps(isVvel)/ifac,            &
     &                                 DTsizeH(isVvel),                 &
     &                                 DTsizeV(isVvel),                 &
     &                                 Kh, Kv,                          &
     &                                 pm, pn, pmon_p, pnom_r,          &
#  ifdef MASKING
     &                                 vmask, pmask,                    &
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  cff=1.0_r8/SQRT(A3d(ic,jc,kc))
                END IF
              ELSE
                cff=0.0_r8
              END IF
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                VnormV(ic,jc,kc)=cff
              END IF
            END DO
          END DO
        END DO
        CALL bc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    VnormV)
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      VnormV)
#  endif
        CALL wrt_norm3d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), idVvel,          &
     &                   ncNRMid(ng), nrmVid(idVvel,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   vmask,                                         &
#  endif
     &                   VnormV)
      END IF
!
!  3D norm at RHO-points.
!
      IF (Master) THEN
        Lsame=.FALSE.
        DO itrc=1,NT(ng)
          is=isTvar(itrc)
          IF (Cnorm(is)) Lsame=.TRUE.
        END DO
        IF (Lsame) THEN
          WRITE (stdout,30) '3D normalization factors at RHO-points'
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
        IF ((Hdecay(isTvar(itrc  ),ng).ne.                              &
     &       Hdecay(isTvar(itrc-1),ng)).or.                             &
     &      (Vdecay(isTvar(itrc  ),ng).ne.                              &
     &       Vdecay(isTvar(itrc-1),ng))) THEN
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
        IF (Cnorm(is)) THEN
          DO kc=1,N(ng)
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
#  ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A3d(ic,jc,kc)=Vscale(ic,jc,kc)
                  END IF
                  CALL ad_conv_r3d_tile (ng, tile, iADM,                &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   NghostPoints,                  &
     &                                   NHsteps(is)/ifac,              &
     &                                   NVsteps(is)/ifac,              &
     &                                   DTsizeH(is),                   &
     &                                   DTsizeV(is),                   &
     &                                   Kh, Kv,                        &
     &                                   pm, pn, pmon_u, pnom_v,        &
#  ifdef MASKING
     &                                   rmask, umask, vmask,           &
#  endif
     &                                   Hz, z_r,                       &
     &                                   A3d)
                  DO k=1,N(ng)
                    DO j=JR_RANGE
                      DO i=IR_RANGE
                        A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                      END DO
                    END DO
                  END DO
                  CALL tl_conv_r3d_tile (ng, tile, iTLM,                &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   NghostPoints,                  &
     &                                   NHsteps(is)/ifac,              &
     &                                   NVsteps(is)/ifac,              &
     &                                   DTsizeH(is),                   &
     &                                   DTsizeV(is),                   &
     &                                   Kh, Kv,                        &
     &                                   pm, pn, pmon_u, pnom_v,        &
#  ifdef MASKING
     &                                   rmask, umask, vmask,           &
#  endif
     &                                   Hz, z_r,                       &
     &                                   A3d)
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    cff=1.0_r8/SQRT(A3d(ic,jc,kc))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (Lsame) THEN
                    DO ntrc=1,NT(ng)
                      VnormR(ic,jc,kc,ntrc)=cff
                    END DO
                  ELSE
                    VnormR(ic,jc,kc,itrc)=cff
                  END IF
                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
      DO itrc=1,NT(ng)
        is=isTvar(itrc)
        IF (Cnorm(is)) THEN
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      VnormR(:,:,:,itrc))
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        VnormR(:,:,:,itrc))
#  endif
          CALL wrt_norm3d (ng, tile, iTLM,                              &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), idTvar(itrc),  &
     &                     ncNRMid(ng), nrmTid(itrc,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                     rmask,                                       &
#  endif
     &                     VnormR(:,:,:,itrc))
        END IF
      END DO

#  ifdef ADJUST_STFLUX
!
!  2D norm at surface tracer fluxes points.
!
      IF (Master) THEN
        Lsame=.FALSE.
        DO itrc=1,NT(ng)
          is=isTsur(itrc)
          IF (Cnorm(is)) Lsame=.TRUE.
        END DO
        IF (Lsame) THEN
          WRITE (stdout,30) '2D normalization factors at STF-points'
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
        IF (Hdecay(isTsur(itrc  ),ng).ne.                               &
     &      Hdecay(isTsur(itrc-1),ng)) THEN
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
        is=isTsur(itrc)
        IF (Cnorm(is)) THEN
          DO jc=Jmin,Jmax
            DO ic=Imin,Imax
#   ifdef MASKING
              compute=0.0_r8
              IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                  &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
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
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  A2d(ic,jc)=Hscale(ic,jc)
                END IF
                CALL ad_conv_r2d_tile (ng, tile, iADM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 NghostPoints,                    &
     &                                 NHsteps(is)/ifac,                &
     &                                 DTsizeH(is),                     &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_u, pnom_v,          &
#   ifdef MASKING
     &                                 rmask, umask, vmask,             &
#   endif
     &                                 A2d)
                DO j=JR_RANGE
                  DO i=IR_RANGE
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
                CALL tl_conv_r2d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 NghostPoints,                    &
     &                                 NHsteps(is)/ifac,                &
     &                                 DTsizeH(is),                     &
     &                                 Kh,                              &
     &                                 pm, pn, pmon_u, pnom_v,          &
#   ifdef MASKING
     &                                 rmask, umask, vmask,             &
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
                IF (Lsame) THEN
                  DO ntrc=1,NT(ng)
                    HnormSTF(ic,jc,ntrc)=cff
                  END DO
                ELSE
                  HnormSTF(ic,jc,itrc)=cff
                END IF
              END IF
            END DO
          END DO
        END IF
      END DO
      DO itrc=1,NT(ng)
        is=isTsur(itrc)
        IF (Cnorm(is)) THEN
          CALL bc_r2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      HnormSTF(:,:,itrc))
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        HnormSTF(:,:,itrc))
#   endif
          CALL wrt_norm2d (ng, tile, iTLM,                              &
     &                     LBi, UBi, LBj, UBj, idTsur(itrc),            &
     &                     ncNRMid(ng), nrmVid(idTsur(itrc),ng),        &
     &                     tNRMindx(ng),                                &
#   ifdef MASKING
     &                     rmask,                                       &
#   endif
     &                     HnormSTF(:,:,itrc))
        END IF
      END DO
#  endif
# endif
!
      IF (Master) THEN
        WRITE (stdout,40)
      END IF

 10   FORMAT (/,' NORMALIZATION - error while writing variable: ',a,/,  &
     &        11x,'into normalization NetCDF file for time record: ',i4)
 20   FORMAT (/,' Background Error Covariance Factors: Exact Method',/)
 30   FORMAT (4x,'Computing ',a)
 40   FORMAT (/)

      RETURN
      END SUBROUTINE normalization_tile

!
!***********************************************************************
      SUBROUTINE randomization_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
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
     &                               VnormR, VnormU, VnormV,            &
#  ifdef ADJUST_STFLUX
     &                               HnormSTF,                          &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
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
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
      USE set_depth_mod
      USE white_noise_mod, ONLY : white_noise2d
# ifdef SOLVE3D
      USE white_noise_mod, ONLY : white_noise3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
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
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:)
#   ifdef ADJUST_STFLUX
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:)
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
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng))
#   ifdef ADJUST_STFLUX
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj)
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
      integer :: i, is, iter, j, status
# ifdef SOLVE3D
      integer :: UBt, itrc, k
# endif
      integer :: start(4), total(4)

      real(r8) :: Aavg, Amax, Amin, FacAvg, FacSqr, cff, val

      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2davg
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
# ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3davg
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
# endif

# include "set_bounds.h"

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
!  Compute background covariance, B, normalization factors using the 
!  randomization approach of Fisher and Courtier (1995). These factors
!  ensure that the diagonal elements of B are equal to unity.
!  Notice that in applications with land/sea masking, the boundary
!  conditions will produce large changes in the covariance structures
!  near the boundary.
!
!  Initialize factors with randon numbers ("white-noise") having an
!  uniform distribution (zero mean and unity variance). Then, scale
!  by the inverse squared root area (2D) or volume (3D) and "color"
!  with the diffusion operator. Iterate this step over a specified
!  number of ensamble members, Nrandom.
!-----------------------------------------------------------------------
!
!  Set randomization summation factors.
!
      FacAvg=1.0_r8/REAL(Nrandom,r8)
      FacSqr=SQRT(REAL(Nrandom-1,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
      tNRMindx(ng)=tNRMindx(ng)+1
      NrecNRM(ng)=NrecNRM(ng)+1
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        start(1)=tNRMindx(ng)
        total(1)=1
        status=nf90_put_var(ncNRMid(ng), nrmVid(idtime,ng),             &
     &                      time(ng:), start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tNRMindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Master) WRITE (stdout,20)
!
!  2D norm at RHO-points.
!
      IF (Cnorm(isFsur)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '2D normalization factors at RHO-points'
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
          CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),            & 
     &                        IstrR, IendR, JstrR, JendR,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Amin, Amax, A2d)
          DO j=JR_RANGE
            DO i=IR_RANGE
              A2d(i,j)=A2d(i,j)*Hscale(i,j)
            END DO
          END DO
          CALL tl_conv_r2d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           NHsteps(isFsur)/ifac,                  &
     &                           DTsizeH(isFsur),                       &
     &                           Kh,                                    &
     &                           pm, pn, pmon_u, pnom_v,                &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
     &                           A2d)
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
!!            HnormR(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormR(i,j)=FacSqr/SQRT(A2dsqr(i,j))
            ELSE
              HnormR(i,j)=0.0_r8
            END IF
# else
!!          HnormR(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
            HnormR(i,j)=FacSqr/SQRT(A2dsqr(i,j))
# endif
          END DO
        END DO
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormR)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormR)
# endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idFsur,                    &
     &                   ncNRMid(ng), nrmVid(idFsur,ng), tNRMindx(ng),  &
# ifdef MASKING
     &                   rmask,                                         &
# endif
     &                   HnormR)
      END IF
!
!  2D norm at U-points.
!
      IF (Cnorm(isUbar)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '2D normalization factors at   U-points'
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
          CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),            & 
     &                        Istr, IendR, JstrR, JendR,                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Amin, Amax, A2d)
          DO j=JU_RANGE
            DO i=IU_RANGE
              A2d(i,j)=A2d(i,j)*Hscale(i,j)
            END DO
          END DO
          CALL tl_conv_u2d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           NHsteps(isUbar)/ifac,                  &
     &                           DTsizeH(isUbar),                       &
     &                           Kh,                                    &
     &                           pm, pn, pmon_r, pnom_p,                &
# ifdef MASKING
     &                           umask, pmask,                          &
# endif
     &                           A2d)
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
!!            HnormU(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormU(i,j)=FacSqr/SQRT(A2dsqr(i,j))
            ELSE
              HnormU(i,j)=0.0_r8
            END IF
# else
!!          HnormU(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
            HnormU(i,j)=FacSqr/SQRT(A2dsqr(i,j))
# endif
          END DO
        END DO
        CALL bc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormU)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormU)
# endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idUbar,                    &
     &                   ncNRMid(ng), nrmVid(idUbar,ng), tNRMindx(ng),  &
# ifdef MASKING
     &                   umask,                                         &
# endif
     &                   HnormU)
      END IF
!
!  2D norm at V-points.
!
      IF (Cnorm(isVbar)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '2D normalization factors at   V-points'
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
          CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),            & 
     &                        IstrR, IendR, Jstr, JendR,                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Amin, Amax, A2d)
          DO j=JV_RANGE
            DO i=IV_RANGE
              A2d(i,j)=A2d(i,j)*Hscale(i,j)
            END DO
          END DO
          CALL tl_conv_v2d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           NHsteps(isVbar)/ifac,                  &
     &                           DTsizeH(isVbar),                       &
     &                           Kh,                                    &
     &                           pm, pn, pmon_p, pnom_r,                &
# ifdef MASKING
     &                           vmask, pmask,                          &
# endif
     &                           A2d)
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
!!            HnormV(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormV(i,j)=FacSqr/SQRT(A2dsqr(i,j))
            ELSE
              HnormV(i,j)=0.0_r8
            END IF
# else
!!          HnormV(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
            HnormV(i,j)=FacSqr/SQRT(A2dsqr(i,j))
# endif
          END DO
        END DO
        CALL bc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormV)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormV)
# endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idVbar,                    &
     &                   ncNRMid(ng), nrmVid(idVbar,ng), tNRMindx(ng),  &
# ifdef MASKING
     &                   vmask,                                         &
# endif
     &                   HnormV)
      END IF

# ifdef ADJUST_WSTRESS
!
!  2D norm at U-stress points.
!
      IF (Cnorm(isUstr)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '2D normalization factors at SUS-points'
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
          CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),            & 
     &                        Istr, IendR, JstrR, JendR,                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Amin, Amax, A2d)
          DO j=JU_RANGE
            DO i=IU_RANGE
              A2d(i,j)=A2d(i,j)*Hscale(i,j)
            END DO
          END DO
          CALL tl_conv_u2d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           NHsteps(isUstr)/ifac,                  &
     &                           DTsizeH(isUstr),                       &
     &                           Kh,                                    &
     &                           pm, pn, pmon_r, pnom_p,                &
#  ifdef MASKING
     &                           umask, pmask,                          &
#  endif
     &                           A2d)
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
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
!!            HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
            ELSE
              HnormSUS(i,j)=0.0_r8
            END IF
#  else
!!          HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
            HnormSUS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
#  endif
          END DO
        END DO
        CALL bc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormSUS)
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormSUS)
#  endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idUsms,                    &
     &                   ncNRMid(ng), nrmVid(idUsms,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   umask,                                         &
#  endif
     &                   HnormSUS)
      END IF
!
!  2D norm at V-stress points.
!
      IF (Cnorm(isVstr)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '2D normalization factors at SVS-points'
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
          CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),            & 
     &                        IstrR, IendR, Jstr, JendR,                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Amin, Amax, A2d)
          DO j=JV_RANGE
            DO i=IV_RANGE
              A2d(i,j)=A2d(i,j)*Hscale(i,j)
            END DO
          END DO
          CALL tl_conv_v2d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           NHsteps(isVstr)/ifac,                  &
     &                           DTsizeH(isVstr),                       &
     &                           Kh,                                    &
     &                           pm, pn, pmon_p, pnom_r,                &
#  ifdef MASKING
     &                           vmask, pmask,                          &
#  endif
     &                           A2d)
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
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
!!            HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
            ELSE
              HnormSVS(i,j)=0.0_r8
            END IF
#  else
!!          HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
            HnormSVS(i,j)=FacSqr/SQRT(A2dsqr(i,j))
#  endif
          END DO
        END DO
        CALL bc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    HnormSVS)
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      HnormSVS)
#  endif
        CALL wrt_norm2d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, idVsms,                    &
     &                   ncNRMid(ng), nrmVid(idVsms,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   vmask,                                         &
#  endif
     &                   HnormSVS)
      END IF
# endif

# ifdef SOLVE3D
!
!  3D norm U-points.
!
      IF (Cnorm(isUvel)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '3D normalization factors at   U-points'
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
          CALL white_noise3d (ng, iTLM, u3dvar, Rscheme(ng),            & 
     &                        Istr, IendR, JstrR, JendR,                &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Amin, Amax, A3d)
          DO k=1,N(ng)
            DO j=JU_RANGE
              DO i=IU_RANGE
                A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
              END DO
            END DO
          END DO
          CALL tl_conv_u3d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints,                          &
     &                           NHsteps(isUvel)/ifac,                  &
     &                           NVsteps(isUvel)/ifac,                  &
     &                           DTsizeH(isUvel),                       &
     &                           DTsizeV(isUvel),                       &
     &                           Kh, Kv,                                &
     &                           pm, pn, pmon_r, pnom_p,                &
#  ifdef MASKING
     &                           umask, pmask,                          &
#  endif
     &                           Hz, z_r,                               &
     &                           A3d)
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
!!              VnormU(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                VnormU(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k))
              ELSE
                VnormU(i,j,k)=0.0_r8
              END IF
#  else
!!            VnormU(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
              VnormU(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k))
#  endif
            END DO
          END DO
        END DO
        CALL bc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    VnormU)
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      VnormU)
#  endif
        CALL wrt_norm3d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), idUvel,          &
     &                   ncNRMid(ng), nrmVid(idUvel,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   umask,                                         &
#  endif
     &                   VnormU)
      END IF
!
!  3D norm at V-points.
!
      IF (Cnorm(isVvel)) THEN
        IF (Master) THEN
          WRITE (stdout,30) '3D normalization factors at   V-points'
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
          CALL white_noise3d (ng, iTLM, v3dvar, Rscheme(ng),            & 
     &                        IstrR, IendR, Jstr, JendR,                &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Amin, Amax, A3d)
          DO k=1,N(ng)
            DO j=JV_RANGE
              DO i=IV_RANGE
                A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
              END DO
            END DO
          END DO
          CALL tl_conv_v3d_tile (ng, tile, iTLM,                        &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints,                          &
     &                           NHsteps(isVvel)/ifac,                  &
     &                           NVsteps(isVvel)/ifac,                  &
     &                           DTsizeH(isVvel),                       &
     &                           DTsizeV(isVvel),                       &
     &                           Kh, Kv,                                &
     &                           pm, pn, pmon_p, pnom_r,                &
#  ifdef MASKING
     &                           vmask, pmask,                          &
#  endif
     &                           Hz, z_r,                               &
     &                           A3d)
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
!!              VnormV(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                VnormV(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k))
              ELSE
                VnormV(i,j,k)=0.0_r8
              END IF
#  else
!!            VnormV(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
              VnormV(i,j,k)=FacSqr/SQRT(A3dsqr(i,j,k))
#  endif
            END DO
          END DO
        END DO
        CALL bc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    VnormV)
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      VnormV)
#  endif
        CALL wrt_norm3d (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), idVvel,          &
     &                   ncNRMid(ng), nrmVid(idVvel,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                   vmask,                                         &
#  endif
     &                   VnormV)
      END IF
!
!  3D norm at RHO-points.
!
      IF (Master) THEN
        Lsame=.FALSE.
        DO itrc=1,NT(ng)
          is=isTvar(itrc)
          IF (Cnorm(is)) Lsame=.TRUE.
        END DO
        IF (Lsame) THEN
          WRITE (stdout,30) '3D normalization factors at RHO-points'
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
        IF ((Hdecay(isTvar(itrc  ),ng).ne.                              &
     &       Hdecay(isTvar(itrc-1),ng)).or.                             &
     &      (Vdecay(isTvar(itrc  ),ng).ne.                              &
     &       Vdecay(isTvar(itrc-1),ng))) THEN
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
        IF (Cnorm(is)) THEN
          DO k=1,N(ng)
            DO j=JR_RANGE
              DO i=IR_RANGE
                A3davg(i,j,k)=0.0_r8
                A3dsqr(i,j,k)=0.0_r8
              END DO
            END DO
          END DO
          DO iter=1,Nrandom
            CALL white_noise3d (ng, iTLM, r3dvar, Rscheme(ng),          & 
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Amin, Amax, A3d)
            DO k=1,N(ng)
              DO j=JR_RANGE
                DO i=IR_RANGE
                  A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                END DO
              END DO
            END DO
            CALL tl_conv_r3d_tile (ng, tile, iTLM,                      &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             NghostPoints,                        &
     &                             NHsteps(is)/ifac,                    &
     &                             NVsteps(is)/ifac,                    &
     &                             DTsizeH(is),                         &
     &                             DTsizeV(is),                         &
     &                             Kh, Kv,                              &
     &                             pm, pn, pmon_u, pnom_v,              &
#  ifdef MASKING
     &                             rmask, umask, vmask,                 &
#  endif
     &                             Hz, z_r,                             &
     &                             A3d)
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
!!                VnormR(i,j,k,itrc)=FacSqr/                            &
!!   &                               SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                  VnormR(i,j,k,itrc)=FacSqr/SQRT(A3dsqr(i,j,k))
                ELSE
                  VnormR(i,j,k,itrc)=0.0_r8
                END IF
#  else
!!              VnormR(i,j,k,itrc)=FacSqr/SQRT(A3dsqr(i,j,k)-Aavg*Aavg)
                VnormR(i,j,k,itrc)=FacSqr/SQRT(A3dsqr(i,j,k))
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
                VnormR(i,j,k,itrc)=VnormR(i,j,k,1)
              END DO
            END DO
          END DO
        END DO
      END IF
      DO itrc=1,NT(ng)
        is=isTvar(itrc)
        IF (Cnorm(is)) THEN
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      VnormR(:,:,:,itrc))
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        VnormR(:,:,:,itrc))
#  endif
          CALL wrt_norm3d (ng, tile, iTLM,                              &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), idTvar(itrc),  &
     &                     ncNRMid(ng), nrmTid(itrc,ng), tNRMindx(ng),  &
#  ifdef MASKING
     &                     rmask,                                       &
#  endif
     &                     VnormR(:,:,:,itrc))
        END IF
      END DO

#  ifdef ADJUST_STFLUX
!
!  2D norm at surface treace flux points.
!
      IF (Master) THEN
        Lsame=.FALSE.
        DO itrc=1,NT(ng)
          is=isTsur(itrc)
          IF (Cnorm(is)) Lsame=.TRUE.
        END DO
        IF (Lsame) THEN
          WRITE (stdout,30) '3D normalization factors at STF-points'
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
        IF (Hdecay(isTvar(itrc  ),ng).ne.                               &
     &      Hdecay(isTvar(itrc-1),ng)) THEN
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
        is=isTsur(itrc)
        IF (Cnorm(is)) THEN
          DO j=JR_RANGE
            DO i=IR_RANGE
              A2davg(i,j)=0.0_r8
              A2dsqr(i,j)=0.0_r8
            END DO
          END DO
          DO iter=1,Nrandom
            CALL white_noise2d (ng, iTLM, r3dvar, Rscheme(ng),          & 
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Amin, Amax, A2d)
            DO j=JR_RANGE
              DO i=IR_RANGE
                A2d(i,j)=A2d(i,j)*Hscale(i,j)
              END DO
            END DO
            CALL tl_conv_r2d_tile (ng, tile, iTLM,                      &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             NghostPoints,                        &
     &                             NHsteps(is)/ifac,                    &
     &                             DTsizeH(is),                         &
     &                             Kh,                                  &
     &                             pm, pn, pmon_u, pnom_v,              &
#   ifdef MASKING
     &                             rmask, umask, vmask,                 &
#   endif
     &                             A2d)
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
!!              HnormSTF(i,j,itrc)=FacSqr/                              &
!!   &                             SQRT(A2dsqr(i,j)-Aavg*Aavg)
                HnormSTF(i,j,itrc)=FacSqr/SQRT(A2dsqr(i,j))
              ELSE
                HnormSTF(i,j,itrc)=0.0_r8
              END IF
#   else
!!            HnormSTF(i,j,itrc)=FacSqr/SQRT(A2dsqr(i,j)-Aavg*Aavg)
              HnormSTF(i,j,itrc)=FacSqr/SQRT(A2dsqr(i,j))
#   endif
            END DO
          END DO
        END IF
      END DO
      IF (Lsame) THEN
        DO itrc=2,NT(ng)
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                VnormR(i,j,k,itrc)=VnormR(i,j,k,1)
              END DO
            END DO
          END DO
        END DO
      END IF
      DO itrc=1,NT(ng)
        is=isTsur(itrc)
        IF (Cnorm(is)) THEN
          CALL bc_r2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      HnormSTF(:,:,itrc))
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, EWperiodic, NSperiodic,     &
     &                        HnormSTF(:,:,itrc))
#   endif
          CALL wrt_norm2d (ng, tile, iTLM,                              &
     &                     LBi, UBi, LBj, UBj, idTsur(itrc),            &
     &                     ncNRMid(ng), nrmVid(idTsur(itrc),ng),        &
     &                     tNRMindx(ng),                                &
#   ifdef MASKING
     &                     rmask,                                       &
#   endif
     &                     HnormSTF(:,:,itrc))
        END IF
      END DO
#  endif
# endif
!
      IF (Master) THEN
        WRITE (stdout,40)
      END IF

 10   FORMAT (/,' RANDOMIZATION - error while writing variable: ',a,/,  &
     &        11x,'into normalization NetCDF file for time record: ',i4)
 20   FORMAT (/,' Error Covariance Factors: Randomization Method',/)
 30   FORMAT (4x,'Computing ',a)
 40   FORMAT (/)

      RETURN
      END SUBROUTINE randomization_tile

!
!***********************************************************************
      SUBROUTINE wrt_norm2d (ng, tile, model,                           &
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, ncid, ncvarid, tindex
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

      integer :: nf_fwrite2d

      real(r8) :: scale

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 2D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize NetCDT to disk.
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
     &                   Amask(LBi,LBj),                                &
# endif
     &                   A(LBi,LBj))
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
      IF (OutThread) THEN
        status=nf90_sync(ncid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,30) TRIM(Vname(1,ifield)), tindex
        CALL my_flush (stdout)
      END IF
!
  10  FORMAT (/,' WRT_NORM2D - error while writing variable: ',a,/,11x, &
     &        'into normalization NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_NORM2D - unable to synchronize normalization ',   &
     &        'NetCDF file disk.')
  30  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i7.7)

      END SUBROUTINE wrt_norm2d

!
!***********************************************************************
      SUBROUTINE wrt_norm3d (ng, tile, model,                           &
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, ncid, ncvarid, tindex
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

      integer :: nf_fwrite3d

      real(r8) :: scale

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 3D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize NetCDT to disk.
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
     &                   Amask(LBi,LBj),                                &
# endif
     &                   A(LBi,LBj,LBk))
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
      IF (OutThread) THEN
        status=nf90_sync(ncid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,30) TRIM(Vname(1,ifield)), tindex
        CALL my_flush (stdout)
      END IF
!
  10  FORMAT (/,' WRT_NORM3D - error while writing variable: ',a,/,11x, &
     &        'into normalization NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_NORM3D - unable to synchronize normalization ',   &
     &        'NetCDF file disk.')
  30  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
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
