#include "cppdefs.h"
      MODULE mod_grid
!
!svn $Id: mod_grid.F 798 2016-05-21 17:04:58Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  IJwaterR   Water points IJ couter for RHO-points masked variables.  !
!  IJwaterU   Water points IJ couter for   U-points masked variables.  !
!  IJwaterV   Water points IJ couter for   V-points masked variables.  !
!  Hz         Thicknesses (m) of vertical RHO-points.                  !
# ifdef ADJUST_BOUNDARY
!  Hz_bry     Thicknesses (m) at the open boundaries; used only for    !
!               4DVAR adjustments.                                     !
# endif
!  Huon       Total U-momentum flux term, Hz*u/pn.                     !
!  Hvom       Total V-momentum flux term, Hz*v/pm.                     !
!  IcePress   Pressure under the ice shelf at RHO-points.              !
!  Rscope     Adjoint sensitivity spatial scope mask at RHO-points.    !
!  Tcline     Width (m) of surface or bottom boundary layer where      !
!               higher vertical resolution is required during          !
!               stretching.                                            !
!  Uscope     Adjoint sensitivity spatial scope mask at U-points.      !
!  Vscope     Adjoint sensitivity spatial scope mask at V-points.      !
!  angler     Angle (radians) between XI-axis and true EAST at         !
!               RHO-points.                                            !
!  CosAngler  Cosine of curvilinear angle, angler.                     !
!  SinAngler  Sine of curvilinear angle, angler.                       !
!  dmde       ETA-derivative of inverse metric factor pm,              !
!               d(1/pm)/d(ETA).                                        !
!  dndx       XI-derivative  of inverse metric factor pn,              !
!               d(1/pn)/d(XI).                                         !
!  f          Coriolis parameter (1/s).                                !
!  fomn       Compound term, f/(pm*pn) at RHO points.                  !
!  grdscl     Grid scale used to adjust horizontal mixing according    !
!               to grid area.                                          !
!  h          Bottom depth (m) at RHO-points.                          !
!  latp       Latitude (degrees_north) at PSI-points.                  !
!  latr       Latitude (degrees_north) at RHO-points.                  !
!  latu       Latitude (degrees_north) at U-points.                    !
!  latv       Latitude (degrees_north) at V-points.                    !
!  lonp       Longitude (degrees_east) at PSI-points.                  !
!  lonr       Longitude (degrees_east) at RHO-points.                  !
!  lonu       Longitude (degrees_east) at U-points.                    !
!  lonv       Longitude (degrees_east) at V-points.                    !
!  omm        RHO-grid area (meters2).                                 !
!  om_p       PSI-grid spacing (meters) in the XI-direction.           !
!  om_r       RHO-grid spacing (meters) in the XI-direction.           !
!  om_u       U-grid spacing (meters) in the XI-direction.             !
!  om_v       V-grid spacing (meters) in the XI-direction.             !
!  on_p       PSI-grid spacing (meters) in the ETA-direction.          !
!  on_r       RHO-grid spacing (meters) in the ETA-direction.          !
!  on_u       U-grid spacing (meters) in the ETA-direction.            !
!  on_v       V-grid spacing (meters) in the ETA-direction.            !
!  pm         Coordinate transformation metric "m" (1/meters)          !
!               associated with the differential distances in XI.      !
!  pmon_p     Compound term, pm/pn at PSI-points.                      !
!  pmon_r     Compound term, pm/pn at RHO-points.                      !
!  pmon_u     Compound term, pm/pn at U-points.                        !
!  pmon_v     Compound term, pm/pn at V-points.                        !
!  pn         Coordinate transformation metric "n" (1/meters)          !
!               associated with the differential distances in ETA.     !
!  pnom_p     Compound term, pn/pm at PSI-points.                      !
!  pnom_r     Compound term, pn/pm at RHO-points.                      !
!  pnom_u     Compound term, pn/pm at U-points.                        !
!  pnom_v     Compound term, pn/pm at V-points.                        !
#ifdef MASKING
!  pmask      Slipperiness time-independent mask at PSI-points:        !
!               (0=Land, 1=Sea, 2=no-slip).                            !
!  rmask      Time-independent mask at RHO-points (0=Land, 1=Sea).     !
!  umask      Time-independent mask at U-points (0=Land, 1=Sea).       !
!  vmask      Time-independent mask at V-points (0=Land, 1=Sea).       !
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!  pmask_avg  Time-averaged full mask at PSI-points (0=dry, 1=wet).    !
!  rmask_avg  Time-averaged full mask at RHO-points (0=dry, 1=wet).    !
!  rmask_avg  Time-averaged full mask at   U-points (0=dry, 1=wet).    !
!  rmask_avg  Time-averaged full mask at   V-points (0=dry, 1=wet).    !
# endif
# if defined AVERAGES2
!  pmask_avg2 Time-average FillValue mask at PSI-points (0=dry, 1=wet) !
!  rmask_avg2 Time-average FillValue mask at RHO-points (0=dry, 1=wet) !
!  rmask_avg2 Time-average FillValue mask at   U-points (0=dry, 1=wet) !
!  rmask_avg2 Time-average FillValue mask at   V-points (0=dry, 1=wet) !
# endif
# ifdef DIAGNOSTICS
!
!  pmask_dia  Diagnostics full mask at PSI-points (0=dry, 1=wet).      !
!  rmask_dia  Diagnostics full mask at RHO-points (0=dry, 1=wet).      !
!  rmask_dia  Diagnostics full mask at   U-points (0=dry, 1=wet).      !
!  rmask_dia  Diagnostics full mask at   V-points (0=dry, 1=wet).      !
# endif
!                                                                      !
!  pmask_full Full mask at PSI-points (0=dry, 1=wet, 2=no-slip).       !
!  rmask_full Full mask at RHO-points (0=dry, 1=wet).                  !
!  rmask_full Full mask at   U-points (0=dry, 1=wet).                  !
!  rmask_full Full mask at   V-points (0=dry, 1=wet).                  !
#endif
#ifdef OUTFLOW_MASK
!  mask_outflow  Mask at RHO-points for ice BC (0=weak, 1=strong)      !
#endif
#ifdef WET_DRY
!  pmask_wet  Wet/Dry mask at PSI-points (0=dry, 1=wet)                !
!  rmask_wet  Wet/Dry mask at RHO-points (0=dry, 1=wet)                !
!  umask_wet  Wet/Dry mask at   U-points (0=dry, 1,2=wet)              !
!  vmask_wet  Wet/Dry mask at   V-points (0=dry, 1,2=wet)              !
!
!  umask_diff Diffusion Wet/Dry mask at U-points (0=dry, 1=wet)        !
!  vmask_diff Diffusion Wet/Dry mask at V-points (0=dry, 1=wet)        !
#endif
#if defined UV_LOGDRAG || defined GLS_MIXING || \
    defined BBL_MODEL  || defined SEDIMENT
!  ZoBot      Bottom roughness length (m).                             !
#endif
#if defined UV_LDRAG
!  rdrag      Linear drag coefficient (m/s).                           !
#elif defined UV_QDRAG
!  rdrag2     Quadratic drag coefficient (nondimensional).             !
#endif
#if defined UV_WAVEDRAG
!  wavedrag   Linear drag coefficient from internal tides (m/s).       !
#endif
!  xp         XI-coordinates (m) at PSI-points.                        !
!  xr         XI-coordinates (m) at RHO-points.                        !
!  xu         XI-coordinates (m) at U-points.                          !
!  xv         XI-coordinates (m) at V-points.                          !
!  yp         ETA-coordinates (m) at PSI-points.                       !
!  yr         ETA-coordinates (m) at RHO-points.                       !
!  yu         ETA-coordinates (m) at U-points.                         !
!  yv         ETA-coordinates (m) at V-points.                         !
!  zice       Depth of ice shelf cavity (m, negative) at               !
!               RHO-points.                                            !
!  z0_r       Time independent depths (m) at horizontal RHO-points and !
!               vertical RHO-points.                                   !
!  z0_w       Time independent depths (m) at horizontal RHO-points and !
!               vertical W-points.                                     !
!  z_r        Actual depths (m) at horizontal RHO-points and           !
!               vertical RHO-points.                                   !
!  z_w        Actual depths (m) at horizontal RHO-points and           !
!               vertical W-points.                                     !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_GRID
!
!  Nonlinear model state.
!
#if defined MASKING && defined PROPAGATOR
          integer, pointer :: IJwaterR(:,:)
          integer, pointer :: IJwaterU(:,:)
          integer, pointer :: IJwaterV(:,:)
#endif

          real(r8), pointer :: angler(:,:)
          real(r8), pointer :: CosAngler(:,:)
          real(r8), pointer :: SinAngler(:,:)

#if defined CURVGRID && defined UV_ADV
          real(r8), pointer :: dmde(:,:)
          real(r8), pointer :: dndx(:,:)
#endif
          real(r8), pointer :: f(:,:)
          real(r8), pointer :: fomn(:,:)
          real(r8), pointer :: grdscl(:,:)
          real(r8), pointer :: h(:,:)
          real(r8), pointer :: latp(:,:)
          real(r8), pointer :: latr(:,:)
          real(r8), pointer :: latu(:,:)
          real(r8), pointer :: latv(:,:)
          real(r8), pointer :: lonp(:,:)
          real(r8), pointer :: lonr(:,:)
          real(r8), pointer :: lonu(:,:)
          real(r8), pointer :: lonv(:,:)
          real(r8), pointer :: omn(:,:)
          real(r8), pointer :: om_p(:,:)
          real(r8), pointer :: om_r(:,:)
          real(r8), pointer :: om_u(:,:)
          real(r8), pointer :: om_v(:,:)
          real(r8), pointer :: on_p(:,:)
          real(r8), pointer :: on_r(:,:)
          real(r8), pointer :: on_u(:,:)
          real(r8), pointer :: on_v(:,:)
          real(r8), pointer :: pm(:,:)
          real(r8), pointer :: pn(:,:)
          real(r8), pointer :: pmon_p(:,:)
          real(r8), pointer :: pmon_r(:,:)
          real(r8), pointer :: pmon_u(:,:)
          real(r8), pointer :: pmon_v(:,:)
          real(r8), pointer :: pnom_p(:,:)
          real(r8), pointer :: pnom_r(:,:)
          real(r8), pointer :: pnom_u(:,:)
          real(r8), pointer :: pnom_v(:,:)
#if defined UV_LOGDRAG || defined GLS_MIXING || \
    defined BBL_MODEL  || defined SEDIMENT
          real(r8), pointer :: ZoBot(:,:)
#endif
#if defined UV_LDRAG
          real(r8), pointer :: rdrag(:,:)
#elif defined UV_QDRAG
          real(r8), pointer :: rdrag2(:,:)
#endif
#if defined UV_WAVEDRAG
          real(r8), pointer :: wavedrag(:,:)
#endif
          real(r8), pointer :: xp(:,:)
          real(r8), pointer :: xr(:,:)
          real(r8), pointer :: xu(:,:)
          real(r8), pointer :: xv(:,:)
          real(r8), pointer :: yp(:,:)
          real(r8), pointer :: yr(:,:)
          real(r8), pointer :: yu(:,:)
          real(r8), pointer :: yv(:,:)
#ifdef SOLVE3D
          real(r8), pointer :: Hz(:,:,:)
# ifdef ADJUST_BOUNDARY
          real(r8), pointer :: Hz_bry(:,:,:)
# endif
          real(r8), pointer :: Huon(:,:,:)
          real(r8), pointer :: Hvom(:,:,:)
          real(r8), pointer :: z0_r(:,:,:)
          real(r8), pointer :: z0_w(:,:,:)
          real(r8), pointer :: z_r(:,:,:)
          real(r8), pointer :: z_v(:,:,:)
          real(r8), pointer :: z_w(:,:,:)
# ifdef ICESHELF
          real(r8), pointer :: IcePress(:,:)
          real(r8), pointer :: zice(:,:)
# endif
#endif
#ifdef MASKING
          real(r8), pointer :: pmask(:,:)
          real(r8), pointer :: rmask(:,:)
          real(r8), pointer :: umask(:,:)
          real(r8), pointer :: vmask(:,:)

# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
          real(r8), pointer :: pmask_avg(:,:)
          real(r8), pointer :: rmask_avg(:,:)
          real(r8), pointer :: umask_avg(:,:)
          real(r8), pointer :: vmask_avg(:,:)
# endif
# if defined AVERAGES2
          real(r8), pointer :: pmask_avg2(:,:)
          real(r8), pointer :: rmask_avg2(:,:)
          real(r8), pointer :: umask_avg2(:,:)
          real(r8), pointer :: vmask_avg2(:,:)
# endif
# ifdef DIAGNOSTICS
          real(r8), pointer :: pmask_dia(:,:)
          real(r8), pointer :: rmask_dia(:,:)
          real(r8), pointer :: umask_dia(:,:)
          real(r8), pointer :: vmask_dia(:,:)
# endif
          real(r8), pointer :: pmask_full(:,:)
          real(r8), pointer :: rmask_full(:,:)
          real(r8), pointer :: umask_full(:,:)
          real(r8), pointer :: vmask_full(:,:)
#endif
#ifdef OUTFLOW_MASK
          real(r8), pointer :: mask_outflow(:,:)
#endif
#ifdef WET_DRY
          real(r8), pointer :: pmask_wet(:,:)
          real(r8), pointer :: rmask_wet(:,:)
          real(r8), pointer :: umask_wet(:,:)
          real(r8), pointer :: vmask_wet(:,:)
          real(r8), pointer :: umask_diff(:,:)
          real(r8), pointer :: vmask_diff(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: rmask_wet_avg(:,:)
# endif
#endif
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
          real(r8), pointer :: Rscope(:,:)
          real(r8), pointer :: Uscope(:,:)
          real(r8), pointer :: Vscope(:,:)
#endif
#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_h(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: tl_Hz(:,:,:)
#  ifdef ADJUST_BOUNDARY
          real(r8), pointer :: tl_Hz_bry(:,:,:)
#  endif
          real(r8), pointer :: tl_Huon(:,:,:)
          real(r8), pointer :: tl_Hvom(:,:,:)
          real(r8), pointer :: tl_z_r(:,:,:)
          real(r8), pointer :: tl_z_w(:,:,:)
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_h(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: ad_Hz(:,:,:)
#  ifdef ADJUST_BOUNDARY
          real(r8), pointer :: ad_Hz_bry(:,:,:)
#  endif
          real(r8), pointer :: ad_Huon(:,:,:)
          real(r8), pointer :: ad_Hvom(:,:,:)
          real(r8), pointer :: ad_z_r(:,:,:)
          real(r8), pointer :: ad_z_w(:,:,:)
# endif
#endif
        END TYPE T_GRID

        TYPE (T_GRID), allocatable :: GRID(:)

      CONTAINS

      SUBROUTINE allocate_grid (ng, LBi, UBi, LBj, UBj, LBij, UBij)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBij, UBij
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( GRID(Ngrids) )
!
!  Nonlinear model state.
!
#if defined MASKING && defined PROPAGATOR
      allocate ( GRID(ng) % IJwaterR(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % IJwaterU(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % IJwaterV(LBi:UBi,LBj:UBj) )
#endif

      allocate ( GRID(ng) % angler(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % CosAngler(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % SinAngler(LBi:UBi,LBj:UBj) )

#if defined CURVGRID && defined UV_ADV
      allocate ( GRID(ng) % dmde(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % dndx(LBi:UBi,LBj:UBj) )
#endif

      allocate ( GRID(ng) % f(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % fomn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % grdscl(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % h(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latu(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latv(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % lonp(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % lonr(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % lonu(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % lonv(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % omn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pm(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_v(LBi:UBi,LBj:UBj) )
#if defined UV_LOGDRAG || defined GLS_MIXING || \
    defined BBL_MODEL  || defined SEDIMENT
      allocate ( GRID(ng) % ZoBot(LBi:UBi,LBj:UBj) )
#endif
#if defined UV_LDRAG
      allocate ( GRID(ng) % rdrag(LBi:UBi,LBj:UBj) )
#elif defined UV_QDRAG
      allocate ( GRID(ng) % rdrag2(LBi:UBi,LBj:UBj) )
#endif
#if defined UV_WAVEDRAG
      allocate ( GRID(ng) % wavedrag(LBi:UBi,LBj:UBj) )
#endif
      allocate ( GRID(ng) % xp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % xr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % xu(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % xv(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yu(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yv(LBi:UBi,LBj:UBj) )

#ifdef SOLVE3D
      allocate ( GRID(ng) % Hz(LBi:UBi,LBj:UBj,N(ng)) )
# ifdef ADJUST_BOUNDARY
      allocate ( GRID(ng) % Hz_bry(LBij:UBij,N(ng),4) )
# endif
      allocate ( GRID(ng) % Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z0_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z0_w(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( GRID(ng) % z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z_v(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z_w(LBi:UBi,LBj:UBj,0:N(ng)) )

# ifdef ICESHELF
      allocate ( GRID(ng) % IcePress(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % zice(LBi:UBi,LBj:UBj) )
# endif

#endif

#ifdef MASKING
      allocate ( GRID(ng) % pmask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask(LBi:UBi,LBj:UBj) )
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
      allocate ( GRID(ng) % pmask_avg(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask_avg(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask_avg(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask_avg(LBi:UBi,LBj:UBj) )
# endif
# if defined AVERAGES2
      allocate ( GRID(ng) % pmask_avg2(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask_avg2(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask_avg2(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask_avg2(LBi:UBi,LBj:UBj) )
# endif
# ifdef DIAGNOSTICS
      allocate ( GRID(ng) % pmask_dia(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask_dia(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask_dia(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask_dia(LBi:UBi,LBj:UBj) )
# endif
      allocate ( GRID(ng) % pmask_full(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask_full(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask_full(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask_full(LBi:UBi,LBj:UBj) )
#endif

#ifdef OUTFLOW_MASK
      allocate ( GRID(ng) % mask_outflow(LBi:UBi,LBj:UBj) )
#endif
#ifdef WET_DRY
      allocate ( GRID(ng) % pmask_wet(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask_wet(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask_wet(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask_wet(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask_diff(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask_diff(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( GRID(ng) % rmask_wet_avg(LBi:UBi,LBj:UBj) )
# endif
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
      allocate ( GRID(ng) % Rscope(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % Uscope(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % Vscope(LBi:UBi,LBj:UBj) )
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( GRID(ng) % tl_h(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( GRID(ng) % tl_Hz(LBi:UBi,LBj:UBj,N(ng)) )
#  ifdef ADJUST_BOUNDARY
      allocate ( GRID(ng) % tl_Hz_bry(LBij:UBij,N(ng),4) )
#  endif
      allocate ( GRID(ng) % tl_Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_z_w(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( GRID(ng) % ad_h(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( GRID(ng) % ad_Hz(LBi:UBi,LBj:UBj,N(ng)) )
#  ifdef ADJUST_BOUNDARY
      allocate ( GRID(ng) % ad_Hz_bry(LBij:UBij,N(ng),4) )
#  endif
      allocate ( GRID(ng) % ad_Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_z_w(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_grid

      SUBROUTINE initialize_grid (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
#ifdef SOLVE3D
      integer :: k
#endif

      real(r8), parameter :: IniVal = 0.0_r8
      real(r8) :: IniMetricVal

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#endif
!
!  Set initialization value that it is special in nexting to just
!  load contact points that have not been initialized from the
!  regular physical grid. This is done to make sure that all these
!  important metric values have been set-up correctly.
!
#ifdef NESTING
      IniMetricVal=spval                   ! very large value
#else
      IniMetricVal=IniVal
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
#if defined MASKING && defined PROPAGATOR
            GRID(ng) % IJwaterR(i,j) = 0
            GRID(ng) % IJwaterU(i,j) = 0
            GRID(ng) % IJwaterV(i,j) = 0
#endif
            GRID(ng) % angler(i,j) = IniMetricVal
            GRID(ng) % CosAngler(i,j) = IniVal
            GRID(ng) % SinAngler(i,j) = IniVal

#if defined CURVGRID && defined UV_ADV
            GRID(ng) % dmde(i,j) = IniMetricVal
            GRID(ng) % dndx(i,j) = IniMetricVal
#endif
            GRID(ng) % f(i,j) = IniMetricVal
            GRID(ng) % fomn(i,j) = IniVal
            GRID(ng) % grdscl(i,j) = IniVal

            GRID(ng) % h(i,j) = IniMetricVal

            GRID(ng) % latp(i,j) = IniVal
            GRID(ng) % latr(i,j) = IniMetricVal
            GRID(ng) % latu(i,j) = IniMetricVal
            GRID(ng) % latv(i,j) = IniMetricVal
            GRID(ng) % lonp(i,j) = IniVal
            GRID(ng) % lonr(i,j) = IniMetricVal
            GRID(ng) % lonu(i,j) = IniMetricVal
            GRID(ng) % lonv(i,j) = IniMetricVal

            GRID(ng) % omn(i,j) = IniVal
            GRID(ng) % om_p(i,j) = IniVal
            GRID(ng) % om_r(i,j) = IniVal
            GRID(ng) % om_u(i,j) = IniVal
            GRID(ng) % om_v(i,j) = IniVal
            GRID(ng) % on_p(i,j) = IniVal
            GRID(ng) % on_r(i,j) = IniVal
            GRID(ng) % on_u(i,j) = IniVal
            GRID(ng) % on_v(i,j) = IniVal

            GRID(ng) % pm(i,j) = IniMetricVal
            GRID(ng) % pn(i,j) = IniMetricVal

            GRID(ng) % pmon_p(i,j) = IniVal
            GRID(ng) % pmon_r(i,j) = IniVal
            GRID(ng) % pmon_u(i,j) = IniVal
            GRID(ng) % pmon_v(i,j) = IniVal
            GRID(ng) % pnom_p(i,j) = IniVal
            GRID(ng) % pnom_r(i,j) = IniVal
            GRID(ng) % pnom_u(i,j) = IniVal
            GRID(ng) % pnom_v(i,j) = IniVal

#if defined UV_LOGDRAG || defined GLS_MIXING || \
    defined BBL_MODEL  || defined SEDIMENT
            GRID(ng) % ZoBot(i,j) = Zob(ng)
#endif
#if defined UV_LDRAG
            GRID(ng) % rdrag(i,j) = rdrg(ng)
#elif defined UV_QDRAG
            GRID(ng) % rdrag2(i,j) = rdrg2(ng)
#endif
#if defined UV_WAVEDRAG
            GRID(ng) % wavedrag(i,j) = Inival
#endif

            GRID(ng) % xp(i,j) = IniVal
            GRID(ng) % xr(i,j) = IniMetricVal
            GRID(ng) % xu(i,j) = IniMetricVal
            GRID(ng) % xv(i,j) = IniMetricVal
            GRID(ng) % yp(i,j) = IniVal
            GRID(ng) % yr(i,j) = IniMetricVal
            GRID(ng) % yu(i,j) = IniMetricVal
            GRID(ng) % yv(i,j) = IniMetricVal

#if defined ICESHELF && defined SOLVE3D
            GRID(ng) % IcePress(i,j) = IniVal
            GRID(ng) % zice(i,j) = IniVal
#endif

#ifdef MASKING
            GRID(ng) % pmask(i,j) = IniVal
            GRID(ng) % rmask(i,j) = IniMetricVal
            GRID(ng) % umask(i,j) = IniMetricVal
            GRID(ng) % vmask(i,j) = IniMetricVal
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
            GRID(ng) % pmask_avg(i,j) = IniVal
            GRID(ng) % rmask_avg(i,j) = IniVal
            GRID(ng) % umask_avg(i,j) = IniVal
            GRID(ng) % vmask_avg(i,j) = IniVal
# endif
# if defined AVERAGES2
            GRID(ng) % pmask_avg2(i,j) = IniVal
            GRID(ng) % rmask_avg2(i,j) = IniVal
            GRID(ng) % umask_avg2(i,j) = IniVal
            GRID(ng) % vmask_avg2(i,j) = IniVal
# endif
# ifdef DIAGNOSTICS
            GRID(ng) % pmask_dia(i,j) = IniVal
            GRID(ng) % rmask_dia(i,j) = IniVal
            GRID(ng) % umask_dia(i,j) = IniVal
            GRID(ng) % vmask_dia(i,j) = IniVal
# endif
            GRID(ng) % pmask_full(i,j) = IniVal
            GRID(ng) % rmask_full(i,j) = IniVal
            GRID(ng) % umask_full(i,j) = IniVal
            GRID(ng) % vmask_full(i,j) = IniVal
#endif

#ifdef OUTFLOW_MASK
            GRID(ng) % mask_outflow(i,j) = IniVal
#endif
#ifdef WET_DRY
            GRID(ng) % pmask_wet(i,j) = IniVal
            GRID(ng) % rmask_wet(i,j) = IniVal
            GRID(ng) % umask_wet(i,j) = IniVal
            GRID(ng) % vmask_wet(i,j) = IniVal
            GRID(ng) % umask_diff(i,j) = IniVal
            GRID(ng) % vmask_diff(i,j) = IniVal
# ifdef SOLVE3D
            GRID(ng) % rmask_wet_avg(i,j) = IniVal
# endif
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
            GRID(ng) % Rscope(i,j) = IniVal
            GRID(ng) % Uscope(i,j) = IniVal
            GRID(ng) % Vscope(i,j) = IniVal
#endif
          END DO

#ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              GRID(ng) % Hz(i,j,k) = IniVal
              GRID(ng) % Huon(i,j,k) = IniVal
              GRID(ng) % Hvom(i,j,k) = IniVal
              GRID(ng) % z0_r(i,j,k) = IniVal
              GRID(ng) % z_r(i,j,k) = IniVal
              GRID(ng) % z_v(i,j,k) = IniVal
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              GRID(ng) % z0_w(i,j,k) = IniVal
              GRID(ng) % z_w(i,j,k) = IniVal
            END DO
          END DO
#endif
        END DO
#if defined ADJUST_BOUNDARY && defined SOLVE3D
        GRID(ng) % Hz_bry = IniVal
#endif
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            GRID(ng) % tl_h(i,j) = IniVal
          END DO
# ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              GRID(ng) % tl_Hz(i,j,k) = IniVal
              GRID(ng) % tl_Huon(i,j,k) = IniVal
              GRID(ng) % tl_Hvom(i,j,k) = IniVal
              GRID(ng) % tl_z_r(i,j,k) = IniVal
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              GRID(ng) % tl_z_w(i,j,k) = IniVal
            END DO
          END DO
# endif
        END DO
# if defined ADJUST_BOUNDARY && defined SOLVE3D
        GRID(ng) % tl_Hz_bry = IniVal
# endif
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            GRID(ng) % ad_h(i,j) = IniVal
          END DO
# ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              GRID(ng) % ad_Hz(i,j,k) = IniVal
              GRID(ng) % ad_Huon(i,j,k) = IniVal
              GRID(ng) % ad_Hvom(i,j,k) = IniVal
              GRID(ng) % ad_z_r(i,j,k) = IniVal
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              GRID(ng) % ad_z_w(i,j,k) = IniVal
            END DO
          END DO
# endif
        END DO
# if defined ADJUST_BOUNDARY && defined SOLVE3D
        GRID(ng) % ad_Hz_bry = IniVal
# endif
      END IF
#endif

      RETURN
      END SUBROUTINE initialize_grid

      END MODULE mod_grid
