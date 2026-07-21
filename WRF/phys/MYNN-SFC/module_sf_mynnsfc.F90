# include "wrfcpp.h"
!>\file module_sf_mynnsfc.F90
!! This file contains
!WRF:MODEL_LAYER:PHYSICS
!
!>\ingroup mynn_sfc
!> This module contain routines to calculate stability parameters, kinematic siscosity
!! in MYNN surface layer scheme
MODULE module_sf_mynnsfc

!-------------------------------------------------------------------
!Modifications implemented by Joseph Olson NOAA/GSL
!The following overviews the current state of this scheme::
!
!   BOTH LAND AND WATER:
!1) Calculation of stability parameter (z/L) taken from Li et al. (2010 BLM)
!   for first iteration of first time step; afterwards, exact calculation
!   using basically the same iterative technique in the module_sf_sfclayrev.F,
!   which leverages Pedro Jimenez's code, and is adapted for MYNN.
!2) Fixed isflux=0 option to turn off scalar fluxes, but keep momentum
!   fluxes for idealized studies (credit: Anna Fitch).
!3) Kinematic viscosity varies with temperature according to Andreas (1989).
!4) Uses the blended Monin-Obukhov flux-profile relationships COARE (Fairall
!   et al 2003) for the unstable regime (a blended mix of Dyer-Hicks 1974 and
!   Grachev et al (2000). Uses Cheng and Brutsaert (2005) for stable conditions.
!5) The following overviews the namelist variables that control the
!   aerodynamic roughness lengths (over water) and the thermal and moisture
!   roughness lengths (defaults are recommended):
!
!   LAND only:
!   "iz0tlnd" namelist option is used to select the following momentum options:
!   (default) =0: Zilitinkevich (1995); Czil now set to 0.095
!             =1: Czil_new (modified according to Chen & Zhang 2008)
!             =2: Modified Yang et al (2002, 2008) - generalized for all landuse
!             =3: constant zt = z0/7.4 (original form; Garratt 1992)
!             =4: GFS - taken from sfc_diff.f, for comparison/testing
!
!   WATER only:
!   "isftcflx" namelist option is used to select the following scalar options:
!   (default) =0: z0, zt, and zq from the COARE algorithm. Set COARE_OPT (below) to
!                 3.0 (Fairall et al. 2003, default)
!                 3.5 (Edson et al 2013)
!             =1: z0 from Davis et al (2008), zt & zq from COARE 3.0/3.5
!             =2: z0 from Davis et al (2008), zt & zq from Garratt (1992)
!             =3: z0 from Taylor and Yelland (2004), zt and zq from COARE 3.0/3.5
!             =4: GFS - taken from sfc_diff.f, for comparison/testing
!
!   SNOW/ICE only:
!   Andreas (2002) snow/ice parameterization for thermal and
!   moisture roughness is used over all gridpoints with snow deeper than
!   0.1 m. This algorithm calculates a z0 for snow (Andreas et al. 2005, BLM),
!   which is only used as part of the thermal and moisture roughness
!   length calculation, not to directly impact the surface winds.
!
! Misc:
!1) Added a more elaborate diagnostic for u10 & V10 for high vertical resolution
!   model configurations but for most model configurations with depth of
!   the lowest half-model level near 10 m, a neutral-log diagnostic is used.
!
!2) Option to activate stochastic parameter perturbations (SPP), which
!   perturb z0, zt, and zq, along with many other parameters in the MYNN-
!   EDMF scheme.
!
!NOTE: This code was primarily tested in combination with the RUC LSM.
!      Performance with the Noah (or other) LSM is relatively unknown.
!-------------------------------------------------------------------
!Include host model constants
      use physcons, only : cp     => con_cp,     & !=7*Rd/2
     &                     grav   => con_g,      & !=9.81
     &                     Rd     => con_rd,     & !=287.
     &                     Rv     => con_rv,     & !=461.6
!     &                     cpv    => con_cvap,   & !=4*Rv
     &                     rovcp  => con_rocp,   & !=Rd/cp
     &                     xlv    => con_hvap,   & !2.5e6
     &                     xlf    => con_hfus,   & !3.5e5
     &                     ep1    => con_fvirt,  & !Rv/Rd - 1
     &                     ep2    => con_eps       !Rd/Rv

!use kind_phys for real-types
    use machine , only : kind_phys

!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!Drive and/or define more constant:
  real(kind_phys), parameter :: ep3           = 1.-ep2
  real(kind_phys), parameter :: g_inv         = 1.0/grav
  real(kind_phys), parameter :: rvovrd        = Rv/Rd
  real(kind_phys), parameter :: wmin          = 0.1    ! Minimum wind speed
  real(kind_phys), parameter :: karman        = 0.4
  real(kind_phys), parameter :: SVP1          = 0.6112
  real(kind_phys), parameter :: SVP2          = 17.67
  real(kind_phys), parameter :: SVP3          = 29.65
  real(kind_phys), parameter :: SVPT0         = 273.15
  real(kind_phys), parameter :: VCONVC        = 1.25
  real(kind_phys), parameter :: onethird      = 1./3.
  real(kind_phys), parameter :: sqrt3         = 1.7320508075688773
  real(kind_phys), parameter :: atan1         = 0.785398163397     !in radians
  real(kind_phys), parameter :: log01         = log(0.01)
  real(kind_phys), parameter :: log05         = log(0.05)
  real(kind_phys), parameter :: log07         = log(0.07)
  real(kind_phys), parameter :: SNOWZ0        = 0.011
  real(kind_phys), parameter :: COARE_OPT     = 3.0  ! 3.0 or 3.5

  !For debugging purposes:
  INTEGER, PARAMETER :: debug_code = 0  !0: no extra ouput
                                        !1: check input
                                        !2: everything - heavy I/O

  REAL(kind_phys), DIMENSION(0:1000 ),SAVE :: psim_stab,psim_unstab, &
                                     psih_stab,psih_unstab
!$acc declare create(psim_stab, psim_unstab, psih_stab, psih_unstab)

CONTAINS

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine
   SUBROUTINE SFCLAY_mynn(                           &
              U3D,V3D,T3D,QV3D,P3D,dz8w,             & !in
              th3d,pi3d,qc3d,                        & !in
              PSFCPA,PBLH,MAVAIL,XLAND,DX,           & !in
              ISFFLX,isftcflx,lsm,lsm_ruc,           & !in
              compute_flux,compute_diag,             & !in
              iz0tlnd,psi_opt,                       & !in
              sigmaf,vegtype,shdmax,ivegsrc,         & !intent(in)
              z0pert,ztpert,                         & !intent(in)
              redrag,sfc_z0_type,                    & !intent(in)
              itimestep,iter,flag_iter,              & !in
              flag_restart,                          & !in
                    wet,       dry,       icy,       & !intent(in)
              tskin_wat, tskin_lnd, tskin_ice,       & !intent(in)
              tsurf_wat, tsurf_lnd, tsurf_ice,       & !intent(in)
               qsfc_wat,  qsfc_lnd,  qsfc_ice,       & !intent(in)
              snowh_wat, snowh_lnd, snowh_ice,       & !intent(in)
                ZNT_wat,   ZNT_lnd,   ZNT_ice,       & !intent(inout)
                UST_wat,   UST_lnd,   UST_ice,       & !intent(inout)
                 cm_wat,    cm_lnd,    cm_ice,       & !intent(inout)
                 ch_wat,    ch_lnd,    ch_ice,       & !intent(inout)
                 rb_wat,    rb_lnd,    rb_ice,       & !intent(inout)
             stress_wat,stress_lnd,stress_ice,       & !intent(inout)
                 fm_wat,    fm_lnd,    fm_ice,       & !intent(inout)
                 fh_wat,    fh_lnd,    fh_ice,       & !intent(inout)
               fm10_wat,  fm10_lnd,  fm10_ice,       & !intent(inout)
                fh2_wat,   fh2_lnd,   fh2_ice,       & !intent(inout)
               HFLX_wat,  HFLX_lnd,  HFLX_ice,       &
               QFLX_wat,  QFLX_lnd,  QFLX_ice,       &
              CH,CHS,CHS2,CQS2,CPM,                  &
              ZNT,USTM,ZOL,MOL,RMOL,                 &
              PSIM,PSIH,                             &
              HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,        &
              QGH,QSFC,                              &
              U10,V10,TH2,T2,Q2,                     &
              GZ1OZ0,WSPD,WSTAR,                     &
              spp_sfc,pattern_spp_sfc,               &
#if defined SWAN_COUPLING || defined WW3_COUPLING
              HWAVE, LWAVEP, DWAVEP,                 &
              PWAVE, Z0_WAV, COSA, SINA,             &
#endif
              ids,ide, jds,jde, kds,kde,             &
              ims,ime, jms,jme, kms,kme,             &
              its,ite, jts,jte, kts,kte,             &
              errmsg, errflg                         )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!-- U3D         3D u-velocity interpolated to theta points (m/s)
!-- V3D         3D v-velocity interpolated to theta points (m/s)
!-- T3D         3D temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- dz8w        3D dz between full levels (m)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- grav        acceleration due to gravity (m/s^2)
!-- ROVCP       R/CP
!-- Rd          gas constant for dry air (J/kg/K)
!-- XLV         latent heat of vaporization for water (J/kg)
!-- PSFCPA      surface pressure (Pa)
!-- ZNT         roughness length (m)
!-- UST         u* in similarity theory (m/s)
!-- USTM        u* in similarity theory (m/s) w* added to WSPD. This is
!               used to couple with TKE scheme but not in MYNN.
!               (as of now, USTM = UST in this version)
!-- PBLH        PBL height from previous time (m)
!-- MAVAIL      surface moisture availability (between 0 and 1)
!-- ZOL         z/L height over Monin-Obukhov length
!-- MOL         T* (similarity theory) (K)
!-- RMOL        Reciprocal of M-O length (/m)
!-- REGIME      flag indicating PBL regime (stable, unstable, etc.)
!-- PSIM        similarity stability function for momentum
!-- PSIH        similarity stability function for heat
!-- XLAND       land mask (1 for land, 2 for water)
!-- HFX         upward heat flux at the surface (W/m^2)
!                  HFX = HFLX * rho * cp
!-- HFLX        upward temperature flux at the surface (K m s^-1)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!                  QFX = QFLX * rho
!-- QFLX        upward moisture flux at the surface (kg kg-1 m s-1)
!-- LH          net upward latent heat flux at surface (W/m^2)
!-- TSK         surface temperature (K)
!-- FLHC        exchange coefficient for heat (W/m^2/K)
!-- FLQC        exchange coefficient for moisture (kg/m^2/s)
!-- CHS         heat/moisture exchange coefficient for LSM (m/s)
!-- QGH         lowest-level saturated mixing ratio
!-- QSFC        qv (specific humidity) at the surface
!-- QSFCMR      qv (mixing ratio) at the surface
!-- U10         diagnostic 10m u wind
!-- V10         diagnostic 10m v wind
!-- TH2         diagnostic 2m theta (K)
!-- T2          diagnostic 2m temperature (K)
!-- Q2          diagnostic 2m mixing ratio (kg/kg)
!-- SNOWH       Snow height (m)
!-- GZ1OZ0      log((z1+ZNT)/ZNT) where ZNT is roughness length
!-- WSPD        wind speed at lowest model level (m/s)
!-- BR          bulk Richardson number in surface layer
!-- ISFFLX      isfflx=1 for surface heat and moisture fluxes
!-- DX          horizontal grid size (m)
!-- SVP1        constant for saturation vapor pressure (=0.6112 kPa)
!-- SVP2        constant for saturation vapor pressure (=17.67 dimensionless)
!-- SVP3        constant for saturation vapor pressure (=29.65 K)
!-- SVPT0       constant for saturation vapor pressure (=273.15 K)
!-- EP1         constant for virtual temperature (Rv/Rd - 1) (dimensionless)
!-- EP2         constant for spec. hum. calc (Rd/Rv = 0.622) (dimensionless)
!-- EP3         constant for spec. hum. calc (1 - Rd/Rv = 0.378 ) (dimensionless)
!-- KARMAN      Von Karman constant
!-- ck          enthalpy exchange coeff at 10 meters
!-- cd          momentum exchange coeff at 10 meters
!-- cka         enthalpy exchange coeff at the lowest model level
!-- cda         momentum exchange coeff at the lowest model level
!-- isftcflx    =0: z0, zt, and zq from COARE3.0/3.5 (Fairall et al 2003/Edson et al 2013)
!   (water      =1: z0 from Davis et al (2008), zt & zq from COARE3.0/3.5
!    only)      =2: z0 from Davis et al (2008), zt & zq from Garratt (1992)
!               =3: z0 from Taylor and Yelland (2004), zt and zq from COARE 3.0/3.5
!-- iz0tlnd     =0: Zilitinkevich (1995) with Czil=0.095,
!   (land       =1: Czil_new (modified according to Chen & Zhang 2008)
!    only)      =2: Modified Yang et al (2002, 2008) - generalized for all landuse
!               =3: constant zt = z0/7.4 (Garratt 1992)
!
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-- errmsg      CCPP error message
!-- errflg      CCPP error code
!=================================================================
! SCALARS
!===================================
      INTEGER,  INTENT(IN)   ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte
      INTEGER,  INTENT(IN)   ::        itimestep,iter
!NAMELIST/CONFIGURATION OPTIONS:
      integer, intent(in)           :: ISFFLX, LSM, LSM_RUC
      INTEGER, OPTIONAL, INTENT(IN) :: ISFTCFLX, IZ0TLND
      INTEGER, OPTIONAL, INTENT(IN) :: spp_sfc, psi_opt
      logical, intent(in) :: compute_flux,compute_diag
      integer, intent(in) :: ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, intent(in) :: flag_restart

!Input data
      integer, dimension(ims:ime), intent(in) :: vegtype
      real(kind_phys), dimension(ims:ime), intent(in) ::           &
     &                    sigmaf,shdmax,z0pert,ztpert
!===================================
! 3D VARIABLES
!===================================
      REAL(kind_phys), DIMENSION( ims:ime, kms:kme )             , &
                            INTENT(IN   )   ::               dz8w, &
                                                             QV3D, &
                                                              P3D, &
                                                              T3D, &
                                                             QC3D, &
                                                          U3D,V3D, &
                                                        th3d,pi3d

      !GJF: This array must be assumed-shape since it is conditionally-allocated
      REAL(kind_phys), DIMENSION( :,: ), OPTIONAL,                 &
                            INTENT(IN) ::         pattern_spp_sfc
!===================================
! 2D VARIABLES
!===================================
      REAL(kind_phys), DIMENSION( ims:ime )                      , &
                            INTENT(IN   )          ::      MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                           PSFCPA, &
                                                               DX

      REAL(kind_phys), DIMENSION( ims:ime )                      , &
                            INTENT(OUT  )          ::     U10,V10, &
                                                        TH2,T2,Q2


      REAL(kind_phys), DIMENSION( ims:ime )                      , &
                            INTENT(INOUT)          ::    HFLX,HFX, &
                                                         QFLX,QFX, &
                                                             RMOL, &
                                                             QSFC, &
                                                              QGH, &
                                                              ZNT, &
                                                              CPM, &
                                                              CHS, &
                                                               CH, &
                                                        FLHC,FLQC, &
                                                      GZ1OZ0,WSPD, &
                                                        PSIM,PSIH, &
                                                        USTM,CHS2, &
                                                        CQS2, WSTAR
      REAL(kind_phys), DIMENSION( ims:ime ),                       &
                            INTENT(INOUT)          ::          LH, &
                                                              ZOL, &
                                                              MOL

      LOGICAL, DIMENSION( ims:ime ), INTENT(IN)    ::              &
&                             wet,  dry,  icy,  flag_iter

      REAL(kind_phys), DIMENSION( ims:ime ), INTENT(IN) ::         &
     &                    tskin_wat, tskin_lnd, tskin_ice,         &
     &                    tsurf_wat, tsurf_lnd, tsurf_ice,         &
     &                    snowh_wat, snowh_lnd, snowh_ice

      REAL(kind_phys), DIMENSION( ims:ime), INTENT(INOUT) ::       &
     &                      ZNT_wat,   ZNT_lnd,   ZNT_ice,         &
     &                      UST_wat,   UST_lnd,   UST_ice,         &
     &                       cm_wat,    cm_lnd,    cm_ice,         &
     &                       ch_wat,    ch_lnd,    ch_ice,         &
     &                       rb_wat,    rb_lnd,    rb_ice,         &
     &                   stress_wat,stress_lnd,stress_ice,         &
     &                       fm_wat,    fm_lnd,    fm_ice,         &
     &                       fh_wat,    fh_lnd,    fh_ice,         &
     &                     fm10_wat,  fm10_lnd,  fm10_ice,         &
     &                      fh2_wat,   fh2_lnd,   fh2_ice,         &
     &                     HFLX_wat,  HFLX_lnd,  HFLX_ice,         &
     &                     QFLX_wat,  QFLX_lnd,  QFLX_ice,         &
     &                     qsfc_wat,  qsfc_lnd,  qsfc_ice

! CCPP error handling
      character(len=*), intent(inout) :: errmsg
      integer,          intent(inout) :: errflg

!ADDITIONAL OUTPUT
!JOE-begin
      REAL(kind_phys), DIMENSION( ims:ime ) :: qstar
!JOE-end
#if defined SWAN_COUPLING || defined WW3_COUPLING
       REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: HWAVE,   &
     &                 LWAVEP, DWAVEP, PWAVE, Z0_WAV, COSA, SINA
#endif
!===================================
! 1D LOCAL ARRAYS
!===================================
      REAL(kind_phys), DIMENSION( its:ite ) ::            U1D,V1D, & !level1 winds
                                                        U1D2,V1D2, & !level2 winds
                                                             QV1D, &
                                                              P1D, &
                                                         T1D,QC1D, &
                                                           dz8w1d, & !level 1 height
                                                           dz2w1d    !level 2 height

      REAL(kind_phys), DIMENSION( its:ite ) ::      rstoch1D

      INTEGER ::  I,J,K,itf,ktf
!-----------------------------------------------------------

      ! Initialize error-handling
      errflg = 0
      errmsg = ''

!$acc enter data copyin( dz8w,U3D,V3D,QV3D,QC3D,P3D,T3D,       &
!$acc                       pattern_spp_sfc)

!$acc enter data copyin( UST_WAT(:), UST_LND(:),  UST_ICE(:),     &
!$acc                    MOL(:),     QFLX(:),     HFLX(:),        &
!$acc                    QSFC(:),    QSFC_WAT(:), QSFC_LND(:),    &
!$acc                    QSFC_ICE(:))

!$acc enter data create( dz8w1d(:),  dz2w1d(:),   U1D(:),         &
!$acc                    V1D(:),     U1D2(:),     V1D2(:),        &
!$acc                    QV1D(:),    QC1D(:),     P1D(:),         &
!$acc                    T1D(:),     rstoch1D(:), qstar(:))


      IF (debug_code >= 1) THEN
        write(*,*)"======= printing of constants:"
        write(*,*)"cp=",    cp," g=",     grav
        write(*,*)"Rd=",    Rd," ep1=",    ep1
        write(*,*)"xlv=",  XLV," xlf=",   XLF
        write(*,*)"ep2=", ep2
      ENDIF

      itf=ite !MIN0(ite,ide-1)
      ktf=kte !MIN0(kte,kde-1)

!$acc parallel loop present(dz8w,U3D,V3D,QV3D,QC3D,P3D,T3D,       &
!$acc                       pattern_spp_sfc,dz8w1d,dz2w1d,U1D,    &
!$acc                       V1D,U1D2,V1D2,QV1D,QC1D,P1D,T1D,      &
!$acc                       rstoch1D,qstar)
      DO i=its,ite
         dz8w1d(I) = dz8w(i,kts)
         dz2w1d(I) = dz8w(i,kts+1)
         U1D(i) =U3D(i,kts)
         V1D(i) =V3D(i,kts)
         !2nd model level winds - for diags with high-res grids
         U1D2(i) =U3D(i,kts+1)
         V1D2(i) =V3D(i,kts+1)
         QV1D(i)=QV3D(i,kts)
         QC1D(i)=QC3D(i,kts)
         P1D(i) =P3D(i,kts)
         T1D(i) =T3D(i,kts)
         if (spp_sfc==1) then
            rstoch1D(i)=pattern_spp_sfc(i,kts)
         else
            rstoch1D(i)=0.0
         endif
         qstar(i)=0.0
      ENDDO

      IF (itimestep==1 .AND. iter==1) THEN
!$acc parallel loop present(U1D,V1D,UST_WAT,UST_LND,UST_ICE,MOL,  &
!$acc                       QFLX,HFLX,QV3D,QSFC,QSFC_WAT,         &
!$acc                       QSFC_LND,QSFC_ICE)
         DO i=its,ite
            IF (.not. flag_restart) THEN
               !Everything here is used before calculated
               if (ust_wat(i) .lt. 1e-4 .or. ust_wat(i) .gt. 3.0) then
                  UST_WAT(i)=MAX(0.04*SQRT(U1D(i)*U1D(i) + V1D(i)*V1D(i)),0.001_kind_phys)
               endif
               if (ust_lnd(i) .lt. 1e-4 .or. ust_lnd(i) .gt. 3.0) then
                  UST_LND(i)=MAX(0.04*SQRT(U1D(i)*U1D(i) + V1D(i)*V1D(i)),0.001_kind_phys)
               endif
               if (ust_ice(i) .lt. 1e-4 .or. ust_ice(i) .gt. 3.0) then
                  UST_ICE(i)=MAX(0.04*SQRT(U1D(i)*U1D(i) + V1D(i)*V1D(i)),0.001_kind_phys)
               endif
               MOL(i)=0.0
            ENDIF ! restart
            QFLX(i)=0.
            HFLX(i)=0.
            if ( LSM == LSM_RUC ) then
            !- qsfc_lnd and qsfc_ice are already available
               QSFC(i)=QV3D(i,kts)/(1.+QV3D(i,kts))
               QSFC_WAT(i)=QSFC(i)
            else
               QSFC(i)=QV3D(i,kts)/(1.+QV3D(i,kts))
               QSFC_WAT(i)=QSFC(i)
               QSFC_LND(i)=QSFC(i)
               QSFC_ICE(i)=QSFC(i)
            endif ! lsm==lsm_ruc
         ENDDO
      ENDIF

!$acc exit data delete( dz8w,U3D,V3D,QV3D,QC3D,P3D,T3D,         &
!$acc                   pattern_spp_sfc, QC1D)

      CALL SFCLAY1D_mynn(flag_iter,                             &
           J,U1D,V1D,T1D,QV1D,P1D,dz8w1d,                       &
           U1D2,V1D2,dz2w1d,                                    &
           PSFCPA,PBLH,MAVAIL,XLAND,DX,                         &
           ISFFLX,isftcflx,iz0tlnd,psi_opt,                     &
           compute_flux,compute_diag,                           &
           sigmaf,vegtype,shdmax,ivegsrc,                       &  !intent(in)
           z0pert,ztpert,                                       &  !intent(in)
           redrag,sfc_z0_type,                                  &  !intent(in)
           itimestep,iter,flag_restart,lsm,lsm_ruc,             &
                  wet,          dry,          icy,              &  !intent(in)
            tskin_wat,    tskin_lnd,    tskin_ice,              &  !intent(in)
            tsurf_wat,    tsurf_lnd,    tsurf_ice,              &  !intent(in)
             qsfc_wat,     qsfc_lnd,     qsfc_ice,              &  !intent(in)
            snowh_wat,    snowh_lnd,    snowh_ice,              &  !intent(in)
              ZNT_wat,      ZNT_lnd,      ZNT_ice,              &  !intent(inout)
              UST_wat,      UST_lnd,      UST_ice,              &  !intent(inout)
               cm_wat,       cm_lnd,       cm_ice,              &  !intent(inout)
               ch_wat,       ch_lnd,       ch_ice,              &  !intent(inout)
               rb_wat,       rb_lnd,       rb_ice,              &  !intent(inout)
           stress_wat,   stress_lnd,   stress_ice,              &  !intent(inout)
               fm_wat,       fm_lnd,       fm_ice,              &  !intent(inout)
               fh_wat,       fh_lnd,       fh_ice,              &  !intent(inout)
             fm10_wat,     fm10_lnd,     fm10_ice,              &  !intent(inout)
              fh2_wat,      fh2_lnd,      fh2_ice,              &
             HFLX_wat,     HFLX_lnd,     HFLX_ice,              &
             QFLX_wat,     QFLX_lnd,     QFLX_ice,              &
           ch,CHS,CHS2,CQS2,CPM,                                &
           ZNT,USTM,ZOL,MOL,RMOL,                               &
           PSIM,PSIH,                                           &
           HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                      &
           QGH,QSFC,U10,V10,TH2,T2,Q2,                          &
           GZ1OZ0,WSPD,wstar,qstar,                             &
           spp_sfc,rstoch1D,                                    &
#if defined SWAN_COUPLING || defined WW3_COUPLING
!          HWAVE(ims,j), LWAVEP(ims,j), DWAVEP(ims,j),          &
!          PWAVE(ims,j), Z0_WAV(ims,j),                         &
           HWAVE, LWAVEP, DWAVEP,                               &
           PWAVE, Z0_WAV, COSA, SINA,                           &
#endif
           ids,ide, jds,jde, kds,kde,                           &
           ims,ime, jms,jme, kms,kme,                           &
           its,ite, jts,jte, kts,kte,                           &
           errmsg, errflg                                       )

!$acc exit data copyout( UST_WAT(:), UST_LND(:),  UST_ICE(:),   &
!$acc                    MOL(:),     QFLX(:),     HFLX(:),      &
!$acc                    QSFC(:),    QSFC_WAT(:), QSFC_LND(:),  &
!$acc                    QSFC_ICE(:))

!$acc exit data delete( dz8w1d(:),   dz2w1d(:),   U1D(:),       &
!$acc                   V1D(:),      U1D2(:),     V1D2(:),      &
!$acc                   QV1D(:),     T1D(:),      P1D(:),       &
!$acc                   rstoch1D(:), qstar(:))

    END SUBROUTINE SFCLAY_MYNN

!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine calculates u*, z/L, and the exchange coefficients
!! which are passed to subsequent scheme to calculate the fluxes.
!! This scheme has options to calculate the fluxes and near-surface
!! diagnostics, as was needed in WRF, but these are skipped for FV3.
   SUBROUTINE SFCLAY1D_mynn(flag_iter,                            &
             J,U1D,V1D,T1D,QV1D,P1D,dz8w1d,U1D2,V1D2,dz2w1d,      &
             PSFCPA,PBLH,MAVAIL,XLAND,DX,                         &
             ISFFLX,isftcflx,iz0tlnd,psi_opt,                     &
             compute_flux,compute_diag,                           &
             sigmaf,vegtype,shdmax,ivegsrc,                       &  !intent(in)
             z0pert,ztpert,                                       &  !intent(in)
             redrag,sfc_z0_type,                                  &  !intent(in)
             itimestep,iter,flag_restart,lsm,lsm_ruc,             &
                    wet,          dry,          icy,              &  !intent(in)
              tskin_wat,    tskin_lnd,    tskin_ice,              &  !intent(in)
              tsurf_wat,    tsurf_lnd,    tsurf_ice,              &  !intent(in)
               qsfc_wat,     qsfc_lnd,     qsfc_ice,              &  !intent(in)
              snowh_wat,    snowh_lnd,    snowh_ice,              &  !intent(in)
                ZNT_wat,      ZNT_lnd,      ZNT_ice,              &  !intent(inout)
                UST_wat,      UST_lnd,      UST_ice,              &  !intent(inout)
                 cm_wat,       cm_lnd,       cm_ice,              &  !intent(inout)
                 ch_wat,       ch_lnd,       ch_ice,              &  !intent(inout)
                 rb_wat,       rb_lnd,       rb_ice,              &  !intent(inout)
             stress_wat,   stress_lnd,   stress_ice,              &  !intent(inout)
               psix_wat,     psix_lnd,     psix_ice,              &  !=fm, intent(inout)
               psit_wat,     psit_lnd,     psit_ice,              &  !=fh, intent(inout)
             psix10_wat,   psix10_lnd,   psix10_ice,              &  !=fm10, intent(inout)
              psit2_wat,    psit2_lnd,    psit2_ice,              &  !=fh2, intent(inout)
               HFLX_wat,     HFLX_lnd,     HFLX_ice,              &
               QFLX_wat,     QFLX_lnd,     QFLX_ice,              &
             ch,CHS,CHS2,CQS2,CPM,                                &
             ZNT,USTM,ZOL,MOL,RMOL,                               &
             PSIM,PSIH,                                           &
             HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                      &
             QGH,QSFC,                                            &
             U10,V10,TH2,T2,Q2,                                   &
             GZ1OZ0,WSPD,wstar,qstar,                             &
             spp_sfc,rstoch1D,                                    &
#if defined SWAN_COUPLING || defined WW3_COUPLING
             HWAVE, LWAVEP, DWAVEP, PWAVE, Z0_WAV, COSA, SINA,    &
#endif
             ids,ide, jds,jde, kds,kde,                           &
             ims,ime, jms,jme, kms,kme,                           &
             its,ite, jts,jte, kts,kte,                           &
             errmsg, errflg                                       )

!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
! SCALARS
!-----------------------------
      INTEGER,  INTENT(IN)        :: ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte, &
                                     J, itimestep, iter, lsm, lsm_ruc
      LOGICAL, INTENT(IN)         :: flag_restart

      REAL(kind_phys), PARAMETER  :: XKA=2.4E-5   !molecular diffusivity
      REAL(kind_phys), PARAMETER  :: PRT=1.       !prandlt number
      REAL(kind_phys), PARAMETER  :: snowh_thresh = 50. !mm

!-----------------------------
! NAMELIST OPTIONS
!-----------------------------
      integer, intent(in) :: ISFFLX
      integer, optional,  intent(in)  :: ISFTCFLX, IZ0TLND
      logical, intent(in) :: compute_flux,compute_diag
      integer, intent(in) :: spp_sfc, psi_opt
      integer, intent(in) :: ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)

!Input data
      integer, dimension(ims:ime), intent(in) :: vegtype
      real(kind_phys), dimension(ims:ime), intent(in) ::           &
      &                    sigmaf,shdmax,z0pert,ztpert

!-----------------------------
! 1D ARRAYS
!-----------------------------
      REAL(kind_phys), DIMENSION( ims:ime ),                       &
                            INTENT(IN)             ::      MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                           PSFCPA, &
                                                               DX

      REAL(kind_phys), DIMENSION( its:ite ),                       &
                            INTENT(IN)             ::     U1D,V1D, &
                                                        U1D2,V1D2, &
                                                         QV1D,P1D, &
                                                              T1D, &
                                                           dz8w1d, &
                                                           dz2w1d

      REAL(kind_phys), DIMENSION( ims:ime ),                       &
                            INTENT(OUT)             ::    QFX,HFX, &
                                                             RMOL
      REAL(kind_phys), DIMENSION( ims:ime ),                       &
                            INTENT(INOUT)           ::  HFLX,QFLX, &
                                                         QGH,QSFC, &
                                                              ZNT, &
                                                              CPM, &
                                                           CHS,CH, &
                                                        FLHC,FLQC, &
                                                           GZ1OZ0, &
                                                             WSPD, &
                                                             PSIM, &
                                                             PSIH, &
                                                             USTM, &
                                                        CHS2,CQS2
      REAL(kind_phys), DIMENSION( ims:ime ),                       &
                            INTENT(INOUT)           ::        MOL, &
                                                              ZOL, &
                                                               LH

      LOGICAL, DIMENSION( ims:ime ), INTENT(IN)    ::              &
     &                wet,     dry,     icy,    flag_iter

      REAL(kind_phys), DIMENSION( ims:ime ), INTENT(in) ::         &
     &                    tskin_wat, tskin_lnd, tskin_ice,         &
     &                    tsurf_wat, tsurf_lnd, tsurf_ice,         &
     &                    snowh_wat, snowh_lnd, snowh_ice

      REAL(kind_phys), DIMENSION( ims:ime ), INTENT(inout) ::      &
     &                      ZNT_wat,   ZNT_lnd,   ZNT_ice,         &
     &                      UST_wat,   UST_lnd,   UST_ice,         &
     &                       cm_wat,    cm_lnd,    cm_ice,         &
     &                       ch_wat,    ch_lnd,    ch_ice,         &
     &                       rb_wat,    rb_lnd,    rb_ice,         &
     &                   stress_wat,stress_lnd,stress_ice,         &
     &                     psix_wat,  psix_lnd,  psix_ice,         &
     &                     psit_wat,  psit_lnd,  psit_ice,         &
     &                   psix10_wat,psix10_lnd,psix10_ice,         &
     &                    psit2_wat, psit2_lnd, psit2_ice,         &
     &                     HFLX_wat,  HFLX_lnd,  HFLX_ice,         &
     &                     QFLX_wat,  QFLX_lnd,  QFLX_ice,         &
     &                     qsfc_wat,  qsfc_lnd,  qsfc_ice

      REAL(kind_phys), DIMENSION( its:ite ),                       &
     &                      INTENT(IN)             ::    rstoch1D

      ! DIAGNOSTIC OUTPUT
      REAL(kind_phys), DIMENSION( ims:ime ),                       &
     &                      INTENT(OUT)            ::    U10, V10, &
     &                                                    TH2, T2, &
     &                                                         Q2

!--------------------------------------------
!JOE-additinal output
      REAL(kind_phys), DIMENSION( ims:ime ),                       &
     &                      INTENT(OUT)            ::       qstar, &
                                                            wstar
!JOE-end

#if defined SWAN_COUPLING || defined WW3_COUPLING
       REAL, DIMENSION( ims:ime ), INTENT(IN) :: HWAVE, LWAVEP,    &
     &                       DWAVEP, PWAVE, Z0_WAV, COSA, SINA
#endif

! CCPP error handling
      character(len=*), intent(inout) :: errmsg
      integer,          intent(inout) :: errflg

! Local fixed-size errmsg character array for error messages on accelerator
! devices distinct from the host (e.g. GPUs). Necessary since OpenACC does
! not support assumed-size (len=*) arrays like errmsg. Additional
! device_errflg integer to denote when device_errmsg needs to be synced
! with errmsg.
      character(len=512) :: device_errmsg
      integer            :: device_errflg

! Special versions of the fixed-size errmsg character array for error messages
! on the device and it's errflag counterpart. These are necessary to ensure
! the return statements at lines 1417 and 2030 are executed only for this
! special case, and not any and all error messages set on the device.
      character(len=512) :: device_special_errmsg
      integer            :: device_special_errflg


!----------------------------------------------------------------
! LOCAL VARS
!----------------------------------------------------------------
      REAL(kind_phys), DIMENSION(its:ite) ::                      &
                 ZA, &    !Height of lowest 1/2 sigma level(m)
                ZA2, &    !Height of 2nd lowest 1/2 sigma level(m)
              THV1D, &    !Theta-v at lowest 1/2 sigma (K)
               TH1D, &    !Theta at lowest 1/2 sigma (K)
               TC1D, &    !T at lowest 1/2 sigma (Celsius)
               TV1D, &    !Tv at lowest 1/2 sigma (K)
              RHO1D, &    !density at lowest 1/2 sigma level
               QVSH, &    !qv at lowest 1/2 sigma (spec humidity)
              PSIH2, &    !M-O stability functions at z=2 m
             PSIM10, &    !M-O stability functions at z=10 m
             PSIH10, &    !M-O stability functions at z=10 m
              WSPDI, &
             GOVRTH, &    !grav/theta
               PSFC, &    !press at surface (Pa/1000)
             QSFCMR, &    !qv at surface (mixing ratio, kg/kg)
              THCON, &    !conversion from temp to theta
          zratio_lnd,   zratio_ice,   zratio_wat, & !z0/zt
             TSK_lnd,      TSK_ice,      TSK_wat, & !absolute temperature
            THSK_lnd,     THSK_ice,     THSK_wat, & !theta
           THVSK_lnd,    THVSK_ice,    THVSK_wat, & !theta-v
          GZ1OZ0_lnd,   GZ1OZ0_ice,   GZ1OZ0_wat, & !LOG((ZA(I)+ZNT(i))/ZNT(i))
          GZ1OZt_lnd,   GZ1OZt_ice,   GZ1OZt_wat, & !LOG((ZA(I)+ZT(i))/ZT(i))
          GZ2OZ0_lnd,   GZ2OZ0_ice,   GZ2OZ0_wat, & !LOG((2.0+ZNT(I))/ZNT(I))
          GZ2OZt_lnd,   GZ2OZt_ice,   GZ2OZt_wat, & !LOG((2.0+ZT(I))/ZT(I))
         GZ10OZ0_lnd,  GZ10OZ0_ice,  GZ10OZ0_wat, & !LOG((10.+ZNT(I))/ZNT(I))
         GZ10OZt_lnd,  GZ10OZt_ice,  GZ10OZt_wat, & !LOG((10.+ZT(I))/ZT(I))
        ZNTstoch_lnd, ZNTstoch_ice, ZNTstoch_wat, &
              ZT_lnd,       ZT_ice,       ZT_wat, &
              ZQ_lnd,       ZQ_ice,       ZQ_wat, &
            PSIQ_lnd,     PSIQ_ice,     PSIQ_wat, &
           PSIQ2_lnd,    PSIQ2_ice,    PSIQ2_wat, &
          QSFCMR_lnd,   QSFCMR_ice,   QSFCMR_wat

      INTEGER ::  N,I,K,L,yesno

      REAL(kind_phys) :: PL,E1,TABS
      REAL(kind_phys) :: WSPD_lnd, WSPD_ice, WSPD_wat
      REAL(kind_phys) :: DTHVDZ,DTHVM,VCONV,ZOL2,ZOL10,ZOLZA,ZOLZ0,ZOLZT
      REAL(kind_phys) :: DTG,DTTHX,PSIQ,PSIQ2,PSIQ10,PSIT10
      REAL(kind_phys) :: FLUXC,VSGD
      REAL(kind_phys) :: restar,VISC,DQG,OLDUST,OLDTST
#if defined SWAN_COUPLING || defined WW3_COUPLING
      REAL(kind_phys) :: CWAVE,UREAL,VREAL,DWIND,R2D,DWAVE1,THWV,Z0WAVE
#endif
      ! Initialize error-handling
      errflg = 0
      errmsg = ''
      device_errflg = errflg
      device_errmsg = errmsg
      device_special_errflg = errflg
      device_special_errmsg = errmsg
!-------------------------------------------------------------------
!$acc update device(psim_stab, psim_unstab, psih_stab, psih_unstab)

!$acc enter data create( ZA,     ZA2,    THV1D,  TH1D,   TC1D,   TV1D,  &
!$acc                    RHO1D,  QVSH,   PSIH2,  PSIM10, PSIH10, WSPDI, &
!$acc                    GOVRTH, PSFC,   THCON,                         &
!$acc                    zratio_lnd,   zratio_ice,   zratio_wat,        &
!$acc                    TSK_lnd,      TSK_ice,      TSK_wat,           &
!$acc                    THSK_lnd,     THSK_ice,     THSK_wat,          &
!$acc                    THVSK_lnd,    THVSK_ice,    THVSK_wat,         &
!$acc                    GZ1OZ0_lnd,   GZ1OZ0_ice,   GZ1OZ0_wat,        &
!$acc                    GZ1OZt_lnd,   GZ1OZt_ice,   GZ1OZt_wat,        &
!$acc                    GZ2OZ0_lnd,   GZ2OZ0_ice,   GZ2OZ0_wat,        &
!$acc                    GZ2OZt_lnd,   GZ2OZt_ice,   GZ2OZt_wat,        &
!$acc                    GZ10OZ0_lnd,  GZ10OZ0_ice,  GZ10OZ0_wat,       &
!$acc                    GZ10OZt_lnd,  GZ10OZt_ice,  GZ10OZt_wat,       &
!$acc                    ZNTstoch_lnd, ZNTstoch_ice, ZNTstoch_wat,      &
!$acc                    ZT_lnd,       ZT_ice,       ZT_wat,            &
!$acc                    ZQ_lnd,       ZQ_ice,       ZQ_wat,            &
!$acc                    PSIQ_lnd,     PSIQ_ice,     PSIQ_wat,          &
!$acc                    PSIQ2_lnd,    PSIQ2_ice,    PSIQ2_wat,         &
!$acc                    QSFCMR_lnd,   QSFCMR_ice,   QSFCMR_wat )

!$acc enter data copyin(flag_iter, dry, wet, icy, CPM, MAVAIL, &
!$acc                   QFX, FLHC, FLQC, CHS, CH, CHS2, CQS2, USTM, &
!$acc                   HFX, LH, wstar, qstar, PBLH, ZOL, MOL, RMOL, &
!$acc                   T2, TH2, Q2, QV1D, PSFCPA, &
!$acc                   WSPD, U10, V10, U1D, V1D, U1D2, V1D2, &
!$acc                   T1D, P1D, rstoch1D, sigmaf, &
!$acc                   shdmax, vegtype, z0pert, ztpert, dx, QGH, &
!$acc                   dz2w1d, dz8w1d, &
!$acc                   stress_wat, stress_lnd, stress_ice, &
!$acc                   rb_wat, rb_lnd, rb_ice, &
!$acc                   tskin_wat, tskin_lnd, tskin_ice, &
!$acc                   tsurf_wat, tsurf_lnd, tsurf_ice, &
!$acc                   psim, psih, &
!$acc                   UST_wat, UST_lnd, UST_ice, &
!$acc                   ZNT_wat, ZNT_lnd, ZNT_ice, &
!$acc                   QSFC, QSFC_lnd, QSFC_wat, QSFC_ice, &
!$acc                   QFLX, QFLX_lnd, QFLX_wat, QFLX_ice, &
!$acc                   HFLX, HFLX_lnd, HFLX_wat, HFLX_ice, &
!$acc                   PSIX_wat, PSIX_lnd, PSIX_ice, &
!$acc                   PSIX10_wat, PSIX10_lnd, PSIX10_ice, &
!$acc                   PSIT2_lnd, PSIT2_wat, PSIT2_ice, &
!$acc                   PSIT_lnd, PSIT_wat, PSIT_ice, &
!$acc                   ch_lnd, ch_wat, ch_ice, &
!$acc                   cm_lnd, cm_wat, cm_ice, &
!$acc                   snowh_lnd, snowh_wat, snowh_ice, &
!$acc                   device_errmsg, device_errflg, &
!$acc                   device_special_errmsg, device_special_errflg)

!$acc parallel loop present(PSFCPA, PSFC, QSFC, T1D, flag_iter, tsurf_lnd,    &
!$acc                       QSFC_wat, QSFCMR_wat, wet, TSK_wat, tskin_wat,    &
!$acc                       QSFC_lnd, QSFCMR_lnd, dry, TSK_lnd, tskin_lnd,    &
!$acc                       QSFC_ice, QSFCMR_ice, icy, TSK_ice, tskin_ice)
      DO I=its,ite

         ! PSFC ( in cmb) is used later in saturation checks
         PSFC(I)=PSFCPA(I)/1000.
         !tgs - do computations if flag_iter(i) = .true.
         if ( flag_iter(i) ) then

         IF (ITIMESTEP == 1) THEN
         !initialize surface specific humidity and mixing ratios for land, ice and water
            IF (wet(i)) THEN
               TSK_wat(I) = tskin_wat(i)
               IF (TSK_wat(I) .LT. 273.15) THEN
                  !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
                  E1=SVP1*EXP(4648*(1./273.15 - 1./TSK_wat(I)) - &
                    & 11.64*LOG(273.15/TSK_wat(I)) + 0.02265*(273.15 - TSK_wat(I)))
               ELSE
                  !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
                  E1=SVP1*EXP(SVP2*(TSK_wat(I)-SVPT0)/(TSK_wat(i)-SVP3))
               ENDIF
               QSFC_wat(I)=EP2*E1/(PSFC(I)-ep3*E1)             !specific humidity
               QSFCMR_wat(I)=EP2*E1/(PSFC(I)-E1)                !mixing ratio
               IF(QSFC_wat(I)>1..or.QSFC_wat(I)<0.) print *,' QSFC_wat(I)',itimestep,i,QSFC_wat(I),TSK_wat(i)
            ENDIF
            IF (dry(i)) THEN
              TSK_lnd(I) = tskin_lnd(i)
              if( lsm == lsm_ruc) then
                QSFCMR_lnd(I)=QSFC_lnd(I)/(1.-QSFC_lnd(I))       !mixing ratio
              else
               TABS = 0.5*(TSK_lnd(I) + T1D(I))
               IF (TABS .LT. 273.15) THEN
                  !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
                  E1=SVP1*EXP(4648*(1./273.15 - 1./TABS) - &
                    & 11.64*LOG(273.15/TABS) + 0.02265*(273.15 - TABS))
               ELSE
                  !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
                  E1=SVP1*EXP(SVP2*(TABS-SVPT0)/(TABS-SVP3))
               ENDIF
                 QSFC_lnd(I)=EP2*E1/(PSFC(I)-ep3*E1)             !specific humidity
                 QSFC_lnd(I)=0.5*(QSFC_lnd(I) + QSFC(I))
                 QSFCMR_lnd(I)=QSFC_lnd(I)/(1.-QSFC_lnd(I))       !mixing ratio
              endif ! lsm
              IF(QSFC_lnd(I)>1..or.QSFC_lnd(I)<0.) print *,' QSFC_lnd(I)',itimestep,i,QSFC_lnd(I),Tskin_lnd(i),tsurf_lnd(i),qsfc(i)
            ENDIF
            IF (icy(i)) THEN
              TSK_ice(I) = tskin_ice(i)
              if( lsm == lsm_ruc) then
                QSFCMR_ice(I)=QSFC_ice(I)/(1.-QSFC_ice(I))        !mixing ratio
              else
               IF (TSK_ice(I) .LT. 273.15) THEN
                  !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
                  E1=SVP1*EXP(4648*(1./273.15 - 1./TSK_ice(I)) - &
                    & 11.64*LOG(273.15/TSK_ice(I)) + 0.02265*(273.15 - TSK_ice(I)))
               ELSE
                  !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
                  E1=SVP1*EXP(SVP2*(TSK_ice(I)-SVPT0)/(TSK_ice(i)-SVP3))
               ENDIF
                 QSFC_ice(I)=EP2*E1/(PSFC(I)-ep3*E1)             !specific humidity
                 QSFCMR_ice(I)=EP2*E1/(PSFC(I)-E1)                !mixing ratio
              endif ! lsm
              IF(QSFC_ice(I)>1..or.QSFC_ice(I)<0.) print *,' QSFC_ice(I)',itimestep,i,QSFC_ice(I),TSK_ice(i)
            ENDIF

         ELSE

            ! Use what comes out of the NST, LSM, SICE after check
            IF (wet(i)) then
               TSK_wat(I) = tskin_wat(i)
               IF (TSK_wat(I) .LT. 273.15) THEN
                  !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
                  E1=SVP1*EXP(4648*(1./273.15 - 1./TSK_wat(I)) - &
                    & 11.64*LOG(273.15/TSK_wat(I)) + 0.02265*(273.15 - TSK_wat(I)))
               ELSE
                  !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
                  E1=SVP1*EXP(SVP2*(TSK_wat(I)-SVPT0)/(TSK_wat(i)-SVP3))
               ENDIF
               QSFC_wat(I)=EP2*E1/(PSFC(I)-ep3*E1)             !specific humidity
            ENDIF
            IF (dry(i).and.(QSFC_lnd(I)>1..or.QSFC_lnd(I)<0.)) then
               !print *,'bad QSFC_lnd(I)',itimestep,iter,i,QSFC_lnd(I),TSKin_lnd(I)
               TABS = 0.5*(TSKin_lnd(I) + T1D(I))
               IF (TABS .LT. 273.15) THEN
                  !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
                  E1=SVP1*EXP(4648*(1./273.15 - 1./TABS) - &
                    & 11.64*LOG(273.15/TABS) + 0.02265*(273.15 - TABS))
               ELSE
                  !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
                  E1=SVP1*EXP(SVP2*(TABS-SVPT0)/(TABS-SVP3))
               ENDIF
                 QSFC_lnd(I)=EP2*E1/(PSFC(I)-ep3*E1)             !specific humidity
                 QSFC_lnd(I)=0.5*(QSFC_lnd(I) + QSFC(I))
            ENDIF
            IF (icy(i).and.(QSFC_ice(I)>1..or.QSFC_ice(I)<0.)) then
               !print *,'bad QSFC_ice(I)',itimestep,iter,i,QSFC_ice(I),TSKin_ice(I)
               IF (TSKin_ice(I) .LT. 273.15) THEN
                  !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
                  E1=SVP1*EXP(4648*(1./273.15 - 1./TSKin_ice(I)) - &
                    & 11.64*LOG(273.15/TSKin_ice(I)) + 0.02265*(273.15 - TSKin_ice(I)))
               ELSE
                  !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
                  E1=SVP1*EXP(SVP2*(TSKin_ice(I)-SVPT0)/(TSKin_ice(i)-SVP3))
               ENDIF
                 QSFC_ice(I)=EP2*E1/(PSFC(I)-ep3*E1)             !specific humidity
            ENDIF

            IF (wet(i)) QSFCMR_wat(I)=QSFC_wat(I)/(1.-QSFC_wat(I))
            IF (dry(i)) QSFCMR_lnd(I)=QSFC_lnd(I)/(1.-QSFC_lnd(I))
            IF (icy(i)) QSFCMR_ice(I)=QSFC_ice(I)/(1.-QSFC_ice(I))

         ENDIF
       endif ! flag_iter
      ENDDO

!$acc serial present(pblh, PSFCPA, dz8w1d, qflx, hflx,                       &
!$acc      dry, tskin_lnd, tsurf_lnd, qsfc_lnd, znt_lnd, ust_lnd, snowh_lnd, &
!$acc      icy, tskin_ice, tsurf_ice, qsfc_ice, znt_ice, ust_ice, snowh_ice, &
!$acc      wet, tskin_wat, tsurf_wat, qsfc_wat, znt_wat, ust_wat, snowh_wat)
      IF (debug_code >= 1) THEN
        write(0,*)"ITIMESTEP=",ITIMESTEP," iter=",iter
        DO I=its,ite
           write(0,*)"=== important input to mynnsfclayer, i:", i
           IF (dry(i)) THEN
             write(0,*)"dry=",dry(i)," pblh=",pblh(i)," tsk=", tskin_lnd(i),&
             " tsurf=", tsurf_lnd(i)," qsfc=", qsfc_lnd(i)," znt=", znt_lnd(i),&
             " ust=", ust_lnd(i)," snowh=", snowh_lnd(i)," psfcpa=",PSFCPA(i),  &
             " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
           ENDIF
           IF (icy(i)) THEN
             write(0,*)"icy=",icy(i)," pblh=",pblh(i)," tsk=", tskin_ice(i),&
             " tsurf=", tsurf_ice(i)," qsfc=", qsfc_ice(i)," znt=", znt_ice(i),&
             " ust=", ust_ice(i)," snowh=", snowh_ice(i),"psfcpa=",PSFCPA(i),  &
             " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
           ENDIF
           IF (wet(i)) THEN
             write(0,*)"wet=",wet(i)," pblh=",pblh(i)," tsk=", tskin_wat(i),&
             " tsurf=", tsurf_wat(i)," qsfc=", qsfc_wat(i)," znt=", znt_wat(i),&
             " ust=", ust_wat(i)," snowh=", snowh_wat(i),"psfcpa=",PSFCPA(i),  &
             " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
           ENDIF
        ENDDO
      ENDIF
!$acc end serial

!$acc parallel loop present(PSFC, PSFCPA, QVSH, QV1D, THCON, flag_iter,        &
!$acc       dry, tskin_lnd, TSK_lnd, tsurf_lnd, THSK_lnd, THVSK_lnd, qsfc_lnd, &
!$acc       icy, tskin_ice, TSK_ice, tsurf_ice, THSK_ice, THVSK_ice, qsfc_ice, &
!$acc       wet, tskin_wat, TSK_wat, tsurf_wat, THSK_wat, THVSK_wat, qsfc_wat)
      DO I=its,ite
         ! PSFC ( in cmb) is used later in saturation checks
         PSFC(I)=PSFCPA(I)/1000.
         QVSH(I)=QV1D(I)/(1.+QV1D(I))        !CONVERT TO SPEC HUM (kg/kg)
         THCON(I)=(100000./PSFCPA(I))**ROVCP
        if( flag_iter(i) ) then
         ! DEFINE SKIN TEMPERATURES FOR LAND/WATER/ICE
         if(dry(i)) then
           TSK_lnd(I) = tskin_lnd(i)
           !TSK_lnd(I) = 0.5 * (tsurf_lnd(i)+tskin_lnd(i))
           ! CONVERT SKIN TEMPERATURES TO POTENTIAL TEMPERATURE:
           THSK_lnd(I) = TSK_lnd(I)*THCON(I)   !(K)
           THVSK_lnd(I) = THSK_lnd(I)*(1.+EP1*qsfc_lnd(I))
           if(THVSK_lnd(I) < 160. .or. THVSK_lnd(I) > 390.) &
           print *,'THVSK_lnd(I)',itimestep,i,THVSK_lnd(I),THSK_lnd(i),tsurf_lnd(i),tskin_lnd(i),qsfc_lnd(i)
         endif
         if(icy(i)) then
           TSK_ice(I) = tskin_ice(i)
           !TSK_ice(I) = 0.5 * (tsurf_ice(i)+tskin_ice(i))
           ! CONVERT SKIN TEMPERATURES TO POTENTIAL TEMPERATURE:
           THSK_ice(I) = TSK_ice(I)*THCON(I)   !(K)
           THVSK_ice(I) = THSK_ice(I)*(1.+EP1*qsfc_ice(I))   !(K)
           if(THVSK_ice(I) < 160. .or. THVSK_ice(I) > 390.) &
           print *,'THVSK_ice(I)',itimestep,i,THVSK_ice(I),THSK_ice(i),tsurf_ice(i),tskin_ice(i),qsfc_ice(i)
         endif
         if(wet(i)) then
           TSK_wat(I) = tskin_wat(i)
           !TSK_wat(I) = 0.5 * (tsurf_wat(i)+tskin_wat(i))
           ! CONVERT SKIN TEMPERATURES TO POTENTIAL TEMPERATURE:
           THSK_wat(I) = TSK_wat(I)*THCON(I)   !(K)
           THVSK_wat(I) = THSK_wat(I)*(1.+EP1*QVSH(I))   !(K)
           if(THVSK_wat(I) < 160. .or. THVSK_wat(I) > 390.) &
           print *,'THVSK_wat(I)',i,THVSK_wat(I),THSK_wat(i),tsurf_wat(i),tskin_wat(i),qsfc_wat(i)
         endif
        endif ! flag_iter
      ENDDO

!$acc parallel loop present(TH1D, T1D, P1D, TC1D)
      DO I=its,ite
         ! CONVERT LOWEST LAYER TEMPERATURE TO POTENTIAL TEMPERATURE:
         TH1D(I)=T1D(I)*(100000./P1D(I))**ROVCP  !(Theta, K)
         TC1D(I)=T1D(I)-273.15                   !(T, Celsius)
      ENDDO

!$acc parallel loop present(THV1D, TH1D, QVSH, TV1D, T1D)
      DO I=its,ite
         ! CONVERT TO VIRTUAL TEMPERATURE
         THV1D(I)=TH1D(I)*(1.+EP1*QVSH(I))             !(K)
         TV1D(I)=T1D(I)*(1.+EP1*QVSH(I))               !(K)
      ENDDO

!$acc parallel loop present(RHO1D, P1D, TV1D, TH1D, ZA, ZA2, dz2w1d, dz8w1d, GOVRTH)
      DO I=its,ite
         RHO1D(I)=P1D(I)/(Rd*TV1D(I))     !now using value calculated in sfc driver
         ZA(I)=0.5*dz8w1d(I)              !height of first half-sigma level
         ZA2(I)=dz8w1d(I) + 0.5*dz2w1d(I) !height of 2nd half-sigma level
         GOVRTH(I)=grav/TH1D(I)
      ENDDO

      !tgs - should QFX and HFX be separate for land, ice and water?
!$acc parallel loop present(QFX, QFLX, RHO1D, HFX, HFLX)
      DO I=its,ite
         QFX(i)=QFLX(i)*RHO1D(I)
         HFX(i)=HFLX(i)*RHO1D(I)*cp
      ENDDO

!$acc serial present(THV1D, TV1D, RHO1D, GOVRTH, &
!$acc                dry, tsk_lnd, thvsk_lnd,    &
!$acc                icy, tsk_ice, thvsk_ice,    &
!$acc                wet, tsk_wat, thvsk_wat)
      IF (debug_code ==2) THEN
        !write(*,*)"ITIMESTEP=",ITIMESTEP
        DO I=its,ite
          write(*,*)"=== derived quantities in mynn sfc layer, i:", i
          write(*,*)" land,      ice,      water"
          write(*,*)"dry=",dry(i)," icy=",icy(i)," wet=",wet(i)
          write(*,*)"tsk=", tsk_lnd(i),tsk_ice(i),tsk_wat(i)
          write(*,*)"thvsk=", thvsk_lnd(i),thvsk_ice(i),thvsk_wat(i)
          write(*,*)"THV1D=", THV1D(i)," TV1D=",TV1D(i)
          write(*,*)"RHO1D=", RHO1D(i)," GOVRTH=",GOVRTH(i)
        ENDDO
      ENDIF
!$acc end serial

!$acc parallel loop present(T1D,P1D,QGH,QV1D,CPM)
      DO I=its,ite
         ! QGH CHANGED TO USE LOWEST-LEVEL AIR TEMP
         ! Q2SAT = QGH IN LSM
         IF (T1D(I) .LT. 273.15) THEN
            !SATURATION VAPOR PRESSURE WRT ICE
            E1=SVP1*EXP(4648.*(1./273.15 - 1./T1D(I)) - &
            &  11.64*LOG(273.15/T1D(I)) + 0.02265*(273.15 - T1D(I)))
         ELSE
            !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
            E1=SVP1*EXP(SVP2*(T1D(I)-SVPT0)/(T1D(I)-SVP3))
         ENDIF
         PL=P1D(I)/1000.
         !QGH(I)=EP2*E1/(PL-ep3*E1)    !specific humidity
         QGH(I)=EP2*E1/(PL-E1)          !mixing ratio
         CPM(I)=CP*(1.+0.84*QV1D(I))
      ENDDO

!$acc serial present(QGH,                       &
!$acc                wet, QSFC_wat, QSFCMR_wat, &
!$acc                dry, QSFC_lnd, QSFCMR_lnd, &
!$acc                icy, QSFC_ice, QSFCMR_ice)
      IF (debug_code == 2) THEN
         write(*,*)"ITIMESTEP=",ITIMESTEP
         DO I=its,ite
            if (wet(i)) then
               write(*,*)"==== q-bombs, i:",i," wet"
               write(*,*)"QSFC_wat=", QSFC_wat(I)," QSFCMR_wat=", QSFCMR_wat(I)," QGH=",QGH(I)
            endif
            if(dry(i)) then
               write(*,*)"==== q-bombs, i:",i," dry"
               write(*,*)"QSFC_lnd=", QSFC_lnd(I)," QSFCMR_lnd=", QSFCMR_lnd(I)," QGH=",QGH(I)
            endif
            if(icy(i)) then
               write(*,*)"==== q-bombs, i:",i," ice"
               write(*,*)"QSFC_ice=", QSFC_ice(I)," QSFCMR_ice=", QSFCMR_ice(I)," QGH=",QGH(I)
            endif
         ENDDO
      ENDIF
!$acc end serial

!$acc parallel loop present(flag_iter,U1D,V1D,WSPD,wet,dry,icy,    &
!$acc                    THV1D,THVSK_wat,THVSK_lnd,THVSK_ice,   &
!$acc                    hfx,RHO1D,qfx,WSTAR,pblh,dx,GOVRTH,ZA, &
!$acc                    TSK_wat,TSK_lnd,TSK_ice,               &
!$acc                    rb_wat,rb_lnd,rb_ice)
      DO I=its,ite
        if( flag_iter(i) ) then
         ! DH* 20200401 - note. A weird bug in Intel 18 on hera prevents using the
         ! normal -O2 optimization in Release mode for this file. Not reproducible
         ! by every user, the bug manifests itself in the resulting wind speed WSPD(I)
         ! being -99.0 despite the assignments in lines 932 and 933. *DH
         WSPD(I)=SQRT(U1D(I)*U1D(I)+V1D(I)*V1D(I))
         WSPD_wat = -99.
         WSPD_ice = -99.
         WSPD_lnd = -99.

         IF (wet(i)) THEN
            DTHVDZ=(THV1D(I)-THVSK_wat(I))
            !--------------------------------------------------------
            ! Calculate the convective velocity scale (WSTAR) and
            ! subgrid-scale velocity (VSGD) following Beljaars (1995, QJRMS)
            ! and Mahrt and Sun (1995, MWR), respectively
            !-------------------------------------------------------
            !tgs - the line below could be used when hflx_wat,qflx_wat are moved from
            !      Interstitial to Sfcprop
            !fluxc = max(hflx_wat(i) + ep1*THVSK_wat(I)*qflx_wat(i),0.)
            fluxc = max(hfx(i)/RHO1D(i)/cp                    &
            &    + ep1*THVSK_wat(I)*qfx(i)/RHO1D(i),0._kind_phys)
            !WSTAR(I) = vconvc*(grav/TSK(i)*pblh(i)*fluxc)**onethird
            WSTAR(I) = vconvc*(grav/TSK_wat(i)*pblh(i)*fluxc)**onethird
            !--------------------------------------------------------
            ! Mahrt and Sun low-res correction - modified for water points (halved)
            ! (for 13 km ~ 0.18 m/s; for 3 km == 0 m/s)
            !--------------------------------------------------------
            VSGD = MIN( 0.25 * (max(dx(i)/5000.-1.,0._kind_phys))**onethird , 0.5_kind_phys)
            WSPD_wat=SQRT(WSPD(I)*WSPD(I)+WSTAR(I)*WSTAR(I)+vsgd*vsgd)
            WSPD_wat=MAX(WSPD_wat,wmin)
            !--------------------------------------------------------
            ! CALCULATE THE BULK RICHARDSON NUMBER OF SURFACE LAYER,
            ! ACCORDING TO AKB(1976), EQ(12).
            !--------------------------------------------------------
            rb_wat(I)=GOVRTH(I)*ZA(I)*DTHVDZ/(WSPD_wat*WSPD_wat)
            rb_wat(I)=MAX(rb_wat(I),-2.0_kind_phys)
            rb_wat(I)=MIN(rb_wat(I), 2.0_kind_phys)
         ENDIF ! end water point

         IF (dry(i)) THEN
            DTHVDZ=(THV1D(I)-THVSK_lnd(I))
            !--------------------------------------------------------
            ! Calculate the convective velocity scale (WSTAR) and
            ! subgrid-scale velocity (VSGD) following Beljaars (1995, QJRMS)
            ! and Mahrt and Sun (1995, MWR), respectively
            !-------------------------------------------------------
            !tgs - the line below could be used when hflx_lnd,qflx_wat are moved from
            !      Interstitial to Sfcprop
            !fluxc = max(hflx_lnd(i) + ep1*THVSK_lnd(I)*qflx_lnd(i),0.)
            fluxc = max(hfx(i)/RHO1D(i)/cp                    &
            &    + ep1*THVSK_lnd(I)*qfx(i)/RHO1D(i),0._kind_phys)
            ! WSTAR(I) = vconvc*(g/TSK(i)*pblh(i)*fluxc)**onethird
            ! increase height scale, assuming that the non-local transoport
            ! from the mass-flux (plume) mixing exceedsd the PBLH.
            WSTAR(I) = vconvc*(grav/TSK_lnd(i)*MIN(1.5*pblh(i),4000._kind_phys)*fluxc)**onethird
            !--------------------------------------------------------
            ! Mahrt and Sun low-res correction
            ! (for 13 km ~ 0.37 m/s; for 3 km == 0 m/s)
            !--------------------------------------------------------
            VSGD = MIN( 0.32 * (max(dx(i)/5000.-1.,0._kind_phys))**onethird , 0.5_kind_phys)
            WSPD_lnd=SQRT(WSPD(I)*WSPD(I)+WSTAR(I)*WSTAR(I)+vsgd*vsgd)
            WSPD_lnd=MAX(WSPD_lnd,wmin)
            !--------------------------------------------------------
            ! CALCULATE THE BULK RICHARDSON NUMBER OF SURFACE LAYER,
            ! ACCORDING TO AKB(1976), EQ(12).
            !--------------------------------------------------------
            rb_lnd(I)=GOVRTH(I)*ZA(I)*DTHVDZ/(WSPD_lnd*WSPD_lnd)
            !From Tilden Meyers:
            !IF (rb_lnd(I) .GE 0.0) THEN
            !   ust_lnd(i)=WSPD_lnd*0.1/(1.0 + 10.0*rb_lnd(I))
            !ELSE
            !   ust_lnd(i)=WSPD_lnd*0.1*(1.0 - 10.0*rb_lnd(I))**onethird
            !ENDIF
            rb_lnd(I)=MAX(rb_lnd(I),-2.0_kind_phys)
            rb_lnd(I)=MIN(rb_lnd(I), 2.0_kind_phys)
         ENDIF ! end land point

         IF (icy(i)) THEN
            DTHVDZ=(THV1D(I)-THVSK_ice(I))
            !--------------------------------------------------------
            ! Calculate the convective velocity scale (WSTAR) and
            ! subgrid-scale velocity (VSGD) following Beljaars (1995, QJRMS)
            ! and Mahrt and Sun (1995, MWR), respectively
            !-------------------------------------------------------
            !tgs - the line below could be used when hflx_ice,qflx_ice are moved from
            !      Interstitial to Sfcprop
            !fluxc = max(hflx_ice(i) + ep1*THVSK_ice(I)*qflx_ice(i)/RHO1D(i),0.)
            fluxc = max(hfx(i)/RHO1D(i)/cp                    &
            &    + ep1*THVSK_ice(I)*qfx(i)/RHO1D(i),0._kind_phys)
            ! WSTAR(I) = vconvc*(g/TSK(i)*pblh(i)*fluxc)**onethird
            ! increase height scale, assuming that the non-local transport
            ! from the mass-flux (plume) mixing exceedsd the PBLH.
            WSTAR(I) = vconvc*(grav/TSK_ice(i)*MIN(1.5*pblh(i),4000._kind_phys)*fluxc)**onethird
            !--------------------------------------------------------
            ! Mahrt and Sun low-res correction
            ! (for 13 km ~ 0.37 m/s; for 3 km == 0 m/s)
            !--------------------------------------------------------
            VSGD = MIN( 0.32 * (max(dx(i)/5000.-1.,0._kind_phys))**onethird , 0.5_kind_phys)
            WSPD_ice=SQRT(WSPD(I)*WSPD(I)+WSTAR(I)*WSTAR(I)+vsgd*vsgd)
            WSPD_ice=MAX(WSPD_ice,wmin)
            !--------------------------------------------------------
            ! CALCULATE THE BULK RICHARDSON NUMBER OF SURFACE LAYER,
            ! ACCORDING TO AKB(1976), EQ(12).
            !--------------------------------------------------------
            rb_ice(I)=GOVRTH(I)*ZA(I)*DTHVDZ/(WSPD_ice*WSPD_ice)
            rb_ice(I)=MAX(rb_ice(I),-2.0_kind_phys)
            rb_ice(I)=MIN(rb_ice(I), 2.0_kind_phys)
         ENDIF ! end ice point

         !NOW CONDENSE THE POSSIBLE WSPD VALUES BY TAKING THE MAXIMUM
         WSPD(I) = MAX(WSPD_ice,WSPD_wat)
         WSPD(I) = MAX(WSPD_lnd,WSPD(I))

         IF (debug_code == 2) THEN
            write(*,*)"===== After rb calc in mynn sfc layer:"
            write(*,*)"ITIMESTEP=",ITIMESTEP
            write(*,*)"WSPD=", WSPD(I)," WSTAR=", WSTAR(I)," vsgd=",vsgd
            IF (icy(i))write(*,*)"rb_ice=", rb_ice(I)," DTHVDZ=",DTHVDZ
            IF (wet(i))write(*,*)"rb_wat=", rb_wat(I)," DTHVDZ=",DTHVDZ
            IF (dry(i))write(*,*)"rb_lnd=", rb_lnd(I)," DTHVDZ=",DTHVDZ
         ENDIF

         ! IF PREVIOUSLY UNSTABLE, DO NOT LET INTO REGIMES 1 AND 2 (STABLE)
         !if (itimestep .GT. 1) THEN
         !    IF(MOL(I).LT.0.)BR(I)=MIN(BR(I),0.0)
         !ENDIF

        endif ! flag_iter
      ENDDO

 1006   format(A,F7.3,A,f9.4,A,f9.5,A,f9.4)
 1007   format(A,F2.0,A,f6.2,A,f7.3,A,f7.2)

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--- BEGIN I-LOOP
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!$acc parallel loop present(flag_iter, PSFCPA, dz8w1d, pblh, &
!$acc                    device_errmsg, device_errflg, &
!$acc                    device_special_errmsg, device_special_errflg, &
!$acc                    wet, dry, icy, &
!$acc                    ZT_wat, ZT_lnd, ZT_ice, &
!$acc                    ZNT_wat, ZNT_lnd, ZNT_ice, &
!$acc                    ZNTstoch_wat, ZNTstoch_lnd, ZNTstoch_ice, &
!$acc                    UST_wat, UST_lnd, UST_ice, &
!$acc                    ZQ_wat, ZQ_lnd, ZQ_ice, &
!$acc                    snowh_wat, snowh_lnd, snowh_ice, &
!$acc                    THVSK_wat, THVSK_lnd, THVSK_ice, &
!$acc                    tskin_wat, tskin_lnd, tskin_ice, &
!$acc                    tsurf_wat, tsurf_lnd, tsurf_ice, &
!$acc                    qsfc_wat, qsfc_lnd, qsfc_ice, &
!$acc                    GZ1OZ0_wat, GZ1OZt_wat, GZ2OZ0_wat, GZ2OZt_wat, GZ10OZ0_wat, GZ10OZt_wat, &
!$acc                    GZ1OZ0_lnd, GZ1OZt_lnd, GZ2OZ0_lnd, GZ2OZt_lnd, GZ10OZ0_lnd, GZ10OZt_lnd, &
!$acc                    GZ1OZ0_ice, GZ1OZt_ice, GZ2OZ0_ice, GZ2OZt_ice, GZ10OZ0_ice, GZ10OZt_ice, &
!$acc                    zratio_wat, zratio_lnd, zratio_ice, &
!$acc                    stress_wat, stress_lnd, stress_ice, &
!$acc                    rb_wat, rb_lnd, rb_ice, &
!$acc                    qflx, qflx_lnd, &
!$acc                    hflx, hflx_lnd, &
!$acc                    psim, psih, psim10, psih10, psih2, &
!$acc                    psix_wat, psix10_wat, psit_wat, psit2_wat, psiq_wat, psiq2_wat, &
!$acc                    psix_lnd, psix10_lnd, psit_lnd, psit2_lnd, psiq_lnd, psiq2_lnd, &
!$acc                    psix_ice, psix10_ice, psit_ice, psit2_ice, psiq_ice, psiq2_ice, &
!$acc                    WSPD, WSPDI, U1D, V1D, TC1D, THV1D, rstoch1D, USTM, ZA, ZOL, QVSH, &
!$acc                    shdmax, vegtype, z0pert, ztpert, mol, rmol, wstar, qstar, sigmaf)

 DO I=its,ite
   if( flag_iter(i) ) then

    !COMPUTE KINEMATIC VISCOSITY (m2/s) Andreas (1989) CRREL Rep. 89-11
    !valid between -173 and 277 degrees C.
    VISC=1.326e-5*(1. + 6.542e-3*TC1D(I) + 8.301e-6*TC1D(I)*TC1D(I) &
                      - 4.84e-9*TC1D(I)*TC1D(I)*TC1D(I))

    IF (wet(i)) THEN
       !--------------------------------------
       ! WATER
       !--------------------------------------
       if (sfc_z0_type >= 0) then ! Avoid calculation is using wave model
          ! CALCULATE z0 (znt)
          !--------------------------------------

          IF (debug_code == 2) THEN
            write(*,*)"=============Input to ZNT over water:"
            write(*,*)"u*:",UST_wat(i)," wspd=",WSPD(i)," visc=",visc," za=",ZA(I)
          ENDIF

#if defined SWAN_COUPLING || defined WW3_COUPLING
# if defined COARE_TAYLOR_YELLAND
          ZNT_wat(I)=MAX(1200.0*HWAVE(I)*                               &
     &           (HWAVE(I)/(LWAVEP(I)+0.001))**4.5+                     &
     &            0.11*VISC/(UST_wat(I)+0.001),1.59E-5)
# elif defined DRENNAN
          CWAVE=MAX(LWAVEP(I)/(PWAVE(I)+0.001),0.1)
          ZNT_wat(I)=MAX(3.35*HWAVE(I)*(MIN(UST_wat(I)/CWAVE,0.1))**3.4+&
     &           0.11*VISC/(UST_wat(I)+0.001),1.59E-5)
# elif defined COARE_OOST
          CWAVE=MAX(LWAVEP(I)/(PWAVE(I)+0.001),0.1)
          ZNT_wat(I)=MAX(25.0/3.141593*LWAVEP(I)*                       &
     &           (MIN(UST_wat(I)/CWAVE,0.1))**4.5+                      &
     &           0.11*VISC/(UST_wat(I)+0.001),1.59E-5)
# elif defined Z0_WAV_SIN
          ZNT_wat(I)=MAX(Z0_WAV(I),1.59E-5)
# elif defined Z0_PORCHETTA
          CWAVE  = MAX(LWAVEP(I)/(PWAVE(I)+0.001),0.1)
          UREAL  = U10(I)*COSA(I)-V10(I)*SINA(I)
          VREAL  = V10(I)*COSA(I)+U10(I)*SINA(I)
          ! radiants to degrees
          R2D  = 45.0/atan(1.0)
          ! wind direction cartesian convention: zero TO north (positive clockwise)
          DWIND=atan2(UREAL,VREAL)*R2D
          ! wind direction nautical convention: zero FROM north (positive clockwise)
          DWIND = DWIND + 180.
          ! wave direction nautical convention: zero FROM north (positive clockwise)
          DWAVE1 = MOD(DWAVEP*R2D + 360.0, 360.0) 
          ! smallest angle between wind and wave direction
          THWV   = MIN(ABS(DWIND-DWAVE1),360.-(ABS(DWIND-DWAVE1)))
          ! Porchetta et al (2019) parameterization
          Z0WAVE = 20.0 * HWAVE(I) * COS(0.45*THWV/R2D) *              &
     &           (MIN(UST_wat(I)/CWAVE,0.1))**(3.8*COS(-0.32*THWV/R2D))
          ! add standard viscous and bounding terms
          ZNT_wat(I)=MAX(Z0WAVE+0.11*VISC/(UST_wat(I)+0.001),1.59E-5)
# else
          ZNT_wat(I)=0.016*UST_wat(I)*UST_wat(I)/G+1.5E-5/ust(i)
# endif
#else
          IF ( PRESENT(ISFTCFLX) ) THEN
             IF ( ISFTCFLX .EQ. 0 ) THEN
                IF (COARE_OPT .EQ. 3.0) THEN
                   !COARE 3.0 (MISLEADING SUBROUTINE NAME)
                   CALL charnock_1955(ZNT_wat(i),UST_wat(i),WSPD(i),visc,ZA(I))
                ELSE
                   !COARE 3.5
                   CALL edson_etal_2013(ZNT_wat(i),UST_wat(i),WSPD(i),visc,ZA(I))
                ENDIF
             ELSEIF ( ISFTCFLX .EQ. 1 .OR. ISFTCFLX .EQ. 2 ) THEN
                CALL davis_etal_2008(ZNT_wat(i),UST_wat(i))
             ELSEIF ( ISFTCFLX .EQ. 3 ) THEN
                CALL Taylor_Yelland_2001(ZNT_wat(i),UST_wat(i),WSPD(i))
             ELSEIF ( ISFTCFLX .EQ. 4 ) THEN
                !GFS surface layer scheme
                CALL GFS_z0_wat(ZNT_wat(i),UST_wat(i),WSPD(i),ZA(I),sfc_z0_type,redrag)
             ENDIF
          ELSE
             !DEFAULT TO COARE 3.0/3.5
             IF (COARE_OPT .EQ. 3.0) THEN
                !COARE 3.0
                CALL charnock_1955(ZNT_wat(i),UST_wat(i),WSPD(i),visc,ZA(I))
             ELSE
                !COARE 3.5
                CALL edson_etal_2013(ZNT_wat(i),UST_wat(i),WSPD(i),visc,ZA(I))
             ENDIF
          ENDIF
       endif !-end wave model check
#endif
#if defined DRAGLIM_DAVIS
          ZNT(I)=MIN(ZNT(I),2.85E-3) !Davis limiting
#endif

       ! add stochastic perturbation of ZNT
       if (spp_sfc==1) then
          ZNTstoch_wat(I)  = MAX(ZNT_wat(I) + ZNT_wat(I)*1.0*rstoch1D(i), 1e-6_kind_phys)
       else
          ZNTstoch_wat(I)  = ZNT_wat(I)
       endif

       IF (debug_code > 1) THEN
          write(*,*)"==========Output ZNT over water:"
          write(*,*)"ZNT:",ZNTstoch_wat(i)
       ENDIF

       !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING NEW ZNT
       ! AHW: Garrattt formula: Calculate roughness Reynolds number
       !      Kinematic viscosity of air (linear approx to
       !      temp dependence at sea level)
       restar=MAX(ust_wat(i)*ZNTstoch_wat(i)/visc, 0.1_kind_phys)

       !--------------------------------------
       !CALCULATE z_t and z_q
       !--------------------------------------
       IF (debug_code > 1) THEN
          write(*,*)"=============Input to ZT over water:"
          write(*,*)"u*:",UST_wat(i)," restar=",restar," visc=",visc
       ENDIF

       IF ( PRESENT(ISFTCFLX) ) THEN
          IF ( ISFTCFLX .EQ. 0 ) THEN
             IF (COARE_OPT .EQ. 3.0) THEN
                CALL fairall_etal_2003(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                       rstoch1D(i),spp_sfc)
             ELSE
                CALL fairall_etal_2014(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                       rstoch1D(i),spp_sfc)
             ENDIF
          ELSEIF ( ISFTCFLX .EQ. 1 ) THEN
             IF (COARE_OPT .EQ. 3.0) THEN
                CALL fairall_etal_2003(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                       rstoch1D(i),spp_sfc)
             ELSE
                CALL fairall_etal_2014(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                       rstoch1D(i),spp_sfc)
             ENDIF
          ELSEIF ( ISFTCFLX .EQ. 2 ) THEN
             CALL garratt_1992(ZT_wat(i),ZQ_wat(i),ZNTstoch_wat(i),restar,2.0_kind_phys)
          ELSEIF ( ISFTCFLX .EQ. 3 ) THEN
             IF (COARE_OPT .EQ. 3.0) THEN
                CALL fairall_etal_2003(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                       rstoch1D(i),spp_sfc)
             ELSE
                CALL fairall_etal_2014(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                       rstoch1D(i),spp_sfc)
             ENDIF
          ELSEIF ( ISFTCFLX .EQ. 4 ) THEN
             !GFS zt formulation
             CALL GFS_zt_wat(ZT_wat(i),ZNTstoch_wat(i),restar,WSPD(i),ZA(i),sfc_z0_type,device_errmsg,device_errflg)
             if(errflg/=0) return

             ZQ_wat(i)=ZT_wat(i)
          ENDIF
       ELSE
          !DEFAULT TO COARE 3.0/3.5
          IF (COARE_OPT .EQ. 3.0) THEN
             CALL fairall_etal_2003(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                    rstoch1D(i),spp_sfc)
          ELSE
             CALL fairall_etal_2014(ZT_wat(i),ZQ_wat(i),restar,UST_wat(i),visc,&
                                    rstoch1D(i),spp_sfc)
          ENDIF
       ENDIF

       IF (debug_code > 1) THEN
         write(*,*)"=============Output ZT & ZQ over water:"
         write(*,*)"ZT:",ZT_wat(i)," ZQ:",ZQ_wat(i)
       ENDIF

       GZ1OZ0_wat(I)= LOG((ZA(I)+ZNTstoch_wat(I))/ZNTstoch_wat(I))
       GZ1OZt_wat(I)= LOG((ZA(I)+ZNTstoch_wat(i))/ZT_wat(i))
       GZ2OZ0_wat(I)= LOG((2.0+ZNTstoch_wat(I))/ZNTstoch_wat(I))
       GZ2OZt_wat(I)= LOG((2.0+ZNTstoch_wat(i))/ZT_wat(i))
       GZ10OZ0_wat(I)=LOG((10.+ZNTstoch_wat(I))/ZNTstoch_wat(I))
       GZ10OZt_wat(I)=LOG((10.+ZNTstoch_wat(i))/ZT_wat(i))
       zratio_wat(i)=ZNTstoch_wat(I)/ZT_wat(I)   !need estimate for Li et al.

    ENDIF !end water point

    IF (dry(I)) THEN

       if ( IZ0TLND .EQ. 4 ) then
          CALL GFS_z0_lnd(ZNT_lnd(i),shdmax(i),ZA(i),vegtype(i),ivegsrc,z0pert(i))
       endif

       ! add stochastic perturbaction of ZNT
       if (spp_sfc==1) then
          ZNTstoch_lnd(I)  = MAX(ZNT_lnd(I) + ZNT_lnd(I)*1.0*rstoch1D(i), 1e-6_kind_phys)
       else
          ZNTstoch_lnd(I)  = ZNT_lnd(I)
       endif
       !add limit to prevent ridiculous values of z0 (more than dz/15)
       ZNTstoch_lnd(I) = min(ZNTstoch_lnd(I), dz8w1d(i)*0.0666_kind_phys)

       !--------------------------------------
       ! LAND
       !--------------------------------------
       !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING DEFAULT ZNT
       restar=MAX(ust_lnd(i)*ZNTstoch_lnd(i)/visc, 0.1_kind_phys)

       !--------------------------------------
       !GET z_t and z_q
       !--------------------------------------
       IF (snowh_lnd(i) > 50.) THEN  ! (mm) Treat as snow cover - use Andreas
          CALL Andreas_2002(ZNTstoch_lnd(i),visc,ust_lnd(i),ZT_lnd(i),ZQ_lnd(i))
       ELSE
          IF ( PRESENT(IZ0TLND) ) THEN
             IF ( IZ0TLND .LE. 1 ) THEN
                CALL zilitinkevich_1995(ZNTstoch_lnd(i),ZT_lnd(i),ZQ_lnd(i),restar,&
                      UST_lnd(I),KARMAN,1.0_kind_phys,IZ0TLND,spp_sfc,rstoch1D(i))
             ELSEIF ( IZ0TLND .EQ. 2 ) THEN
                ! DH note - at this point, qstar is either not initialized
                ! or initialized to zero, but certainly not set correctly
                device_special_errmsg = 'Logic error: qstar is not set correctly when calling Yang_2008'
                device_special_errflg = 1
#ifndef _OPENACC
! Necessary since OpenACC does not support branching in parallel code
! Must sync errmsg and errflg with device_errmsg and device_errflg, respectively
! so that proper error message and error flag codes are returned.
                errmsg = device_special_errmsg
                errflg = device_special_errflg
                return
#endif
                CALL Yang_2008(ZNTSTOCH_lnd(i),ZT_lnd(i),ZQ_lnd(i),UST_lnd(i),MOL(I),&
                              qstar(I),restar,visc)
             ELSEIF ( IZ0TLND .EQ. 3 ) THEN
                !Original MYNN in WRF-ARW used this form:
                CALL garratt_1992(ZT_lnd(i),ZQ_lnd(i),ZNTSTOCH_lnd(i),restar,1.0_kind_phys)
             ELSEIF ( IZ0TLND .EQ. 4 ) THEN
                !GFS:
                CALL GFS_zt_lnd(ZT_lnd(i),ZNTSTOCH_lnd(i),sigmaf(i),ztpert(i),UST_lnd(i))
                ZQ_lnd(i)=ZT_lnd(i)
             ENDIF
          ELSE
             !DEFAULT TO ZILITINKEVICH
             CALL zilitinkevich_1995(ZNTSTOCH_lnd(i),ZT_lnd(i),ZQ_lnd(i),restar,&
                         UST_lnd(I),KARMAN,1.0_kind_phys,0,spp_sfc,rstoch1D(i))
          ENDIF
       ENDIF

       IF (ZNTstoch_lnd(i) < 1E-8 .OR. Zt_lnd(i) < 1E-10) THEN
         write(0,*)"===(land) capture bad input in mynn sfc layer, i=:",i
         write(0,*)" ZNT=", ZNTstoch_lnd(i)," ZT=",Zt_lnd(i)
         write(0,*)" tsk=", tskin_lnd(i)," restar=",restar,&
         " tsurf=", tsurf_lnd(i)," qsfc=", qsfc_lnd(i)," znt=", znt_lnd(i),&
         " ust=", ust_lnd(i)," snowh=", snowh_lnd(i),"psfcpa=",PSFCPA(i), &
         " dz=",dz8w1d(i)," qflx=",qflx_lnd(i)," hflx=",hflx_lnd(i)," hpbl=",pblh(i)
       ENDIF

       GZ1OZ0_lnd(I)= LOG((ZA(I)+ZNTstoch_lnd(I))/ZNTstoch_lnd(I))
       GZ1OZt_lnd(I)= LOG((ZA(I)+ZNTstoch_lnd(i))/ZT_lnd(i))
       GZ2OZ0_lnd(I)= LOG((2.0+ZNTstoch_lnd(I))/ZNTstoch_lnd(I))
       GZ2OZt_lnd(I)= LOG((2.0+ZNTstoch_lnd(i))/ZT_lnd(i))
       GZ10OZ0_lnd(I)=LOG((10.+ZNTstoch_lnd(I))/ZNTstoch_lnd(I))
       GZ10OZt_lnd(I)=LOG((10.+ZNTstoch_lnd(i))/ZT_lnd(i))
       zratio_lnd(i)=ZNTstoch_lnd(I)/ZT_lnd(I)   !need estimate for Li et al.

    ENDIF !end land point

    IF (icy(I)) THEN

       ! add stochastic perturbaction of ZNT
       if (spp_sfc==1) then
          ZNTstoch_ice(I)  = MAX(ZNT_ice(I) + ZNT_ice(I)*1.0*rstoch1D(i), 1e-6_kind_phys)
       else
          ZNTstoch_ice(I)  = ZNT_ice(I)
       endif

       !--------------------------------------
       ! ICE
       !--------------------------------------
       !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING DEFAULT ZNT
       restar=MAX(ust_ice(i)*ZNTstoch_ice(i)/visc, 0.1_kind_phys)
       !--------------------------------------
       !GET z_t and z_q
       !--------------------------------------
       CALL Andreas_2002(ZNTstoch_ice(i),visc,ust_ice(i),ZT_ice(i),ZQ_ice(i))

       GZ1OZ0_ice(I)= LOG((ZA(I)+ZNTstoch_ice(I))/ZNTstoch_ice(I))
       GZ1OZt_ice(I)= LOG((ZA(I)+ZNTstoch_ice(i))/ZT_ice(i))
       GZ2OZ0_ice(I)= LOG((2.0+ZNTstoch_ice(I))/ZNTstoch_ice(I))
       GZ2OZt_ice(I)= LOG((2.0+ZNTstoch_ice(i))/ZT_ice(i))
       GZ10OZ0_ice(I)=LOG((10.+ZNTstoch_ice(I))/ZNTstoch_ice(I))
       GZ10OZt_ice(I)=LOG((10.+ZNTstoch_ice(i))/ZT_ice(i))
       zratio_ice(i)=ZNTstoch_ice(I)/ZT_ice(I)   !need estimate for Li et al.

    ENDIF  !end ice point

    !Capture a representative ZNT
    !tgs - should this be changed for fractional grid or fractional sea ice?
    IF (dry(i)) THEN
       ZNT(i)=ZNTstoch_lnd(I)
    ELSEIF (wet(i)) THEN
       ZNT(i)=ZNTstoch_wat(I)
    ELSEIF (icy(i)) THEN
       ZNT(i)=ZNTstoch_ice(I)
    ENDIF

    !--------------------------------------------------------------------
    !--- DIAGNOSE STABILITY FUNCTIONS FOR THE APPROPRIATE STABILITY CLASS:
    !    THE STABILITY CLASSES ARE DETERMINED BY THE BULK RICHARDSON NUMBER.
    !--------------------------------------------------------------------

    IF (wet(i)) THEN
       IF (rb_wat(I) .GT. 0.0) THEN

          IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
             !COMPUTE z/L first guess:
             CALL Li_etal_2010(ZOL(I),rb_wat(I),ZA(I)/ZNTstoch_wat(I),zratio_wat(I))
             !ZOL(I)=ZA(I)*KARMAN*grav*MOL(I)/(TH1D(I)*MAX(UST_wat(I)*UST_wat(I),0.0001))
             ZOL(I)=MAX(ZOL(I),0.0_kind_phys)
             ZOL(I)=MIN(ZOL(I),20._kind_phys)

             IF (debug_code >= 1) THEN
               IF (ZNTstoch_wat(i) < 1E-8 .OR. Zt_wat(i) < 1E-10) THEN
                 write(0,*)"===(wet) capture bad input in mynn sfc layer, i=:",i
                 write(0,*)"rb=", rb_wat(I)," ZNT=", ZNTstoch_wat(i)," ZT=",Zt_wat(i)
                 write(0,*)" tsk=", tskin_wat(i)," prev z/L=",ZOL(I),&
                 " tsurf=", tsurf_wat(i)," qsfc=", qsfc_wat(i)," znt=", znt_wat(i),&
                 " ust=", ust_wat(i)," snowh=", snowh_wat(i),"psfcpa=",PSFCPA(i),  &
                 " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
               ENDIF
             ENDIF

             !Use Pedros iterative function to find z/L
             !zol(I)=zolri(rb_wat(I),ZA(I),ZNTstoch_wat(I),ZT_wat(I),ZOL(I),psi_opt)
             !Use brute-force method
             zol(I)=zolrib(rb_wat(I),ZA(I),ZNTstoch_wat(I),zt_wat(I),GZ1OZ0_wat(I),GZ1OZt_wat(I),ZOL(I),psi_opt)
          ENDIF ! restart
          ZOL(I)=MAX(ZOL(I),0.0_kind_phys)
          ZOL(I)=MIN(ZOL(I),20._kind_phys)

          zolzt = zol(I)*zt_wat(I)/ZA(I)                ! zt/L
          zolz0 = zol(I)*ZNTstoch_wat(I)/ZA(I)          ! z0/L
          zolza = zol(I)*(za(I)+ZNTstoch_wat(I))/za(I)  ! (z+z0/L
          zol10 = zol(I)*(10.+ZNTstoch_wat(I))/za(I)    ! (10+z0)/L
          zol2  = zol(I)*(2.+ZNTstoch_wat(I))/za(I)     ! (2+z0)/L

          !COMPUTE PSIM and PSIH
          !CALL PSI_Suselj_Sood_2010(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Beljaars_Holtslag_1991(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),ZT_wat(I),ZNTstoch_wat(I),ZA(I))
          !CALL PSI_CB2005(PSIM(I),PSIH(I),zolza,zolz0)
          ! or use tables
          psim(I)=psim_stable(zolza,psi_opt)-psim_stable(zolz0,psi_opt)
          psih(I)=psih_stable(zolza,psi_opt)-psih_stable(zolzt,psi_opt)
          psim10(I)=psim_stable(zol10,psi_opt)-psim_stable(zolz0,psi_opt)
          psih10(I)=psih_stable(zol10,psi_opt)-psih_stable(zolz0,psi_opt)
          psih2(I)=psih_stable(zol2,psi_opt)-psih_stable(zolzt,psi_opt)

          ! 1.0 over Monin-Obukhov length
          RMOL(I)= ZOL(I)/ZA(I)

       ELSEIF(rb_wat(I) .EQ. 0.) THEN
          !=========================================================
          !-----CLASS 3; FORCED CONVECTION/NEUTRAL:
          !=========================================================

          PSIM(I)=0.0
          PSIH(I)=PSIM(I)
          PSIM10(I)=0.
          PSIH10(I)=0.
          PSIH2(I)=0.

          ZOL(I)  =0.
          RMOL(I) =0.

       ELSEIF(rb_wat(I) .LT. 0.)THEN
          !==========================================================
          !-----CLASS 4; FREE CONVECTION:
          !==========================================================

          !COMPUTE z/L first guess:
          IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
             CALL Li_etal_2010(ZOL(I),rb_wat(I),ZA(I)/ZNTstoch_wat(I),zratio_wat(I))
             !ZOL(I)=ZA(I)*KARMAN*grav*MOL(I)/(TH1D(I)*MAX(UST_wat(I)*UST_wat(I),0.001))
             ZOL(I)=MAX(ZOL(I),-20.0_kind_phys)
             ZOL(I)=MIN(ZOL(I),0.0_kind_phys)

             IF (debug_code >= 1) THEN
               IF (ZNTstoch_wat(i) < 1E-8 .OR. Zt_wat(i) < 1E-10) THEN
                 write(0,*)"===(wet) capture bad input in mynn sfc layer, i=:",i
                 write(0,*)"rb=", rb_wat(I)," ZNT=", ZNTstoch_wat(i)," ZT=",Zt_wat(i)
                 write(0,*)" tsk=", tskin_wat(i)," wstar=",wstar(i)," prev z/L=",ZOL(I),&
                 " tsurf=", tsurf_wat(i)," qsfc=", qsfc_wat(i)," znt=", znt_wat(i),&
                 " ust=", ust_wat(i)," snowh=", snowh_wat(i),"psfcpa=",PSFCPA(i),  &
                 " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
               ENDIF
             ENDIF

             !Use Pedros iterative function to find z/L
             !zol(I)=zolri(rb_wat(I),ZA(I),ZNTstoch_wat(I),ZT_wat(I),ZOL(I),psi_opt)
             !Use brute-force method
             zol(I)=zolrib(rb_wat(I),ZA(I),ZNTstoch_wat(I),zt_wat(I),GZ1OZ0_wat(I),GZ1OZt_wat(I),ZOL(I),psi_opt)
          ENDIF ! restart
          ZOL(I)=MAX(ZOL(I),-20.0_kind_phys)
          ZOL(I)=MIN(ZOL(I),0.0_kind_phys)

          zolzt = zol(I)*zt_wat(I)/ZA(I)                 ! zt/L
          zolz0 = zol(I)*ZNTstoch_wat(I)/ZA(I)           ! z0/L
          zolza = zol(I)*(za(I)+ZNTstoch_wat(I))/za(I)   ! (z+z0/L
          zol10 = zol(I)*(10.+ZNTstoch_wat(I))/za(I)     ! (10+z0)/L
          zol2  = zol(I)*(2.+ZNTstoch_wat(I))/za(I)      ! (2+z0)/L

          !COMPUTE PSIM and PSIH
          !CALL PSI_Suselj_Sood_2010(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Hogstrom_1996(PSIM(I),PSIH(I),ZOL(I), ZT_wat(I), ZNTstoch_wat(I), ZA(I))
          !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),ZT_wat(I),ZNTstoch_wat(I),ZA(I))
          ! use tables
          psim(I)=psim_unstable(zolza,psi_opt)-psim_unstable(zolz0,psi_opt)
          psih(I)=psih_unstable(zolza,psi_opt)-psih_unstable(zolzt,psi_opt)
          psim10(I)=psim_unstable(zol10,psi_opt)-psim_unstable(zolz0,psi_opt)
          psih10(I)=psih_unstable(zol10,psi_opt)-psih_unstable(zolz0,psi_opt)
          psih2(I)=psih_unstable(zol2,psi_opt)-psih_unstable(zolzt,psi_opt)

          !---LIMIT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND
          !---HIGH ROUGHNESS.  THIS PREVENTS DENOMINATOR IN FLUXES
          !---FROM GETTING TOO SMALL
          PSIH(I)=MIN(PSIH(I),0.9*GZ1OZt_wat(I))
          PSIM(I)=MIN(PSIM(I),0.9*GZ1OZ0_wat(I))
          PSIH2(I)=MIN(PSIH2(I),0.9*GZ2OZt_wat(I))
          PSIM10(I)=MIN(PSIM10(I),0.9*GZ10OZ0_wat(I))
          PSIH10(I)=MIN(PSIH10(I),0.9*GZ10OZt_wat(I))

          RMOL(I) = ZOL(I)/ZA(I)

       ENDIF

       ! CALCULATE THE RESISTANCE:
       PSIX_wat(I)  =MAX(GZ1OZ0_wat(I)-PSIM(I)   , 1.0_kind_phys)    ! = fm
       PSIX10_wat(I)=MAX(GZ10OZ0_wat(I)-PSIM10(I), 1.0_kind_phys)    ! = fm10
       PSIT_wat(I)  =MAX(GZ1OZt_wat(I)-PSIH(I)   , 1.0_kind_phys)    ! = fh
       PSIT2_wat(I) =MAX(GZ2OZt_wat(I)-PSIH2(I)  , 1.0_kind_phys)    ! = fh2
       PSIQ_wat(I)  =MAX(LOG((ZA(I)+ZQ_wat(i))/ZQ_wat(I))-PSIH(I) ,1.0_kind_phys)
       PSIQ2_wat(I) =MAX(LOG((2.0+ZQ_wat(i))/ZQ_wat(I))-PSIH2(I) ,1.0_kind_phys)

    ENDIF ! end water points

    IF (dry(i)) THEN
       IF (rb_lnd(I) .GT. 0.0) THEN

          IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
             !COMPUTE z/L first guess:
             CALL Li_etal_2010(ZOL(I),rb_lnd(I),ZA(I)/ZNTstoch_lnd(I),zratio_lnd(I))
             !ZOL(I)=ZA(I)*KARMAN*grav*MOL(I)/(TH1D(I)*MAX(UST_lnd(I)*UST_lnd(I),0.0001))
             ZOL(I)=MAX(ZOL(I),0.0_kind_phys)
             ZOL(I)=MIN(ZOL(I),20._kind_phys)

             IF (debug_code >= 1) THEN
               IF (ZNTstoch_lnd(i) < 1E-8 .OR. Zt_lnd(i) < 1E-10) THEN
                 write(0,*)"===(land) capture bad input in mynn sfc layer, i=:",i
                 write(0,*)"rb=", rb_lnd(I)," ZNT=", ZNTstoch_lnd(i)," ZT=",Zt_lnd(i)
                 write(0,*)" tsk=", tskin_lnd(i)," prev z/L=",ZOL(I),&
                 " tsurf=", tsurf_lnd(i)," qsfc=", qsfc_lnd(i)," znt=", znt_lnd(i),&
                 " ust=", ust_lnd(i)," snowh=", snowh_lnd(i),"psfcpa=",PSFCPA(i),  &
                 " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
               ENDIF
             ENDIF

             !Use Pedros iterative function to find z/L
             !zol(I)=zolri(rb_lnd(I),ZA(I),ZNTstoch_lnd(I),ZT_lnd(I),ZOL(I),psi_opt)
             !Use brute-force method
             zol(I)=zolrib(rb_lnd(I),ZA(I),ZNTstoch_lnd(I),zt_lnd(I),GZ1OZ0_lnd(I),GZ1OZt_lnd(I),ZOL(I),psi_opt)
          ENDIF ! restart
          ZOL(I)=MAX(ZOL(I),0.0_kind_phys)
          ZOL(I)=MIN(ZOL(I),20._kind_phys)

          zolzt = zol(I)*zt_lnd(I)/ZA(I)                ! zt/L
          zolz0 = zol(I)*ZNTstoch_lnd(I)/ZA(I)          ! z0/L
          zolza = zol(I)*(za(I)+ZNTstoch_lnd(I))/za(I)  ! (z+z0/L
          zol10 = zol(I)*(10.+ZNTstoch_lnd(I))/za(I)    ! (10+z0)/L
          zol2  = zol(I)*(2.+ZNTstoch_lnd(I))/za(I)     ! (2+z0)/L

          !COMPUTE PSIM and PSIH
          !CALL PSI_Beljaars_Holtslag_1991(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Zilitinkevich_Esau_2007(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),ZT_lnd(I),ZNTstoch_lnd(I),ZA(I))
          !CALL PSI_CB2005(PSIM(I),PSIH(I),zolza,zolz0)
          psim(I)=psim_stable(zolza,psi_opt)-psim_stable(zolz0,psi_opt)
          psih(I)=psih_stable(zolza,psi_opt)-psih_stable(zolzt,psi_opt)
          psim10(I)=psim_stable(zol10,psi_opt)-psim_stable(zolz0,psi_opt)
          psih10(I)=psih_stable(zol10,psi_opt)-psih_stable(zolz0,psi_opt)
          psih2(I)=psih_stable(zol2,psi_opt)-psih_stable(zolzt,psi_opt)

          ! 1.0 over Monin-Obukhov length
          RMOL(I)= ZOL(I)/ZA(I)

       ELSEIF(rb_lnd(I) .EQ. 0.) THEN
          !=========================================================
          !-----CLASS 3; FORCED CONVECTION/NEUTRAL:
          !=========================================================

          PSIM(I)=0.0
          PSIH(I)=PSIM(I)
          PSIM10(I)=0.
          PSIH10(I)=0.
          PSIH2(I)=0.

          ZOL(I)  =0.
          RMOL(I) =0.

       ELSEIF(rb_lnd(I) .LT. 0.)THEN
          !==========================================================
          !-----CLASS 4; FREE CONVECTION:
          !==========================================================

          IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
             !COMPUTE z/L first guess:
             CALL Li_etal_2010(ZOL(I),rb_lnd(I),ZA(I)/ZNTstoch_lnd(I),zratio_lnd(I))
             !ZOL(I)=ZA(I)*KARMAN*grav*MOL(I)/(TH1D(I)*MAX(UST_lnd(I)*UST_lnd(I),0.001))
             ZOL(I)=MAX(ZOL(I),-20.0_kind_phys)
             ZOL(I)=MIN(ZOL(I),0.0_kind_phys)

             IF (debug_code >= 1) THEN
               IF (ZNTstoch_lnd(i) < 1E-8 .OR. Zt_lnd(i) < 1E-10) THEN
                 write(0,*)"===(land) capture bad input in mynn sfc layer, i=:",i
                 write(0,*)"rb=", rb_lnd(I)," ZNT=", ZNTstoch_lnd(i)," ZT=",Zt_lnd(i)
                 write(0,*)" tsk=", tskin_lnd(i)," wstar=",wstar(i)," prev z/L=",ZOL(I),&
                 " tsurf=", tsurf_lnd(i)," qsfc=", qsfc_lnd(i)," znt=", znt_lnd(i),&
                 " ust=", ust_lnd(i)," snowh=", snowh_lnd(i),"psfcpa=",PSFCPA(i),  &
                 " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
               ENDIF
             ENDIF

             !Use Pedros iterative function to find z/L
             !zol(I)=zolri(rb_lnd(I),ZA(I),ZNTstoch_lnd(I),ZT_lnd(I),ZOL(I),psi_opt)
             !Use brute-force method
             zol(I)=zolrib(rb_lnd(I),ZA(I),ZNTstoch_lnd(I),zt_lnd(I),GZ1OZ0_lnd(I),GZ1OZt_lnd(I),ZOL(I),psi_opt)
          ENDIF ! restart
          ZOL(I)=MAX(ZOL(I),-20.0_kind_phys)
          ZOL(I)=MIN(ZOL(I),0.0_kind_phys)

          zolzt = zol(I)*zt_lnd(I)/ZA(I)                 ! zt/L
          zolz0 = zol(I)*ZNTstoch_lnd(I)/ZA(I)           ! z0/L
          zolza = zol(I)*(za(I)+ZNTstoch_lnd(I))/za(I)   ! (z+z0/L
          zol10 = zol(I)*(10.+ZNTstoch_lnd(I))/za(I)     ! (10+z0)/L
          zol2  = zol(I)*(2.+ZNTstoch_lnd(I))/za(I)      ! (2+z0)/L

          !COMPUTE PSIM and PSIH
          !CALL PSI_Hogstrom_1996(PSIM(I),PSIH(I),ZOL(I), ZT_lnd(I), ZNTstoch_lnd(I), ZA(I))
          !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),ZT_lnd(I),ZNTstoch_lnd(I),ZA(I))
          ! use tables
          psim(I)=psim_unstable(zolza,psi_opt)-psim_unstable(zolz0,psi_opt)
          psih(I)=psih_unstable(zolza,psi_opt)-psih_unstable(zolzt,psi_opt)
          psim10(I)=psim_unstable(zol10,psi_opt)-psim_unstable(zolz0,psi_opt)
          psih10(I)=psih_unstable(zol10,psi_opt)-psih_unstable(zolz0,psi_opt)
          psih2(I)=psih_unstable(zol2,psi_opt)-psih_unstable(zolzt,psi_opt)

          !---LIMIT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND
          !---HIGH ROUGHNESS.  THIS PREVENTS DENOMINATOR IN FLUXES
          !---FROM GETTING TOO SMALL
          PSIH(I)=MIN(PSIH(I),0.9*GZ1OZt_lnd(I))
          PSIM(I)=MIN(PSIM(I),0.9*GZ1OZ0_lnd(I))
          PSIH2(I)=MIN(PSIH2(I),0.9*GZ2OZt_lnd(I))
          PSIM10(I)=MIN(PSIM10(I),0.9*GZ10OZ0_lnd(I))
          PSIH10(I)=MIN(PSIH10(I),0.9*GZ10OZt_lnd(I))

          RMOL(I) = ZOL(I)/ZA(I)

       ENDIF

       ! CALCULATE THE RESISTANCE:
       PSIX_lnd(I)  =MAX(GZ1OZ0_lnd(I)-PSIM(I), 1.0_kind_phys)
       PSIX10_lnd(I)=MAX(GZ10OZ0_lnd(I)-PSIM10(I), 1.0_kind_phys)
       PSIT_lnd(I)  =MAX(GZ1OZt_lnd(I)-PSIH(I) , 1.0_kind_phys)
       PSIT2_lnd(I) =MAX(GZ2OZt_lnd(I)-PSIH2(I), 1.0_kind_phys)
       PSIQ_lnd(I)  =MAX(LOG((ZA(I)+ZQ_lnd(i))/ZQ_lnd(I))-PSIH(I) ,1.0_kind_phys)
       PSIQ2_lnd(I) =MAX(LOG((2.0+ZQ_lnd(i))/ZQ_lnd(I))-PSIH2(I) ,1.0_kind_phys)

    ENDIF ! end land points

    IF (icy(i)) THEN
       IF (rb_ice(I) .GT. 0.0) THEN

          IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
             !COMPUTE z/L first guess:
             CALL Li_etal_2010(ZOL(I),rb_ice(I),ZA(I)/ZNTstoch_ice(I),zratio_ice(I))
             !ZOL(I)=ZA(I)*KARMAN*grav*MOL(I)/(TH1D(I)*MAX(UST_ice(I)*UST_ice(I),0.0001))
             ZOL(I)=MAX(ZOL(I),0.0_kind_phys)
             ZOL(I)=MIN(ZOL(I),20._kind_phys)

             IF (debug_code >= 1) THEN
               IF (ZNTstoch_ice(i) < 1E-8 .OR. Zt_ice(i) < 1E-10) THEN
                 write(0,*)"===(ice) capture bad input in mynn sfc layer, i=:",i
                 write(0,*)"rb=", rb_ice(I)," ZNT=", ZNTstoch_ice(i)," ZT=",Zt_ice(i)
                 write(0,*)" tsk=", tskin_ice(i)," prev z/L=",ZOL(I),&
                 " tsurf=", tsurf_ice(i)," qsfc=", qsfc_ice(i)," znt=", znt_ice(i),&
                 " ust=", ust_ice(i)," snowh=", snowh_ice(i),"psfcpa=",PSFCPA(i),  &
                 " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
               ENDIF
             ENDIF

             !Use Pedros iterative function to find z/L
             !zol(I)=zolri(rb_ice(I),ZA(I),ZNTstoch_ice(I),ZT_ice(I),ZOL(I),psi_opt)
             !Use brute-force method
             zol(I)=zolrib(rb_ice(I),ZA(I),ZNTstoch_ice(I),zt_ice(I),GZ1OZ0_ice(I),GZ1OZt_ice(I),ZOL(I),psi_opt)
          ENDIF ! restart
          ZOL(I)=MAX(ZOL(I),0.0_kind_phys)
          ZOL(I)=MIN(ZOL(I),20._kind_phys)

          zolzt = zol(I)*zt_ice(I)/ZA(I)                ! zt/L
          zolz0 = zol(I)*ZNTstoch_ice(I)/ZA(I)          ! z0/L
          zolza = zol(I)*(za(I)+ZNTstoch_ice(I))/za(I)  ! (z+z0/L
          zol10 = zol(I)*(10.+ZNTstoch_ice(I))/za(I)    ! (10+z0)/L
          zol2  = zol(I)*(2.+ZNTstoch_ice(I))/za(I)     ! (2+z0)/L

          !COMPUTE PSIM and PSIH
          !CALL PSI_Beljaars_Holtslag_1991(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_Zilitinkevich_Esau_2007(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),ZT_ice(I),ZNTstoch_ice(I),ZA(I))
          !CALL PSI_CB2005(PSIM(I),PSIH(I),zolza,zolz0)
          psim(I)=psim_stable(zolza,psi_opt)-psim_stable(zolz0,psi_opt)
          psih(I)=psih_stable(zolza,psi_opt)-psih_stable(zolzt,psi_opt)
          psim10(I)=psim_stable(zol10,psi_opt)-psim_stable(zolz0,psi_opt)
          psih10(I)=psih_stable(zol10,psi_opt)-psih_stable(zolz0,psi_opt)
          psih2(I)=psih_stable(zol2,psi_opt)-psih_stable(zolzt,psi_opt)

          ! 1.0 over Monin-Obukhov length
          RMOL(I)= ZOL(I)/ZA(I)

       ELSEIF(rb_ice(I) .EQ. 0.) THEN
          !=========================================================
          !-----CLASS 3; FORCED CONVECTION/NEUTRAL:
          !=========================================================

          PSIM(I)=0.0
          PSIH(I)=PSIM(I)
          PSIM10(I)=0.
          PSIH10(I)=0.
          PSIH2(I)=0.

          ZOL(I)  =0.
          RMOL(I) =0.

       ELSEIF(rb_ice(I) .LT. 0.)THEN
          !==========================================================
          !-----CLASS 4; FREE CONVECTION:
          !==========================================================

          IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
             !COMPUTE z/L first guess:
             CALL Li_etal_2010(ZOL(I),rb_ice(I),ZA(I)/ZNTstoch_ice(I),zratio_ice(I))
             !ZOL(I)=ZA(I)*KARMAN*grav*MOL(I)/(TH1D(I)*MAX(UST_ice(I)*UST_ice(I),0.001))
             ZOL(I)=MAX(ZOL(I),-20.0_kind_phys)
             ZOL(I)=MIN(ZOL(I),0.0_kind_phys)

             IF (debug_code >= 1) THEN
               IF (ZNTstoch_ice(i) < 1E-8 .OR. Zt_ice(i) < 1E-10) THEN
                 write(0,*)"===(ice) capture bad input in mynn sfc layer, i=:",i
                 write(0,*)"rb=", rb_ice(I)," ZNT=", ZNTstoch_ice(i)," ZT=",Zt_ice(i)
                 write(0,*)" tsk=", tskin_ice(i)," wstar=",wstar(i)," prev z/L=",ZOL(I),&
                 " tsurf=", tsurf_ice(i)," qsfc=", qsfc_ice(i)," znt=", znt_ice(i),&
                 " ust=", ust_ice(i)," snowh=", snowh_ice(i),"psfcpa=",PSFCPA(i),  &
                 " dz=",dz8w1d(i)," qflx=",qflx(i)," hflx=",hflx(i)," hpbl=",pblh(i)
               ENDIF
             ENDIF

             !Use Pedros iterative function to find z/L
             !zol(I)=zolri(rb_ice(I),ZA(I),ZNTstoch_ice(I),ZT_ice(I),ZOL(I),psi_opt)
             !Use brute-force method
             zol(I)=zolrib(rb_ice(I),ZA(I),ZNTstoch_ice(I),zt_ice(I),GZ1OZ0_ice(I),GZ1OZt_ice(I),ZOL(I),psi_opt)
          ENDIF ! restart
          ZOL(I)=MAX(ZOL(I),-20.0_kind_phys)
          ZOL(I)=MIN(ZOL(I),0.0_kind_phys)

          zolzt = zol(I)*zt_ice(I)/ZA(I)                 ! zt/L
          zolz0 = zol(I)*ZNTstoch_ice(I)/ZA(I)           ! z0/L
          zolza = zol(I)*(za(I)+ZNTstoch_ice(I))/za(I)   ! (z+z0/L
          zol10 = zol(I)*(10.+ZNTstoch_ice(I))/za(I)     ! (10+z0)/L
          zol2  = zol(I)*(2.+ZNTstoch_ice(I))/za(I)      ! (2+z0)/L

          !COMPUTE PSIM and PSIH
          !CALL PSI_Hogstrom_1996(PSIM(I),PSIH(I),ZOL(I), ZT_ice(I), ZNTstoch_ice(I), ZA(I))
          !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
          !CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),ZT_ice(I),ZNTstoch_ice(I),ZA(I))
          ! use tables
          psim(I)=psim_unstable(zolza,psi_opt)-psim_unstable(zolz0,psi_opt)
          psih(I)=psih_unstable(zolza,psi_opt)-psih_unstable(zolzt,psi_opt)
          psim10(I)=psim_unstable(zol10,psi_opt)-psim_unstable(zolz0,psi_opt)
          psih10(I)=psih_unstable(zol10,psi_opt)-psih_unstable(zolz0,psi_opt)
          psih2(I)=psih_unstable(zol2,psi_opt)-psih_unstable(zolzt,psi_opt)

          !---LIMIT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND
          !---HIGH ROUGHNESS.  THIS PREVENTS DENOMINATOR IN FLUXES
          !---FROM GETTING TOO SMALL
          PSIH(I)=MIN(PSIH(I),0.9*GZ1OZt_ice(I))
          PSIM(I)=MIN(PSIM(I),0.9*GZ1OZ0_ice(I))
          PSIH2(I)=MIN(PSIH2(I),0.9*GZ2OZt_ice(I))
          PSIM10(I)=MIN(PSIM10(I),0.9*GZ10OZ0_ice(I))
          PSIH10(I)=MIN(PSIH10(I),0.9*GZ10OZt_ice(I))

          RMOL(I) = ZOL(I)/ZA(I)

       ENDIF

       ! CALCULATE THE RESISTANCE:
       PSIX_ice(I)  =MAX(GZ1OZ0_ice(I)-PSIM(I)   , 1.0_kind_phys)
       PSIX10_ice(I)=MAX(GZ10OZ0_ice(I)-PSIM10(I), 1.0_kind_phys)
       PSIT_ice(I)  =MAX(GZ1OZt_ice(I)-PSIH(I)   , 1.0_kind_phys)
       PSIT2_ice(I) =MAX(GZ2OZt_ice(I)-PSIH2(I)  , 1.0_kind_phys)
       PSIQ_ice(I)  =MAX(LOG((ZA(I)+ZQ_ice(i))/ZQ_ice(I))-PSIH(I) ,1.0_kind_phys)
       PSIQ2_ice(I) =MAX(LOG((2.0+ZQ_ice(i))/ZQ_ice(I))-PSIH2(I)  ,1.0_kind_phys)

    ENDIF ! end ice points

    !------------------------------------------------------------
    !-----COMPUTE THE FRICTIONAL VELOCITY:
    !------------------------------------------------------------

    IF (wet(I)) THEN
       ! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE
       OLDUST = UST_wat(I)
       UST_wat(I)=0.5*UST_wat(I)+0.5*KARMAN*WSPD(I)/PSIX_wat(I)
       !NON-AVERAGED:
       !UST_wat(I)=KARMAN*WSPD(I)/PSIX_wat(I)
       stress_wat(i)=ust_wat(i)**2

       ! Compute u* without vconv for use in HFX calc when isftcflx > 0
       WSPDI(I)=MAX(SQRT(U1D(I)*U1D(I)+V1D(I)*V1D(I)), wmin)
       !tgs - should USTM be separater for dry, icy, wet?
       USTM(I)=0.5*USTM(I)+0.5*KARMAN*WSPDI(I)/PSIX_wat(I)

       ! for possible future changes in sea-ice fraction from 0 to >0:
       if (.not. icy(i)) ust_ice(i)=ust_wat(i)
    ENDIF ! end water points

    IF (dry(I)) THEN
       ! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE
       OLDUST = UST_lnd(I)
       UST_lnd(I)=0.5*UST_lnd(I)+0.5*KARMAN*WSPD(I)/PSIX_lnd(I)
       !NON-AVERAGED:
       !UST_lnd(I)=KARMAN*WSPD(I)/PSIX_lnd(I)
       !From Tilden Meyers:
       !IF (rb_lnd(I) .GE. 0.0) THEN
       !   ust_lnd(i)=WSPD_lnd*0.1/(1.0 + 10.0*rb_lnd(I))
       !ELSE
       !   ust_lnd(i)=WSPD_lnd*0.1*(1.0 - 10.0*rb_lnd(I))**onethird
       !ENDIF
       UST_lnd(I)=MAX(UST_lnd(I),0.005_kind_phys)
       stress_lnd(i)=ust_lnd(i)**2

       !set ustm = ust over land.
       !tgs - should USTM be separater for dry, icy, wet?
       USTM(I)=UST_lnd(I)
    ENDIF ! end water points

    IF (icy(I)) THEN
       ! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE
       OLDUST = UST_ice(I)
       UST_ice(I)=0.5*UST_ice(I)+0.5*KARMAN*WSPD(I)/PSIX_ice(I)
       !NON-AVERAGED:
       !UST_ice(I)=KARMAN*WSPD(I)/PSIX_ice(I)
       UST_ice(I)=MAX(UST_ice(I),0.005_kind_phys)
       stress_ice(i)=ust_ice(i)**2

       !Set ustm = ust over ice.
       !tgs - should USTM be separate for for dry, icy, wet?
       USTM(I)=UST_ice(I)

       ! for possible future changes in sea-ice fraction from 1 to <1:
       !tgs - sea ice can be <1 now
       if (.not. wet(i)) ust_wat(i)=ust_ice(i)
    ENDIF ! end ice points

    !----------------------------------------------------
    !----COMPUTE THE TEMPERATURE SCALE (a.k.a. FRICTION TEMPERATURE, T*, or MOL)
    !----AND COMPUTE THE MOISTURE SCALE (or q*)
    !----------------------------------------------------

    !tgs - should we have MOL and qstar separate for dry, icy and wet?
    IF (wet(I)) THEN
       DTG=THV1D(I)-THVSK_wat(I)
       OLDTST=MOL(I)
       MOL(I)=KARMAN*DTG/PSIT_wat(I)/PRT
       !t_star(I) = -HFX(I)/(UST(I)*CPM(I)*RHO1D(I))
       !t_star(I) = MOL(I)
       !----------------------------------------------------
       DQG=(QVSH(i)-qsfc_wat(i))*1000.   !(kg/kg -> g/kg)
       qstar(I)=KARMAN*DQG/PSIQ_wat(I)/PRT
    ENDIF

    IF (dry(I)) THEN
       DTG=THV1D(I)-THVSK_lnd(I)
       OLDTST=MOL(I)
       MOL(I)=KARMAN*DTG/PSIT_lnd(I)/PRT
       !t_star(I) = -HFX(I)/(UST(I)*CPM(I)*RHO1D(I))
       !t_star(I) = MOL(I)
       !----------------------------------------------------
       DQG=(QVSH(i)-qsfc_lnd(i))*1000.   !(kg/kg -> g/kg)
       qstar(I)=KARMAN*DQG/PSIQ_lnd(I)/PRT
    ENDIF

    IF (icy(I)) THEN
       DTG=THV1D(I)-THVSK_ice(I)
       OLDTST=MOL(I)
       MOL(I)=KARMAN*DTG/PSIT_ice(I)/PRT
       !t_star(I) = -HFX(I)/(UST(I)*CPM(I)*RHO1D(I))
       !t_star(I) = MOL(I)
       !----------------------------------------------------
       DQG=(QVSH(i)-qsfc_ice(i))*1000.   !(kg/kg -> g/kg)
       qstar(I)=KARMAN*DQG/PSIQ_ice(I)/PRT
    ENDIF

   endif ! flag_iter
 ENDDO   ! end i-loop

#ifdef _OPENACC
! Necessary since OpenACC does not support branching in parallel code.
! Must sync host errflg, errmsg to determine if return must be triggered
! and correct error message and error flag code returned.
! This code is being executed on the HOST side only, pulling data from DEVICE.
!$acc exit data copyout(device_special_errflg, device_special_errmsg)
 IF (device_special_errflg /= 0) THEN
    errflg = device_special_errflg
    errmsg = device_special_errmsg
    return
 ENDIF
#endif

!$acc serial present(wet, dry, icy,                  &
!$acc        PSIM, PSIH, CPM, RHO1D, ZOL, wspd, MOL, &
!$acc        wstar, qstar, THV1D, HFX, MAVAIL, QVSH, &
!$acc         THVSK_wat,    THVSK_lnd,    THVSK_ice, &
!$acc           UST_wat,      UST_lnd,      UST_ice, &
!$acc      ZNTstoch_wat, ZNTstoch_lnd, ZNTstoch_ice, &
!$acc            zt_wat,       zt_lnd,       zt_ice)
 IF (debug_code == 2) THEN
    DO I=its,ite
       IF(wet(i))write(*,*)"==== AT END OF MAIN LOOP, i=",i, "(wet)"
       IF(dry(i))write(*,*)"==== AT END OF MAIN LOOP, i=",i, "(land)"
       IF(icy(i))write(*,*)"==== AT END OF MAIN LOOP, i=",i, "(ice)"
       write(*,*)"z/L:",ZOL(I)," wspd:",wspd(I)," Tstar:",MOL(I)
       IF(wet(i))write(*,*)"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
                           " DTHV:",THV1D(I)-THVSK_wat(I)
       IF(dry(i))write(*,*)"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
                           " DTHV:",THV1D(I)-THVSK_lnd(I)
       IF(icy(i))write(*,*)"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
                           " DTHV:",THV1D(I)-THVSK_ice(i)
       write(*,*)"CPM:",CPM(I)," RHO1D:",RHO1D(I)," q*:",qstar(I)," T*:",MOL(I)
       IF(wet(i))write(*,*)"U*:",UST_wat(I)," Z0:",ZNTstoch_wat(I)," Zt:",zt_wat(I)
       IF(dry(i))write(*,*)"U*:",UST_lnd(I)," Z0:",ZNTstoch_lnd(I)," Zt:",zt_lnd(I)
       IF(icy(i))write(*,*)"U*:",UST_ice(I)," Z0:",ZNTstoch_ice(I)," Zt:",zt_ice(I)
       write(*,*)"hfx:",HFX(I)," MAVAIL:",MAVAIL(I)," QVSH(I):",QVSH(I)
       write(*,*)"============================================="
    ENDDO ! end i-loop
 ENDIF
!$acc end serial

   !----------------------------------------------------------
   !  COMPUTE SURFACE HEAT AND MOISTURE FLUXES
   !----------------------------------------------------------
!$acc parallel loop present(flag_iter, dry, wet, icy, &
!$acc                    QFX, HFX, FLHC, FLQC, LH, CHS, CH, CHS2, CQS2, &
!$acc                    RHO1D, MAVAIL, USTM, &
!$acc                    UST_lnd, UST_wat, UST_ice, &
!$acc                    PSIQ_lnd, PSIT_lnd, PSIX_lnd, &
!$acc                    PSIQ_wat, PSIT_wat, PSIX_wat, &
!$acc                    PSIQ_ice, PSIT_ice, PSIX_ice, &
!$acc                    PSIQ2_lnd, PSIT2_lnd, &
!$acc                    PSIQ2_wat, PSIT2_wat, &
!$acc                    PSIQ2_ice, PSIT2_ice, &
!$acc                    QSFC, QSFC_lnd, QSFC_wat, QSFC_ice, &
!$acc                    QFLX, QFLX_lnd, QFLX_wat, QFLX_ice, &
!$acc                    HFLX, HFLX_lnd, HFLX_wat, HFLX_ice, &
!$acc                    QSFCMR_lnd, QSFCMR_wat, QSFCMR_ice, &
!$acc                    QV1D, WSPD, WSPDI, CPM, TH1D, &
!$acc                    THSK_lnd, THSK_wat, THSK_ice, &
!$acc                    ch_lnd, ch_wat, ch_ice, &
!$acc                    cm_lnd, cm_wat, cm_ice)
 DO I=its,ite
  if( flag_iter(i) ) then

    IF (ISFFLX .LT. 1) THEN

       QFX(i)  = 0.
       HFX(i)  = 0.
       HFLX(i) = 0.
       FLHC(I) = 0.
       FLQC(I) = 0.
       LH(I)   = 0.
       CHS(I)  = 0.
       CH(I)   = 0.
       CHS2(i) = 0.
       CQS2(i) = 0.
       ch_wat(I)= 0.
       cm_wat(I)= 0.
       ch_lnd(I)= 0.
       cm_lnd(I)= 0.
       ch_ice(I)= 0.
       cm_ice(I)= 0.

    ELSE

      IF (dry(i)) THEN

         !------------------------------------------
         ! CALCULATE THE EXCHANGE COEFFICIENTS FOR HEAT (FLHC)
         ! AND MOISTURE (FLQC)
         !------------------------------------------
         !tgs - should FLQC, FLHC be separate for dry, icy and wet?
         FLQC(I)=RHO1D(I)*MAVAIL(I)*UST_lnd(I)*KARMAN/PSIQ_lnd(i)
         FLHC(I)=RHO1D(I)*CPM(I)*UST_lnd(I)*KARMAN/PSIT_lnd(I)

         IF (compute_flux) THEN
            !----------------------------------
            ! COMPUTE SURFACE MOISTURE FLUX:
            !----------------------------------
            QFX(I)=FLQC(I)*(QSFCMR_lnd(I)-QV1D(I))
            !QFX(I)=FLQC(I)*(QSFC_lnd(I)-QV1D(I))
            QFX(I)=MAX(QFX(I),-0.02_kind_phys)      !allows small neg QFX
            LH(i)=XLV*QFX(i)
            ! BWG, 2020-06-17: Mod next 2 lines for fractional
            QFLX_lnd(i)=QFX(i)/RHO1D(i)
            QFLX(i)=QFLX_lnd(i)

            !----------------------------------
            ! COMPUTE SURFACE HEAT FLUX:
            !----------------------------------
            !HFX(I)=FLHC(I)*(THSK_lnd(I)-TH1D(I))
            HFX(I)=RHO1D(I)*CPM(I)*KARMAN*WSPD(i)/PSIX_lnd(I)*KARMAN/PSIT_lnd(I)*(THSK_lnd(I)-TH1D(i))
            HFX(I)=MAX(HFX(I),-250._kind_phys)
            ! BWG, 2020-06-17: Mod next 2 lines for fractional
            HFLX_lnd(I)=HFX(I)/(RHO1D(I)*cpm(I))
            HFLX(I)=HFLX_lnd(I)
         ENDIF

         !TRANSFER COEFF FOR SOME LSMs:
         !CHS(I)=UST(I)*KARMAN/(ALOG(KARMAN*UST(I)*ZA(I) &
         !       /XKA+ZA(I)/ZL)-PSIH(I))

         !tgs - should QSFC, CHS, CHS2 and CQS2 be separate for dry, icy and wet?
         CHS(I)=UST_lnd(I)*KARMAN/PSIT_lnd(I)

         !THESE ARE USED FOR 2-M DIAGNOSTICS ONLY
         CQS2(I)=UST_lnd(I)*KARMAN/PSIQ2_lnd(i)
         CHS2(I)=UST_lnd(I)*KARMAN/PSIT2_lnd(I)

         QSFC(I)=QSFC_lnd(I)

      ELSEIF (wet(i)) THEN

         !------------------------------------------
         ! CALCULATE THE EXCHANGE COEFFICIENTS FOR HEAT (FLHC)
         ! AND MOISTURE (FLQC)
         !------------------------------------------
         FLQC(I)=RHO1D(I)*MAVAIL(I)*UST_wat(I)*KARMAN/PSIQ_wat(i)
         FLHC(I)=RHO1D(I)*CPM(I)*UST_wat(I)*KARMAN/PSIT_wat(I)

         IF (compute_flux) THEN
            !----------------------------------
            ! COMPUTE SURFACE MOISTURE FLUX:
            !----------------------------------
            QFX(I)=FLQC(I)*(QSFCMR_wat(I)-QV1D(I))
            !QFX(I)=FLQC(I)*(QSFC_wat(I)-QV1D(I))
            QFX(I)=MAX(QFX(I),-0.02_kind_phys)      !allows small neg QFX
            LH(I)=XLV*QFX(I)
            ! BWG, 2020-06-17: Mod next 2 lines for fractional
            QFLX_wat(i)=QFX(i)/RHO1D(i)
            QFLX(i)=QFLX_wat(i)

            !----------------------------------
            ! COMPUTE SURFACE HEAT FLUX:
            !----------------------------------
            !HFX(I)=FLHC(I)*(THSK_wat(I)-TH1D(I))
            HFX(I)=RHO1D(I)*CPM(I)*KARMAN*WSPD(i)/PSIX_wat(I)*KARMAN/PSIT_wat(I)*(THSK_wat(I)-TH1D(i))
            IF ( PRESENT(ISFTCFLX) ) THEN
               IF ( ISFTCFLX.NE.0 ) THEN
                  ! AHW: add dissipative heating term
                  HFX(I)=HFX(I)+RHO1D(I)*USTM(I)*USTM(I)*WSPDI(I)
               ENDIF
            ENDIF
            ! BWG, 2020-06-17: Mod next 2 lines for fractional
            HFLX_wat(I)=HFX(I)/(RHO1D(I)*cpm(I))
            HFLX(I)=HFLX_wat(I)
         ENDIF

         !TRANSFER COEFF FOR SOME LSMs:
         !CHS(I)=UST(I)*KARMAN/(ALOG(KARMAN*UST(I)*ZA(I) &
         !       /XKA+ZA(I)/ZL)-PSIH(I))
         CHS(I)=UST_wat(I)*KARMAN/PSIT_wat(I)

         !THESE ARE USED FOR 2-M DIAGNOSTICS ONLY
         CQS2(I)=UST_wat(I)*KARMAN/PSIQ2_wat(i)
         CHS2(I)=UST_wat(I)*KARMAN/PSIT2_wat(I)

         QSFC(I)=QSFC_wat(I)

      ELSEIF (icy(i)) THEN

         !------------------------------------------
         ! CALCULATE THE EXCHANGE COEFFICIENTS FOR HEAT (FLHC)
         ! AND MOISTURE (FLQC)
         !------------------------------------------
         FLQC(I)=RHO1D(I)*MAVAIL(I)*UST_ice(I)*KARMAN/PSIQ_ice(i)
         FLHC(I)=RHO1D(I)*CPM(I)*UST_ice(I)*KARMAN/PSIT_ice(I)

         IF (compute_flux) THEN
            !----------------------------------
            ! COMPUTE SURFACE MOISTURE FLUX:
            !----------------------------------
            QFX(I)=FLQC(I)*(QSFCMR_ice(I)-QV1D(I))
            !QFX(I)=FLQC(I)*(QSFC_ice(I)-QV1D(I))
            QFX(I)=MAX(QFX(I),-0.02_kind_phys)      !allows small neg QFX
            LH(I)=XLF*QFX(I)
            ! BWG, 2020-06-17: Mod next 2 lines for fractional
            QFLX_ice(i)=QFX(i)/RHO1D(i)
            QFLX(i)=QFLX_ice(i)

            !----------------------------------
            ! COMPUTE SURFACE HEAT FLUX:
            !----------------------------------
            !HFX(I)=FLHC(I)*(THSK_ice(I)-TH1D(I))
            HFX(I)=RHO1D(I)*CPM(I)*KARMAN*WSPD(i)/PSIX_ice(I)*KARMAN/PSIT_ice(I)*(THSK_ice(I)-TH1D(i))
            HFX(I)=MAX(HFX(I),-250._kind_phys)
            ! BWG, 2020-06-17: Mod next 2 lines for fractional
            HFLX_ice(I)=HFX(I)/(RHO1D(I)*cpm(I))
            HFLX(I)=HFLX_ice(I)
         ENDIF

         !TRANSFER COEFF FOR SOME LSMs:
         !CHS(I)=UST(I)*KARMAN/(ALOG(KARMAN*UST(I)*ZA(I) &
         !       /XKA+ZA(I)/ZL)-PSIH(I))
         CHS(I)=UST_ice(I)*KARMAN/PSIT_ice(I)

         !THESE ARE USED FOR 2-M DIAGNOSTICS ONLY
         CQS2(I)=UST_ice(I)*KARMAN/PSIQ2_ice(i)
         CHS2(I)=UST_ice(I)*KARMAN/PSIT2_ice(I)

         QSFC(I)=QSFC_ice(I)

      ENDIF

      IF (debug_code > 1) THEN
         write(*,*)"QFX=",QFX(I),"FLQC=",FLQC(I)
         if(icy(i))write(*,*)"ice, MAVAIL:",MAVAIL(I)," u*=",UST_ice(I)," psiq=",PSIQ_ice(i)
         if(dry(i))write(*,*)"lnd, MAVAIL:",MAVAIL(I)," u*=",UST_lnd(I)," psiq=",PSIQ_lnd(i)
         if(wet(i))write(*,*)"ocn, MAVAIL:",MAVAIL(I)," u*=",UST_wat(I)," psiq=",PSIQ_wat(i)
      ENDIF

      ! The exchange coefficient for cloud water is assumed to be the
      ! same as that for heat. CH is multiplied by WSPD.
      ch(i)=flhc(i)/( cpm(i)*RHO1D(i) )

      !-----------------------------------------
      !--- COMPUTE EXCHANGE COEFFICIENTS FOR FV3
      !-----------------------------------------
      IF (wet(i)) THEN
         ch_wat(I)=(karman/psix_wat(I))*(karman/psit_wat(i))
         cm_wat(I)=(karman/psix_wat(I))*(karman/psix_wat(I))
      ENDIF
      IF (dry(i)) THEN
         ch_lnd(I)=(karman/psix_lnd(I))*(karman/psit_lnd(i))
         cm_lnd(I)=(karman/psix_lnd(I))*(karman/psix_lnd(I))
      ENDIF
      IF (icy(i)) THEN
         ch_ice(I)=(karman/psix_ice(I))*(karman/psit_ice(i))
         cm_ice(I)=(karman/psix_ice(I))*(karman/psix_ice(I))
      ENDIF

   ENDIF !end ISFFLX option
  endif ! flag_iter
ENDDO ! end i-loop

IF (compute_diag) then
   !$acc parallel loop present(flag_iter, dry, wet, icy, &
   !$acc                    ZA, ZA2, T2, TH2, TH1D, Q2, QV1D, PSFCPA, &
   !$acc                    THSK_lnd, THSK_wat, THSK_ice, &
   !$acc                    QSFC_lnd, QSFC_wat, QSFC_ice, &
   !$acc                    U10, V10, U1D, V1D, U1D2, V1D2, &
   !$acc                    ZNTstoch_lnd, ZNTstoch_lnd, ZNTstoch_ice, &
   !$acc                    PSIX_lnd, PSIX_wat, PSIX_ice, &
   !$acc                    PSIX10_lnd, PSIX10_wat, PSIX10_ice, &
   !$acc                    PSIT2_lnd, PSIT2_wat, PSIT2_ice, &
   !$acc                    PSIT_lnd, PSIT_wat, PSIT_ice, &
   !$acc                    PSIQ2_lnd, PSIQ2_wat, PSIQ2_ice, &
   !$acc                    PSIQ_lnd, PSIQ_wat, PSIQ_ice)
   DO I=its,ite
     if( flag_iter(i) ) then
      !-----------------------------------------------------
      !COMPUTE DIAGNOSTICS
      !-----------------------------------------------------
      !COMPUTE 10 M WNDS
      !-----------------------------------------------------
      ! If the lowest model level is close to 10-m, use it
      ! instead of the flux-based diagnostic formula.
      if (ZA(i) .le. 7.0) then
         ! high vertical resolution
         if(ZA2(i) .gt. 7.0 .and. ZA2(i) .lt. 13.0) then
            !use 2nd model level
            U10(I)=U1D2(I)
            V10(I)=V1D2(I)
         else
            IF (dry(i)) THEN
              !U10(I)=U1D(I)*PSIX10_lnd(I)/PSIX_lnd(I)
              !V10(I)=V1D(I)*PSIX10_lnd(I)/PSIX_lnd(I)
              !use neutral-log:
              U10(I)=U1D(I)*log(10./ZNTstoch_lnd(I))/log(ZA(I)/ZNTstoch_lnd(I))
              V10(I)=V1D(I)*log(10./ZNTstoch_lnd(I))/log(ZA(I)/ZNTstoch_lnd(I))
            ELSEIF (wet(i)) THEN
              U10(I)=U1D(I)*log(10./ZNTstoch_wat(I))/log(ZA(I)/ZNTstoch_wat(I))
              V10(I)=V1D(I)*log(10./ZNTstoch_wat(I))/log(ZA(I)/ZNTstoch_wat(I))
            ELSEIF (icy(i)) THEN
              U10(I)=U1D(I)*log(10./ZNTstoch_ice(I))/log(ZA(I)/ZNTstoch_ice(I))
              V10(I)=V1D(I)*log(10./ZNTstoch_ice(I))/log(ZA(I)/ZNTstoch_ice(I))
            ENDIF
         endif
      elseif (ZA(i) .gt. 7.0 .and. ZA(i) .lt. 13.0) then
         !moderate vertical resolution
         IF (dry(i)) THEN
            !U10(I)=U1D(I)*PSIX10_lnd(I)/PSIX_lnd(I)
            !V10(I)=V1D(I)*PSIX10_lnd(I)/PSIX_lnd(I)
            !use neutral-log:
            U10(I)=U1D(I)*log(10./ZNTstoch_lnd(I))/log(ZA(I)/ZNTstoch_lnd(I))
            V10(I)=V1D(I)*log(10./ZNTstoch_lnd(I))/log(ZA(I)/ZNTstoch_lnd(I))
         ELSEIF (wet(i)) THEN
            U10(I)=U1D(I)*log(10./ZNTstoch_wat(I))/log(ZA(I)/ZNTstoch_wat(I))
            V10(I)=V1D(I)*log(10./ZNTstoch_wat(I))/log(ZA(I)/ZNTstoch_wat(I))
         ELSEIF (icy(i)) THEN
            U10(I)=U1D(I)*log(10./ZNTstoch_ice(I))/log(ZA(I)/ZNTstoch_ice(I))
            V10(I)=V1D(I)*log(10./ZNTstoch_ice(I))/log(ZA(I)/ZNTstoch_ice(I))
         ENDIF
      else
         ! very coarse vertical resolution
         IF (dry(i)) THEN
            U10(I)=U1D(I)*PSIX10_lnd(I)/PSIX_lnd(I)
            V10(I)=V1D(I)*PSIX10_lnd(I)/PSIX_lnd(I)
         ELSEIF (wet(i)) THEN
            U10(I)=U1D(I)*PSIX10_wat(I)/PSIX_wat(I)
            V10(I)=V1D(I)*PSIX10_wat(I)/PSIX_wat(I)
         ELSEIF (icy(i)) THEN
            U10(I)=U1D(I)*PSIX10_ice(I)/PSIX_ice(I)
            V10(I)=V1D(I)*PSIX10_ice(I)/PSIX_ice(I)
         ENDIF
      endif

      !-----------------------------------------------------
      !COMPUTE 2m T, TH, AND Q
      !THESE WILL BE OVERWRITTEN FOR LAND POINTS IN THE LSM
      !-----------------------------------------------------
      IF (dry(i)) THEN
         DTG=TH1D(I)-THSK_lnd(I)
         TH2(I)=THSK_lnd(I)+DTG*PSIT2_lnd(I)/PSIT_lnd(I)
         !***  BE CERTAIN THAT THE 2-M THETA IS BRACKETED BY
         !***  THE VALUES AT THE SURFACE AND LOWEST MODEL LEVEL.
         IF ((TH1D(I)>THSK_lnd(I) .AND. (TH2(I)<THSK_lnd(I) .OR. TH2(I)>TH1D(I))) .OR. &
            (TH1D(I)<THSK_lnd(I) .AND. (TH2(I)>THSK_lnd(I) .OR. TH2(I)<TH1D(I)))) THEN
            TH2(I)=THSK_lnd(I) + 2.*(TH1D(I)-THSK_lnd(I))/ZA(I)
         ENDIF
         T2(I)=TH2(I)*(PSFCPA(I)/100000.)**ROVCP

         Q2(I)=QSFC_lnd(I)+(QV1D(I)-QSFC_lnd(I))*PSIQ2_lnd(i)/PSIQ_lnd(i)
         Q2(I)= MAX(Q2(I), MIN(QSFC_lnd(I), QV1D(I)))
         Q2(I)= MIN(Q2(I), 1.05*QV1D(I))
      ELSEIF (wet(i)) THEN
         DTG=TH1D(I)-THSK_wat(I)
         TH2(I)=THSK_wat(I)+DTG*PSIT2_wat(I)/PSIT_wat(I)
         !***  BE CERTAIN THAT THE 2-M THETA IS BRACKETED BY
         !***  THE VALUES AT THE SURFACE AND LOWEST MODEL LEVEL.
         IF ((TH1D(I)>THSK_wat(I) .AND. (TH2(I)<THSK_wat(I) .OR. TH2(I)>TH1D(I))) .OR. &
            (TH1D(I)<THSK_wat(I) .AND. (TH2(I)>THSK_wat(I) .OR. TH2(I)<TH1D(I)))) THEN
            TH2(I)=THSK_wat(I) + 2.*(TH1D(I)-THSK_wat(I))/ZA(I)
         ENDIF
         T2(I)=TH2(I)*(PSFCPA(I)/100000.)**ROVCP

         Q2(I)=QSFC_wat(I)+(QV1D(I)-QSFC_wat(I))*PSIQ2_wat(i)/PSIQ_wat(i)
         Q2(I)= MAX(Q2(I), MIN(QSFC_wat(I), QV1D(I)))
         Q2(I)= MIN(Q2(I), 1.05*QV1D(I))
      ELSEIF (icy(i)) THEN
         DTG=TH1D(I)-THSK_ice(I)
         TH2(I)=THSK_ice(I)+DTG*PSIT2_ice(I)/PSIT_ice(I)
         !***  BE CERTAIN THAT THE 2-M THETA IS BRACKETED BY
         !***  THE VALUES AT THE SURFACE AND LOWEST MODEL LEVEL.
         IF ((TH1D(I)>THSK_ice(I) .AND. (TH2(I)<THSK_ice(I) .OR. TH2(I)>TH1D(I))) .OR. &
            (TH1D(I)<THSK_ice(I) .AND. (TH2(I)>THSK_ice(I) .OR. TH2(I)<TH1D(I)))) THEN
            TH2(I)=THSK_ice(I) + 2.*(TH1D(I)-THSK_ice(I))/ZA(I)
         ENDIF
         T2(I)=TH2(I)*(PSFCPA(I)/100000.)**ROVCP

         Q2(I)=QSFC_ice(I)+(QV1D(I)-QSFC_ice(I))*PSIQ2_ice(i)/PSIQ_ice(i)
         Q2(I)= MAX(Q2(I), MIN(QSFC_ice(I), QV1D(I)))
         Q2(I)= MIN(Q2(I), 1.05*QV1D(I))
      ENDIF
     endif ! flag_iter
   ENDDO
ENDIF ! end compute_diag

!-----------------------------------------------------
! DEBUG - SUSPICIOUS VALUES
!-----------------------------------------------------
!$acc serial present(dry, wet, icy, CPM, MAVAIL, &
!$acc       HFX, LH, wstar, RHO1D, PBLH, ZOL, ZA, MOL, &
!$acc       PSIM, PSIH, WSTAR, T1D, TH1D, THV1D, QVSH, &
!$acc       UST_wat, UST_lnd, UST_ice, &
!$acc       THSK_wat, THSK_lnd, THSK_ice, &
!$acc       THVSK_wat, THVSK_lnd, THVSK_ice, &
!$acc       ZNTstoch_wat, ZNTstoch_lnd, ZNTstoch_ice, &
!$acc       ZT_wat, ZT_lnd, ZT_ice, &
!$acc       QSFC_wat, QSFC_lnd, QSFC_ice, &
!$acc       PSIX_wat, PSIX_lnd, PSIX_ice)
IF ( debug_code == 2) THEN
   DO I=its,ite
      yesno = 0
      IF (compute_flux) THEN
        IF (HFX(I) > 1200. .OR. HFX(I) < -700.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "HFX: ",HFX(I)
            yesno = 1
        ENDIF
        IF (LH(I)  > 1200. .OR. LH(I)  < -700.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "LH: ",LH(I)
            yesno = 1
        ENDIF
      ENDIF
      IF (wet(i)) THEN
         IF (UST_wat(I) < 0.0 .OR. UST_wat(I) > 4.0 )THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "UST_wat: ",UST_wat(I)
            yesno = 1
         ENDIF
      ENDIF
      IF (dry(i)) THEN
         IF (UST_lnd(I) < 0.0 .OR. UST_lnd(I) > 4.0 )THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "UST_lnd: ",UST_lnd(I)
            yesno = 1
         ENDIF
      ENDIF
      IF (icy(i)) THEN
         IF (UST_ice(I) < 0.0 .OR. UST_ice(I) > 4.0 )THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "UST_ice: ",UST_ice(I)
            yesno = 1
         ENDIF
      ENDIF
      IF (WSTAR(I)<0.0 .OR. WSTAR(I) > 6.0)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "WSTAR: ",WSTAR(I)
            yesno = 1
      ENDIF
      IF (RHO1D(I)<0.0 .OR. RHO1D(I) > 1.6 )THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "rho: ",RHO1D(I)
            yesno = 1
      ENDIF
      IF (dry(i)) THEN
         IF (QSFC_lnd(I)*1000. <0.0 .OR. QSFC_lnd(I)*1000. >40.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "QSFC_lnd: ",QSFC_lnd(I)
            yesno = 1
         ENDIF
      ENDIF
      IF (PBLH(I)<0. .OR. PBLH(I)>6000.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            I,J, "PBLH: ",PBLH(I)
            yesno = 1
      ENDIF

      IF (yesno == 1) THEN
        IF (wet(i)) THEN
           print*," OTHER INFO over water:"
           print*,"z/L:",ZOL(I)," U*:",UST_wat(I)," Tstar:",MOL(I)
           print*,"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
              " DTHV:",THV1D(I)-THVSK_wat(I)
           print*,"CPM:",CPM(I)," RHO1D:",RHO1D(I)," L:",&
              ZOL(I)/ZA(I)," DTH:",TH1D(I)-THSK_wat(I)
           print*," Z0:",ZNTstoch_wat(I)," Zt:",ZT_wat(I)," za:",za(I)
           print*,"MAVAIL:",MAVAIL(I)," QSFC_wat(I):",&
              QSFC_wat(I)," QVSH(I):",QVSH(I)
           print*,"PSIX=",PSIX_wat(I)," T1D(i):",T1D(i)
           write(*,*)"============================================="
        ENDIF
        IF (dry(i)) THEN
           print*," OTHER INFO over land:"
           print*,"z/L:",ZOL(I)," U*:",UST_lnd(I),&
              " Tstar:",MOL(I)
           print*,"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
              " DTHV:",THV1D(I)-THVSK_lnd(I)
           print*,"CPM:",CPM(I)," RHO1D:",RHO1D(I)," L:",&
              ZOL(I)/ZA(I)," DTH:",TH1D(I)-THSK_lnd(I)
           print*," Z0:",ZNTstoch_lnd(I)," Zt:",ZT_lnd(I)," za:",za(I)
           print*," MAVAIL:",MAVAIL(I)," QSFC_lnd(I):",&
              QSFC_lnd(I)," QVSH(I):",QVSH(I)
           print*,"PSIX=",PSIX_lnd(I)," T1D(i):",T1D(i)
           write(*,*)"============================================="
        ENDIF
        IF (icy(i)) THEN
           print*," OTHER INFO:"
           print*,"z/L:",ZOL(I)," U*:",UST_ice(I),&
              " Tstar:",MOL(I)
           print*,"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
              " DTHV:",THV1D(I)-THVSK_ice(I)
           print*,"CPM:",CPM(I)," RHO1D:",RHO1D(I)," L:",&
              ZOL(I)/ZA(I)," DTH:",TH1D(I)-THSK_ice(I)
           print*," Z0:",ZNTstoch_ice(I)," Zt:",ZT_ice(I)," za:",za(I)
           print*," MAVAIL:",MAVAIL(I)," QSFC_ice(I):",&
              QSFC_ice(I)," QVSH(I):",QVSH(I)
           print*,"PSIX=",PSIX_ice(I)," T1D(i):",T1D(i)
           write(*,*)"============================================="
        ENDIF
      ENDIF
    ENDDO ! end i-loop
 ENDIF ! end debug option
!$acc end serial

!$acc exit data copyout(CPM, FLHC, FLQC, CHS, CH, CHS2, CQS2,&
!$acc                   USTM, wstar, qstar, ZOL, MOL, RMOL,  &
!$acc                   HFX, QFX, LH, QSFC, QFLX, HFLX,      &
!$acc                   T2, TH2, Q2, WSPD, U10, V10,         &
!$acc                   QGH, psim, psih,                     &
!$acc                   stress_wat, stress_lnd, stress_ice,  &
!$acc                   rb_wat,     rb_lnd,     rb_ice,      &
!$acc                   UST_wat,    UST_lnd,    UST_ice,     &
!$acc                   ZNT_wat,    ZNT_lnd,    ZNT_ice,     &
!$acc                   QSFC_lnd,   QSFC_wat,   QSFC_ice,    &
!$acc                   QFLX_lnd,   QFLX_wat,   QFLX_ice,    &
!$acc                   HFLX_lnd,   HFLX_wat,   HFLX_ice,    &
!$acc                   PSIX_wat,   PSIX_lnd,   PSIX_ice,    &
!$acc                   PSIX10_wat, PSIX10_lnd, PSIX10_ice,  &
!$acc                   PSIT2_lnd,  PSIT2_wat,  PSIT2_ice,   &
!$acc                   PSIT_lnd,   PSIT_wat,   PSIT_ice,    &
!$acc                   ch_lnd,     ch_wat,     ch_ice,      &
!$acc                   cm_lnd,     cm_wat,     cm_ice,      &
!$acc                   device_errmsg, device_errflg)

! Final sync of device and host error flags and messages
IF (device_errflg /= 0) THEN
    errflg = device_errflg
    errmsg = device_errmsg
ENDIF

!$acc exit data delete( flag_iter, dry, wet, icy, dx,         &
!$acc                   MAVAIL, PBLH, PSFCPA, z0pert, ztpert, &
!$acc                   QV1D, U1D, V1D, U1D2, V1D2, T1D, P1D, &
!$acc                   rstoch1D, sigmaf, shdmax, vegtype,    &
!$acc                   dz2w1d, dz8w1d,                       &
!$acc                   snowh_wat, snowh_lnd, snowh_ice,      &
!$acc                   tskin_wat, tskin_lnd, tskin_ice,      &
!$acc                   tsurf_wat, tsurf_lnd, tsurf_ice)

!$acc exit data delete(  ZA,     ZA2,    THV1D,  TH1D,   TC1D,   TV1D,  &
!$acc                    RHO1D,  QVSH,   PSIH2,  PSIM10, PSIH10, WSPDI, &
!$acc                    GOVRTH, PSFC,   THCON,                         &
!$acc                    zratio_lnd,   zratio_ice,   zratio_wat,        &
!$acc                    TSK_lnd,      TSK_ice,      TSK_wat,           &
!$acc                    THSK_lnd,     THSK_ice,     THSK_wat,          &
!$acc                    THVSK_lnd,    THVSK_ice,    THVSK_wat,         &
!$acc                    GZ1OZ0_lnd,   GZ1OZ0_ice,   GZ1OZ0_wat,        &
!$acc                    GZ1OZt_lnd,   GZ1OZt_ice,   GZ1OZt_wat,        &
!$acc                    GZ2OZ0_lnd,   GZ2OZ0_ice,   GZ2OZ0_wat,        &
!$acc                    GZ2OZt_lnd,   GZ2OZt_ice,   GZ2OZt_wat,        &
!$acc                    GZ10OZ0_lnd,  GZ10OZ0_ice,  GZ10OZ0_wat,       &
!$acc                    GZ10OZt_lnd,  GZ10OZt_ice,  GZ10OZt_wat,       &
!$acc                    ZNTstoch_lnd, ZNTstoch_ice, ZNTstoch_wat,      &
!$acc                    ZT_lnd,       ZT_ice,       ZT_wat,            &
!$acc                    ZQ_lnd,       ZQ_ice,       ZQ_wat,            &
!$acc                    PSIQ_lnd,     PSIQ_ice,     PSIQ_wat,          &
!$acc                    PSIQ2_lnd,    PSIQ2_ice,    PSIQ2_wat,         &
!$acc                    QSFCMR_lnd,   QSFCMR_ice,   QSFCMR_wat )

END SUBROUTINE SFCLAY1D_mynn
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the thermal and moisture roughness lengths
!! from Zilitinkevich (1995) and Zilitinkevich et al. (2001) over
!! land and water, respectively.
!!
!! MODS:
!! 20120705 : added IZ0TLND option. Note: This option was designed
!!            to work with the Noah LSM and may be specific for that
!!            LSM only. Tests with RUC LSM showed no improvements.
  SUBROUTINE zilitinkevich_1995(Z_0,Zt,Zq,restar,ustar,KARMAN,&
        & landsea,IZ0TLND2,spp_sfc,rstoch)

       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)       :: Z_0,restar,ustar,KARMAN,landsea
       INTEGER,  OPTIONAL,   INTENT(IN)  :: IZ0TLND2
       REAL(kind_phys), INTENT(OUT)      :: Zt,Zq
       REAL(kind_phys) :: CZIL       !=0.100 in Chen et al. (1997)
                                     !=0.075 in Zilitinkevich (1995)
                                     !=0.500 in Lemone et al. (2008)
       INTEGER,  INTENT(IN)  ::    spp_sfc
       REAL(kind_phys), INTENT(IN)  :: rstoch


       IF (landsea-1.5 .GT. 0) THEN    !WATER

          !THIS IS BASED ON Zilitinkevich, Grachev, and Fairall (2001;
          !Their equations 15 and 16).
          IF (restar .LT. 0.1) THEN
             Zt = Z_0*EXP(KARMAN*2.0)
             Zt = MIN( Zt, 6.0e-5_kind_phys)
             Zt = MAX( Zt, 2.0e-9_kind_phys)
             Zq = Z_0*EXP(KARMAN*3.0)
             Zq = MIN( Zq, 6.0e-5_kind_phys)
             Zq = MAX( Zq, 2.0e-9_kind_phys)
          ELSE
             Zt = Z_0*EXP(-KARMAN*(4.0*SQRT(restar)-3.2))
             Zt = MIN( Zt, 6.0e-5_kind_phys)
             Zt = MAX( Zt, 2.0e-9_kind_phys)
             Zq = Z_0*EXP(-KARMAN*(4.0*SQRT(restar)-4.2))
             Zq = MIN( Zt, 6.0e-5_kind_phys)
             Zq = MAX( Zt, 2.0e-9_kind_phys)
          ENDIF

       ELSE                             !LAND

          !Option to modify CZIL according to Chen & Zhang, 2009
          IF ( IZ0TLND2 .EQ. 1 ) THEN
             CZIL = 10.0 ** ( -0.40 * ( Z_0 / 0.07 ) )
          ELSE
             CZIL = 0.095 !0.075 !0.10
          END IF

          Zt = Z_0*EXP(-KARMAN*CZIL*SQRT(restar))
          Zt = MIN( Zt, 0.75*Z_0)

          Zq = Z_0*EXP(-KARMAN*CZIL*SQRT(restar))
          Zq = MIN( Zq, 0.75*Z_0)

! stochastically perturb thermal and moisture roughness length.
! currently set to half the amplitude:
          if (spp_sfc==1) then
             Zt = Zt + Zt * 0.5 * rstoch
             Zt = MAX(Zt, 0.0001_kind_phys)
             Zq = Zt
          endif

       ENDIF

   END SUBROUTINE zilitinkevich_1995
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
   SUBROUTINE davis_etal_2008(Z_0,ustar)

    !a.k.a. : Donelan et al. (2004)
    !This formulation for roughness length was designed to match
    !the labratory experiments of Donelan et al. (2004).
    !This is an update version from Davis et al. 2008, which
    !corrects a small-bias in Z_0 (AHW real-time 2012).

       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: ustar
       REAL(kind_phys), INTENT(OUT) :: Z_0
       REAL(kind_phys) :: ZW, ZN1, ZN2
       REAL(kind_phys), PARAMETER :: OZO=1.59E-5

       !OLD FORM: Z_0 = 10.*EXP(-10./(ustar**onethird))
       !NEW FORM:

       ZW  = MIN((ustar/1.06)**(0.3),1.0_kind_phys)
       ZN1 = 0.011*ustar*ustar*g_inv + OZO
       ZN2 = 10.*exp(-9.5*ustar**(-onethird)) + &
             0.11*1.5E-5/MAX(ustar,0.01_kind_phys)
             !0.11*1.5E-5/AMAX1(ustar,0.01)
       Z_0 = (1.0-ZW) * ZN1 + ZW * ZN2

       Z_0 = MAX( Z_0, 1.27e-7_kind_phys)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3_kind_phys)  !Davis et al. (2008)

   END SUBROUTINE davis_etal_2008
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!>This formulation for roughness length was designed account for.
!!wave steepness.
   SUBROUTINE Taylor_Yelland_2001(Z_0,ustar,wsp10)
       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: ustar,wsp10
       REAL(kind_phys), INTENT(OUT) :: Z_0
       REAL(kind_phys), parameter   :: pi=3.14159265
       REAL(kind_phys) :: hs, Tp, Lp

       !hs is the significant wave height
        hs = 0.0248*(wsp10**2.)
       !Tp dominant wave period
        Tp = 0.729*MAX(wsp10,0.1_kind_phys)
       !Lp is the wavelength of the dominant wave
        Lp = grav*Tp**2/(2*pi)

       Z_0 = 1200.*hs*(hs/Lp)**4.5
       Z_0 = MAX( Z_0, 1.27e-7_kind_phys)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3_kind_phys)  !Davis et al. (2008)

   END SUBROUTINE Taylor_Yelland_2001
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!>This version of Charnock's relation employs a varying
!! Charnock parameter, similar to COARE3.0 [Fairall et al. (2003)].
!! The Charnock parameter CZC is varied from .011 to .018.
!! between 10-m wsp = 10 and 18..
   SUBROUTINE charnock_1955(Z_0,ustar,wsp10,visc,zu)
       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: ustar, visc, wsp10, zu
       REAL(kind_phys), INTENT(OUT) :: Z_0
       REAL(kind_phys), PARAMETER   :: CZO2=0.011
       REAL(kind_phys)              :: CZC    ! variable charnock "constant"
       REAL(kind_phys)              :: wsp10m ! logarithmically calculated 10 m

       wsp10m = wsp10*log(10./1e-4)/log(zu/1e-4)
       CZC = CZO2 + 0.007*MIN(MAX((wsp10m-10.)/8., 0._kind_phys), 1.0_kind_phys)

       Z_0 = CZC*ustar*ustar*g_inv + (0.11*visc/MAX(ustar,0.05_kind_phys))
       Z_0 = MAX( Z_0, 1.27e-7_kind_phys)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3_kind_phys)  !Davis et al. (2008)

   END SUBROUTINE charnock_1955
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This version of Charnock's relation employs a varying
!!Charnock parameter, taken from COARE 3.5 [Edson et al. (2001, JPO)].
!!The Charnock parameter CZC is varied from about .005 to .028
!!between 10-m wind speeds of 6 and 19 m/s.
   SUBROUTINE edson_etal_2013(Z_0,ustar,wsp10,visc,zu)
       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: ustar, visc, wsp10, zu
       REAL(kind_phys), INTENT(OUT) :: Z_0
       REAL(kind_phys), PARAMETER   :: m=0.0017, b=-0.005
       REAL(kind_phys)              :: CZC    ! variable charnock "constant"
       REAL(kind_phys)              :: wsp10m ! logarithmically calculated 10 m

       wsp10m = wsp10*log(10/1e-4)/log(zu/1e-4)
       wsp10m = MIN(19._kind_phys, wsp10m)
       CZC    = m*wsp10m + b
       CZC    = MAX(CZC, 0.0_kind_phys)

       Z_0 = CZC*ustar*ustar*g_inv + (0.11*visc/MAX(ustar,0.07_kind_phys))
       Z_0 = MAX( Z_0, 1.27e-7_kind_phys)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3_kind_phys)  !Davis et al. (2008)

   END SUBROUTINE edson_etal_2013
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This formulation for the thermal and moisture roughness lengths
!! (Zt and Zq) relates them to Z0 via the roughness Reynolds number (Ren).
!!This formula comes from Fairall et al. (2003). It is modified from
!!the original Garratt-Brutsaert model to better fit the COARE/HEXMAX
!!data. The formula for land uses a constant ratio (Z_0/7.4) taken
!!from Garratt (1992).
   SUBROUTINE garratt_1992(Zt,Zq,Z_0,Ren,landsea)
       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: Ren, Z_0,landsea
       REAL(kind_phys), INTENT(OUT) :: Zt,Zq
       REAL(kind_phys)              :: Rq
       REAL(kind_phys), PARAMETER   :: e=2.71828183

       IF (landsea-1.5 .GT. 0) THEN    !WATER

          Zt = Z_0*EXP(2.0 - (2.48*(Ren**0.25)))
          Zq = Z_0*EXP(2.0 - (2.28*(Ren**0.25)))

          Zq = MIN( Zq, 5.5e-5_kind_phys)
          Zq = MAX( Zq, 2.0e-9_kind_phys)
          Zt = MIN( Zt, 5.5e-5_kind_phys)
          Zt = MAX( Zt, 2.0e-9_kind_phys) !same lower limit as ECMWF
       ELSE                            !LAND
          Zq = Z_0/(e**2.)      !taken from Garratt (1980,1992)
          Zt = Zq
       ENDIF

    END SUBROUTINE garratt_1992
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!>This formulation for thermal and moisture roughness length (Zt and Zq)
!! as a function of the roughness Reynolds number (Ren) comes from the
!! COARE3.0 formulation, empirically derived from COARE and HEXMAX data
!! [Fairall et al. (2003)]. Edson et al. (2004; JGR) suspected that this
!!relationship overestimated the scalar roughness lengths for low Reynolds
!!number flows, so an optional smooth flow relationship, taken from Garratt
!!(1992, p. 102), is available for flows with Ren < 2.
!!
!!This is for use over water only.
    SUBROUTINE fairall_etal_2003(Zt,Zq,Ren,ustar,visc,rstoch,spp_sfc)
       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)   :: Ren,ustar,visc,rstoch
       INTEGER, INTENT(IN)           :: spp_sfc
       REAL(kind_phys), INTENT(OUT)  :: Zt,Zq

       IF (Ren .le. 2.) then

          Zt = (5.5e-5)*(Ren**(-0.60))
          Zq = Zt
          !FOR SMOOTH SEAS, CAN USE GARRATT
          !Zq = 0.2*visc/MAX(ustar,0.1)
          !Zq = 0.3*visc/MAX(ustar,0.1)

       ELSE

          !FOR ROUGH SEAS, USE COARE
          Zt = (5.5e-5)*(Ren**(-0.60))
          Zq = Zt

       ENDIF

       if (spp_sfc==1) then
          Zt = Zt + Zt * 0.5 * rstoch
          Zq = Zt
       endif

       Zt = MIN(Zt,1.0e-4_kind_phys)
       Zt = MAX(Zt,2.0e-9_kind_phys)

       Zq = MIN(Zt,1.0e-4_kind_phys)
       Zq = MAX(Zt,2.0e-9_kind_phys)

    END SUBROUTINE fairall_etal_2003
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This formulation for thermal and moisture roughness length (Zt and Zq)
!! as a function of the roughness Reynolds number (Ren) comes from the
!! COARE 3.5/4.0 formulation, empirically derived from COARE and HEXMAX data
!! The actual reference is unknown. This was passed along by Jim Edson (personal communication).
!! This is for use over water only, preferably open ocean.
    SUBROUTINE fairall_etal_2014(Zt,Zq,Ren,ustar,visc,rstoch,spp_sfc)
       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: Ren,ustar,visc,rstoch
       INTEGER, INTENT(IN)          :: spp_sfc
       REAL(kind_phys), INTENT(OUT) :: Zt,Zq

       !Zt = (5.5e-5)*(Ren**(-0.60))
       Zt = MIN(1.6E-4_kind_phys, 5.8E-5/(Ren**0.72))
       Zq = Zt

       IF (spp_sfc ==1) THEN
          Zt = MAX(Zt + Zt*0.5*rstoch,2.0e-9_kind_phys)
          Zq = MAX(Zt + Zt*0.5*rstoch,2.0e-9_kind_phys)
       ELSE
          Zt = MAX(Zt,2.0e-9_kind_phys)
          Zq = MAX(Zt,2.0e-9_kind_phys)
       ENDIF

    END SUBROUTINE fairall_etal_2014
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This is a modified version of Yang et al (2002 QJRMS, 2008 JAMC)
!! and Chen et al (2010, J of Hydromet). Although it was originally
!! designed for arid regions with bare soil, it is modified
!! here to perform over a broader spectrum of vegetation.
!!
!!The original formulation relates the thermal roughness length (Zt)
!!to u* and T*:
!!
!! Zt = ht * EXP(-beta*(ustar**0.5)*(ABS(tstar)**0.25))
!!
!!where ht = Renc*visc/ustar and the critical Reynolds number
!!(Renc) = 70. Beta was originally = 10 (2002 paper) but was revised
!!to 7.2 (in 2008 paper). Their form typically varies the
!!ratio Z0/Zt by a few orders of magnitude (1-1E4).
!!
!!This modified form uses beta = 1.5 and a variable Renc (function of Z_0),
!!so zt generally varies similarly to the Zilitinkevich form (with Czil = 0.1)
!!for very small or negative surface heat fluxes but can become close to the
!!Zilitinkevich with Czil = 0.2 for very large HFX (large negative T*).
!!Also, the exponent (0.25) on tstar was changed to 1.0, since we found
!!Zt was reduced too much for low-moderate positive heat fluxes.
!!
!!This should only be used over land!
       SUBROUTINE Yang_2008(Z_0,Zt,Zq,ustar,tstar,qst,Ren,visc)

       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: Z_0, Ren, ustar, tstar, qst, visc
       REAL(kind_phys) ::      ht,     &! roughness height at critical Reynolds number
                               tstar2, &! bounded T*, forced to be non-positive
                               qstar2, &! bounded q*, forced to be non-positive
                               Z_02,   &! bounded Z_0 for variable Renc2 calc
                               Renc2    ! variable Renc, function of Z_0
       REAL(kind_phys), INTENT(OUT) :: Zt,Zq
       REAL(kind_phys), PARAMETER   :: Renc=300., & !old constant Renc
                                       beta=1.5,  & !important for diurnal variation
                                       m=170.,    & !slope for Renc2 function
                                       b=691.       !y-intercept for Renc2 function

       Z_02 = MIN(Z_0,0.5_kind_phys)
       Z_02 = MAX(Z_02,0.04_kind_phys)
       Renc2= b + m*log(Z_02)
       ht     = Renc2*visc/MAX(ustar,0.01_kind_phys)
       tstar2 = MIN(tstar, 0.0_kind_phys)
       qstar2 = MIN(qst,0.0_kind_phys)

       Zt     = ht * EXP(-beta*(ustar**0.5)*(ABS(tstar2)**1.0))
       Zq     = ht * EXP(-beta*(ustar**0.5)*(ABS(qstar2)**1.0))
       !Zq     = Zt

       Zt = MIN(Zt, Z_0/2.0)
       Zq = MIN(Zq, Z_0/2.0)

    END SUBROUTINE Yang_2008
!--------------------------------------------------------------------
!  Taken from the GFS (sfc_diff.f) for comparison
!>\ingroup mynn_sfc
    SUBROUTINE GFS_z0_lnd(z0max,shdmax,z1,vegtype,ivegsrc,z0pert)

        !$acc routine seq
        REAL(kind_phys), INTENT(OUT)  :: z0max
        REAL(kind_phys), INTENT(IN)   :: shdmax,z1,z0pert
        INTEGER, INTENT(IN)           :: vegtype,ivegsrc
        REAL(kind_phys)               :: tem1, tem2

!            z0max = max(1.0e-6, min(0.01 * z0max, z1))
!already converted into meters in the wrapper
            z0max = max(1.0e-6_kind_phys, min(z0max, z1))
!** xubin's new z0  over land
            tem1  = 1.0 - shdmax
            tem2  = tem1 * tem1
            tem1  = 1.0  - tem2

            if( ivegsrc == 1 ) then

              if (vegtype == 10) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype == 6) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              elseif (vegtype == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            elseif (ivegsrc == 2 ) then

              if (vegtype == 7) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype == 8) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype == 9) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              elseif (vegtype == 11) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            endif

! mg, sfc-perts: add surface perturbations to z0max over land
            if (z0pert /= 0.0 ) then
              z0max = z0max * (10.**z0pert)
            endif

            z0max = max(z0max, 1.0e-6_kind_phys)

    END SUBROUTINE GFS_z0_lnd
!--------------------------------------------------------------------
!  Taken from the GFS (sfc_diff.f) for comparison
!>\ingroup mynn_sfc
    SUBROUTINE GFS_zt_lnd(ztmax,z0max,sigmaf,ztpert,ustar_lnd)

        !$acc routine seq
        REAL(kind_phys), INTENT(OUT)  :: ztmax
        REAL(kind_phys), INTENT(IN)   :: z0max,sigmaf,ztpert,ustar_lnd
        REAL(kind_phys)               :: czilc, tem1, tem2
        REAL(kind_phys), PARAMETER    :: ca = 0.4

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
           czilc = 0.8

           tem1  = 1.0 - sigmaf
           ztmax = z0max*exp( - tem1*tem1                 &
    &                     * czilc*ca*sqrt(ustar_lnd*(0.01/1.5e-05)))
!
!            czilc = 10.0 ** (- 4. * z0max) ! Trier et al. (2011, WAF)
!            ztmax = z0max * exp( - czilc * ca      &
!     &            * 258.2 * sqrt(ustar_lnd*z0max) )


! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if (ztpert /= 0.0) then
              ztmax = ztmax * (10.**ztpert)
            endif
            ztmax = max(ztmax, 1.0e-6_kind_phys)

    END SUBROUTINE GFS_zt_lnd
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
    SUBROUTINE GFS_z0_wat(z0rl_wat,ustar_wat,WSPD,z1,sfc_z0_type,redrag)

        !$acc routine seq
        REAL(kind_phys), INTENT(OUT)  :: z0rl_wat
        REAL(kind_phys), INTENT(INOUT):: ustar_wat
        REAL(kind_phys), INTENT(IN)   :: wspd,z1
        LOGICAL, INTENT(IN)           :: redrag
        INTEGER, INTENT(IN)           :: sfc_z0_type
        REAL(kind_phys)               :: z0,z0max,wind10m
        REAL(kind_phys), PARAMETER    :: charnock = 0.014, z0s_max=.317e-2

!            z0           = 0.01 * z0rl_wat
!Already converted to meters in the wrapper
            z0           = z0rl_wat
            z0max        = max(1.0e-6_kind_phys, min(z0,z1))
            ustar_wat    = sqrt(grav * z0 / charnock)
            wind10m      = wspd*log(10./1e-4)/log(z1/1e-4)
            !wind10m      = sqrt(u10m(i)*u10m(i)+v10m(i)*v10m(i))
!
            if (sfc_z0_type >= 0) then
              if (sfc_z0_type == 0) then
                z0 = (charnock / grav) * ustar_wat * ustar_wat

! mbek -- toga-coare flux algorithm
!               z0 = (charnock / g) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!               cc = ustar(i) * z0 / rnu
!               pp = cc / (1. + cc)
!               ff = g * arnu / (charnock * ustar(i) ** 3)
!               z0 = arnu / (ustar(i) * ff ** pp)

                if (redrag) then
                  !z0rl_wat = 100.0 * max(min(z0, z0s_max), 1.e-7)
                  z0rl_wat = max(min(z0, z0s_max), 1.e-7_kind_phys)
                else
                  !z0rl_wat = 100.0 * max(min(z0,.1), 1.e-7)
                  z0rl_wat = max(min(z0,.1_kind_phys), 1.e-7_kind_phys)
                endif

              elseif (sfc_z0_type == 6) then   ! wang
                 call znot_m_v6(wind10m, z0)  ! wind, m/s, z0, m
                 !z0rl_wat = 100.0 * z0          ! cm
              elseif (sfc_z0_type == 7) then   ! wang
                 call znot_m_v7(wind10m, z0)  ! wind, m/s, z0, m
                 !z0rl_wat = 100.0 * z0          ! cm
              else
                 z0rl_wat = 1.0e-6
              endif

            endif

    END SUBROUTINE GFS_z0_wat
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
    SUBROUTINE GFS_zt_wat(ztmax,z0rl_wat,restar,WSPD,z1,sfc_z0_type,device_errmsg,device_errflg)
        !$acc routine seq
        real(kind_phys), INTENT(OUT)  :: ztmax
        real(kind_phys), INTENT(IN)   :: wspd,z1,z0rl_wat,restar
        INTEGER, INTENT(IN)           :: sfc_z0_type

! Using device_errmsg and device_errflg rather than the CCPP errmsg and errflg
! so that this subroutine can be run on an accelerator device with OpenACC.
!        character(len=*), intent(out) :: errmsg
!        integer,          intent(out) :: errflg
        character(len=512), intent(out) :: device_errmsg
        integer,            intent(out) :: device_errflg

        real(kind_phys)               :: z0,z0max,wind10m,rat,ustar_wat
        real(kind_phys), PARAMETER    :: charnock = 0.014, z0s_max=.317e-2

        ! Initialize error-handling
!        errflg = 0
!        errmsg = ''
        device_errflg = 0
        device_errmsg = ''

!            z0           = 0.01 * z0rl_wat
!Already converted to meters in the wrapper
            z0           = z0rl_wat
            z0max        = max(1.0e-6_kind_phys, min(z0,z1))
            ustar_wat    = sqrt(grav * z0 / charnock)
            wind10m      = wspd*log(10./1e-4)/log(z1/1e-4)

!**  test xubin's new z0

!           ztmax  = z0max

!input            restar = max(ustar_wat(i)*z0max*visi, 0.000001)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat   = min(7.0_kind_phys, 2.67 * sqrt(sqrt(restar)) - 2.57)
            ztmax = max(z0max * exp(-rat), 1.0e-6_kind_phys)
!
            if (sfc_z0_type == 6) then
              call znot_t_v6(wind10m, ztmax)   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type == 7) then
              call znot_t_v7(wind10m, ztmax)   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type > 0) then
              write(0,*)'not a valid option for sfc_z0_type=',sfc_z0_type
!              errflg = 1
!              errmsg = 'ERROR(GFS_zt_wat): sfc_z0_type not valid.'
              device_errflg = 1
              device_errmsg = 'ERROR(GFS_zt_wat): sfc_z0_type not valid.'
              return

            endif

    END SUBROUTINE GFS_zt_wat
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! add fitted z0,zt curves for hurricane application (used in HWRF/HMON)
!! Weiguo Wang, 2019-0425

      SUBROUTINE znot_m_v6(uref, znotm)
      !$acc routine seq
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate areodynamical roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
! For high winds, try to fit available observational data
!
! Bin Liu, NOAA/NCEP/EMC 2017
!
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      REAL(kind_phys), INTENT(IN) :: uref
      REAL(kind_phys), INTENT(OUT):: znotm
      REAL(kind_phys), PARAMETER  :: p13 = -1.296521881682694e-02,     &
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,&
     &      p10 = -8.396975715683501e+00,                              &

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,&
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,&
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,&

     &      p35 =  1.840430200185075e-07, p34 = -2.793849676757154e-05,&
     &      p33 =  1.735308193700643e-03, p32 = -6.139315534216305e-02,&
     &      p31 =  1.255457892775006e+00, p30 = -1.663993561652530e+01,&

     &      p40 =  4.579369142033410e-04


       if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp(p10 + uref * (p11 + uref * (p12 + uref*p13)))
       elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23       &
     &              + uref * (p24 + uref * p25))))
       elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33  &
     &                   + uref * (p34 + uref * p35)))))
       elseif ( uref > 53.0) then
         znotm = p40
       else
          print*, 'Wrong input uref value:',uref
       endif

      END SUBROUTINE znot_m_v6
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> Calculate scalar roughness over water with input 10-m wind 
!! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
!! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!!
!!\author Bin Liu, NOAA/NCEP/EMC 2017
!   
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
      SUBROUTINE znot_t_v6(uref, znott)

      !$acc routine seq
      IMPLICIT NONE
!
      REAL(kind_phys), INTENT(IN) :: uref
      REAL(kind_phys), INTENT(OUT):: znott
      REAL(kind_phys), PARAMETER  ::      p00 =  1.100000000000000e-04,&
     &      p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,&
     &      p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,&
     &      p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,&

     &      p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,&
     &      p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,&
     &      p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,&

     &      p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,&
     &      p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,&
     &      p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,&

     &      p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,&
     &      p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,&
     &      p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,&

     &      p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,&
     &      p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,&
     &      p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,&
     &      p50 = -1.036679430885215e-02,                              &

     &      p60 =  4.751256171799112e-05

      if (uref >= 0.0 .and. uref < 5.9 ) then
         znott = p00
      elseif (uref >= 5.9 .and. uref <= 15.4) then
         znott = p10 + uref * (p11 + uref * (p12 + uref * (p13  &
     &               + uref * (p14 + uref * p15))))
      elseif (uref > 15.4 .and. uref <= 21.6) then
         znott = p20 + uref * (p21 + uref * (p22 + uref * (p23  &
     &               + uref * (p24 + uref * p25))))
      elseif (uref > 21.6 .and. uref <= 42.2) then
         znott = p30 + uref * (p31 + uref * (p32 + uref * (p33  &
     &               + uref * (p34 + uref * p35))))
      elseif ( uref > 42.2 .and. uref <= 53.3) then
         znott = p40 + uref * (p41 + uref * (p42 + uref * (p43  &
     &               + uref * (p44 + uref * p45))))
      elseif ( uref > 53.3 .and. uref <= 80.0) then
         znott = p50 + uref * (p51 + uref * (p52 + uref * (p53  &
     &               + uref * (p54 + uref * (p55 + uref * p56)))))
      elseif ( uref > 80.0) then
         znott = p60
      else
         print*, 'Wrong input uref value:',uref
      endif

      END SUBROUTINE znot_t_v6

!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!> Calculate areodynamical roughness over water with input 10-m wind
!! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
!! For high winds, try to fit available observational data
!! Comparing to znot_t_v6, slightly decrease Cd for higher wind speed
!!   
!!\author Bin Liu, NOAA/NCEP/EMC 2018
      SUBROUTINE znot_m_v7(uref, znotm)

      !$acc routine seq
      IMPLICIT NONE
!
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      REAL(kind_phys), INTENT(IN) :: uref
      REAL(kind_phys), INTENT(OUT):: znotm

      REAL(kind_phys), PARAMETER  ::      p13 = -1.296521881682694e-02,&
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,&
     &      p10 = -8.396975715683501e+00,                              &

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,&
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,&
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,&

     &      p35 =  1.897534489606422e-07, p34 = -3.019495980684978e-05,&
     &      p33 =  1.931392924987349e-03, p32 = -6.797293095862357e-02,&
     &      p31 =  1.346757797103756e+00, p30 = -1.707846930193362e+01,&

     &      p40 =  3.371427455376717e-04

      if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp( p10 + uref * (p11 + uref * (p12 + uref * p13)))
      elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23        &
     &              + uref * (p24 + uref * p25))))
      elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33   &
     &                   + uref * (p34 + uref * p35)))))
      elseif ( uref > 53.0) then
        znotm = p40
      else
        print*, 'Wrong input uref value:',uref
      endif

      END SUBROUTINE znot_m_v7
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> Calculate scalar roughness over water with input 10-m wind
!! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
!! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!! To be compatible with the slightly decreased Cd for higher wind speed
!!    
!!\author Bin Liu, NOAA/NCEP/EMC 2018
      SUBROUTINE znot_t_v7(uref, znott)

      !$acc routine seq
      IMPLICIT NONE
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

      REAL(kind_phys), INTENT(IN) :: uref
      REAL(kind_phys), INTENT(OUT):: znott
      REAL(kind_phys), PARAMETER  ::      p00 =  1.100000000000000e-04,&
     &      p15 = -9.193764479895316e-10, p14 =  7.052217518653943e-08,&
     &      p13 = -2.163419217747114e-06, p12 =  3.342963077911962e-05,&
     &      p11 = -2.633566691328004e-04, p10 =  8.644979973037803e-04,&

     &      p25 = -9.402722450219142e-12, p24 =  1.325396583616614e-09,&
     &      p23 = -7.299148051141852e-08, p22 =  1.982901461144764e-06,&
     &      p21 = -2.680293455916390e-05, p20 =  1.484341646128200e-04,&

     &      p35 =  7.921446674311864e-12, p34 = -1.019028029546602e-09,&
     &      p33 =  5.251986927351103e-08, p32 = -1.337841892062716e-06,&
     &      p31 =  1.659454106237737e-05, p30 = -7.558911792344770e-05,&

     &      p45 = -2.694370426850801e-10, p44 =  5.817362913967911e-08,&
     &      p43 = -5.000813324746342e-06, p42 =  2.143803523428029e-04,&
     &      p41 = -4.588070983722060e-03, p40 =  3.924356617245624e-02,&

     &      p56 = -1.663918773476178e-13, p55 =  6.724854483077447e-11,&
     &      p54 = -1.127030176632823e-08, p53 =  1.003683177025925e-06,&
     &      p52 = -5.012618091180904e-05, p51 =  1.329762020689302e-03,&
     &      p50 = -1.450062148367566e-02, p60 =  6.840803042788488e-05

        if (uref >= 0.0 .and. uref < 5.9 ) then
           znott = p00
        elseif (uref >= 5.9 .and. uref <= 15.4) then
           znott = p10 + uref * (p11 + uref * (p12 + uref * (p13     &
     &                 + uref * (p14 + uref * p15))))
        elseif (uref > 15.4 .and. uref <= 21.6) then
           znott = p20 + uref * (p21 + uref * (p22 + uref * (p23     &
     &                 + uref * (p24 + uref * p25))))
        elseif (uref > 21.6 .and. uref <= 42.6) then
           znott = p30 + uref * (p31 + uref * (p32 + uref * (p33     &
     &                 + uref * (p34 + uref * p35))))
        elseif ( uref > 42.6 .and. uref <= 53.0) then
           znott = p40 + uref * (p41 + uref * (p42 + uref * (p43     &
     &                 + uref * (p44 + uref * p45))))
        elseif ( uref > 53.0 .and. uref <= 80.0) then
           znott = p50 + uref * (p51 + uref * (p52 + uref * (p53     &
     &                 + uref * (p54 + uref * (p55 + uref * p56)))))
        elseif ( uref > 80.0) then
           znott = p60
        else
           print*, 'Wrong input uref value:',uref
        endif

        END SUBROUTINE znot_t_v7

!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This is taken from Andreas (2002; J. of Hydromet) and
!! Andreas et al. (2005; BLM).
!!
!! This should only be used over snow/ice!
    SUBROUTINE Andreas_2002(Z_0,bvisc,ustar,Zt,Zq)

       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: Z_0, bvisc, ustar
       REAL(kind_phys), INTENT(OUT) :: Zt, Zq
       REAL(kind_phys) :: Ren2, zntsno

       REAL(kind_phys), PARAMETER  ::                               &
                           bt0_s=1.25,  bt0_t=0.149,  bt0_r=0.317,  &
                           bt1_s=0.0,   bt1_t=-0.55,  bt1_r=-0.565, &
                           bt2_s=0.0,   bt2_t=0.0,    bt2_r=-0.183

       REAL(kind_phys), PARAMETER  ::                               &
                           bq0_s=1.61,  bq0_t=0.351,  bq0_r=0.396,  &
                           bq1_s=0.0,   bq1_t=-0.628, bq1_r=-0.512, &
                           bq2_s=0.0,   bq2_t=0.0,    bq2_r=-0.180

      !Calculate zo for snow (Andreas et al. 2005, BLM)
       zntsno = 0.135*bvisc/ustar +                                 &
               (0.035*(ustar*ustar)*g_inv) *                        &
               (5.*exp(-1.*(((ustar - 0.18)/0.1)*((ustar - 0.18)/0.1))) + 1.)
       Ren2 = ustar*zntsno/bvisc

       ! Make sure that Re is not outside of the range of validity
       ! for using their equations
       IF (Ren2 .gt. 1000.) Ren2 = 1000.

       IF (Ren2 .le. 0.135) then

          Zt = zntsno*EXP(bt0_s + bt1_s*LOG(Ren2) + bt2_s*LOG(Ren2)**2)
          Zq = zntsno*EXP(bq0_s + bq1_s*LOG(Ren2) + bq2_s*LOG(Ren2)**2)

       ELSE IF (Ren2 .gt. 0.135 .AND. Ren2 .lt. 2.5) then

          Zt = zntsno*EXP(bt0_t + bt1_t*LOG(Ren2) + bt2_t*LOG(Ren2)**2)
          Zq = zntsno*EXP(bq0_t + bq1_t*LOG(Ren2) + bq2_t*LOG(Ren2)**2)

       ELSE

          Zt = zntsno*EXP(bt0_r + bt1_r*LOG(Ren2) + bt2_r*LOG(Ren2)**2)
          Zq = zntsno*EXP(bq0_r + bq1_r*LOG(Ren2) + bq2_r*LOG(Ren2)**2)

       ENDIF

    END SUBROUTINE Andreas_2002
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the stability functions based off
!! of Hogstrom (1996).
    SUBROUTINE PSI_Hogstrom_1996(psi_m, psi_h, zL, Zt, Z_0, Za)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL, Zt, Z_0, Za
       REAL(kind_phys), INTENT(OUT) :: psi_m, psi_h
       REAL(kind_phys)  :: x, x0, y, y0, zmL, zhL

       zmL = Z_0*zL/Za
       zhL = Zt*zL/Za

       IF (zL .gt. 0.) THEN  !STABLE (not well tested - seem large)

          psi_m = -5.3*(zL - zmL)
          psi_h = -8.0*(zL - zhL)

       ELSE                 !UNSTABLE

          x = (1.-19.0*zL)**0.25
          x0= (1.-19.0*zmL)**0.25
          y = (1.-11.6*zL)**0.5
          y0= (1.-11.6*zhL)**0.5

          psi_m = 2.*LOG((1.+x)/(1.+x0)) +         &
                    &LOG((1.+x**2.)/(1.+x0**2.)) - &
                    &2.0*ATAN(x) + 2.0*ATAN(x0)
          psi_h = 2.*LOG((1.+y)/(1.+y0))

       ENDIF

    END SUBROUTINE PSI_Hogstrom_1996
!--------------------------------------------------------------------
!> \ingroup mynn_sfc
!> This subroutine returns the stability functions based off
!! of Hogstrom (1996), but with different constants compatible
!! with Dyer and Hicks (1970/74?). This formulation is used for
!! testing/development by Nakanishi (personal communication).
    SUBROUTINE PSI_DyerHicks(psi_m, psi_h, zL, Zt, Z_0, Za)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL, Zt, Z_0, Za
       REAL(kind_phys), INTENT(OUT) :: psi_m, psi_h
       REAL(kind_phys)  :: x, x0, y, y0, zmL, zhL

       zmL = Z_0*zL/Za  !Zo/L
       zhL = Zt*zL/Za   !Zt/L

       IF (zL .gt. 0.) THEN  !STABLE

          psi_m = -5.0*(zL - zmL)
          psi_h = -5.0*(zL - zhL)

       ELSE                 !UNSTABLE

          x = (1.-16.*zL)**0.25
          x0= (1.-16.*zmL)**0.25

          y = (1.-16.*zL)**0.5
          y0= (1.-16.*zhL)**0.5

          psi_m = 2.*LOG((1.+x)/(1.+x0)) +         &
                    &LOG((1.+x**2.)/(1.+x0**2.)) - &
                    &2.0*ATAN(x) + 2.0*ATAN(x0)
          psi_h = 2.*LOG((1.+y)/(1.+y0))

       ENDIF

    END SUBROUTINE PSI_DyerHicks
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the stability functions based off
!! of Beljaar and Holtslag 1991, which is an extension of Holtslag
!! and Debruin 1989.
    SUBROUTINE PSI_Beljaars_Holtslag_1991(psi_m, psi_h, zL)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL
       REAL(kind_phys), INTENT(OUT) :: psi_m, psi_h
       REAL(kind_phys), PARAMETER   :: a=1., b=0.666, c=5., d=0.35

       IF (zL .lt. 0.) THEN  !UNSTABLE

          WRITE(*,*)"WARNING: Universal stability functions from"
          WRITE(*,*)"        Beljaars and Holtslag (1991) should only"
          WRITE(*,*)"        be used in the stable regime!"
          psi_m = 0.
          psi_h = 0.

       ELSE                 !STABLE

          psi_m = -(a*zL + b*(zL -(c/d))*exp(-d*zL) + (b*c/d))
          psi_h = -((1.+.666*a*zL)**1.5 + &
                  b*(zL - (c/d))*exp(-d*zL) + (b*c/d) -1.)

       ENDIF

    END SUBROUTINE PSI_Beljaars_Holtslag_1991
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the stability functions come from
!! Zilitinkevich and Esau (2007, BM), which are formulatioed from the
!! "generalized similarity theory" and tuned to the LES DATABASE64
!! to determine their dependence on z/L.
    SUBROUTINE PSI_Zilitinkevich_Esau_2007(psi_m, psi_h, zL)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL
       REAL(kind_phys), INTENT(OUT) :: psi_m, psi_h
       REAL(kind_phys), PARAMETER   :: Cm=3.0, Ct=2.5

       IF (zL .lt. 0.) THEN  !UNSTABLE

          WRITE(*,*)"WARNING: Universal stability function from"
          WRITE(*,*)"        Zilitinkevich and Esau (2007) should only"
          WRITE(*,*)"        be used in the stable regime!"
          psi_m = 0.
          psi_h = 0.

       ELSE                 !STABLE

          psi_m = -Cm*(zL**(5./6.))
          psi_h = -Ct*(zL**(4./5.))

       ENDIF

    END SUBROUTINE PSI_Zilitinkevich_Esau_2007
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the flux-profile relationships
!! of Businger el al. 1971.
    SUBROUTINE PSI_Businger_1971(psi_m, psi_h, zL)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL
       REAL(kind_phys), INTENT(OUT) :: psi_m, psi_h
       REAL(kind_phys)  :: x, y
       REAL(kind_phys), PARAMETER  ::  Pi180 = 3.14159265/180.

       IF (zL .lt. 0.) THEN  !UNSTABLE

          x = (1. - 15.0*zL)**0.25
          y = (1. - 9.0*zL)**0.5

          psi_m = LOG(((1.+x)/2.)**2.) + &
                 &LOG((1.+x**2.)/2.) - &
                 &2.0*ATAN(x) + Pi180*90.
          psi_h = 2.*LOG((1.+y)/2.)

       ELSE                 !STABLE

          psi_m = -4.7*zL
          psi_h = -(4.7/0.74)*zL

       ENDIF

    END SUBROUTINE PSI_Businger_1971
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns flux-profile relatioships based off
!!of Lobocki (1993), which is derived from the MY-level 2 model.
!!Suselj and Sood (2010) applied the surface layer length scales
!!from Nakanishi (2001) to get this new relationship. These functions
!!are more agressive (larger magnitude) than most formulations. They
!!showed improvement over water, but untested over land.
    SUBROUTINE PSI_Suselj_Sood_2010(psi_m, psi_h, zL)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL
       REAL(kind_phys), INTENT(OUT) :: psi_m, psi_h
       REAL(kind_phys), PARAMETER   :: Rfc=0.19, Ric=0.183, PHIT=0.8

       IF (zL .gt. 0.) THEN  !STABLE

          psi_m = -(zL/Rfc + 1.1223*EXP(1.-1.6666/zL))
          !psi_h = -zL*Ric/((Rfc**2.)*PHIT) + 8.209*(zL**1.1091)
          !THEIR EQ FOR PSI_H CRASHES THE MODEL AND DOES NOT MATCH
          !THEIR FIG 1. THIS EQ (BELOW) MATCHES THEIR FIG 1 BETTER:
          psi_h = -(zL*Ric/((Rfc**2.)*5.) + 7.09*(zL**1.1091))

       ELSE                 !UNSTABLE

          psi_m = 0.9904*LOG(1. - 14.264*zL)
          psi_h = 1.0103*LOG(1. - 16.3066*zL)

       ENDIF

    END SUBROUTINE PSI_Suselj_Sood_2010
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine returns the stability functions based off
!! of Cheng and Brutseart (2005, BLM), for use in stable conditions only.
!! The returned values are the combination of psi((za+zo)/L) - psi(z0/L)
    SUBROUTINE PSI_CB2005(psim1,psih1,zL,z0L)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(IN)  :: zL,z0L
       REAL(kind_phys), INTENT(OUT) :: psim1,psih1

       psim1 = -6.1*LOG(zL + (1.+ zL**2.5)**0.4)            &
               -6.1*LOG(z0L + (1.+ z0L**2.5)**0.4)
       psih1 = -5.5*log(zL + (1.+ zL**1.1)**0.90909090909)  &
               -5.5*log(z0L + (1.+ z0L**1.1)**0.90909090909)

    END SUBROUTINE PSI_CB2005
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine returns a more robust z/L that best matches
!! the z/L from Hogstrom (1996) for unstable conditions and Beljaars
!! and Holtslag (1991) for stable conditions.
    SUBROUTINE Li_etal_2010(zL, Rib, zaz0, z0zt)

       !$acc routine seq
       IMPLICIT NONE
       REAL(kind_phys), INTENT(OUT)  :: zL
       REAL(kind_phys), INTENT(IN) :: Rib, zaz0, z0zt
       REAL(kind_phys) :: alfa, beta, zaz02, z0zt2
       REAL(kind_phys), PARAMETER  ::                                &
                          & au11=0.045,   bu11=0.003,   bu12=0.0059, &
                          & bu21=-0.0828, bu22=0.8845,  bu31=0.1739, &
                          & bu32=-0.9213, bu33=-0.1057
       REAL(kind_phys), PARAMETER  ::                                &
                          & aw11=0.5738,  aw12=-0.4399, aw21=-4.901, &
                          & aw22=52.50,   bw11=-0.0539, bw12=1.540,  &
                          & bw21=-0.669,  bw22=-3.282
       REAL(kind_phys), PARAMETER  ::                                &
                          & as11=0.7529,  as21=14.94,   bs11=0.1569, &
                          & bs21=-0.3091, bs22=-1.303

       !set limits according to Li et al (2010), p 157.
       zaz02=zaz0
       IF (zaz0 .lt. 100.0) zaz02=100.
       IF (zaz0 .gt. 100000.0) zaz02=100000.

       !set more limits according to Li et al (2010)
       z0zt2=z0zt
       IF (z0zt .lt. 0.5) z0zt2=0.5
       IF (z0zt .gt. 100.0) z0zt2=100.

       alfa = LOG(zaz02)
       beta = LOG(z0zt2)

       IF (Rib .le. 0.0) THEN
          zL = au11*alfa*Rib**2 + (                   &
               &  (bu11*beta + bu12)*alfa**2 +        &
               &  (bu21*beta + bu22)*alfa    +        &
               &  (bu31*beta**2 + bu32*beta + bu33))*Rib
          !if(zL .LT. -15 .OR. zl .GT. 0.)print*,"VIOLATION Rib<0:",zL
          zL = MAX(zL,-15._kind_phys) !LIMITS SET ACCORDING TO Li et al (2010)
          zL = MIN(zL,0._kind_phys)   !Figure 1.
       ELSEIF (Rib .gt. 0.0 .AND. Rib .le. 0.2) THEN
          zL = ((aw11*beta + aw12)*alfa +             &
             &  (aw21*beta + aw22))*Rib**2 +          &
             & ((bw11*beta + bw12)*alfa +             &
             &  (bw21*beta + bw22))*Rib
          !if(zL .LT. 0 .OR. zl .GT. 4)print*,"VIOLATION 0<Rib<0.2:",zL
          zL = MIN(zL,4._kind_phys) !LIMITS APPROX SET ACCORDING TO Li et al (2010)
          zL = MAX(zL,0._kind_phys) !THEIR FIGURE 1B.
       ELSE
          zL = (as11*alfa + as21)*Rib + bs11*alfa +   &
             &  bs21*beta + bs22
          !if(zL .LE. 1 .OR. zl .GT. 23)print*,"VIOLATION Rib>0.2:",zL
          zL = MIN(zL,20._kind_phys) !LIMITS ACCORDING TO Li et al (2010), THIER
                           !FIGUE 1C.
          zL = MAX(zL,1._kind_phys)
       ENDIF

    END SUBROUTINE Li_etal_2010
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
      REAL(kind_phys) function zolri(ri,za,z0,zt,zol1,psi_opt)

      !> This iterative algorithm was taken from the revised surface layer
      !! scheme in WRF-ARW, written by Pedro Jimenez and Jimy Dudhia and
      !! summarized in Jimenez et al. (2012, MWR). This function was adapted
      !! to input the thermal roughness length, zt, (as well as z0) and use initial
      !! estimate of z/L.

      IMPLICIT NONE
      REAL(kind_phys), INTENT(IN) :: ri,za,z0,zt,zol1
      INTEGER, INTENT(IN) :: psi_opt
      REAL(kind_phys) :: x1,x2,fx1,fx2
      INTEGER :: n
      INTEGER, PARAMETER :: nmax = 20
      !REAL(kind_phys), DIMENSION(nmax):: zLhux

      if (ri.lt.0.)then
         x1=zol1 - 0.02  !-5.
         x2=0.
      else
         x1=0.
         x2=zol1 + 0.02 !5.
      endif

      n=1
      fx1=zolri2(x1,ri,za,z0,zt,psi_opt)
      fx2=zolri2(x2,ri,za,z0,zt,psi_opt)

      Do While (abs(x1 - x2) > 0.01 .and. n < nmax)
        if(abs(fx2).lt.abs(fx1))then
          x1=x1-fx1/(fx2-fx1)*(x2-x1)
          fx1=zolri2(x1,ri,za,z0,zt,psi_opt)
          zolri=x1
        else
          x2=x2-fx2/(fx2-fx1)*(x2-x1)
          fx2=zolri2(x2,ri,za,z0,zt,psi_opt)
          zolri=x2
        endif
        n=n+1
        !print*," n=",n," x1=",x1," x2=",x2
        !zLhux(n)=zolri
      enddo

      if (n==nmax .and. abs(x1 - x2) >= 0.01) then
         !if convergence fails, use approximate values:
         CALL Li_etal_2010(zolri, ri, za/z0, z0/zt)
         !zLhux(n)=zolri
         !print*,"iter FAIL, n=",n," Ri=",ri," z0=",z0
      else
         !print*,"SUCCESS,n=",n," Ri=",ri," z0=",z0
      endif

      end function
!-------------------------------------------------------------------
      REAL(kind_phys) function zolri2(zol2,ri2,za,z0,zt,psi_opt)

      ! INPUT: =================================
      ! zol2 - estimated z/L
      ! ri2  - calculated bulk Richardson number
      ! za   - 1/2 depth of first model layer
      ! z0   - aerodynamic roughness length
      ! zt   - thermal roughness length
      ! OUTPUT: ================================
      ! zolri2 - delta Ri

      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: psi_opt
      REAL(kind_phys), INTENT(IN)    :: ri2,za,z0,zt
      REAL(kind_phys), INTENT(INOUT) :: zol2
      REAL(kind_phys) :: zol20,zol3,psim1,psih1,psix2,psit2,zolt

      if(zol2*ri2 .lt. 0.)zol2=0.  ! limit zol2 - must be same sign as ri2

      zol20=zol2*z0/za ! z0/L
      zol3=zol2+zol20  ! (z+z0)/L
      zolt=zol2*zt/za  ! zt/L

      if (ri2.lt.0) then
         !psix2=log((za+z0)/z0)-(psim_unstable(zol3)-psim_unstable(zol20))
         !psit2=log((za+zt)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
         psit2=MAX(log((za+z0)/zt)-(psih_unstable(zol3,psi_opt)-psih_unstable(zolt,psi_opt)), 1.0_kind_phys)
         psix2=MAX(log((za+z0)/z0)-(psim_unstable(zol3,psi_opt)-psim_unstable(zol20,psi_opt)), 1.0_kind_phys)
      else
         !psix2=log((za+z0)/z0)-(psim_stable(zol3)-psim_stable(zol20))
         !psit2=log((za+zt)/zt)-(psih_stable(zol3)-psih_stable(zol20))
         psit2=MAX(log((za+z0)/zt)-(psih_stable(zol3,psi_opt)-psih_stable(zolt,psi_opt)), 1.0_kind_phys)
         psix2=MAX(log((za+z0)/z0)-(psim_stable(zol3,psi_opt)-psim_stable(zol20,psi_opt)),1.0_kind_phys)
      endif

      zolri2=zol2*psit2/psix2**2 - ri2
      !print*,"  target ri=",ri2," est ri=",zol2*psit2/psix2**2

      end function
!====================================================================

      REAL(kind_phys) function zolrib(ri,za,z0,zt,logz0,logzt,zol1,psi_opt)

      !$acc routine seq
      ! This iterative algorithm to compute z/L from bulk-Ri

      IMPLICIT NONE
      REAL(kind_phys), INTENT(IN) :: ri,za,z0,zt,logz0,logzt
      INTEGER, INTENT(IN)         :: psi_opt
      REAL(kind_phys), INTENT(INOUT) :: zol1
      REAL(kind_phys) :: zol20,zol3,zolt,zolold
      INTEGER :: n
      INTEGER, PARAMETER :: nmax = 20
      !REAL(kind_phys), DIMENSION(nmax):: zLhux
      REAL(kind_phys) :: psit2,psix2

      !print*,"+++++++INCOMING: z/L=",zol1," ri=",ri
      if (zol1*ri .lt. 0.) THEN
         !print*,"begin: WRONG QUADRANTS: z/L=",zol1," ri=",ri
         zol1=0.
      endif

      if (ri .lt. 0.) then
        zolold=-99999.
        zolrib=-66666.
      else
        zolold=99999.
        zolrib=66666.
      endif
      n=1

      DO While (abs(zolold - zolrib) > 0.01 .and. n < nmax)

        if(n==1)then
          zolold=zol1
        else
          zolold=zolrib
        endif
        zol20=zolold*z0/za ! z0/L
        zol3=zolold+zol20  ! (z+z0)/L
        zolt=zolold*zt/za  ! zt/L
        !print*,"z0/L=",zol20," (z+z0)/L=",zol3," zt/L=",zolt
        if (ri.lt.0) then
           !psit2=log((za+zt)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
           !psit2=log((za+z0)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
           psit2=MAX(logzt-(psih_unstable(zol3,psi_opt)-psih_unstable(zolt,psi_opt)), 1.0_kind_phys)
           psix2=MAX(logz0-(psim_unstable(zol3,psi_opt)-psim_unstable(zol20,psi_opt)), 1.0_kind_phys)
        else
           !psit2=log((za+zt)/zt)-(psih_stable(zol3)-psih_stable(zol20))
           !psit2=log((za+z0)/zt)-(psih_stable(zol3)-psih_stable(zol20))
           psit2=MAX(logzt-(psih_stable(zol3,psi_opt)-psih_stable(zolt,psi_opt)), 1.0_kind_phys)
           psix2=MAX(logz0-(psim_stable(zol3,psi_opt)-psim_stable(zol20,psi_opt)), 1.0_kind_phys)
        endif
        !print*,"n=",n," psit2=",psit2," psix2=",psix2
        zolrib=ri*psix2**2/psit2
        !zLhux(n)=zolrib
        n=n+1
      enddo

      if (n==nmax .and. abs(zolold - zolrib) > 0.01 ) then
         !print*,"iter FAIL, n=",n," Ri=",ri," z/L=",zolri
         !if convergence fails, use approximate values:
         CALL Li_etal_2010(zolrib, ri, za/z0, z0/zt)
         !zLhux(n)=zolrib
         !print*,"FAILED, n=",n," Ri=",ri," z0=",z0
         !print*,"z/L=",zLhux(1:nmax)
      else
         !if(zolrib*ri .lt. 0.) THEN
         !   !print*,"end: WRONG QUADRANTS: z/L=",zolrib," ri=",ri
         !   !CALL Li_etal_2010(zolrib, ri, za/z0, z0/zt)
         !endif
         !print*,"SUCCESS,n=",n," Ri=",ri," z0=",z0
      endif

      end function
!====================================================================
!>\ingroup mynn_sfc
!!
   SUBROUTINE psi_init(psi_opt,errmsg,errflg)

    integer                       :: N,psi_opt
    real(kind_phys)               :: zolf
    character(len=*), intent(out) :: errmsg
    integer, intent(out)          :: errflg

    if (psi_opt == 0) then
       DO N=0,1000
          ! stable function tables
          zolf = float(n)*0.01
          psim_stab(n)=psim_stable_full(zolf)
          psih_stab(n)=psih_stable_full(zolf)

          ! unstable function tables
          zolf = -float(n)*0.01
          psim_unstab(n)=psim_unstable_full(zolf)
          psih_unstab(n)=psih_unstable_full(zolf)
       ENDDO
    else
       DO N=0,1000
          ! stable function tables
          zolf = float(n)*0.01
          psim_stab(n)=psim_stable_full_gfs(zolf)
          psih_stab(n)=psih_stable_full_gfs(zolf)

          ! unstable function tables
          zolf = -float(n)*0.01
          psim_unstab(n)=psim_unstable_full_gfs(zolf)
          psih_unstab(n)=psih_unstable_full_gfs(zolf)
       ENDDO
    endif

    !Simple test to see if initialization worked:
    if (psim_stab(1) < 0. .AND. psih_stab(1) < 0. .AND. &
        psim_unstab(1) > 0. .AND. psih_unstab(1) > 0.) then
       errmsg = 'In MYNN SFC, Psi tables have been initialized'
       errflg = 0
    else
       errmsg = 'Error in MYNN SFC: Problem initializing psi tables'
       errflg = 1
    endif

   END SUBROUTINE psi_init
! ==================================================================
! ... integrated similarity functions from MYNN...
!
!>\ingroup mynn_sfc
   real(kind_phys) function psim_stable_full(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf

        !psim_stable_full=-6.1*log(zolf+(1+zolf**2.5)**(1./2.5))
        psim_stable_full=-6.1*log(zolf+(1+zolf**2.5)**0.4)

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_stable_full(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf

        !psih_stable_full=-5.3*log(zolf+(1+zolf**1.1)**(1./1.1))
        psih_stable_full=-5.3*log(zolf+(1+zolf**1.1)**0.9090909090909090909)

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psim_unstable_full(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf,x,ym,psimc,psimk

        x=(1.-16.*zolf)**.25
        !psimk=2*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*ATAN(1.)
        !psimk=2.*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*atan1
        psimk=2.*LOG(0.5*(1+X))+LOG(0.5*(1+X*X))-2.*ATAN(X)+2.*atan1

        ym=(1.-10.*zolf)**onethird
        !psimc=(3./2.)*log((ym**2.+ym+1.)/3.)-sqrt(3.)*ATAN((2.*ym+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
        psimc=1.5*log((ym**2 + ym+1.)*onethird)-sqrt3*ATAN((2.*ym+1)/sqrt3)+4.*atan1/sqrt3

        psim_unstable_full=(psimk+zolf**2*(psimc))/(1+zolf**2.)

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_unstable_full(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf,y,yh,psihc,psihk

        y=(1.-16.*zolf)**.5
        !psihk=2.*log((1+y)/2.)
        psihk=2.*log((1+y)*0.5)

        yh=(1.-34.*zolf)**onethird
        !psihc=(3./2.)*log((yh**2.+yh+1.)/3.)-sqrt(3.)*ATAN((2.*yh+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
        psihc=1.5*log((yh**2.+yh+1.)*onethird)-sqrt3*ATAN((2.*yh+1)/sqrt3)+4.*atan1/sqrt3

        psih_unstable_full=(psihk+zolf**2*(psihc))/(1+zolf**2)

   end function

! ==================================================================
! ... integrated similarity functions from GFS...
!
!>\ingroup mynn_sfc
!!
   REAL(kind_phys) function psim_stable_full_gfs(zolf)
        !$acc routine seq
        REAL(kind_phys) :: zolf
        REAL(kind_phys), PARAMETER :: alpha4 = 20.
        REAL(kind_phys) :: aa

        aa     = sqrt(1. + alpha4 * zolf)
        psim_stable_full_gfs  = -1.*aa + log(aa + 1.)

   end function

!>\ingroup mynn_sfc
!!
   real(kind_phys) function psih_stable_full_gfs(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf
        real(kind_phys), PARAMETER :: alpha4 = 20.
        real(kind_phys) :: bb

        bb     = sqrt(1. + alpha4 * zolf)
        psih_stable_full_gfs  = -1.*bb + log(bb + 1.)

   end function

!>\ingroup mynn_sfc
!!
   real(kind_phys) function psim_unstable_full_gfs(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf
        real(kind_phys) :: hl1,tem1
        real(kind_phys), PARAMETER :: a0=-3.975,  a1=12.32,  &
                           b1=-7.755,  b2=6.041

        if (zolf .ge. -0.5) then
           hl1   = zolf
           psim_unstable_full_gfs  = (a0  + a1*hl1)  * hl1   / (1.+ (b1+b2*hl1)  *hl1)
        else
           hl1   = -zolf
           tem1  = 1.0 / sqrt(hl1)
           psim_unstable_full_gfs  = log(hl1) + 2. * sqrt(tem1) - .8776
        end if

   end function

!>\ingroup mynn_sfc
!!
   real(kind_phys) function psih_unstable_full_gfs(zolf)
        !$acc routine seq
        real(kind_phys) :: zolf
        real(kind_phys) :: hl1,tem1
        real(kind_phys), PARAMETER :: a0p=-7.941, a1p=24.75, &
                           b1p=-8.705, b2p=7.899

        if (zolf .ge. -0.5) then
           hl1   = zolf
           psih_unstable_full_gfs  = (a0p + a1p*hl1) * hl1   / (1.+ (b1p+b2p*hl1)*hl1)
        else
           hl1   = -zolf
           tem1  = 1.0 / sqrt(hl1)
           psih_unstable_full_gfs  = log(hl1) + .5 * tem1 + 1.386
        end if

   end function

!>\ingroup mynn_sfc
!! look-up table functions - or, if beyond -10 < z/L < 10, recalculate
   real(kind_phys) function psim_stable(zolf,psi_opt)
        !$acc routine seq
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psim_stable = psim_stab(nzol) + rzol*(psim_stab(nzol+1)-psim_stab(nzol))
        else
           if (psi_opt == 0) then
              psim_stable = psim_stable_full(zolf)
           else
              psim_stable = psim_stable_full_gfs(zolf)
           endif
        endif

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_stable(zolf,psi_opt)
        !$acc routine seq
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psih_stable = psih_stab(nzol) + rzol*(psih_stab(nzol+1)-psih_stab(nzol))
        else
           if (psi_opt == 0) then
              psih_stable = psih_stable_full(zolf)
           else
              psih_stable = psih_stable_full_gfs(zolf)
           endif
        endif

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psim_unstable(zolf,psi_opt)
        !$acc routine seq
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psim_unstable = psim_unstab(nzol) + rzol*(psim_unstab(nzol+1)-psim_unstab(nzol))
        else
           if (psi_opt == 0) then
              psim_unstable = psim_unstable_full(zolf)
           else
              psim_unstable = psim_unstable_full_gfs(zolf)
           endif
        endif

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_unstable(zolf,psi_opt)
        !$acc routine seq
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psih_unstable = psih_unstab(nzol) + rzol*(psih_unstab(nzol+1)-psih_unstab(nzol))
        else
           if (psi_opt == 0) then
              psih_unstable = psih_unstable_full(zolf)
           else
              psih_unstable = psih_unstable_full_gfs(zolf)
           endif
        endif

   end function
!========================================================================

END MODULE module_sf_mynnsfc
