!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group             John C. Warner   !
!    Licensed under a MIT/X style license              Neil K. Ganju   !
!    See License_ROMS.md                               Alexis Beudin   !
!================================================ Tarandeep S. Kalra ===
!                                                                      !
!  Sumerged Aquiatic Vegetation Model:                                 !
!                                                                      !
!  This module defines the model state variables, parameters, and I/O  !
!  indices. It includes its allocation and initialization routines.    !
!                                                                      !
!  Dimension Parameters:                                               !
!                                                                      !
!    NVEG               Number of kernel aquatic vegetation types      !
!    NVEGP              Number of vegetation properties                !
!                                                                      !
!  VEG(ng) Structure Variables: derived type T_VEG                     !
!                                                                      !
!    plant              Vegetation variable properties:                !
!                         plant(:,:,:,isDens) => density               !
!                         plant(:,:,:,isDiam) => diameter              !
!                         plant(:,:,:,isHght) => height                !
!                         plant(:,:,:,isThck) => thickness             !
!    ru_veg             3D U-momentum RHS term (accounting for all     !
!                         vegetation types)                            !
!    rv_veg             3D V-momentum RHS term (accounting for all     !
!                         vegetation types)                            !
!    ru_veg_loc         3D U-momentum RHS term (accounting for only    !
!                         local vegetation types)                      !
!    rv_veg_loc         3D V-momentum RHS term (accounting for only    !
!                         local vegetation types)                      !
!    step2d_uveg        2D U-momentum RHS term                         !
!    step2d_vveg        2D V-momentum RHS term                         !
#ifdef VEG_FLEX
!    bend               Bending for each vegetation                    !
!    Lveg               Effective blade length                         !
# endif
# ifdef VEG_TURB
!    tke_veg            Turbulent kinetic energy from vegetation       !
!    gls_veg            Length scale change from vegetation            !
#endif
#if defined VEG_SWAN_COUPLING && defined VEG_STREAMING
!    dissip_veg         Dissipation from Wave model due to vegetation  !
!    Cdwave_veg         Spectral Cd from Wave model due to vegetation  !
!    BWDXL_veg          Wave X-streaming effect due to vegetation      !
!    BWDYL_veg          Wave Y-streaming effect due to vegetation      !
#endif
#ifdef MARSH_WAVE_THRUST
!    marsh_mask         Marsh mask at cell center (RHO-points)         !
!    umask_marsh        March mask at u-face cell boundary             !
!    vmask_marsh        March mask at v-face cell boundary             !
!    Thrust_xi          Wave thrust on marsh xi-faces                  !
!    Thrust_eta         Wave thrust on marsh eta-faces                 !
!    Thrust_total       Total wave thrust magnitue on marsh cells      !
#  if defined MARSH_SED_EROSION
!    marsh_flux_out     Total sediment flux out from marsh cells       !
#  endif
#  if defined MARSH_RETREAT
!    marsh_retreat      Amount of marsh retreat                        !
#  endif
#  if defined MARSH_TIDAL_RANGE
!    zeta_max_rec       Marsh mean higher high water (MHHW)            !
!    zeta_min_rec       Marsh mean lower low water (MLLW)              !
!    marsh_tidal_range  Marsh mean tidal range (MHHW minus MLLW)       !
#  endif
#  if defined MARSH_VERT_GROWTH
!    marsh_high_water   Read or record mean high water                 !
!    marsh_low_water    Read mean low water                            !
!    marsh_biomass_peak Peak biomass on marsh                          !
!    marsh_vert_rate    Vertical rate of marsh growth (m/yr)           !
!    marsh_accret       Total accretion in marsh elevation (m)         !
#  endif
#endif
!                                                                      !
!  Vegetation Model Input Configuration Variables:                     !
!                                                                      !
!  CD_VEG               Drag coefficient from each vegtation type      !
!  E_VEG                Youngs modulus from each vegetation type       !
!  VEG_MASSDEN          Mass density from each vegetation type         !
!  VEGHMIXCOEF          Viscosity coefficient from vegetation boundary !
!                                                                      !
!  KFAC_MARSH           Marsh sediment erodibility coefficient         !
!  SCARP_HGHT           Absolute change in scarp height to convert     !
!                         marsh to open water cell                     !
!                                                                      !
!  PAR_FAC1             Marsh parabolic curve growth parameter 1       !
!  PAR_FAC2             Marsh parabolic curve growth parameter 2       !
!  TDAYS_MARSH_GROWTH   Growing number of days for marsh               !
!  MARSH_BULK_DENS      Bulk density for marsh organic sediment        !
!  NUGP                 Fraction of below ground biomass               !
!  BMAX                 Marsh peak biomass production                  !
!  CHIREF               Fraction of recalcitrant Carbon                !
!  ALPHA_PDENS          Growth parameter 1 for marsh plant density     !
!  BETA_PDENS           Growth parameter 2 for marsh plant density     !
!  ALPHA_PHGHT          Growth parameter 1 for marsh plant height      !
!  BETA_PHGHT           Growth parameter 2 for marsh plant height      !
!  ALPHA_PDIAM          Growth parameter 1 for marsh plant diameter    !
!  BETA_PDIAM           Growth parameter 2 for marsh plant diameter    !
!                                                                      !
!  Submerged Acuatic Plant Property Indices:                           !
!                                                                      !
!  isDens               Plant density, individuals per unit area (m2)  !
!  isDiam               Dominant plant mean diameter (m) per cell      !
!  isHght               Dominant plant mean height (m) per cell        !
!  isThck               Dominant plant mean thickness (m) per cell     !
!                                                                      !
!  Vegetation Model I/O Metadata Indices:                              !
!                                                                      !
!  idvprp(:)            Aquatic vegetation properties, 1:NVEGP         !
!  idTmfo(:)            Sediment flux out marsh cells, 1:NST classes   !
!                         (defined and allocated in mod_sediment)      !
!                                                                      !
!  idWdvg               Wave dissipation due to vegetation             !
!  idCdvg               Spectral drag froms waves and vegetation       !
!  idTims               Marsh mask: [1] marsh cell, [0] non-marsh cell !
!  idTtot               Total lateral wave thrust on marsh cells       !
!  idTmmr               Amount of marshland retreat                    !
!  idTmtr               Mean tidal range (MHHW-MLLW) for marsh growth  !
!  idTmhw               Mean higher high water (MHHW) in marsh cells   !
!  idTmlw               Mean lower low water (MLLW) in marsh cells     !
!  idTmbp               Below ground biomass for marsh growth          !
!  idTmvg               Marsh vertical growth rate                     !
!  idTmvt               Amount of marsh vertical growth, accretion     !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Beudin, A.,  Kalra, T.S., Ganju, N.K., Warner, J.C., 2017:          !
!    Development of a coupled wave-flow-vegetation interaction         !
!    model, Computers & Geosciences, Vol 100, 76-86,                   !
!    doi:10.1016/j.cageo.2016.12.010.                                  !
!                                                                      !
!  Kalra, T.S., Ganju, N.K., Aretxabaleta, A.L., Carr, J.A., Zafer,    !
!    D., Moriarty, J.M., 2021: Modeling Marsh Dynamics Using a 3-D     !
!    Coupled Wave-Flow-Sediment Model, Front. Mar. Sci., Vol 8,        !
!    doi:10.3389/fmars.2021.740921.                                    !
!                                                                      !
!======================================================================!
!                                                                      !
      USE mod_kinds
!
      implicit none
!
!-----------------------------------------------------------------------
!  Submerged aquatic vegetation model parameters and I/O identification
!  indices.
!-----------------------------------------------------------------------
!
      integer :: idWdvg    ! wave dissipation due to vegetation
      integer :: idCdvg    ! spectral drag due to waves and vegetation
      integer :: idTims    ! marsk mask: [1] marsh cell, [0] non-marsh
      integer :: idTtot    ! lateral wave thrust on marsh cells
      integer :: idTmmr    ! amount of marshland retreat
      integer :: idTmtr    ! mean tidal range leading marsh growth
      integer :: idTmhw    ! mean high water in marsh cells (MHWW)
      integer :: idTmlw    ! mean low  water in marsh cells
      integer :: idTmbp    ! below ground biomass for marsh growth
      integer :: idTmvg    ! marsh vertical growth rate
      integer :: idTmvt    ! amount of marsh vertical growth, accretion
!
!  State properties sub-indices for each aquatic dominant plant type.
!
      integer, parameter :: isDens = 1   ! density, individuals per area
      integer, parameter :: isDiam = 2   ! mean diameter
      integer, parameter :: isHght = 3   ! mean height
      integer, parameter :: isThck = 4   ! mean thickness
!
      integer, parameter :: NVEGP  = 4    ! number of plant properties
!                                           (set larger to sub-index)
      integer, dimension(NVEGP) :: idvprp ! aquatic plant properties
!
!  Number of vegetation types or groups. For example, seagrasses,
!  salt marshes, mangroves, and other herbaceous plants.
!
      integer :: NVEG
!
!-----------------------------------------------------------------------
!  Submerged Aquatic Vegetation Model Configuration Variables.
!-----------------------------------------------------------------------
!
# if defined MARSH_TIDAL_RANGE
      integer ::  NTIMES_MARSH
!
# endif
#if defined VEG_DRAG || defined VEG_BIOMASS
      real(r8), allocatable :: E_VEG(:,:)
      real(r8), allocatable :: CD_VEG(:,:)
      real(r8), allocatable :: VEG_MASSDENS(:,:)
      real(r8), allocatable :: VEGHMIXCOEF(:,:)
#endif
!
#ifdef MARSH_DYNAMICS
# if defined MARSH_SED_EROSION
      real(r8), allocatable :: KFAC_MARSH(:)
#  if defined MARSH_RETREAT
      real(r8), allocatable :: SCARP_HGHT(:)
#  endif
# endif
# if defined MARSH_VERT_GROWTH
      real(r8), allocatable :: PAR_FAC1(:), PAR_FAC2(:)
      real(r8), allocatable :: TDAYS_MARSH_GROWTH(:)
!     real(r8), allocatable :: MARSH_BULK_DENS(:)
      real(r8), allocatable :: NUGP(:)
      real(r8), allocatable :: BMAX(:)
      real(r8), allocatable :: CHIREF(:)
#  if defined MARSH_BIOMASS_VEG
      real(r8), allocatable :: ALPHA_PDENS(:), BETA_PDENS(:)
      real(r8), allocatable :: ALPHA_PHGHT(:), BETA_PHGHT(:)
      real(r8), allocatable :: ALPHA_PDIAM(:), BETA_PDIAM(:)
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Submerged Squatic Vegetation Model State Variables, TYPE T_VEG.
!-----------------------------------------------------------------------
!
      TYPE T_VEG
!
#if defined VEG_DRAG || defined VEG_BIOMASS
        real(r8), pointer :: plant(:,:,:,:)
#endif
#ifdef VEG_DRAG
        real(r8), pointer :: ru_veg(:,:,:)
        real(r8), pointer :: rv_veg(:,:,:)
        real(r8), pointer :: ru_loc_veg(:,:,:,:)
        real(r8), pointer :: rv_loc_veg(:,:,:,:)
        real(r8), pointer :: step2d_uveg(:,:)
        real(r8), pointer :: step2d_vveg(:,:)
        real(r8), pointer :: Lveg(:,:,:)
#endif
#ifdef VEG_FLEX
        real(r8), pointer :: bend(:,:,:)
#endif
#ifdef VEG_TURB
        real(r8), pointer :: tke_veg(:,:,:)
        real(r8), pointer :: gls_veg(:,:,:)
#endif
#ifdef VEG_HMIXING
        real(r8), pointer :: visc2d_r_veg(:,:)
        real(r8), pointer :: visc3d_r_veg(:,:,:)
#endif
#if defined VEG_SWAN_COUPLING && defined VEG_STREAMING
        real(r8), pointer :: dissip_veg(:,:)
        real(r8), pointer :: Cdwave_veg(:,:)
        real(r8), pointer :: BWDXL_veg(:,:,:)
        real(r8), pointer :: BWDYL_veg(:,:,:)
#endif
#ifdef MARSH_DYNAMICS
        real(r8), pointer :: marsh_mask(:,:)
# ifdef MARSH_WAVE_THRUST
        real(r8), pointer :: umask_marsh(:,:)
        real(r8), pointer :: vmask_marsh(:,:)
        real(r8), pointer :: Thrust_xi(:,:)
        real(r8), pointer :: Thrust_eta(:,:)
        real(r8), pointer :: Thrust_total(:,:)
# endif
# if defined MARSH_SED_EROSION
        real(r8), pointer :: marsh_flux_out(:,:,:)
# endif
# if defined MARSH_RETREAT
        real(r8), pointer :: marsh_retreat(:,:)
# endif
# if defined MARSH_STOCH
        real(r8), pointer :: marsh_stoch(:,:)
# endif
# if defined MARSH_TIDAL_RANGE
        real(r8), pointer :: zeta_max1(:,:)
        real(r8), pointer :: zeta_min1(:,:)
        real(r8), pointer :: zeta_max_rec(:,:,:)
        real(r8), pointer :: zeta_min_rec(:,:,:)
        real(r8), pointer :: marsh_tidal_range(:,:)
        real(r8) :: counter_loc_rl
# endif
# if defined MARSH_VERT_GROWTH
        real(r8), pointer :: marsh_high_water(:,:)
        real(r8), pointer :: marsh_low_water(:,:)
        real(r8), pointer :: marsh_biomass_peak(:,:)
        real(r8), pointer :: marsh_vert_rate(:,:)
        real(r8), pointer :: marsh_accret(:,:)
# endif
#endif
!
      END TYPE T_VEG
!
      TYPE (T_VEG), allocatable :: VEG(:)
!
      PUBLIC  :: allocate_vegetation
      PUBLIC  :: deallocate_vegetation
      PUBLIC  :: initialize_vegetation
!
      CONTAINS
!
      SUBROUTINE allocate_vegetation (ng, model, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  It allocates module variables.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param, ONLY : Dmem, iNLM, N, NST, Ngrids
!
!  Imported variable declarations.
!
      integer,           intent(in) :: ng, model
      integer, optional, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate Submerged Aquatic Vegetation Model Configuration Variables.
!-----------------------------------------------------------------------
!
      CONF_VARS : IF (model.eq.0) THEN

#if defined VEG_DRAG || defined VEG_BIOMASS
        IF (.not.allocated(E_VEG)) THEN
          allocate ( E_VEG(NVEG,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(NVEG*Ngrids,r8)
        END IF
!
        IF (.not.allocated(CD_VEG)) THEN
          allocate ( CD_VEG(NVEG,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(NVEG*Ngrids,r8)
        END IF
!
        IF (.not.allocated(VEG_MASSDENS)) THEN
          allocate ( VEG_MASSDENS(NVEG,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(NVEG*Ngrids,r8)
        END IF
!
        IF (.not.allocated(VEGHMIXCOEF)) THEN
          allocate ( VEGHMIXCOEF(NVEG,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(NVEG*Ngrids,r8)
        END IF
#endif

#ifdef MARSH_DYNAMICS
# if defined MARSH_SED_EROSION
!
        IF (.not.allocated(KFAC_MARSH)) THEN
          allocate ( KFAC_MARSH(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF
#  if defined MARSH_RETREAT
!
        IF (.not.allocated(SCARP_HGHT)) THEN
          allocate ( SCARP_HGHT(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
#  endif
# endif

# if defined MARSH_VERT_GROWTH
!
        IF (.not.allocated(PAR_FAC1)) THEN
          allocate ( PAR_FAC1(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(PAR_FAC2)) THEN
          allocate ( PAR_FAC2(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(TDAYS_MARSH_GROWTH)) THEN
          allocate ( TDAYS_MARSH_GROWTH(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      

!!      IF (.not.allocated(MARSH_BULK_DENS)) THEN
!!        allocate ( MARSH_BULK_DENS(Ngrids) )
!!        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
!!      END IF      
!
        IF (.not.allocated(NUGP)) THEN
          allocate ( NUGP(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(BMAX)) THEN
          allocate ( BMAX(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(CHIREF)) THEN
          allocate ( CHIREF(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      

#  if defined MARSH_BIOMASS_VEG
!
        IF (.not.allocated(ALPHA_PDENS)) THEN
          allocate ( ALPHA_PDENS(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(ALPHA_PHGHT)) THEN
          allocate ( ALPHA_PHGHT(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(ALPHA_PDIAM)) THEN
          allocate ( ALPHA_PDIAM(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(BETA_PDENS)) THEN
          allocate ( BETA_PDENS(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(BETA_PHGHT)) THEN
          allocate ( BETA_PHGHT(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
!
        IF (.not.allocated(BETA_PDIAM)) THEN
          allocate ( BETA_PDIAM(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        END IF      
#  endif
# endif
#endif
      END IF CONF_VARS
!
!-----------------------------------------------------------------------
!  Allocate Submerged Aquatic Vegetation Model State Variables.
!-----------------------------------------------------------------------
!
!  Allocate module structure.
!
      IF ((ng.eq.1).and.(model.eq.iNLM)) THEN
        allocate ( VEG(Ngrids) )!
      END IF
!
!  Allocate nonlinear state variables.
!
      NLM_VARS : IF (model.eq.iNLM) THEN
!
        size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)  ! horizontal array size
!
#if defined VEG_DRAG || defined VEG_BIOMASS
        allocate ( VEG(ng) % plant(LBi:UBi,LBj:UBj,NVEG,NVEGP) )
        Dmem(ng)=Dmem(ng)+REAL(NVEG*NVEGP,r8)*size2d
!
#endif
#ifdef VEG_DRAG
        allocate ( VEG(ng) % ru_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
!
        allocate ( VEG(ng) % rv_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
!
        allocate ( VEG(ng) % ru_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NVEG,r8)*size2d
!
        allocate ( VEG(ng) % rv_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NVEG,r8)*size2d
!
        allocate ( VEG(ng) % step2d_uveg(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % step2d_vveg(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % Lveg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

# ifdef VEG_FLEX
!
        allocate ( VEG(ng) % bend(LBi:UBi,LBj:UBj,NVEG) )
        Dmem(ng)=Dmem(ng)+REAL(NVEG,r8)*size2d
# endif
# ifdef VEG_HMIXING
!
        allocate ( VEG(ng) % visc2d_r_veg(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % visc3d_r_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
# endif
# ifdef VEG_TURB
!
        allocate ( VEG(ng) % tke_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
!
        allocate ( VEG(ng) % gls_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
# endif
#endif
#if defined VEG_SWAN_COUPLING && defined VEG_STREAMING
!
        allocate ( VEG(ng) % dissip_veg(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % Cdwave_veg(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % BWDXL_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
!
        allocate ( VEG(ng) % BWDYL_veg(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
#endif
#ifdef MARSH_DYNAMICS
!
        allocate ( VEG(ng) % marsh_mask(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
# ifdef MARSH_WAVE_THRUST
!
        allocate ( VEG(ng) % umask_marsh(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % vmask_marsh(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % Thrust_xi (LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % Thrust_eta(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % Thrust_total(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
# endif
# if defined MARSH_SED_EROSION
!
        allocate ( VEG(ng) % marsh_flux_out(LBi:UBi,LBj:UBj,NST) )
        Dmem(ng)=Dmem(ng)+REAL(NST,r8)*size2d
# endif
# if defined MARSH_RETREAT
!
        allocate ( VEG(ng) % marsh_retreat(LBi:UBi,LBj:UBj ) )
        Dmem(ng)=Dmem(ng)+size2d
# endif
# if defined MARSH_STOCH
!
        allocate ( VEG(ng) % marsh_stoch(LBi:UBi,LBj:UBj ) )
        Dmem(ng)=Dmem(ng)+size2d
# endif
# if defined MARSH_TIDAL_RANGE
!
        allocate ( VEG(ng) % zeta_max1(LBi:UBi,LBj:UBj ) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % zeta_min1(LBi:UBi,LBj:UBj ) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % zeta_max_rec(LBi:UBi,LBj:UBj,NTIMES_MARSH ) )
        Dmem(ng)=Dmem(ng)+REAL(NTIMES_MARSH,r8)*size2d
!
        allocate ( VEG(ng) % zeta_min_rec(LBi:UBi,LBj:UBj,NTIMES_MARSH ) )
        Dmem(ng)=Dmem(ng)+REAL(NTIMES_MARSH,r8)*size2d
!
        allocate ( VEG(ng) % marsh_tidal_range(LBi:UBi,LBj:UBj))
        Dmem(ng)=Dmem(ng)+size2d
# endif
# if defined MARSH_VERT_GROWTH
!
        allocate ( VEG(ng) % marsh_high_water(LBi:UBi,LBj:UBj))
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % marsh_low_water(LBi:UBi,LBj:UBj))
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % marsh_biomass_peak(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % marsh_vert_rate(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( VEG(ng) % marsh_accret(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
# endif
#endif
      END IF NLM_VARS
!
      RETURN
      END SUBROUTINE allocate_vegetation
!
      SUBROUTINE deallocate_vegetation (ng)
!
!=======================================================================
!                                                                      !
!  It deallocates module variables.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param, ONLY : Ngrids
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!-----------------------------------------------------------------------
!  Deallocate variables and structure
!-----------------------------------------------------------------------
!
!  Deallocate configuration variables.
!
      IF (ng.eq.Ngrids) THEN

#if defined VEG_DRAG || defined VEG_BIOMASS
        IF (allocated(E_VEG))         deallocate ( E_VEG )
        IF (allocated(CD_VEG))        deallocate ( CD_VEG )
        IF (allocated(VEG_MASSDENS))  deallocate ( VEG_MASSDENS )
        IF (allocated(VEGHMIXCOEF))   deallocate ( VEGHMIXCOEF )
#endif

#ifdef MARSH_DYNAMICS
# if defined MARSH_SED_EROSION
        IF (allocated(KFAC_MARSH))    deallocate ( KFAC_MARSH )
#  if defined MARSH_RETREAT
        IF (allocated(SCARP_HGHT))    deallocate ( SCARP_HGHT )
#  endif
# endif

# if defined MARSH_VERT_GROWTH
        IF (allocated(PAR_FAC1))      deallocate ( PAR_FAC1 )
        IF (allocated(PAR_FAC2))      deallocate ( PAR_FAC2 )
        IF (allocated(NUGP))          deallocate ( NUGP )
        IF (allocated(BMAX))          deallocate ( BMAX )
        IF (allocated(CHIREF))        deallocate ( CHIREF )
#  if defined MARSH_BIOMASS_VEG
        IF (allocated(ALPHA_PDENS))   deallocate ( ALPHA_PDENS )
        IF (allocated(ALPHA_PHGHT))   deallocate ( ALPHA_PHGHT )
        IF (allocated(ALPHA_PDIAM))   deallocate ( ALPHA_PDIAM )
        IF (allocated(BETA_PDENS))    deallocate ( BETA_PDENS )
        IF (allocated(BETA_PHGHT))    deallocate ( BETA_PHGHT )
        IF (allocated(BETA_PDIAM))    deallocate ( BETA_PDIAM )
#  endif
# endif
#endif
      END IF
!
!  Deallocate derived-type VEG structure.
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(VEG)) deallocate ( VEG )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_vegetation
!
      SUBROUTINE initialize_vegetation (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param,  ONLY : BOUNDS, Dmem, iNLM, N, NST
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ip, it, iveg, j, k
!
      real(r8), parameter :: IniVal = 0.0_r8
!
#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
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
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize vegetation structure variables.
!-----------------------------------------------------------------------
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
#if defined VEG_DRAG || defined VEG_BIOMASS
        DO ip=1,NVEGP
          DO iveg=1,NVEG
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                VEG(ng) % plant(i,j,iveg,ip) = IniVal
              END DO
            END DO
          END DO
        END DO
!
#endif
#ifdef VEG_DRAG
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            DO k=1,N(ng)
              VEG(ng) % ru_veg(i,j,k)  = IniVal
              VEG(ng) % rv_veg(i,j,k)  = IniVal
              VEG(ng) % Lveg(i,j,k)    = IniVal
              DO iveg=1,NVEG
                VEG(ng) % ru_loc_veg(i,j,k,iveg) = IniVal
                VEG(ng) % rv_loc_veg(i,j,k,iveg) = IniVal
              END DO               
            END DO
            VEG(ng) % step2d_uveg(i,j) = IniVal
            VEG(ng) % step2d_vveg(i,j) = IniVal
          END DO
        END DO
#endif

#ifdef VEG_FLEX
        DO iveg=1,NVEG
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % bend(i,j,iveg) = IniVal
            END DO
          END DO
        END DO
#endif

#ifdef VEG_TURB
!
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % tke_veg(i,j,k) = IniVal
              VEG(ng) % gls_veg(i,j,k) = IniVal
            END DO
          END DO
        END DO
#endif

#if defined VEG_SWAN_COUPLING && defined VEG_STREAMING
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            VEG(ng) % dissip_veg(i,j)    = IniVal
            VEG(ng) % Cdwave_veg(i,j)    = IniVal
            DO k=1,N(ng)
              VEG(ng) % BWDXL_veg(i,j,k) = IniVal
              VEG(ng) % BWDYL_veg(i,j,k) = IniVal
            END DO
          END DO
        END DO
#endif

#ifdef MARSH_DYNAMICS
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            VEG(ng) % marsh_mask(i,j)         = IniVal
# ifdef MARSH_WAVE_THRUST
            VEG(ng) % umask_marsh(i,j)        = IniVal
            VEG(ng) % vmask_marsh(i,j)        = IniVal
            VEG(ng) % Thrust_xi (i,j)         = IniVal
            VEG(ng) % Thrust_eta(i,j)         = IniVal
            VEG(ng) % Thrust_total(i,j)       = IniVal
#  ifdef MARSH_SED_EROSION
            DO k=1,NST
              VEG(ng) % marsh_flux_out(i,j,k) = IniVal
            END DO
#  endif
# endif
# if defined MARSH_RETREAT
            VEG(ng) % marsh_retreat(i,j)      = IniVal
# endif
# if defined MARSH_STOCH
            VEG(ng) % marsh_stoch(i,j)        = IniVal
# endif 
# if defined MARSH_TIDAL_RANGE
            DO it=1,NTIMES_MARSH
              VEG(ng) % zeta_max_rec(i,j,it)  = IniVal
              VEG(ng) % zeta_min_rec(i,j,it)  = IniVal
            END DO
            VEG(ng) % zeta_max1(i,j)          = -10.0_r8
            VEG(ng) % zeta_min1(i,j)          =  10.0_r8
            VEG(ng) % marsh_tidal_range(i,j)  = IniVal
# endif
# if defined MARSH_VERT_GROWTH
            VEG(ng) % marsh_high_water(i,j)   = IniVal
            VEG(ng) % marsh_low_water(i,j)    = IniVal
            VEG(ng) % marsh_biomass_peak(i,j) = IniVal
            VEG(ng) % marsh_vert_rate(i,j)    = IniVal
            VEG(ng) % marsh_accret(i,j)       = IniVal
# endif
          END DO
        END DO
# if defined MARSH_TIDAL_RANGE
        VEG(ng) % counter_loc_rl = 1.0_r8
# endif
#endif
      END IF
!
      RETURN
      END SUBROUTINE initialize_vegetation
