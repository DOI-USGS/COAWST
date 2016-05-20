!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Nemuro ecosystem model:                              !
!                                                                      !
!  AlphaPS     Small Phytoplankton photochemical reaction coefficient: !
!                initial slope (low light) of the P-I curve,           !
!                [1/(W/m2) 1/day].                                     !
!  AlphaPL     Large Phytoplankton photochemical reaction coefficient: !
!                initial slope (low light) of the P-I curve,           !
!                [1/(W/m2) 1/day].                                     !
!  AlphaZL     Large Zooplankton assimilation efficiency,              !
!                [nondimemsional].                                     !
!  AlphaZP     Predator Zooplankton assimilation efficiency,           !
!                [nondimemsional].                                     !
!  AlphaZS     Small Zooplankton assimilation efficiency,              !
!                [nondimemsional].                                     !
!  AttPL       Light attenuation due to Large Phytoplankton, self-     !
!                shading coefficient, [m2/millimole_N].                !
!  AttPS       Light attenuation due to Small Phytoplankton, self-     !
!                shading coefficient, [m2/millimole_N].                !
!  AttSW       Light attenuation due to sea water, [1/m].              !
!  BetaPL      Large Phytoplankton photoinhibition coefficient,        !
!                [1/(W/m2) 1/day].                                     !
!  BetaPS      Small Phytoplankton photoinhibition coefficient,        !
!                [1/(W/m2) 1/day].                                     !
!  BetaZL      Large Zooplankton growth efficiency [nondimensional].   !
!  BetaZP      Predator Zooplankton growth efficiency [nondimensional].!
!  BetaZS      Small Zooplankton growth efficiency [nondimensional].   !
!  BioIter     Maximum number of iterations to achieve convergence of  !
!                the nonlinear solution.                               !
!  GammaL      Large Phytoplankton ratio of extracellular excretion to !
!                photosynthesis [nondimensional].                      !
!  GammaS      Small Phytoplankton ratio of extracellular excretion to !
!                photosynthesis [nondimensional].                      !
!  GRmaxLpl    Large Zooplankton maximum grazing rate on Large         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxLps    Large Zooplankton maximum grazing rate on Small         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxLzs    Small Zooplankton maximum grazing rate on Small         !
!                Zooplankton at 0 Celsius, [1/day].                    !
!  GRmaxPpl    Predator Zooplankton maximum grazing rate on Large      !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxPzl    Predator Zooplankton maximum grazing rate on Large      !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxPzs    Predator Zooplankton maximum grazing rate on Small      !
!                Zooplankton at 0 Celsius, [1/day].                    !
!  GRmaxSps    Small Zooplankton maximum grazing rate on Small         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxSpl    Small Zooplankton maximum grazing rate on Large         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  KD2N        Temperature coefficient for DON to NH4 decomposition,   !
!                [1/Celsius].                                          !
!  KGppL       Large Phytoplankton temperature coefficient for         !
!                photosynthetic rate, [1/Celsius].                     !
!  KGppS       Small Phytoplankton temperature coefficient for         !
!                photosynthetic rate, [1/Celsius].                     !
!  KGraL       Large Zooplankton temperature coefficient for grazing,  !
!                [1/Celsius].                                          !
!  KGraP       Predator Zooplankton temperature coefficient for        !
!                grazing,[1/Celsius].                                  !
!  KGraS       Small Zooplankton temperature coefficient for grazing,  !
!                [1/Celsius].                                          !
!  KMorPL      Large Phytoplankton temperature coefficient for         !
!                mortality, [1/Celsius].                               !
!  KMorPS      Small Phytoplankton temperature coefficient for         !
!                mortality, [1/Celsius].                               !
!  KMorZL      Large Zooplankton temperature coefficient for           !
!                mortality, [1/Celsius].                               !
!  KMorZP      Predator Zooplankton temperature coefficient for        !
!                mortality, [1/Celsius].                               !
!  KMorZS      Small Zooplankton temperature coefficient for           !
!                mortality, [1/Celsius].                               !
!  KNit        Temperature coefficient for nitrification (NH4 to NO3)  !
!                decomposition, [1/Celsius].                           !
!  KNH4L       Large Phytoplankton half satuation constant for NH4,    !
!                [millimole_N/m3].                                     !
!  KNH4S       Small Phytoplankton half satuation constant for NH4,    !
!                [millimole_N/m3].                                     !
!  KNO3L       Large Phytoplankton half satuation constant for NO3,    !
!                [millimole_N/m3].                                     !
!  KNO3S       Small Phytoplankton half satuation constant for NO3,    !
!                [millimole_N/m3].                                     !
!  KO2S        Temperature coefficient for Opal to SiOH4 decomposition,!
!                [1/Celsius].                                          !
!  KP2D        Temperature coefficient for PON to DON decomposition,   !
!                [1/Celsius].                                          !
!  KP2N        Temperature coefficient for PON to NH4 decomposition,   !
!                [1/Celsius].                                          !
!  KPL2ZS      Small Zooplankton half-saturation coefficient for       !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPL2ZL      Large Zooplankton half-saturation coefficient for       !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPL2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPS2ZL      Larg Zooplankton half-saturation coefficient for        !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KPS2ZS      Small Zooplankton half-saturation coefficient for       !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KResPL      Large Phytoplankton temperature coefficient for         !
!                respiration, [1/Celsius].                             !
!  KResPS      Small Phytoplankton temperature coefficient for         !
!                respiration, [1/Celsius].                             !
!  KSiL        Large Phytoplankton half satuation constant for SiOH4,  !
!                [millimole_Si/m3].                                    !
!  KZL2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Large Zooplankton [millimole_N/m3]^2.    !
!  KZS2ZL      Large Zooplankton half-saturation coefficient for       !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KZS2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Small Zooplankton [millimole_N/m3]^2.    !
!  LamL        Large Zooplankton Ivlev constant, [m3/millimole_N].     !
!  LamP        Predator Zooplankton Ivlev constant, [m3/millimole_N].  !
!  LamS        Small Zooplankton Ivlev constant, [m3/millimole_N].     !
!  MorPL0      Large Phytoplankton mortality rate at 0 Celsius,        !
!                [m3/millimole_N 1/day].                               !
!  MorPS0      Small Phytoplankton mortality rate at 0 Celsius,        !
!                [m3/millimole_N 1/day].                               !
!  MorZL0      Large Zooplankton mortality rate at 0 Celsius,          !
!                [m3/millimole_N 1/day].                               !
!  MorZP0      Predator Zooplankton mortality rate at 0 Celsius,       !
!                [m3/millimole_N 1/day].                               !
!  MorZS0      Small Zooplankton mortality rate at 0 Celsius,          !
!                [m3/millimole_N 1/day].                               !
!  Nit0        Nitrification (NH4 to NO3) rate at 0 Celsius, [1/day].  !
!  PARfrac     Fraction of shortwave radiation that is available for   !
!                photosyntesis [nondimensional].                       !
!  PL2ZSstar   Small Zooplankton threshold value for grazing on        !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PL2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PL2ZPstar   Predator Zooplankton threshold value for grazing on     !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PS2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  PS2ZSstar   Small Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  PusaiL      Large Phytoplankton Ammonium inhibition coefficient,    !
!                [m3/millimole_N].                                     !
!  PusaiPL     Predator Zooplankton grazing on Large Phytoplankton     !
!                inhibition coefficient, [m3/millimole_N].             !
!  PusaiS      Small Phytoplankton Ammonium inhibition coefficient,    !
!                [m3/millimole_N].                                     !
!  PusaiZS     Predator Zooplankton grazing on Small Zooplankton       !
!                inhibition coefficient, [m3/millimole_N].             !
!  ResPL0      Large Phytoplankton respiration rate at 0 Celsius,      !
!                [1/day].                                              !
!  ResPS0      Small Phytoplankton respiration rate at 0 Celsius,      !
!                [1/day].                                              !
!  RSiN        Si:N ratio [millimole_Si/millimole_N].                  !
!  setVOpal    Opal Settling (sinking) velocity [m/day].               !
!  setVPON     PON Settling (sinking) velocity [m/day].                !
!  VD2N0       DON to NH4 decomposition rate at 0 Celsius, [1/day].    !
!  VmaxL       Maximum Large Phytoplankton photosynthetic rate [1/day] !
!                in the absence of photoinhibition under optimal light.!
!  VmaxS       Maximum Small Phytoplankton photosynthetic rate [1/day] !
!                in the absence of photoinhibition under optimal light.!
!  VO2S0       Opal to Silicate decomposition rate at 0 Celsius,       !
!                [1/day].                                              !
!  VP2D0       PON to DON decomposition rate at 0 Celsius, [1/day].    !
!  VP2N0       PON to NH4 decomposition rate at 0 Celsius, [1/day].    !
!  ZL2ZPstar   Small Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  ZS2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Small Zooplankton, [millimole_N/m3].                  !
!  ZS2ZPstar   Predator Zooplankton threshold value for grazing on     !
!                Small Zooplankton, [millimole_N/m3].                  !
!                                                                      !
! Parameters for iron limitation                                       !
!  T_Fe      Iron uptake timescale, [day].                             !
!  A_Fe      Empirical Fe:C power, [-].                                !
!  B_Fe      Empirical Fe:C coefficient, [1/M-C].                      !
!  SK_FeC    Small phytoplankton Fe:C at F=0.5, [muM-Fe/M-C].          !
!  LK_FeC    Large phytoplankton Fe:C at F=0.5, [muM-Fe/M-C].          !
!  FeRR      Fe remineralization rate, [1/day].                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iLphy                  ! Large Phytoplankton biomass
      integer :: iSphy                  ! Small Phytoplankton biomass
      integer :: iLzoo                  ! Large Zooplankton biomass
      integer :: iSzoo                  ! Small Zooplankton biomass
      integer :: iPzoo                  ! Predator Zooplankton biomass
      integer :: iNO3_                  ! Nitrate concentration
      integer :: iNH4_                  ! Ammonium concentration
      integer :: iPON_                  ! Particulate Organic Nitrogen
      integer :: iDON_                  ! Dissolved Organic Nitrogen
      integer :: iSiOH                  ! Silicate concentration
      integer :: iopal                  ! Particulate organic silica
# ifdef IRON_LIMIT
      integer :: iFeSp                  ! Small phytoplankton iron
      integer :: iFeLp                  ! Large phytoplankton iron
      integer :: iFeD_                  ! Available dissolved iron
# endif
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: AlphaPL(:)            ! 1/(W/m2) 1/day
      real(r8), allocatable :: AlphaPS(:)            ! 1/(W/m2) 1/day
      real(r8), allocatable :: AlphaZL(:)            ! nondimensional
      real(r8), allocatable :: AlphaZP(:)            ! nondimensional
      real(r8), allocatable :: AlphaZS(:)            ! nondimensional
      real(r8), allocatable :: AttPL(:)              ! m2/mmole_N
      real(r8), allocatable :: AttPS(:)              ! m2/mmole_N
      real(r8), allocatable :: AttSW(:)              ! 1/m
      real(r8), allocatable :: BetaPL(:)             ! 1/(W/m2) 1/day
      real(r8), allocatable :: BetaPS(:)             ! 1/(W/m2) 1/day
      real(r8), allocatable :: BetaZS(:)             ! nondimensional
      real(r8), allocatable :: BetaZL(:)             ! nondimensional
      real(r8), allocatable :: BetaZP(:)             ! nondimensional
      real(r8), allocatable :: GammaL(:)             ! nondimensional
      real(r8), allocatable :: GammaS(:)             ! nondimensional
      real(r8), allocatable :: GRmaxLpl(:)           ! 1/day
      real(r8), allocatable :: GRmaxLps(:)           ! 1/day
      real(r8), allocatable :: GRmaxLzs(:)           ! 1/day
      real(r8), allocatable :: GRmaxPpl(:)           ! 1/day
      real(r8), allocatable :: GRmaxPzl(:)           ! 1/day
      real(r8), allocatable :: GRmaxPzs(:)           ! 1/day
      real(r8), allocatable :: GRmaxSpl(:)           ! 1/day
      real(r8), allocatable :: GRmaxSps(:)           ! 1/day
      real(r8), allocatable :: KD2N(:)               ! 1/Celsius
      real(r8), allocatable :: KGppL(:)              ! 1/Celsius
      real(r8), allocatable :: KGppS(:)              ! 1/Celsius
      real(r8), allocatable :: KGraL(:)              ! 1/Celsius
      real(r8), allocatable :: KGraP(:)              ! 1/Celsius
      real(r8), allocatable :: KGraS(:)              ! 1/Celsius
      real(r8), allocatable :: KMorPL(:)             ! 1/Celsius
      real(r8), allocatable :: KMorPS(:)             ! 1/Celsius
      real(r8), allocatable :: KMorZL(:)             ! 1/Celsius
      real(r8), allocatable :: KMorZP(:)             ! 1/Celsius
      real(r8), allocatable :: KMorZS(:)             ! 1/Celsius
      real(r8), allocatable :: KNH4L(:)              ! mmole_N/m3
      real(r8), allocatable :: KNH4S(:)              ! mmole_N/m3
      real(r8), allocatable :: KNit(:)               ! 1/Celsius
      real(r8), allocatable :: KNO3L(:)              ! mmole_N/m3
      real(r8), allocatable :: KNO3S(:)              ! mmole_N/m3
      real(r8), allocatable :: KO2S(:)               ! 1/Celsius
      real(r8), allocatable :: KP2D(:)               ! 1/Celsius
      real(r8), allocatable :: KP2N(:)               ! 1/Celsius
      real(r8), allocatable :: KPL2ZL(:)             ! mmole_N/m3
      real(r8), allocatable :: KPL2ZS(:)             ! mmole_N/m3
      real(r8), allocatable :: KPS2ZL(:)             ! mmole_N/m3
      real(r8), allocatable :: KPS2ZS(:)             ! mmole_N/m3
      real(r8), allocatable :: KPL2ZP(:)             ! mmole_N/m3
      real(r8), allocatable :: KResPL(:)             ! 1/Celsius
      real(r8), allocatable :: KResPS(:)             ! 1/Celsius
      real(r8), allocatable :: KSiL(:)               ! mmole_Si/m3
      real(r8), allocatable :: KZL2ZP(:)             ! mmole_N/m3
      real(r8), allocatable :: KZS2ZL(:)             ! mmole_N/m3
      real(r8), allocatable :: KZS2ZP(:)             ! mmole_N/m3
      real(r8), allocatable :: LamL(:)               ! m3/mmole_N
      real(r8), allocatable :: LamP(:)               ! m3/mmole_N
      real(r8), allocatable :: LamS(:)               ! m3/mmole_N
      real(r8), allocatable :: MorPL0(:)             ! m3/mmole_N/day
      real(r8), allocatable :: MorPS0(:)             ! m3/mmole_N/day
      real(r8), allocatable :: MorZL0(:)             ! m3/mmole_N 1/day
      real(r8), allocatable :: MorZP0(:)             ! m3/mmole_N 1/day
      real(r8), allocatable :: MorZS0(:)             ! m3/mmole_N 1/day
      real(r8), allocatable :: Nit0(:)               ! 1/day
      real(r8), allocatable :: PARfrac(:)            ! nondimensional
      real(r8), allocatable :: PusaiL(:)             ! m3/mmole_N
      real(r8), allocatable :: PusaiPL(:)            ! m3/mmole_N
      real(r8), allocatable :: PusaiS(:)             ! m3/mmole_N
      real(r8), allocatable :: PusaiZS(:)            ! m3/mmole_N
      real(r8), allocatable :: PL2ZLstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PL2ZPstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PL2ZSstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PS2ZLstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PS2ZSstar(:)          ! mmole_N/m3
      real(r8), allocatable :: ResPL0(:)             ! 1/day
      real(r8), allocatable :: ResPS0(:)             ! 1/day
      real(r8), allocatable :: RSiN(:)               ! mmole_Si/mmole_N
      real(r8), allocatable :: setVOpal(:)           ! m/day
      real(r8), allocatable :: setVPON(:)            ! m/day
      real(r8), allocatable :: VD2N0(:)              ! 1/day
      real(r8), allocatable :: VmaxL(:)              ! 1/day
      real(r8), allocatable :: VmaxS(:)              ! 1/day
      real(r8), allocatable :: VO2S0(:)              ! 1/day
      real(r8), allocatable :: VP2D0(:)              ! 1/day
      real(r8), allocatable :: VP2N0(:)              ! 1/day
      real(r8), allocatable :: ZL2ZPstar(:)          ! mmole_N/m3
      real(r8), allocatable :: ZS2ZLstar(:)          ! mmole_N/m3
      real(r8), allocatable :: ZS2ZPstar(:)          ! mmole_N/m3
#ifdef IRON_LIMIT
      real(r8), allocatable :: T_Fe(:)               ! day
      real(r8), allocatable :: A_Fe(:)               ! nondimensional
      real(r8), allocatable :: B_Fe(:)               ! 1/M-C
      real(r8), allocatable :: SK_FeC(:)             ! muM-Fe/M-C
      real(r8), allocatable :: LK_FeC(:)             ! muM-Fe/M-C
      real(r8), allocatable :: FeRR(:)               ! 1/day
#endif

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
#  ifdef IRON_LIMIT
      NBT = 14
#  else
      NBT = 11
#  endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AlphaPL)) THEN
        allocate ( AlphaPL(Ngrids) )
      END IF
      IF (.not.allocated(AlphaPS)) THEN
        allocate ( AlphaPS(Ngrids) )
      END IF
      IF (.not.allocated(AlphaZL)) THEN
        allocate ( AlphaZL(Ngrids) )
      END IF
      IF (.not.allocated(AlphaZP)) THEN
        allocate ( AlphaZP(Ngrids) )
      END IF
      IF (.not.allocated(AlphaZS)) THEN
        allocate ( AlphaZS(Ngrids) )
      END IF
      IF (.not.allocated(AttPL)) THEN
        allocate ( AttPL(Ngrids) )
      END IF
      IF (.not.allocated(AttPS)) THEN
        allocate ( AttPS(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(BetaPL)) THEN
        allocate ( BetaPL(Ngrids) )
      END IF
      IF (.not.allocated(BetaPS)) THEN
        allocate ( BetaPS(Ngrids) )
      END IF
      IF (.not.allocated(BetaZS)) THEN
        allocate ( BetaZS(Ngrids) )
      END IF
      IF (.not.allocated(BetaZL)) THEN
        allocate ( BetaZL(Ngrids) )
      END IF
      IF (.not.allocated(BetaZP)) THEN
        allocate ( BetaZP(Ngrids) )
      END IF
      IF (.not.allocated(GammaL)) THEN
        allocate ( GammaL(Ngrids) )
      END IF
      IF (.not.allocated(GammaS)) THEN
        allocate ( GammaS(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxLpl)) THEN
        allocate ( GRmaxLpl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxLps)) THEN
        allocate ( GRmaxLps(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxLzs)) THEN
        allocate ( GRmaxLzs(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxPpl)) THEN
        allocate ( GRmaxPpl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxPzl)) THEN
        allocate ( GRmaxPzl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxPzs)) THEN
        allocate ( GRmaxPzs(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxSpl)) THEN
        allocate ( GRmaxSpl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxSps)) THEN
        allocate ( GRmaxSps(Ngrids) )
      END IF
      IF (.not.allocated(KD2N)) THEN
        allocate ( KD2N(Ngrids) )
      END IF
      IF (.not.allocated(KGppL)) THEN
        allocate ( KGppL(Ngrids) )
      END IF
      IF (.not.allocated(KGppS)) THEN
        allocate ( KGppS(Ngrids) )
      END IF
      IF (.not.allocated(KGraL)) THEN
        allocate ( KGraL(Ngrids) )
      END IF
      IF (.not.allocated(KGraP)) THEN
        allocate ( KGraP(Ngrids) )
      END IF
      IF (.not.allocated(KGraS)) THEN
        allocate ( KGraS(Ngrids) )
      END IF
      IF (.not.allocated(KMorPL)) THEN
        allocate ( KMorPL(Ngrids) )
      END IF
      IF (.not.allocated(KMorPS)) THEN
        allocate ( KMorPS(Ngrids) )
      END IF
      IF (.not.allocated(KMorZL)) THEN
        allocate ( KMorZL(Ngrids) )
      END IF
      IF (.not.allocated(KMorZP)) THEN
        allocate ( KMorZP(Ngrids) )
      END IF
      IF (.not.allocated(KMorZS)) THEN
        allocate ( KMorZS(Ngrids) )
      END IF
      IF (.not.allocated(KNH4L)) THEN
        allocate ( KNH4L(Ngrids) )
      END IF
      IF (.not.allocated(KNH4S)) THEN
        allocate ( KNH4S(Ngrids) )
      END IF
      IF (.not.allocated(KNit)) THEN
        allocate ( KNit(Ngrids) )
      END IF
      IF (.not.allocated(KNO3L)) THEN
        allocate ( KNO3L(Ngrids) )
      END IF
      IF (.not.allocated(KNO3S)) THEN
        allocate ( KNO3S(Ngrids) )
      END IF
      IF (.not.allocated(KO2S)) THEN
        allocate ( KO2S(Ngrids) )
      END IF
      IF (.not.allocated(KP2D)) THEN
        allocate ( KP2D(Ngrids) )
      END IF
      IF (.not.allocated(KP2N)) THEN
        allocate ( KP2N(Ngrids) )
      END IF
      IF (.not.allocated(KPL2ZL)) THEN
        allocate ( KPL2ZL(Ngrids) )
      END IF
      IF (.not.allocated(KPL2ZS)) THEN
        allocate ( KPL2ZS(Ngrids) )
      END IF
      IF (.not.allocated(KPS2ZL)) THEN
        allocate ( KPS2ZL(Ngrids) )
      END IF
      IF (.not.allocated(KPS2ZS)) THEN
        allocate ( KPS2ZS(Ngrids) )
      END IF
      IF (.not.allocated(KPL2ZP)) THEN
        allocate ( KPL2ZP(Ngrids) )
      END IF
      IF (.not.allocated(KResPL)) THEN
        allocate ( KResPL(Ngrids) )
      END IF
      IF (.not.allocated(KResPS)) THEN
        allocate ( KResPS(Ngrids) )
      END IF
      IF (.not.allocated(KSiL)) THEN
        allocate ( KSiL(Ngrids) )
      END IF
      IF (.not.allocated(KZL2ZP)) THEN
        allocate ( KZL2ZP(Ngrids) )
      END IF
      IF (.not.allocated(KZS2ZL)) THEN
        allocate ( KZS2ZL(Ngrids) )
      END IF
      IF (.not.allocated(KZS2ZP)) THEN
        allocate ( KZS2ZP(Ngrids) )
      END IF
      IF (.not.allocated(LamL)) THEN
        allocate ( LamL(Ngrids) )
      END IF
      IF (.not.allocated(LamP)) THEN
        allocate ( LamP(Ngrids) )
      END IF
      IF (.not.allocated(LamS)) THEN
        allocate ( LamS(Ngrids) )
      END IF
      IF (.not.allocated(MorPL0)) THEN
        allocate ( MorPL0(Ngrids) )
      END IF
      IF (.not.allocated(MorPS0)) THEN
        allocate ( MorPS0(Ngrids) )
      END IF
      IF (.not.allocated(MorZL0)) THEN
        allocate ( MorZL0(Ngrids) )
      END IF
      IF (.not.allocated(MorZP0)) THEN
        allocate ( MorZP0(Ngrids) )
      END IF
      IF (.not.allocated(MorZS0)) THEN
        allocate ( MorZS0(Ngrids) )
      END IF
      IF (.not.allocated(Nit0)) THEN
        allocate ( Nit0(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(PusaiL)) THEN
        allocate ( PusaiL(Ngrids) )
      END IF
      IF (.not.allocated(PusaiPL)) THEN
        allocate ( PusaiPL(Ngrids) )
      END IF
      IF (.not.allocated(PusaiS)) THEN
        allocate ( PusaiS(Ngrids) )
      END IF
      IF (.not.allocated(PusaiZS)) THEN
        allocate ( PusaiZS(Ngrids) )
      END IF
      IF (.not.allocated(PL2ZLstar)) THEN
        allocate ( PL2ZLstar(Ngrids) )
      END IF
      IF (.not.allocated(PL2ZPstar)) THEN
        allocate ( PL2ZPstar(Ngrids) )
      END IF
      IF (.not.allocated(PL2ZSstar)) THEN
        allocate ( PL2ZSstar(Ngrids) )
      END IF
      IF (.not.allocated(PS2ZLstar)) THEN
        allocate ( PS2ZLstar(Ngrids) )
      END IF
      IF (.not.allocated(PS2ZSstar)) THEN
        allocate ( PS2ZSstar(Ngrids) )
      END IF
      IF (.not.allocated(ResPL0)) THEN
        allocate ( ResPL0(Ngrids) )
      END IF
      IF (.not.allocated(ResPS0)) THEN
        allocate ( ResPS0(Ngrids) )
      END IF
      IF (.not.allocated(RSiN)) THEN
        allocate ( RSiN(Ngrids) )
      END IF
      IF (.not.allocated(setVOpal)) THEN
        allocate ( setVOpal(Ngrids) )
      END IF
      IF (.not.allocated(setVPON)) THEN
        allocate ( setVPON(Ngrids) )
      END IF
      IF (.not.allocated(VD2N0)) THEN
        allocate ( VD2N0(Ngrids) )
      END IF
      IF (.not.allocated(VmaxL)) THEN
        allocate ( VmaxL(Ngrids) )
      END IF
      IF (.not.allocated(VmaxS)) THEN
        allocate ( VmaxS(Ngrids) )
      END IF
      IF (.not.allocated(VO2S0)) THEN
        allocate ( VO2S0(Ngrids) )
      END IF
      IF (.not.allocated(VP2D0)) THEN
        allocate ( VP2D0(Ngrids) )
      END IF
      IF (.not.allocated(VP2N0)) THEN
        allocate ( VP2N0(Ngrids) )
      END IF
      IF (.not.allocated(ZL2ZPstar)) THEN
        allocate ( ZL2ZPstar(Ngrids) )
      END IF
      IF (.not.allocated(ZS2ZLstar)) THEN
        allocate ( ZS2ZLstar(Ngrids) )
      END IF
      IF (.not.allocated(ZS2ZPstar)) THEN
        allocate ( ZS2ZPstar(Ngrids) )
      END IF
#ifdef IRON_LIMIT
      IF (.not.allocated(T_Fe)) THEN
        allocate ( T_Fe(Ngrids) )
      END IF
      IF (.not.allocated(A_Fe)) THEN
        allocate ( A_Fe(Ngrids) )
      END IF
      IF (.not.allocated(B_Fe)) THEN
        allocate ( B_Fe(Ngrids) )
      END IF
      IF (.not.allocated(SK_FeC)) THEN
        allocate ( SK_FeC(Ngrids) )
      END IF
      IF (.not.allocated(LK_FeC)) THEN
        allocate ( LK_FeC(Ngrids) )
      END IF
      IF (.not.allocated(FeRR)) THEN
        allocate ( FeRR(Ngrids) )
      END IF
#endif
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iSphy=ic+1
      iLphy=ic+2
      iSzoo=ic+3
      iLzoo=ic+4
      iPzoo=ic+5
      iNO3_=ic+6
      iNH4_=ic+7
      iPON_=ic+8
      iDON_=ic+9
      iSiOH=ic+10
      iopal=ic+11
# ifdef IRON_LIMIT
      iFeSp=ic+12
      iFeLp=ic+13
      iFeD_=ic+14
# endif

      RETURN
      END SUBROUTINE initialize_biology
