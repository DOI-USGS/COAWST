#include "cppdefs.h"
      MODULE mod_biology
#ifdef BIOLOGY
!
!svn $Id: mod_biology.F 1589 2008-07-24 04:02:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================

# ifdef BIO_FASHAM
!                                                                      !
!  Parameters for Fasham et al. (1990) type model:                     !
!                                                                      !
!   AttSW    Light attenuation due to sea water [1/m].                 !
!   AttChl   Light attenuation by Chlorophyll [1/(mg_Chl m2)].         !
!   BioIter  Maximum number of iterations to achieve convergence       !
!              of the nonlinear solution.                              !
!   Chl2C_m  Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].        !
!   ChlMin   Chlorophill minimum threshold value [mg_Chl/m3].          !
!   CoagR    Coagulation rate: agregation rate of SDeN + Phyt ==> LDeN !
!              [1/day].                                                !
!   D_p5NH4  Half-saturation radiation for nitrification inhibition    !
!              [Watts/m2].                                             !
!   I_thNH4  Radiation threshold for nitrification inhibition          !
!              [Watts/m2].                                             !
!   K_NH4    Inverse half-saturation for Phytoplankton NH4 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_NO3    Inverse half-saturation for Phytoplankton NO3 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_Phy    Zooplankton half-saturation, squared constant for         !
!              ingestion [mmol_N/m3]^2.                                !
!   LDeRR    Large Detrital re-mineralization rate [1/day].            !
!   NitriR   Nitrification rate: oxidation of NH4 to NO3 [1/day].      !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   PhyCN    Phytoplankton Carbon:Nitrogen ratio [mol_C/mol_N].        !
!   PhyIP    Phytoplankton NH4 inhibition parameter [1/(mmol_N)].      !
!   PhyIS    Phytoplankton, initial slope of the P-I curve             !
!              [mg_C/(mg_Chl W m-2 day)].                              !
!   ZooMin   Phytoplankton minimum threshold value [mmol_N/m3].        !
!   PhyMR    Phytoplankton mortality rate [1/day] to small detritus.   !
!   SDeAR    Small detritus aggregation rate into Large detritus       !
!              [1/day].                                                !
!   SDeBR    Small Detrital breakdown to NH4 rate [1/day].             !
!   SDeRR    Large Detrital re-mineralization rate [1/day].            !
!   Vp0      Eppley temperature-limited and light-limited growth       !
!              tuning parameter [nondimensional].                      !
!   wLDet    Vertical sinking velocities for Large Detritus            !
!              fraction [m/day].                                       !
!   wPhy     Vertical sinking velocity for Phytoplankton               !
!              fraction [m/day].                                       !
!   wSDet    Vertical sinking velocities for Small Detritus            !
!              fraction [m/day].                                       !
!   ZooAE_N  Zooplankton nitrogen assimilation efficiency fraction     !
!              [nondimensional].                                       !
!   ZooBM    Zooplankton basal metabolism [1/day].                     !
!   ZooCN    Zooplankton Carbon:Nitrogen ratio [mol_C/mol_N].          !
!   ZooER    Zooplankton specific excretion rate [1/day].              !
!   ZooGR    Zooplankton maximum growth rate [1/day].                  !
!   ZooMin   Zooplankton minimum threshold value [mmol_N/m3].          !
!   ZooMR    Zooplankton mortality to Detritus [1/day].                !
!   pCO2air  CO2 partial pressure in the air [ppmv].                   !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: AttChl        ! 1/(mg_Chl m2)
        real(r8), dimension(Ngrids) :: Chl2C_m       ! mg_Chl/mg_C
        real(r8), dimension(Ngrids) :: ChlMin        ! mg_Chl/m3
        real(r8), dimension(Ngrids) :: CoagR         ! 1/day
        real(r8), dimension(Ngrids) :: D_p5NH4       ! Watts/m2
        real(r8), dimension(Ngrids) :: I_thNH4       ! Watts/m2
        real(r8), dimension(Ngrids) :: K_NH4         ! m3/mmol_N
        real(r8), dimension(Ngrids) :: K_NO3         ! m3/mmol_N
        real(r8), dimension(Ngrids) :: K_Phy         ! (mmol_N/m3)^2
        real(r8), dimension(Ngrids) :: LDeRRN        ! 1/day
        real(r8), dimension(Ngrids) :: LDeRRC        ! 1/day
        real(r8), dimension(Ngrids) :: NitriR        ! 1/day
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: PhyCN         ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: PhyIP         ! 1/mmol_N
        real(r8), dimension(Ngrids) :: PhyIS         ! 1/(Watts m-2 day)
        real(r8), dimension(Ngrids) :: PhyMin        ! mmol_N/m3
        real(r8), dimension(Ngrids) :: PhyMR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeAR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeBR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeRRN        ! 1/day
        real(r8), dimension(Ngrids) :: SDeRRC        ! 1/day
        real(r8), dimension(Ngrids) :: Vp0           ! nondimensional
        real(r8), dimension(Ngrids) :: wLDet         ! m/day 
        real(r8), dimension(Ngrids) :: wPhy          ! m/day
        real(r8), dimension(Ngrids) :: wSDet         ! m/day
        real(r8), dimension(Ngrids) :: ZooAE_N       ! nondimensional
        real(r8), dimension(Ngrids) :: ZooBM         ! 1/day
        real(r8), dimension(Ngrids) :: ZooCN         ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: ZooER         ! 1/day
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMin        ! mmol_N/m3
        real(r8), dimension(Ngrids) :: ZooMR         ! 1/day
        real(r8), dimension(Ngrids) :: pCO2air       ! ppmv

# endif
# ifdef BIO_LIMADONEY
!                                                                      !
!  Parameters for Lima and Doney (2004) model:                         !
!                                                                      !
!   AttSW     Light attenuation due to seawater [1/m].                 !
!   AttChl    Light attenuation by chlorophyll [1/(mg_Chl m2)].        !
!   PARfrac   Photosynthetically active fraction of shortwave          !
!               radiation.                                             !
!   alphaPI   Initial slope of P_I curve (in carbon).                  !
!   thetaN0   Maximum Chlorophyll to phytoplankton Nitrogen ratio.     !
!   Tref      Reference temperature.                                   !
!   e_a       Temperature coefficient.                                 !
!   p_cref    Maximum phytoplankton C-spec. growth rate at Tref.       !
!   z_umax1   Maximum grazing rate on small phytoplankton at Tref.     !
!   z_umax2   Maximum grazing rate on large phytoplankton at Tref.     !
!   p1_excr   Small phytoplankton excretion rate.                      !
!   p2_excr   Large phytoplankton excretion rate.                      !
!   p1_mort2  Small phytoplankton quadratic mortality rate.            !
!   p2_mort2  Large phytoplankton quadratic mortality rate.            !
!   z_mort    Zooplankton mortality rate.                              !
!   z_mort2   zooplankton quadratic mortality rate.                    !
!   d1_rem    Small detrital remineralization rate.                    !
!   d2_rem    Large detrital remineralization rate.                    !
!   wDet2     Large detritial sinking rate.                            !
!   p1_kno3   Small phytoplankton NO3 uptake coefficient.              !
!   p2_kno3   Large phytoplankton NO3 uptake coefficient.              !
!   p1_knh4   Small phytoplankton NH4 uptake coefficient.              !
!   p2_knh4   Large phytoplankton NH4 uptake coefficient.              !
!   p2_si     Large phytoplankton Si uptake coefficient.               !
!   p_nitr    Nitrification rate.                                      !
!   falloc0   Minimum recycling allocation factor.                     !
!   si_alloc  Silica recycling allocation factor.                      !
!   z_grz     Grazing coefficient.                                     !
!   zegest0   Zooplankton egestion coefficient.                        !
!   r_ref     Respiration rate at temperature Tref.                    !
!   q_min     Minimum phytoplankton N:C ratio (mol/mol).               !
!   q_max     Maximum phytoplankton N:C ratio (mol/mol).               !
!   q_simin   Minimum large phytoplankton Si:C ratio (mol/mol).        !
!   q_simax   Maximum large phytoplankton Si:C ratio (mol/mol).        !
!   vncref    Value of VNC at Tref.                                    !
!   biolambda cost of biosynthesis (mmolC/mmolN) GD98.                 !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: AttChl        ! 1(mg_Chl m2)
        real(r8), dimension(Ngrids) :: PARfrac
        real(r8), dimension(Ngrids) :: alphaPI
        real(r8), dimension(Ngrids) :: thetaN0
        real(r8), dimension(Ngrids) :: tref
        real(r8), dimension(Ngrids) :: e_a
        real(r8), dimension(Ngrids) :: p_cref
        real(r8), dimension(Ngrids) :: z_umax1
        real(r8), dimension(Ngrids) :: z_umax2
        real(r8), dimension(Ngrids) :: p1_excr
        real(r8), dimension(Ngrids) :: p2_excr
        real(r8), dimension(Ngrids) :: p1_mort2
        real(r8), dimension(Ngrids) :: p2_mort2
        real(r8), dimension(Ngrids) :: z_mort
        real(r8), dimension(Ngrids) :: z_mort2
        real(r8), dimension(Ngrids) :: d1_rem
        real(r8), dimension(Ngrids) :: d2_rem
        real(r8), dimension(Ngrids) :: wDet2         ! m/day    
        real(r8), dimension(Ngrids) :: p1_kno3
        real(r8), dimension(Ngrids) :: p2_kno3
        real(r8), dimension(Ngrids) :: p1_knh4
        real(r8), dimension(Ngrids) :: p2_knh4
        real(r8), dimension(Ngrids) :: p2_si
        real(r8), dimension(Ngrids) :: p_nitr
        real(r8), dimension(Ngrids) :: falloc0
        real(r8), dimension(Ngrids) :: si_alloc
        real(r8), dimension(Ngrids) :: z_grz
        real(r8), dimension(Ngrids) :: zegest0
        real(r8), dimension(Ngrids) :: r_ref
        real(r8), dimension(Ngrids) :: q_min
        real(r8), dimension(Ngrids) :: q_max
        real(r8), dimension(Ngrids) :: q_simin
        real(r8), dimension(Ngrids) :: q_simax
        real(r8), dimension(Ngrids) :: vncref
        real(r8), dimension(Ngrids) :: biolambda

# endif
# ifdef NEMURO
!
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
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AlphaPL       ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: AlphaPS       ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: AlphaZL       ! nondimensional
        real(r8), dimension(Ngrids) :: AlphaZP       ! nondimensional
        real(r8), dimension(Ngrids) :: AlphaZS       ! nondimensional
        real(r8), dimension(Ngrids) :: AttPL         ! m2/mmole_N
        real(r8), dimension(Ngrids) :: AttPS         ! m2/mmole_N
        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: BetaPL        ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: BetaPS        ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: BetaZS        ! nondimensional
        real(r8), dimension(Ngrids) :: BetaZL        ! nondimensional
        real(r8), dimension(Ngrids) :: BetaZP        ! nondimensional
        real(r8), dimension(Ngrids) :: GammaL        ! nondimensional
        real(r8), dimension(Ngrids) :: GammaS        ! nondimensional
        real(r8), dimension(Ngrids) :: GRmaxLpl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxLps      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxLzs      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxPpl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxPzl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxPzs      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxSps      ! 1/day
        real(r8), dimension(Ngrids) :: KD2N          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGppL         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGppS         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGraL         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGraP         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGraS         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorPL        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorPS        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorZL        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorZP        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorZS        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KNH4L         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KNH4S         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KNit          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KNO3L         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KNO3S         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KO2S          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KP2D          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KP2N          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KPL2ZL        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPS2ZL        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPS2ZS        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPL2ZP        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KResPL        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KResPS        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KSiL          ! mmole_Si/m3
        real(r8), dimension(Ngrids) :: KZL2ZP        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KZS2ZL        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KZS2ZP        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: LamL          ! m3/mmole_N
        real(r8), dimension(Ngrids) :: LamP          ! m3/mmole_N
        real(r8), dimension(Ngrids) :: LamS          ! m3/mmole_N
        real(r8), dimension(Ngrids) :: MorPL0        ! m3/mmole_N/day
        real(r8), dimension(Ngrids) :: MorPS0        ! m3/mmole_N/day
        real(r8), dimension(Ngrids) :: MorZL0        ! m3/mmole_N 1/day
        real(r8), dimension(Ngrids) :: MorZP0        ! m3/mmole_N 1/day
        real(r8), dimension(Ngrids) :: MorZS0        ! m3/mmole_N 1/day
        real(r8), dimension(Ngrids) :: Nit0          ! 1/day
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: PusaiL        ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PusaiPL       ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PusaiS        ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PusaiZS       ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PL2ZLstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PL2ZPstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PS2ZLstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PS2ZSstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: ResPL0        ! 1/day
        real(r8), dimension(Ngrids) :: ResPS0        ! 1/day
        real(r8), dimension(Ngrids) :: RSiN          ! mmole_Si/mmole_N
        real(r8), dimension(Ngrids) :: setVOpal      ! m/day
        real(r8), dimension(Ngrids) :: setVPON       ! m/day
        real(r8), dimension(Ngrids) :: VD2N0         ! 1/day
        real(r8), dimension(Ngrids) :: VmaxL         ! 1/day
        real(r8), dimension(Ngrids) :: VmaxS         ! 1/day
        real(r8), dimension(Ngrids) :: VO2S0         ! 1/day
        real(r8), dimension(Ngrids) :: VP2D0         ! 1/day
        real(r8), dimension(Ngrids) :: VP2N0         ! 1/day
        real(r8), dimension(Ngrids) :: ZL2ZPstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: ZS2ZLstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: ZS2ZPstar     ! mmole_N/m3
# endif

# ifdef NPZD_FRANKS
!
!  Parameters for Franks et al. (1986) type model:                     !
!                                                                      !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  DetRR     Detritus remineraliztion rate, [1/day].                   !
!  K_ext     Light extinction coefficient, [1/m].                      !
!  K_NO3     Inverse half-saturation for phytoplankton nitrate uptake  !
!              [1/(millimole_N m-3)].                                  !
!  K_Phy     Phytoplankton saturation coefficient, [millimole_N m-3].  !
!  PhyMR     Phytoplankton senescence/mortality rate, [1/day].         !
!  Vm_NO3    Nitrate uptake rate, [1/day].                             !
!  wDet      Detrital sinking rate, [m/day].                           !
!  ZooGR     Zooplankton maximum growth rate, [1/day].                 !
!  ZooMR     Zooplankton mortality rate, [1/day].                      !
!  ZooMD     Zooplankton death bits rate, [1/day].                     !
!  ZooGA     Zooplankton grazing inefficiency, [nondimensional].       !
!  ZooEC     Zooplankton excreted fraction, [nondimensional].          !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

#  ifdef ANA_BIOLOGY
        real(r8), allocatable :: BioIni(:,:)
#  endif
        real(r8), dimension(Ngrids) :: DetRR         ! 1/day
        real(r8), dimension(Ngrids) :: K_ext         ! 1/m
        real(r8), dimension(Ngrids) :: K_NO3         ! 1/(mmol/m3)
        real(r8), dimension(Ngrids) :: K_Phy         ! mmol/m3
        real(r8), dimension(Ngrids) :: PhyMR         ! 1/day
        real(r8), dimension(Ngrids) :: Vm_NO3        ! 1/day
        real(r8), dimension(Ngrids) :: wDet          ! m/day
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wDet
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wDet
#  endif
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMD         ! 1/day
        real(r8), dimension(Ngrids) :: ZooGA         ! nondimensional
        real(r8), dimension(Ngrids) :: ZooEC         ! nondimensional
# endif
# ifdef NPZD_POWELL
!
!  Parameters for Franks et al. (1986) type model:                     !
!                                                                      !
!  AttPhy    Light attenuation due to phytoplankton (self-shading      !
!              coefficient), [m2/millimole_N].                         !
!  AttSW     Light attenuation due to sea water, [1/m].                !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  DetRR     Detritus remineraliztion rate, [1/day].                   !
!  K_NO3     Half-saturation for phytoplankton nitrate uptake          !
!              [millimole_N m-3].                                      !
!  Ivlev     Ivlev constant for zooplankton grazing parameterization,  !
!              [nondimensional].                                       !
!  PARfrac   Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!  PhyIS     Phytoplankton, initial slope of the P-I curve [m2/W].     !
!  PhyMRD    Phytoplankton mortality rate to the Detritus pool,        !
!              [1/day].                                                !
!  PhyMRN    Phytoplankton mortality rate to the Nitrogen pool,        !
!              [1/day].                                                !
!  Vm_NO3    Nitrate uptake rate, [1/day].                             !
!  wDet      Detrital sinking rate, [m/day].                           !
!  wPhy      Phytoplankton sinking rate, [m/day].                      !
!  ZooEED    Zooplankton excretion efficiency to Detritus pool,        !
!              {nondimensional].                                       !
!  ZooEEN    Zooplankton excretion efficiency to Nitrogen pool,        !
!              {nondimensional].                                       !
!  ZooGR     Zooplankton grazing rate, [1/day].                        !
!  ZooMRD    Zooplankton mortality rate to Detritus pool, [1/day].     !
!  ZooMRN    Zooplankton mortality rate to Nitrogen pool, [1/day].     !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

#  ifdef ANA_BIOLOGY
        real(r8), allocatable :: BioIni(:,:)
#  endif
        real(r8), dimension(Ngrids) :: AttPhy        ! m2/mmole
        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: DetRR         ! 1/day
        real(r8), dimension(Ngrids) :: K_NO3         ! mmol/m3
        real(r8), dimension(Ngrids) :: Ivlev         ! nondimensional
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_PARfrac
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_PARfrac
#  endif
        real(r8), dimension(Ngrids) :: PhyIS         ! m2/W
        real(r8), dimension(Ngrids) :: PhyMRD        ! 1/day
        real(r8), dimension(Ngrids) :: PhyMRN        ! 1/day
        real(r8), dimension(Ngrids) :: Vm_NO3        ! 1/day
        real(r8), dimension(Ngrids) :: wDet          ! m/day
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wDet
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wDet
#  endif
        real(r8), dimension(Ngrids) :: wPhy          ! m/day
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wPhy
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wPhy
#  endif
        real(r8), dimension(Ngrids) :: ZooEED        ! nondimensional
        real(r8), dimension(Ngrids) :: ZooEEN        ! nondimensional
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMRD        ! 1/day
        real(r8), dimension(Ngrids) :: ZooMRN        ! 1/day
# endif
# ifdef ECOSIM
!
!  EcoSim Model Phytoplaknton Parameters:                              !
!                                                                      !
!                                                                      !
!  HsNO3          Half-saturation for phytoplankton NO3 uptake         !
!                   (micromole_NO3/liter).                             !
!  HsNH4          Half-saturation for phytoplankton NH4 uptake         !
!                   (micromole_NH4/liter).                             !
!  HsSiO          Half-saturation for phytoplankton SiO uptake         !
!                   (micromole_SiO/liter).                             !
!  HsPO4          Half-saturation for phytoplankton PO4 uptake         !
!                   (micromole_PO4/liter).                             !
!  HsFe           Half-saturation for phytoplankton Fe uptake          !
!                  (micromole_Fe/liter).                               !
!  GtALG_max      Maximum 24 hour growth rate (1/day).                 !
!  PhyTbase       Phytoplankton temperature base for exponential       !
!                   response to temperature (Celsius).                 !
!  PhyTfac        Phytoplankton exponential temperature factor         !
!                   (1/Celsius).                                       !
!  BET_           Nitrate uptake inhibition for NH4 (l/micromole).     !
!  maxC2nALG      Maximum phytoplankton C:N ratio                      !
!                   (micromole_C/micromole_N).                         !
!  minC2nALG      Balanced phytoplankton C:N ratio                     !
!                   (micromole_C/micromole_N).                         !
!  C2nALGminABS   Absolute minimum phytoplankton C:N ratio             !
!                   (micromole_C/micromole_N).                         !
!  maxC2SiALG     Maximum phytoplankton C:Si ratio                     !
!                   (micromole_C/micromole_Si).                        !
!  minC2SiALG     Balanced phytoplankton C:Si ratio                    !
!                   (micromole_C/micromole_Si).                        !
!  C2SiALGminABS  Absolute minimum phytoplankton C:Si ratio            !
!                  (micromole_C/micromole_Si).                         !
!  maxC2pALG      Maximum phytoplankton C:P ratio                      !
!                   (micromole_C/micromole_P).                         !
!  minC2pALG      Balanced phytoplankton C:P ratio                     !
!                   (micromole_C/micromole_P).                         !
!  C2pALGminABS   Absolute minimum phytoplankton C:P ratio             !
!                   (micromole_C/micromole_P).                         !
!  maxC2FeALG     Maximum phytoplankton C:Fe ratio                     !
!                   (micromole_C/micromole_Fe).                        !
!  minC2FeALG     Balanced phytoplankton C:Fe ratio                    !
!                   (micromole_C/micromole_Fe).                        !
!  C2FeALGminABS  Absolute minimum phytoplankton C:Fe ratio            !
!                   (micromole_C/micromole_Fe).                        !
!  qu_yld         Maximum quantum yield                                !
!                   (micromole_C/micromole_quanta).                    !
!  E0_comp        Compensation light level (micromole_quanta).         !
!  E0_inhib       Light level for onset of photoinhibition             !
!                   (micromole_quanta).                                !
!  inhib_fac      Exponential decay factor for light limited growth    !
!                   (1/micromole_quanta).                              !
!  C2Chl_max      Maximum lighted limited (nutrient replete) C:Chl     !
!                   ratio (microgram_C/microgram_Chl).                 !
!  mxC2Cl         Rate of change in the lighted limited C:Chl ratio    !
!                   (microgram_C/microgram_Chl/micromole_quanta).      !
!  b_C2Cl         Minimum lighted limited (nutrient replete) C:Chl     !
!                   ratio (microgram_C/microgram_Chl).                 !
!  mxC2Cn         Rate of change in the nutrient limited C:Chl ratio   !
!                   [(ug_C/ug_Chl)/(umole_C/umole_N)].                 !
!  b_C2Cn         Minimum nutrient limited C:Chl ratio                 !
!                   (microgram_C/microgram_Chl).                       !
!  mxPacEff       Rate of change in package effect                     !
!                   [1/(microgram_C/microgram_Chl)].                   !
!  b_PacEff       Maximum package effect                               !
!                   [1/(microgram_C/microgram_Chl)].                   !
!  mxChlB         Rate of change in the Chl_b:Chl_a ratio              !
!                   [(ug_Chl_b/ug_Chl_a)/(ug_C/ug_Chl_a)].             !
!  b_ChlB         Maximum Chl_b:Chl_a ratio                            !
!                   (microgram_Chl_b/microgram_Chl_a).                 !
!  mxChlC         Rate of change in the Chl_c:Chl_a ratio              !             
!                   [(ug_Chl_c/ug_Chl_a)/(ug_C/ug_Chl_a)].             !
!  b_ChlC         Maximum Chl_c:Chl_a ratio                            !
!                   (microgram_Chl_c/microgram_Chl_a).                 !
!  mxPSC          Rate of change in the PSC:Chl_a ratio                !
!                   [(ug_PSC/ug_Chl_a)/(ug_C/ug_Chl_a)].               !
!  b_PSC          Maximum PSC:Chl_a ratio                              !
!                  (microgram_Chl_c/microgram_Chl_a).                  !
!  mxPPC          Rate of change in the PPC:Chl_a ratio                !
!                   [(ug_PPC/ug_Chl_a)/(ug_C/ug_Chl_a)].               !
!  b_PPC          Maximum PPC:Chl_a ratio                              !
!                  (microgram_Chl_c/microgram_Chl_a).                  !
!  mxLPUb         Rate of change in the LPUb:Chl_a ratio               !
!                   [(ug_LPUb/ug_Chl_a)/(ug_C/ug_Chl_a)].              !
!  b_LPUb         Maximum LPUb:Chl_a ratio                             !
!                   (microgram_HPUb/microgram_Chl_a).                  !
!  mxHPUb         Rate of change in the HPUb:Chl_a ratio               !
!                   [(ug_HPUb/ug_Chl_a)/(ug_C/ug_Chl_a)].              !
!  b_HPUb         Maximum HPUb:Chl_a ratio                             !
!                   (microgram_HPUb/microgram_Chl_a).                  !
!  FecDOC         Proportion of grazing stress which is apportioned    !
!                   to DOM (nondimensional).                           !
!  FecPEL         Proportion of grazing stress which is apportioned    !
!                   to fecal pellets (nondimesional).                  !
!  FecCYC         Proportion of grazing stress which is apportioned    !
!                   to direct remineralization (nondimensional).       !
!  ExALG          Proportion of daily production that is lost to       !
!                   excretion (nondimensional).                        !
!  WS             Phytoplankton sinking speed (meters/day).            !
!  HsGRZ          Phytoplankton grazing parameter (nondimensional).    !
!  MinRefuge      Refuge Phytoplankton population (micromole_C/liter). !
!  RefugeDep      Maximum Refuge Phytoplankton depth (meters).         !
!  Norm_Vol       Normalized Volume factor (nondimensional).           !
!  Norm_Surf      Normalized surface area factor (nondimensional).     !
!  HsDOP          Half Saturation Constant for DOP uptake              !
!                   (micromole_DOP/liter).                             !
!  C2pALKPHOS     C:P ratio where DOP uptake begins                    !
!                   (micromole_C/micromole_P).                         !
!  HsDON          Half Saturation Constant for DON uptake              !
!                   (micromole_DON/liter).                             !
!  C2nNupDON      C:N ratio where DON uptake begins                    !
!                   (micromole_C/micromole_N).                         !
!                                                                      !
! Bacteria Parameters:                                                 !
!                                                                      !
!  HsDOC_ba       Half saturation constant for bacteria DOC uptake     !
!                   (micromole_DOC/liter).                             !
!  GtBAC_max      Maximum 24 hour bacterial growth rate (1/day).       !
!  BacTbase       Phytoplankton temperature base for exponential       !
!                   response to temperature, (Celsius).                !
!  BacTfac        Phytoplankton exponential temperature factor         !
!                   (1/Celsius).                                       !
!  C2nBAC         Carbon to Nitrogen ratio of Bacteria                 !
!                   (micromole_C/micromole_N).                         !
!  C2pBAC         Carbon to Phosphorus ratio of Bacteria               !
!                   (micromole_C/micromole_P).                         !
!  C2FeBAC        Carbon to Iron ratio of Bacteria                     !
!                   (micromole_C/micromole_Fe)                         !
!  BacDOC         Proportion of bacterial grazing stress which is      !
!                   apportioned to DOM (nondimensional).               !
!  BacPEL         Proportion of bacterial grazing stress which is      !
!                   apportioned to fecal pellets (nondimensional).     !
!  BacCYC         Proportion of bacterial grazing stress which is      !
!                   apportioned to direct remineralization             !
!                   (nondimensional).                                  !
!  ExBAC_c        Bacterial recalcitrant carbon excretion as a         !
!                   proportion of uptake (nondimensional)              !
!  ExBacC2N       Bacterial recalcitrant excretion carbon to nitrogen  !
!                   ratio (micromole_C/micromole_N).                   !
!  Bac_Ceff       Bacterial gross growth carbon efficiency             !
!                   (nondimensional).                                  !
!  RtNIT          Maximum bacterial nitrification rate (1/day).        !
!  HsNIT          Half saturation constant for bacterial nitrification !
!                   (micromole NH4/liter)                              !
!                                                                      !
! Dissolved Organic Matter Parameters:                                 !
!                                                                      !
!  cDOCfrac_c     Colored fraction of DOC production from              !
!                   phytoplankton and bacterial losses                 !
!                   (nondimensional).                                  !
!  RtUVR_DIC      UV degradation of DOC into DIC at 410 nanometers     !
!                   (micromole/meter/liter/hour).                      !
!  RtUVR_DIC      UV degradation of DOC into colorless labile DOC at   !
!                   410 nanometers (micromole/meter/liter/hour).       !
!                                                                      !
! Fecal and detritus Parameters:                                       !
!                                                                      !
!  WF             Fecal sinking flux (meters/day).                     !
!  RegTbase       Fecal regeneration temperature base for exponential  !
!                   response to temperature (Celsius).                 !
!  RegTfac        Fecal regeneration exponential temperature factor    !
!                   (1/Celsius).                                       !
!  RegCR          Fecal carbon regeneration rate (1/day).              !
!  RegNR          Fecal nitrogen regeneration rate (1/day).            !
!  RegSR          Fecal silica regeneration rate (1/day).              !
!  RegPR          Fecal phosphorus regeneration rate (1/day).          !
!  RegFR          Fecal iron regeneration rate (1/day).                !
!                                                                      !
!======================================================================!
!
      USE mod_param
      USE mod_eclight
      USE mod_scalars

      implicit none
!
!-----------------------------------------------------------------------
!  Standard input parameters.
!-----------------------------------------------------------------------
!
!  Number of biological iterations.
!
      integer, dimension(Ngrids) :: BioIter
!
!  Control flags.
!
      logical, dimension(Ngrids) :: RtUVR_flag
      logical, dimension(Ngrids) :: NFIX_flag
      logical, dimension(Ngrids) :: Regen_flag
!
!  Phytoplankton parameters.
!
      real(r8), dimension(Nphy,Ngrids) :: HsNO3
      real(r8), dimension(Nphy,Ngrids) :: HsNH4
      real(r8), dimension(Nphy,Ngrids) :: HsSiO
      real(r8), dimension(Nphy,Ngrids) :: HsPO4
      real(r8), dimension(Nphy,Ngrids) :: HsFe
      real(r8), dimension(Nphy,Ngrids) :: GtALG_max
      real(r8), dimension(Nphy,Ngrids) :: PhyTbase
      real(r8), dimension(Nphy,Ngrids) :: PhyTfac
      real(r8), dimension(Nphy,Ngrids) :: BET_
      real(r8), dimension(Nphy,Ngrids) :: maxC2nALG
      real(r8), dimension(Nphy,Ngrids) :: minC2nALG
      real(r8), dimension(Nphy,Ngrids) :: C2nALGminABS
      real(r8), dimension(Nphy,Ngrids) :: maxC2SiALG
      real(r8), dimension(Nphy,Ngrids) :: minC2SiALG
      real(r8), dimension(Nphy,Ngrids) :: C2SiALGminABS
      real(r8), dimension(Nphy,Ngrids) :: maxC2pALG
      real(r8), dimension(Nphy,Ngrids) :: minC2pALG
      real(r8), dimension(Nphy,Ngrids) :: C2pALGminABS
      real(r8), dimension(Nphy,Ngrids) :: maxC2FeALG
      real(r8), dimension(Nphy,Ngrids) :: minC2FeALG
      real(r8), dimension(Nphy,Ngrids) :: C2FeALGminABS
      real(r8), dimension(Nphy,Ngrids) :: qu_yld
      real(r8), dimension(Nphy,Ngrids) :: E0_comp
      real(r8), dimension(Nphy,Ngrids) :: E0_inhib
      real(r8), dimension(Nphy,Ngrids) :: inhib_fac
      real(r8), dimension(Nphy,Ngrids) :: C2CHL_max
      real(r8), dimension(Nphy,Ngrids) :: mxC2Cl
      real(r8), dimension(Nphy,Ngrids) :: b_C2Cl
      real(r8), dimension(Nphy,Ngrids) :: mxC2Cn
      real(r8), dimension(Nphy,Ngrids) :: b_C2Cn
      real(r8), dimension(Nphy,Ngrids) :: mxPacEff
      real(r8), dimension(Nphy,Ngrids) :: b_PacEff
      real(r8), dimension(Nphy,Ngrids) :: mxChlB
      real(r8), dimension(Nphy,Ngrids) :: b_ChlB
      real(r8), dimension(Nphy,Ngrids) :: mxChlC
      real(r8), dimension(Nphy,Ngrids) :: b_ChlC
      real(r8), dimension(Nphy,Ngrids) :: mxPSC
      real(r8), dimension(Nphy,Ngrids) :: b_PSC
      real(r8), dimension(Nphy,Ngrids) :: mxPPC
      real(r8), dimension(Nphy,Ngrids) :: b_PPC
      real(r8), dimension(Nphy,Ngrids) :: mxLPUb
      real(r8), dimension(Nphy,Ngrids) :: b_LPUb
      real(r8), dimension(Nphy,Ngrids) :: mxHPUb
      real(r8), dimension(Nphy,Ngrids) :: b_HPUb
      real(r8), dimension(Nphy,Ngrids) :: FecDOC
      real(r8), dimension(Nphy,Nfec,Ngrids) :: FecPEL
      real(r8), dimension(Nphy,Ngrids) :: FecCYC
      real(r8), dimension(Nphy,Ngrids) :: ExALG
      real(r8), dimension(Nphy,Ngrids) :: WS
      real(r8), dimension(Nphy,Ngrids) :: HsGRZ
      real(r8), dimension(Nphy,Ngrids) :: MinRefuge
      real(r8), dimension(Nphy,Ngrids) :: RefugeDep
      real(r8), dimension(Nphy,Ngrids) :: Norm_Vol
      real(r8), dimension(Nphy,Ngrids) :: Norm_Surf
      real(r8), dimension(Nphy,Ngrids) :: HsDOP
      real(r8), dimension(Nphy,Ngrids) :: C2pALKPHOS
      real(r8), dimension(Nphy,Ngrids) :: HsDON
      real(r8), dimension(Nphy,Ngrids) :: C2nNupDON
!
!  Bacteria parameters.
!
      real(r8), dimension(Nbac,Ngrids) :: HsDOC_ba
      real(r8), dimension(Nbac,Ngrids) :: GtBAC_max
      real(r8), dimension(Nbac,Ngrids) :: BacTbase
      real(r8), dimension(Nbac,Ngrids) :: BacTfac
      real(r8), dimension(Ngrids) :: C2nBAC
      real(r8), dimension(Ngrids) :: C2pBAC
      real(r8), dimension(Ngrids) :: C2FeBAC
      real(r8), dimension(Ngrids) :: BacDOC
      real(r8), dimension(Ngrids) :: BacPEL
      real(r8), dimension(Ngrids) :: BacCYC
      real(r8), dimension(Ngrids) :: ExBAC_c
      real(r8), dimension(Ngrids) :: ExBacC2N
      real(r8), dimension(Ngrids) :: Bac_Ceff
      real(r8), dimension(Ngrids) :: RtNIT
      real(r8), dimension(Ngrids) :: HsNIT  
!
!  DOM parameters.
!
      real(r8), dimension(Ndom,Ngrids) :: cDOCfrac_c
      real(r8), dimension(Ngrids) :: RtUVR_DIC
      real(r8), dimension(Ngrids) :: RtUVR_DOC
!
!  Fecal parameters.
!
      real(r8), dimension(Nfec,Ngrids) :: WF
      real(r8), dimension(Nfec,Ngrids) :: RegTbase
      real(r8), dimension(Nfec,Ngrids) :: RegTfac
      real(r8), dimension(Nfec,Ngrids) :: RegCR
      real(r8), dimension(Nfec,Ngrids) :: RegNR
      real(r8), dimension(Nfec,Ngrids) :: RegSR
      real(r8), dimension(Nfec,Ngrids) :: RegPR
      real(r8), dimension(Nfec,Ngrids) :: RegFR
!
!-----------------------------------------------------------------------
!  Internal parameters.
!-----------------------------------------------------------------------
!    
!  Spectral band width used in light calculations.

      real(r8), parameter :: DLAM  = 5.0_r8
!
!  Flags used for testing purposes.
!
      real(r8), parameter :: SMALL  = 1.0e-6_r8
      real(r8), parameter :: VSMALL = 1.0e-14_r8
      real(r8), parameter :: LARGE  = 1.0e+10_r8
      real(r8), parameter :: VLARGE = 1.0e+50_r8
!
!  Array indexes for frequently used constituents.
!
      integer, parameter :: ilab=1    ! labile index for DOC.
      integer, parameter :: irct=2    ! relict index for DOC.
      integer, parameter :: ichl=1    ! pigment index for chlorophyll-a
      integer, parameter :: isfc=1    ! slow fecal pellet index
      integer, parameter :: iffc=2    ! fast fecal pellet index
!
!  Phytoplankton calculated paramters.
!
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2nALG   ! inverse C2nALG
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2SiALG  ! inverse C2SiALG
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2pALG   ! inverse C2pALG
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2FeALG  ! inverse C2FeALG
!
!  Bacteria calculated parameters.
!
      real(r8), dimension(Ngrids) :: N2cBAC
      real(r8), dimension(Ngrids) :: P2cBAC
      real(r8), dimension(Ngrids) :: Fe2cBAC
      real(r8), dimension(Nbac,Ngrids) :: HsNH4_ba
      real(r8), dimension(Nbac,Ngrids) :: HsPO4_ba
      real(r8), dimension(Nbac,Ngrids) :: HsFe_ba
      real(r8), dimension(Ngrids) :: R_ExBAC_c
      real(r8), dimension(Ngrids) :: ExBAC_n
      real(r8), dimension(Ngrids) :: Frac_ExBAC_n
      real(r8), dimension(Ngrids) :: I_Bac_Ceff
!
!  Absorption parameters.
!
      real(r8), dimension(NBands) :: wavedp   ! a and b factor
      real(r8), dimension(Ndom) :: aDOC410    ! CDM absorption at 410
      real(r8), dimension(Ndom) :: aDOC300    ! CDM absorption at 300

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine initializes several parameters in module "mod_biology" !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: ibac, iband, ifec, iphy, ng
!
!-----------------------------------------------------------------------
!  Derived parameters.
!-----------------------------------------------------------------------
!
!  Convert rates from day-1 to second-1.
!
      DO ng=1,Ngrids
        DO iphy=1,Nphy
          GtALG_max(iphy,ng)=GtALG_max(iphy,ng)*sec2day
          ExALG(iphy,ng)=ExALG(iphy,ng)*sec2day
          HsGRZ(iphy,ng)=HsGRZ(iphy,ng)*sec2day
          WS(iphy,ng)=WS(iphy,ng)*sec2day
        END DO
        DO ibac=1,Nbac
          GtBAC_max(ibac,ng)=GtBAC_max(ibac,ng)*sec2day
        END DO
        DO ifec=1,Nfec
          WF(ifec,ng)=WF(ifec,ng)*sec2day
        END DO
        RtNIT(ng)=RtNIT(ng)*sec2day
      END DO
!
!  Calculated reciprocal phytoplankton parameters.
!
      DO ng=1,Ngrids
        DO iphy=1,Nphy
          IF (maxC2nALG(iphy,ng).gt.SMALL) THEN
            ImaxC2nALG(iphy,ng)=1.0_r8/maxC2nALG(iphy,ng)
          ELSE
            ImaxC2nALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2SiALG(iphy,ng).gt.SMALL) THEN
            ImaxC2SiALG(iphy,ng)=1.0_r8/maxC2SiALG(iphy,ng)
          ELSE
            ImaxC2SiALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2pALG(iphy,ng).gt.SMALL) THEN
            ImaxC2pALG(iphy,ng)=1.0_r8/maxC2pALG(iphy,ng)
          ELSE
            ImaxC2pALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2FeALG(iphy,ng).gt.SMALL) THEN
            ImaxC2FeALG(iphy,ng)=1.0_r8/maxC2FeALG(iphy,ng)
          ELSE
            ImaxC2FeALG(iphy,ng)=0.0_r8
          END IF
        END DO
      END DO
!
!  Calculated bacterial parameters.
!
      DO ng=1,Ngrids
        DO ibac=1,Nbac
          HsNH4_ba(ibac,ng)=HsDOC_ba(ibac,ng)/C2nBAC(ng)
          HsPO4_ba(ibac,ng)=HsDOC_ba(ibac,ng)/C2pBAC(ng)
          HsFe_ba (ibac,ng)=HsDOC_ba(ibac,ng)/C2FeBAC(ng)
        END DO
      END DO
!
!  Inverse parameters for computational efficiency.
!
      DO ng=1,Ngrids
        N2cBAC(ng)=1.0_r8/C2nBAC(ng)
        P2cBAC(ng)=1.0_r8/C2pBAC(ng)
        Fe2cBAC(ng)=1.0_r8/C2FeBAC(ng)
        I_Bac_Ceff(ng)=1.0_r8/Bac_Ceff(ng) 
      END DO
!
!  Reciprocal of non baterial recalcitran carbon excretion.
!
      DO ng=1,Ngrids
        R_ExBAC_c(ng)=1.0_r8/(1.0_r8-ExBAC_c(ng))
      END DO
!
!  Bacterial recalcitrant nitrogen excretion as a function of uptake.
!
      DO ng=1,Ngrids
        ExBAC_n(ng)=ExBAC_c(ng)*C2nBAC(ng)/ExBacC2N(ng)
        Frac_ExBAC_n(ng)=1.0_r8-ExBAC_n(ng)
      END DO
!
!  Scale UV degradation parameters.
!
      DO ng=1,Ngrids
        RtUVR_DIC(ng)=RtUVR_DIC(ng)/3600.0_r8
        RtUVR_DOC(ng)=RtUVR_DOC(ng)/3600.0_r8 
      END DO
!
!  If applicable, zero-out fecal regeneration parameters.
!
      DO ng=1,Ngrids
        IF (Regen_flag(ng)) THEN
          DO ifec=1,Nfec
            RegCR(ifec,ng)=RegCR(ifec,ng)*sec2day
            RegNR(ifec,ng)=RegNR(ifec,ng)*sec2day
            RegPR(ifec,ng)=RegPR(ifec,ng)*sec2day
            RegFR(ifec,ng)=RegFR(ifec,ng)*sec2day
            RegSR(ifec,ng)=RegSR(ifec,ng)*sec2day
          END DO
        ELSE
          DO ifec=1,Nfec
            RegCR(ifec,ng)=0.0_r8
            RegNR(ifec,ng)=0.0_r8
            RegPR(ifec,ng)=0.0_r8
            RegFR(ifec,ng)=0.0_r8
            RegSR(ifec,ng)=0.0_r8
          END DO
        END IF
      END DO
!
!  Spectral dependency for scattering and backscattering.
!
      DO iband=1,NBands
        wavedp(iband)=(550.0_r8/(397.0_r8+REAL(iband,r8)*DLAM))
      END DO
!
!  Calculated IOP parameter values.
!
      aDOC410(ilab)=aDOC(ilab,1)*EXP(0.014_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC410(irct)=aDOC(irct,1)*EXP(0.025_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC300(ilab)=EXP(0.0145_r8*(410.0_r8-300.0_r8))
      aDOC300(irct)=EXP(0.0145_r8*(410.0_r8-300.0_r8))

      RETURN
      END SUBROUTINE initialize_biology

# endif
#endif

      END MODULE mod_biology
