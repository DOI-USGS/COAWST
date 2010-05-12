!
!svn $Id: ecosim_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Ecosim tracer parameters:                                           !
!                                                                      !
!  NBands         Number of spectral irradiance bands.                 !
!  Nbac           Number of bacteria constituents.                     !
!  Ndom           Number of dissolved matter constituents.             !
!  Nfec           Number of fecal matter constituents.                 !
!  Nphy           Number of phytoplankton constituents.                !
!  Npig           Number of pigment constituents.                      !
!  PHY            Indices of phytoplankton species considered.         !
!  PIG            Phytoplankton-pigment matrix.                        !
!                                                                      !
!  EcoSim Model Phytoplaknton Parameters:                              !
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
      USE mod_scalars
!
      implicit none
!
!-----------------------------------------------------------------------
!  Bio-optical EcoSim parameters.
!-----------------------------------------------------------------------
!
      integer, parameter :: NBands = 60      ! spectral bands
      integer, parameter :: Nbac = 1         ! bacteria constituents
      integer, parameter :: Ndom = 2         ! DOM constituents
      integer, parameter :: Nfec = 2         ! Fecal constituents
      integer, parameter :: Nphy = 4         ! Phytoplankton groups
      integer, parameter :: Npig = 7         ! Pigments
!
!  Determine number of EcoSim biological tracer. Currently, there is a
!  maximum of seven phytoplankton species and seven different pigments:
!
! [1] small diatom           [1] chlorophyll-a
! [2] large diatom           [2] chlorophyll-b
! [3] small dinoflagellate   [3] chlorophyll-c
! [4] large dinoflagellate   [4] photosythetic carotenoids
! [5] synechococcus          [5] photoprotective carotenoids
! [6] small prochlorococcus  [6] low  urobilin phycoeurythin carotenoids
! [7] large prochlorococcus  [7] high urobilin phycoeurythin carotenoids
!
!  The phytoplankton/pigment matrix is as follows:
!
!               P h y t o p l a n k t o n
!              [1]   [2]   [3]   [4]   [5]   [6]   [7]
!
!       t [7]   0     0     0     0     1     0     0
!       n [6]   0     0     0     0     0     0     0
!       e [5]   1     1     1     1     1     1     1
!       m [4]   1     1     1     1     0     0     0
!       g [3]   1     1     1     1     0     0     0
!       i [2]   0     0     0     0     0     1     1
!       P [1]   1     1     1     1     1     1     1
!
      integer, parameter, dimension(7,7) :: PIG = reshape (             &
     &                                      (/ 1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 1, 1,     &
     &                                         1, 1, 1, 1, 0, 0, 0,     &
     &                                         1, 1, 1, 1, 0, 0, 0,     &
     &                                         1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 1, 0, 0 /),  &
     &                                      (/ 7, 7 /) )
!
!  Set phytoplankton species to consider (see above classification):
!
      integer, parameter, dimension(Nphy) :: PHY = (/ 1, 2, 4, 5 /)
!
!-----------------------------------------------------------------------
!  Set biological tracer identification indices.
!-----------------------------------------------------------------------
!
      integer, allocatable :: idbio(:)  ! Biological tracers

      integer :: iBacC(Nbac)        ! Bacteria, Carbon group
      integer :: iBacN(Nbac)        ! Bacteria, Nitrogen group
      integer :: iBacP(Nbac)        ! Bacteria, Phosphorous group
      integer :: iBacF(Nbac)        ! Bacteria, Iron group
      integer :: iCDMC(Ndom)        ! Color degradational matter
      integer :: iDOMC(Ndom)        ! DOM, Carbon group
      integer :: iDOMN(Ndom)        ! DOM, Nitrogen group
      integer :: iDOMP(Ndom)        ! DOM, Phosphorous group
      integer :: iFecC(Nfec)        ! Fecal matter, Carbon group
      integer :: iFecN(Nfec)        ! Fecal matter, Nitrogen group
      integer :: iFecP(Nfec)        ! Fecal matter, Phosphorous group
      integer :: iFecF(Nfec)        ! Fecal matter, Iron group
      integer :: iFecS(Nfec)        ! Fecal matter, Silica group
      integer :: iPhyC(Nphy)        ! Phytoplankton, Carbon group
      integer :: iPhyN(Nphy)        ! Phytoplankton, Nitrogen group
      integer :: iPhyP(Nphy)        ! Phytoplankton, Phosphorous group
      integer :: iPhyF(Nphy)        ! Phytoplankton, Iron group
      integer :: iPhyS(Nphy)        ! Phytoplankton, Silica group
      integer :: iPigs(Nphy,Npig)   ! Phytoplankton, pigment group
      integer :: iNO3_              ! Nitrate concentration
      integer :: iNH4_              ! Ammonium concentration
      integer :: iPO4_              ! Phosphate concentration
      integer :: iFeO_              ! Iron concentration
      integer :: iSiO_              ! Silica concentration
      integer :: iDIC_              ! Dissolved inorganic Carbon
      integer :: FirstPig           ! Index of first tracer pigment
!
!-----------------------------------------------------------------------
!  EcoSim group names used on standard output.
!-----------------------------------------------------------------------
!
      character (len=16), dimension(Nbac) :: BacName
      character (len=11), dimension(Ndom) :: DomName
      character (len=13), dimension(Nfec) :: FecName
      character (len=21), dimension(Nphy) :: PhyName

      character (len=39), dimension(7) :: PigName =                     &
     &          (/ 'chlorophyll-a                          ',           &
     &             'chlorophyll-b                          ',           &
     &             'chlorophyll-c                          ',           &
     &             'photosynthetic carotenoids             ',           &
     &             'photoprotective carotenoids            ',           &
     &             'low urobilin phycoeurythin carotenoids ',           &
     &             'high urobilin phycoeurythin carotenoids' /)
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
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic, j, ng

      character (len=21), dimension(7) :: PhyGroups =                   &
     &                                 (/ 'small diatom         ',      &
     &                                    'large diatom         ',      &
     &                                    'small dinoflagellate ',      &
     &                                    'large dinoflagellate ',      &
     &                                    'synechococcus        ',      &
     &                                    'small prochlorococcus',      &
     &                                    'large prochlorococcus' /)
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
!       Nutrients: NO3, NO4, PO4, FeO, SiO, DIC  (6)
!        Bacteria: C, Fe, N, P                   (Nbac*4)
!             DOM: CDM, C, N, P                  (Ndom*4)
!           Fecal: C, Fe, N, P, Si               (Nfec*5)
!    Phytoplakton: C, Fe, N, P                   (Nfec*4 + Si)
!        Pigments: look table
!
      NBT=6+(Nbac*4)+(Ndom*4)+(Nfec*5)+(Nphy*4)
!
!  Add phytoplankton silica constituents.
!
      DO i=1,Nphy
        IF (PHY(i).le.2) NBT=NBT+1
      END DO
!
!  Add pigments.  Check phytoplankton-pigment table for values greater
!  than zero.
!
      DO j=1,Npig
        DO i=1,Nphy
          IF (PIG(PHY(i),j).eq.1) NBT=NBT+1
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
!
!  Set identification indices.
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iDIC_=ic+1
      iFeO_=ic+2
      iNH4_=ic+3
      iNO3_=ic+4
      iPO4_=ic+5
      iSiO_=ic+6
      ic=ic+6
      DO i=1,Nbac
        iBacC(i)=ic+1
        iBacF(i)=ic+2
        iBacN(i)=ic+3
        iBacP(i)=ic+4
        ic=ic+4
      END DO
      DO i=1,Ndom
        iCDMC(i)=ic+1
        iDOMC(i)=ic+2
        iDOMN(i)=ic+3
        iDOMP(i)=ic+4
        ic=ic+4
      END DO
      DO i=1,Nfec
        iFecC(i)=ic+1
        iFecF(i)=ic+2
        iFecN(i)=ic+3
        iFecP(i)=ic+4
        iFecS(i)=ic+5
        ic=ic+5
      END DO
      DO i=1,Nphy
        iPhyC(i)=ic+1
        iPhyF(i)=ic+2
        iPhyN(i)=ic+3
        iPhyP(i)=ic+4
        ic=ic+4
      END DO
      DO i=1,Nphy
        IF (PHY(i).le.2) THEN
          ic=ic+1
          iPhyS(i)=ic
        ELSE
          iPhyS(i)=0
        END IF
      END DO
      FirstPig=ic+1
      DO j=1,Npig
        DO i=1,Nphy
          iPigs(i,j)=0
          IF (PIG(PHY(i),j).eq.1) THEN
            ic=ic+1
            iPigs(i,j)=ic
          END IF
        END DO
      END DO
!
!  Set EcoSim group variable names.
!
      DO i=1,Nphy
        PhyName(i)=PhyGroups(PHY(i))
      END DO
      DO i=1,Nbac
        WRITE (BacName(i),'(a,1x,i1)') 'Bacteria Group', i
      END DO
      DO i=1,Ndom
        WRITE (DomName(i),'(a,1x,i1)') 'DOM Group', i
      END DO
      DO i=1,Nfec
        WRITE (FecName(i),'(a,1x,i1)') 'Fecal Group', i
      END DO

      RETURN
      END SUBROUTINE initialize_biology
