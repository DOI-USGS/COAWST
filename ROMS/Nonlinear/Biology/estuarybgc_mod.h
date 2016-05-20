!
!svn $Id: estuarybgc_mod.h 2232 2012-01-03 18:55:20Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Fennel et al. (2006) model:                          !
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
#ifdef SPECTRAL_LIGHT
!   SIGATRB  NAP absorb. cross section at 440 nm [meter2 gram-1].      !
!   STRB     Spectral slope of NAP absorption                          !
!   BLTRB    Baseline NAP absorption [meter2 gram-1].                  !
!   SIGBTRB  Scattering cross section of turbidity [meter2 gram-1].    !
!   ETASPEC  Scattering spectral exponent.                             !
!   BB2B     Particulate backscattering ratio.                         !
!   Ndom     Number of dissolved matter constituents.                  !
#endif
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iNO3_                  ! Nitrate concentration
      integer :: iNH4_                  ! Ammonium concentration
      integer :: iChlo                  ! Chlorophyll concentration
      integer :: iPhyt                  ! Phytoplankton concentration
      integer :: iZoop                  ! Zooplankton concentration
      integer :: iLDeN                  ! Large detritus N-concentration
      integer :: iSDeN                  ! Small detritus N-concentration
#ifdef CARBON
      integer :: iLDeC                  ! Large detritus C-concentration
      integer :: iSDeC                  ! Small detritus C-concentration
      integer :: iTIC_                  ! Total inorganic carbon
      integer :: iTAlk                  ! Total alkalinity
#endif
#ifdef OXYGEN
      integer :: iOxyg                  ! Dissolved oxygen concentration
#endif
#ifdef SPECTRAL_LIGHT
      integer :: idPARo                 ! PAR    Photosynthetically Available Radiation (PAR)
      integer :: idPARs                 ! PARs            Spectral PAR
      integer :: idSpKd                 ! SpKd            Spectral Attenuation
#endif
#ifdef SAV_BIOMASS
      integer :: iddinw                 ! DINwcr  Dissolved Inorganic Nitrogen (water column)
      integer :: iddins                 ! DINsed  Dissolved Inorganic Nitrogen (sediment column) 
      integer :: iddowc                 ! DOwcr   Dissolved Oxygen (water column) 
      integer :: idwsvl                 ! DINwcr_sav Dissolved Inorganic N in water column due to SAV
      integer :: idsagb                 ! AGB     Above ground biomass 
      integer :: idsbgb                 ! BGB     Below ground biomass 
      integer :: idsvpp                 ! PP
      integer :: idsvam                 ! AGM
      integer :: idsgar                 ! AGAR
      integer :: idsvbr                 ! AGBR
      integer :: idsvrs                 ! SEARS
      integer :: idsvbg                 ! AGBG
      integer :: idsvag                 ! BGAG
      integer :: idsbgr                 ! BGR
      integer :: idsbgm                 ! BGM

#endif
#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Biological 2D diagnostic variable IDs.
!
      integer, allocatable :: iDbio2(:)       ! 2D biological terms

      integer  :: iCOfx                       ! air-sea CO2 flux
      integer  :: iDNIT                       ! denitrification flux
      integer  :: ipCO2                       ! partial pressure of CO2
      integer  :: iO2fx                       ! air-sea O2 flux
!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms

      integer  :: iPPro = 1                   ! primary productivity
      integer  :: iNO3u = 2                   ! NO3 uptake
#endif
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: AttSW(:)              ! 1/m
      real(r8), allocatable :: AttChl(:)             ! 1/(mg_Chl m2)
      real(r8), allocatable :: Chl2C_m(:)            ! mg_Chl/mg_C
      real(r8), allocatable :: ChlMin(:)             ! mg_Chl/m3
      real(r8), allocatable :: CoagR(:)              ! 1/day
      real(r8), allocatable :: D_p5NH4(:)            ! Watts/m2
      real(r8), allocatable :: I_thNH4(:)            ! Watts/m2
      real(r8), allocatable :: K_NH4(:)              ! m3/mmol_N
      real(r8), allocatable :: K_NO3(:)              ! m3/mmol_N
      real(r8), allocatable :: K_Phy(:)              ! (mmol_N/m3)^2
      real(r8), allocatable :: LDeRRN(:)             ! 1/day
      real(r8), allocatable :: LDeRRC(:)             ! 1/day
      real(r8), allocatable :: NitriR(:)             ! 1/day
      real(r8), allocatable :: PARfrac(:)            ! nondimensional
      real(r8), allocatable :: PhyCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: PhyIP(:)              ! 1/mmol_N
      real(r8), allocatable :: PhyIS(:)              ! 1/(Watts m-2 day)
      real(r8), allocatable :: PhyMin(:)             ! mmol_N/m3
      real(r8), allocatable :: PhyMR(:)              ! 1/day
      real(r8), allocatable :: SDeAR(:)              ! 1/day
      real(r8), allocatable :: SDeBR(:)              ! 1/day
      real(r8), allocatable :: SDeRRN(:)             ! 1/day
      real(r8), allocatable :: SDeRRC(:)             ! 1/day
      real(r8), allocatable :: Vp0(:)                ! nondimensional
      real(r8), allocatable :: wLDet(:)              ! m/day
      real(r8), allocatable :: wPhy(:)               ! m/day
      real(r8), allocatable :: wSDet(:)              ! m/day
      real(r8), allocatable :: ZooAE_N(:)            ! nondimensional
      real(r8), allocatable :: ZooBM(:)              ! 1/day
      real(r8), allocatable :: ZooCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: ZooER(:)              ! 1/day
      real(r8), allocatable :: ZooGR(:)              ! 1/day
      real(r8), allocatable :: ZooMin(:)             ! mmol_N/m3
      real(r8), allocatable :: ZooMR(:)              ! 1/day
      real(r8), allocatable :: pCO2air(:)            ! ppmv
      
#ifdef SPECTRAL_LIGHT
      integer, parameter :: NBands = 60      ! spectral bands
!  Spectral band width used in light calculations.
      real(r8), parameter :: DLAM  = 5.0_r8
      integer, parameter :: Nphy = 1         ! Phytoplankton groups
      integer, parameter :: Npig = 1         ! Pigments
      real(r8), dimension(NBands) :: wavedp   ! a and b factor
# ifdef CDOM_VARIABLE
      integer, parameter :: Ndom = 2         ! DOM constituents
      integer, allocatable :: iCDMC(:)  ! Color degradational matter
      integer, parameter :: ilab=1    ! labile index for DOC.
      integer, parameter :: irct=2    ! relict index for DOC.
      character (len=11), dimension(Ndom) :: DomName
# endif
      real(r8), allocatable :: SIGATRB(:)             
      real(r8), allocatable :: STRB(:)                
      real(r8), allocatable :: BLTRB(:)             
      real(r8), allocatable :: SIGBTRB(:)             
      real(r8), allocatable :: ETASPEC(:)             
      real(r8), allocatable :: BB2B(:)             
#endif
#ifdef SAV_BIOMASS
      integer,  allocatable :: GMODopt(:)
      real(r8), allocatable :: KNSED(:)             
      real(r8), allocatable :: KNWC(:)             
      real(r8), allocatable :: TOPT(:)             
      real(r8), allocatable :: THTA(:)             
      real(r8), allocatable :: THTA2(:)             
      real(r8), allocatable :: SCL(:)             
      real(r8), allocatable :: SCL2(:)             
      real(r8), allocatable :: KI(:)             
      real(r8), allocatable :: SR(:)             
      real(r8), allocatable :: LMBAMX(:)             
      real(r8), allocatable :: KMAG(:)             
      real(r8), allocatable :: ARSC(:)             
      real(r8), allocatable :: ARC(:)             
      real(r8), allocatable :: BSRC(:)             
      real(r8), allocatable :: RC(:)             
      real(r8), allocatable :: RTStTL(:)             
      real(r8), allocatable :: DOWNt(:)             
      real(r8), allocatable :: TRNS(:)             
      real(r8), allocatable :: TCRIT(:)             
      real(r8), allocatable :: KM(:)             
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
#ifdef SPECTRAL_LIGHT
      integer :: iband
#endif
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
#ifdef CARBON
# ifdef OXYGEN
      NBT=12
# else
      NBT=11
# endif
#else
# ifdef OXYGEN
      NBT=8
# else
      NBT=7
# endif
#endif
# ifdef CDOM_VARIABLE
      NBT=NBT+Ndom
# endif

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio3d=2
      NDbio2d=0
# ifdef DENITRIFICATION
      NDbio2d=NDbio2d+1
# endif
# ifdef CARBON
      NDbio2d=NDbio2d+2
# endif
# ifdef OXYGEN
      NDbio2d=NDbio2d+1
# endif
!
!  Initialize biology diagnostic indices.
!
      ic=0
# ifdef DENITRIFICATION
      iDNIT=ic+1
      ic=ic+1
# endif
# ifdef CARBON
      iCOfx=ic+1
      ipCO2=ic+2
      ic=ic+2
# endif
# ifdef OXYGEN
      iO2fx=ic+1
# endif
#endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(AttChl)) THEN
        allocate ( AttChl(Ngrids) )
      END IF
      IF (.not.allocated(Chl2C_m)) THEN
        allocate ( Chl2C_m(Ngrids) )
      END IF
      IF (.not.allocated(ChlMin)) THEN
        allocate ( ChlMin(Ngrids) )
      END IF
      IF (.not.allocated(CoagR)) THEN
        allocate ( CoagR(Ngrids) )
      END IF
      IF (.not.allocated(D_p5NH4)) THEN
        allocate ( D_p5NH4(Ngrids) )
      END IF
      IF (.not.allocated(I_thNH4)) THEN
        allocate ( I_thNH4(Ngrids) )
      END IF
      IF (.not.allocated(K_NH4)) THEN
        allocate ( K_NH4(Ngrids) )
      END IF
      IF (.not.allocated(K_NO3)) THEN
        allocate ( K_NO3(Ngrids) )
      END IF
      IF (.not.allocated(K_Phy)) THEN
        allocate ( K_Phy(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRN)) THEN
        allocate ( LDeRRN(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRC)) THEN
        allocate ( LDeRRC(Ngrids) )
      END IF
      IF (.not.allocated(NitriR)) THEN
        allocate ( NitriR(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(PhyCN)) THEN
        allocate ( PhyCN(Ngrids) )
      END IF
      IF (.not.allocated(PhyIP)) THEN
        allocate ( PhyIP(Ngrids) )
      END IF
      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
      END IF
      IF (.not.allocated(PhyMin)) THEN
        allocate ( PhyMin(Ngrids) )
      END IF
      IF (.not.allocated(PhyMR)) THEN
        allocate ( PhyMR(Ngrids) )
      END IF
      IF (.not.allocated(SDeAR)) THEN
        allocate ( SDeAR(Ngrids) )
      END IF
      IF (.not.allocated(SDeBR)) THEN
        allocate ( SDeBR(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRN)) THEN
        allocate ( SDeRRN(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRC)) THEN
        allocate ( SDeRRC(Ngrids) )
      END IF
      IF (.not.allocated(Vp0)) THEN
        allocate ( Vp0(Ngrids) )
      END IF
      IF (.not.allocated(wLDet)) THEN
        allocate ( wLDet(Ngrids) )
      END IF
      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
      END IF
      IF (.not.allocated(wSDet)) THEN
        allocate ( wSDet(Ngrids) )
      END IF
      IF (.not.allocated(ZooAE_N)) THEN
        allocate ( ZooAE_N(Ngrids) )
      END IF
      IF (.not.allocated(ZooBM)) THEN
        allocate ( ZooBM(Ngrids) )
      END IF
      IF (.not.allocated(ZooCN)) THEN
        allocate ( ZooCN(Ngrids) )
      END IF
      IF (.not.allocated(ZooER)) THEN
        allocate ( ZooER(Ngrids) )
      END IF
      IF (.not.allocated(ZooGR)) THEN
        allocate ( ZooGR(Ngrids) )
      END IF
      IF (.not.allocated(ZooMin)) THEN
        allocate ( ZooMin(Ngrids) )
      END IF
      IF (.not.allocated(ZooMR)) THEN
        allocate ( ZooMR(Ngrids) )
      END IF
      IF (.not.allocated(pCO2air)) THEN
        allocate ( pCO2air(Ngrids) )
      END IF
#ifdef SPECTRAL_LIGHT
      IF (.not.allocated(SIGATRB)) THEN
        allocate ( SIGATRB(Ngrids) )
      END IF
      IF (.not.allocated(STRB)) THEN
        allocate ( STRB(Ngrids) )
      END IF
      IF (.not.allocated(BLTRB)) THEN
        allocate ( BLTRB(Ngrids) )
      END IF
      IF (.not.allocated(SIGBTRB)) THEN
        allocate ( SIGBTRB(Ngrids) )
      END IF
      IF (.not.allocated(ETASPEC)) THEN
        allocate ( ETASPEC(Ngrids) )
      END IF
      IF (.not.allocated(BB2B)) THEN
        allocate ( BB2B(Ngrids) )
      END IF
#endif
#ifdef SAV_BIOMASS
      IF (.not.allocated(GMODopt)) THEN
        allocate ( GMODopt(Ngrids) )
      END IF
      IF (.not.allocated(KNSED)) THEN
        allocate ( KNSED(Ngrids) )
      END IF
      IF (.not.allocated(KNWC)) THEN
        allocate ( KNWC(Ngrids) )
      END IF
      IF (.not.allocated(TOPT)) THEN
        allocate ( TOPT(Ngrids) )
      END IF
      IF (.not.allocated(THTA)) THEN
        allocate ( THTA(Ngrids) )
      END IF
      IF (.not.allocated(THTA2)) THEN
        allocate ( THTA2(Ngrids) )
      END IF
      IF (.not.allocated(SCL)) THEN
        allocate ( SCL(Ngrids) )
      END IF
      IF (.not.allocated(SCL2)) THEN
        allocate ( SCL2(Ngrids) )
      END IF
      IF (.not.allocated(KI)) THEN
        allocate ( KI(Ngrids) )
      END IF
      IF (.not.allocated(SR)) THEN
        allocate ( SR(Ngrids) )
      END IF
      IF (.not.allocated(LMBAMX)) THEN
        allocate ( LMBAMX(Ngrids) )
      END IF
      IF (.not.allocated(KMAG)) THEN
        allocate ( KMAG(Ngrids) )
      END IF
      IF (.not.allocated(ARSC)) THEN
        allocate ( ARSC(Ngrids) )
      END IF
      IF (.not.allocated(ARC)) THEN
        allocate ( ARC(Ngrids) )
      END IF
      IF (.not.allocated(BSRC)) THEN
        allocate ( BSRC(Ngrids) )
      END IF
      IF (.not.allocated(RC)) THEN
        allocate ( RC(Ngrids) )
      END IF
      IF (.not.allocated(RtStTL)) THEN
        allocate (RtStTL(Ngrids) )
      END IF
      IF (.not.allocated(DOWNt)) THEN
        allocate ( DOWNt(Ngrids) )
      END IF
      IF (.not.allocated(TRNS)) THEN
        allocate ( TRNS(Ngrids) )
      END IF
      IF (.not.allocated(TCRIT)) THEN
        allocate ( TCRIT(Ngrids) )
      END IF
      IF (.not.allocated(KM)) THEN
        allocate ( KM(Ngrids) )
      END IF
#endif
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
# ifdef CDOM_VARIABLE
      IF (.not.allocated(iCDMC)) THEN
        allocate ( iCDMC(Ndom) )
      END IF
# endif

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
      END IF
#endif
#ifdef SPECTRAL_LIGHT
!
!  Spectral dependency for scattering and backscattering.
!
      DO iband=1,NBands
        wavedp(iband)=(550.0_r8/(397.0_r8+REAL(iband,r8)*DLAM))
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iNO3_=ic+1
      iNH4_=ic+2
      iChlo=ic+3
      iPhyt=ic+4
      iZoop=ic+5
      iLDeN=ic+6
      iSDeN=ic+7
      ic=ic+7
# ifdef CDOM_VARIABLE
      DO i=1,Ndom
        iCDMC(i)=ic+1
        ic=ic+1
      END DO
# endif
# ifdef CARBON
      iLDeC=ic+1
      iSDeC=ic+2
      iTIC_=ic+3
      iTAlk=ic+4
      ic=ic+4
# endif
# ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
# endif

      RETURN
      END SUBROUTINE initialize_biology
