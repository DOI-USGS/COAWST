!
!svn $Id: fennel_mod.h 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!   K_PO4    Inverse half-saturation for Phytoplankton PO4 uptake      !
!              [m3/(mmol_P)].                                          !
!   K_Phy    Zooplankton half-saturation, squared constant for         !
!              ingestion [mmol_N/m3]^2.                                !
!   LDeRR    Large Detrital re-mineralization rate [1/day].            !
!   NitriR   Nitrification rate: oxidation of NH4 to NO3 [1/day].      !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   PhyCN    Phytoplankton Carbon:Nitrogen ratio [mol_C/mol_N].        !
!   R_P2N    Phytoplankton Phosphorus:Nitrogen ratio [mol_P/mol_N].    !
!   PhyIP    Phytoplankton NH4 inhibition parameter [1/(mmol_N)].      !
!   PhyIS    Phytoplankton, initial slope of the P-I curve             !
!              [1/(W m-2 day)].                                        !
!   ZooMin   Phytoplankton minimum threshold value [mmol_N/m3].        !
!   PhyMR    Phytoplankton mortality rate [1/day] to small detritus.   !
!   SDeAR    Small detritus aggregation rate into Large detritus       !
!              [1/day].                                                !
!   SDeBR    Small Detrital breakdown to NH4 rate [1/day].             !
!   SDeRR    Large Detrital re-mineralization rate [1/day].            !
!   RDeRR    River Detrital re-mineralization rate [1/day].            !
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
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iNO3_                  ! Nitrate concentration
      integer :: iNH4_                  ! Ammonium concentration
#ifdef PO4
      integer :: iPO4_                  ! Phosphate concentration
#endif
      integer :: iChlo                  ! Chlorophyll concentration
      integer :: iPhyt                  ! Phytoplankton concentration
      integer :: iZoop                  ! Zooplankton concentration
      integer :: iLDeN                  ! Large detritus N-concentration
      integer :: iSDeN                  ! Small detritus N-concentration
#ifdef RIVER_DON
      integer :: iRDeN                  ! River detritus N-concentration
#endif
#ifdef CARBON
      integer :: iLDeC                  ! Large detritus C-concentration
      integer :: iSDeC                  ! Small detritus C-concentration
      integer :: iTIC_                  ! Total inorganic carbon
      integer :: iTAlk                  ! Total alkalinity
# ifdef RIVER_DON
      integer :: iRDeC                  ! River detritus C-concentration
# endif
#endif
#ifdef OXYGEN
      integer :: iOxyg                  ! Dissolved oxygen concentration
#endif
#ifdef ODU
      integer :: iODU_                  ! Dissolved oxygen demand units
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
# if defined SEDBIO_COUP
      integer  :: isdO2                       !Seabed-sea diff. O2 flx
      integer  :: iseO2                       !Seabed-sea eros. O2 flx
      integer  :: isdNO                       !Seabed-sea diff. NO3 flx
      integer  :: iseNO                       !Seabed-sea eros. NO3 flx
      integer  :: isdNH                       !Seabed-sea diff. NH4 flx
      integer  :: iseNH                       !Seabed-sea eros. NH4 flx
      integer  :: isdOD                       !Seabed-sea diff. ODU flx
      integer  :: iseOD                       !Seabed-sea eros. ODU flx
# endif
!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms

      integer  :: iPPro = 1                   ! primary productivity
      integer  :: iNO3u = 2                   ! NO3 uptake
      integer  :: iNifx = 3                   ! Nitrification flux
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
      real(r8), allocatable :: K_PO4(:)              ! m3/mmol_P
      real(r8), allocatable :: K_Phy(:)              ! (mmol_N/m3)^2
      real(r8), allocatable :: LDeRRN(:)             ! 1/day
      real(r8), allocatable :: LDeRRC(:)             ! 1/day
      real(r8), allocatable :: NitriR(:)             ! 1/day
      real(r8), allocatable :: PARfrac(:)            ! nondimensional
      real(r8), allocatable :: PhyCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: R_P2N(:)              ! mol_P/mol_N
      real(r8), allocatable :: PhyIP(:)              ! 1/mmol_N
      real(r8), allocatable :: PhyIS(:)              ! 1/(Watts m-2 day)
      real(r8), allocatable :: PhyMin(:)             ! mmol_N/m3
      real(r8), allocatable :: PhyMR(:)              ! 1/day
      real(r8), allocatable :: SDeAR(:)              ! 1/day
      real(r8), allocatable :: SDeBR(:)              ! 1/day
      real(r8), allocatable :: SDeRRN(:)             ! 1/day
      real(r8), allocatable :: SDeRRC(:)             ! 1/day
      real(r8), allocatable :: RDeRRN(:)             ! 1/day
      real(r8), allocatable :: RDeRRC(:)             ! 1/day
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
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
#ifdef CARBON
# ifdef OXYGEN
#  if defined PO4 && defined RIVER_DON
      NBT=15
#  elif defined RIVER_DON && !defined PO4
      NBT=14
#  elif defined PO4 && !defined RIVER_DON
      NBT=13
#  else
      NBT=12
#  endif
# else
#  if defined PO4 && defined RIVER_DON
      NBT=14
#  elif defined RIVER_DON && !defined PO4
      NBT=13
#  elif defined PO4 && !defined RIVER_DON
      NBT=12
#  else
      NBT=11
#  endif
# endif
#else
# ifdef OXYGEN
#  if defined PO4 && defined RIVER_DON
      NBT=10
#  elif defined PO4 || defined RIVER_DON
      NBT=9
#  else
      NBT=8
#  endif
# else
#  if defined PO4 && defined RIVER_DON
      NBT=9
#  elif defined PO4 || defined RIVER_DON
      NBT=8
#  else
      NBT=7
#  endif
# endif
#endif
#ifdef ODU
      NBT=NBT+1
#endif

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio3d=3
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
# ifdef SEDBIO_COUP
      NDbio2d=NDbio2d+8
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
      ic=ic+1
# endif
# ifdef SEDBIO_COUP
      isdO2=ic+1
      iseO2=ic+2
      isdNO=ic+3
      iseNO=ic+4
      isdNH=ic+5
      iseNH=ic+6
      isdOD=ic+7
      iseOD=ic+8
      ic=ic+8
# endif
#endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(AttChl)) THEN
        allocate ( AttChl(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Chl2C_m)) THEN
        allocate ( Chl2C_m(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ChlMin)) THEN
        allocate ( ChlMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(CoagR)) THEN
        allocate ( CoagR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(D_p5NH4)) THEN
        allocate ( D_p5NH4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(I_thNH4)) THEN
        allocate ( I_thNH4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(K_NH4)) THEN
        allocate ( K_NH4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(K_NO3)) THEN
        allocate ( K_NO3(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(K_PO4)) THEN
        allocate ( K_PO4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(K_Phy)) THEN
        allocate ( K_Phy(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LDeRRN)) THEN
        allocate ( LDeRRN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(LDeRRC)) THEN
        allocate ( LDeRRC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(NitriR)) THEN
        allocate ( NitriR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyCN)) THEN
        allocate ( PhyCN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(R_P2N)) THEN
        allocate ( R_P2N(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyIP)) THEN
        allocate ( PhyIP(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyMin)) THEN
        allocate ( PhyMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyMR)) THEN
        allocate ( PhyMR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(SDeAR)) THEN
        allocate ( SDeAR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(SDeBR)) THEN
        allocate ( SDeBR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(SDeRRN)) THEN
        allocate ( SDeRRN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(SDeRRC)) THEN
        allocate ( SDeRRC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(RDeRRN)) THEN
        allocate ( RDeRRN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(RDeRRC)) THEN
        allocate ( RDeRRC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Vp0)) THEN
        allocate ( Vp0(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(wLDet)) THEN
        allocate ( wLDet(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(wSDet)) THEN
        allocate ( wSDet(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooAE_N)) THEN
        allocate ( ZooAE_N(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooBM)) THEN
        allocate ( ZooBM(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooCN)) THEN
        allocate ( ZooCN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooER)) THEN
        allocate ( ZooER(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooGR)) THEN
        allocate ( ZooGR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooMin)) THEN
        allocate ( ZooMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooMR)) THEN
        allocate ( ZooMR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(pCO2air)) THEN
        allocate ( pCO2air(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
        Dmem(1)=Dmem(1)+REAL(NBT,r8)
      END IF

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
        Dmem(1)=Dmem(1)+REAL(NDbio2d,r8)
      END IF

      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
        Dmem(1)=Dmem(1)+REAL(NDbio3d,r8)
      END IF
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
# ifdef RIVER_DON
      iRDeN=ic+1
      ic=ic+1
# endif
# ifdef PO4
      iPO4_=ic+1
      ic=ic+1
# endif
# ifdef CARBON
      iLDeC=ic+1
      iSDeC=ic+2
      iTIC_=ic+3
      iTAlk=ic+4
      ic=ic+4
#  ifdef RIVER_DON
      iRDeC=ic+1
      ic=ic+1
#  endif
# endif
# ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
# endif
# ifdef ODU
      iODU_=ic+1
      ic=ic+1
# endif

      RETURN
      END SUBROUTINE initialize_biology
