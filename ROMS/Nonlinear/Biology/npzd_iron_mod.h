!
!svn $Id: npzd_iron_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Powell et al. (2006) ecosystem model with iron       !
!  limitation:                                                         !
!                                                                      !
!  AttPhy     Light attenuation due to phytoplankton (self-shading     !
!               coefficient), [m2/millimole_N].                        !
!  AttSW      Light attenuation due to sea water, [1/m].               !
!  A_Fe       Empirical Fe:C power, [nondimensional].                  !
!  B_Fe       Empirical Fe:C coefficient, [1/M-C].                     !
!  BioIter    Maximum number of iterations to achieve convergence of   !
!               the nonlinear solution.                                !
!  BioIni     Initial concentration for analytical initial (uniform)   !
!               conditions.                                            !
!  DetRR      Detritus remineraliztion rate, [1/day].                  !
!  FeRR       Fe remineralization rate, [1/day].                       !
!  FeHmin     Minimum bathymetry value (meter; positive) considered to !
!               nudge dissolved iron over the shelf (h <= FeHmin).     !
!  FeMax      Dissolved iron value to nudge over the shelf to simulate !
!               Fe coastal source.                                     !
!  FeNudgTime Dissolved iron nudging time scale (days) over the shelf. !
!               Inverse scale will be computed internally.             !
!  K_FeC      Fe:C at F=0.5, [muM-Fe/M-C].                             !
!  K_NO3      Inverse half-saturation for phytoplankton nitrate uptake !
!               [1/(millimole_N m-3)].                                 !
!  Ivlev      Ivlev constant for zooplankton grazin parameterization,  !
!               [nondimensional].                                      !
!  PARfrac    Fraction of shortwave radiation that is available for    !
!               photosyntesis [non-dimensional].                       !
!  PhyIS      Phytoplankton, initial slope of the P-I curve [m2/W].    !
!  PhyMRD     Phytoplankton mortality rate to the Detritus pool,       !
!               [1/day].                                               !
!  PhyMRN     Phytoplankton mortality rate to the Nitrogen pool,       !
!              [1/day].                                                !
!  T_Fe       Iron uptake timescale, [day].                            !
!  Vm_NO3     Nitrate uptake rate, [1/day].                            !
!  wDet       Detrital sinking rate, [m/day].                          !
!  wPhy       Phytoplankton sinking rate, [m/day].                     !
!  ZooEED     Zooplankton excretion efficiency to Detritus pool,       !
!               [nondimensional].                                      !
!  ZooEEN     Zooplankton excretion efficiency to Nitrogen pool,       !
!               [nondimensional].                                      !
!  ZooGR      Zooplankton grazing rate, [1/day].                       !
!  ZooMRD     Zooplankton mortality rate to Detritus pool, [1/day].    !
!  ZooMRN     Zooplankton mortality rate to Nitrogen pool, [1/day].    !
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
      integer :: iPhyt                  ! Phytoplankton concentration
      integer :: iZoop                  ! Zooplankton concentration
      integer :: iSDet                  ! Small detritus concentration
#ifdef IRON_LIMIT
      integer :: iFphy                  ! Phytoplankton-associated iron
      integer :: iFdis                  ! Available disolved iron
#endif
!
!  Biological parameters.
!
      integer, dimension(Ngrids) :: BioIter

      real(r8), allocatable :: BioIni(:,:)
      real(r8), dimension(Ngrids) :: AttPhy          ! m2/mmole
      real(r8), dimension(Ngrids) :: AttSW           ! 1/m
      real(r8), dimension(Ngrids) :: DetRR           ! 1/day
      real(r8), dimension(Ngrids) :: K_NO3           ! 1/(mmol/m3)
      real(r8), dimension(Ngrids) :: Ivlev           ! nondimensional
      real(r8), dimension(Ngrids) :: PARfrac         ! nondimensional
#ifdef IRON_LIMIT
      real(r8), dimension(Ngrids) :: A_Fe            ! nondimensional
      real(r8), dimension(Ngrids) :: B_Fe            ! 1/M-C
      real(r8), dimension(Ngrids) :: FeRR            ! 1/day
# ifdef IRON_RELAX
      real(r8), dimension(Ngrids) :: FeHmin          ! m
      real(r8), dimension(Ngrids) :: FeMax           ! mmole/m3
      real(r8), dimension(Ngrids) :: FeNudgTime      ! day
# endif
      real(r8), dimension(Ngrids) :: K_FeC           ! muM-Fe/M-C
      real(r8), dimension(Ngrids) :: T_Fe            ! day
#endif
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_PARfrac      ! nondimensional
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_PARfrac      ! nondimensional
#endif
      real(r8), dimension(Ngrids) :: PhyIS           ! m2/W
      real(r8), dimension(Ngrids) :: PhyMRD          ! 1/day
      real(r8), dimension(Ngrids) :: PhyMRN          ! 1/day
      real(r8), dimension(Ngrids) :: Vm_NO3          ! 1/day
      real(r8), dimension(Ngrids) :: wDet            ! m/day
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_wDet         ! m/day
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_wDet         ! m/day
#endif
      real(r8), dimension(Ngrids) :: wPhy            ! m/day
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_wPhy         ! m/day
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_wPhy         ! m/day
#endif
      real(r8), dimension(Ngrids) :: ZooEED          ! nondimensional
      real(r8), dimension(Ngrids) :: ZooEEN          ! nondimensional
      real(r8), dimension(Ngrids) :: ZooGR           ! 1/day
      real(r8), dimension(Ngrids) :: ZooMRD          ! 1/day
      real(r8), dimension(Ngrids) :: ZooMRN          ! 1/day

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
#ifdef IRON_LIMIT
      NBT=6
#else
      NBT=4
#endif
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
      iNO3_=ic+1
      iPhyt=ic+2
      iZoop=ic+3
      iSDet=ic+4
#ifdef IRON_LIMIT
      iFphy=ic+5
      iFdis=ic+6
#endif

      RETURN
      END SUBROUTINE initialize_biology
