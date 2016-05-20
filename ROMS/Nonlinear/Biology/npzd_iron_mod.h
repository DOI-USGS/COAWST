!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      integer :: iFdis                  ! Available dissolved iron
#endif
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: BioIni(:,:)
      real(r8), allocatable :: AttPhy(:)       ! m2/mmole
      real(r8), allocatable :: AttSW(:)        ! 1/m
      real(r8), allocatable :: DetRR(:)        ! 1/day
      real(r8), allocatable :: K_NO3(:)        ! 1/(mmol/m3)
      real(r8), allocatable :: Ivlev(:)        ! nondimensional
      real(r8), allocatable :: PARfrac(:)      ! nondimensional
#ifdef IRON_LIMIT
      real(r8), allocatable :: A_Fe(:)         ! nondimensional
      real(r8), allocatable :: B_Fe(:)         ! 1/M-C
      real(r8), allocatable :: FeRR(:)         ! 1/day
# ifdef IRON_RELAX
      real(r8), allocatable :: FeHmin(:)       ! m
      real(r8), allocatable :: FeMax(:)        ! mmole/m3
      real(r8), allocatable :: FeNudgTime(:)   ! day
# endif
      real(r8), allocatable :: K_FeC(:)        ! muM-Fe/M-C
      real(r8), allocatable :: T_Fe(:)         ! day
#endif
#ifdef TANGENT
      real(r8), allocatable :: tl_PARfrac(:)   ! nondimensional
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_PARfrac(:)   ! nondimensional
#endif
      real(r8), allocatable :: PhyIS(:)        ! m2/W
      real(r8), allocatable :: PhyMRD(:)       ! 1/day
      real(r8), allocatable :: PhyMRN(:)       ! 1/day
      real(r8), allocatable :: Vm_NO3(:)       ! 1/day
      real(r8), allocatable :: wDet(:)         ! m/day
#ifdef TANGENT
      real(r8), allocatable :: tl_wDet(:)      ! m/day
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_wDet(:)      ! m/day
#endif
      real(r8), allocatable :: wPhy(:)         ! m/day
#ifdef TANGENT
      real(r8), allocatable :: tl_wPhy(:)      ! m/day
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_wPhy(:)      ! m/day
#endif
      real(r8), allocatable :: ZooEED(:)       ! nondimensional
      real(r8), allocatable :: ZooEEN(:)       ! nondimensional
      real(r8), allocatable :: ZooGR(:)        ! 1/day
      real(r8), allocatable :: ZooMRD(:)       ! 1/day
      real(r8), allocatable :: ZooMRN(:)       ! 1/day

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
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AttPhy)) THEN
        allocate ( AttPhy(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(DetRR)) THEN
        allocate ( DetRR(Ngrids) )
      END IF
      IF (.not.allocated(K_NO3)) THEN
        allocate ( K_NO3(Ngrids) )
      END IF
      IF (.not.allocated(Ivlev)) THEN
        allocate ( Ivlev(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
#ifdef IRON_LIMIT
      IF (.not.allocated(A_Fe)) THEN
        allocate ( A_Fe(Ngrids) )
      END IF
      IF (.not.allocated(B_Fe)) THEN
        allocate ( B_Fe(Ngrids) )
      END IF
      IF (.not.allocated(FeRR)) THEN
        allocate ( FeRR(Ngrids) )
      END IF
# ifdef IRON_RELAX
      IF (.not.allocated(FeHmin)) THEN
        allocate ( FeHmin(Ngrids) )
      END IF
      IF (.not.allocated(FeMax)) THEN
        allocate ( FeMax(Ngrids) )
      END IF
      IF (.not.allocated(FeNudgTime)) THEN
        allocate ( FeNudgTime(Ngrids) )
      END IF
# endif
      IF (.not.allocated(K_FeC)) THEN
        allocate ( K_FeC(Ngrids) )
      END IF
      IF (.not.allocated(T_Fe)) THEN
        allocate ( T_Fe(Ngrids) )
      END IF
#endif
#ifdef TANGENT
      IF (.not.allocated(tl_PARfrac)) THEN
        allocate ( tl_PARfrac(Ngrids) )
      END IF
#endif
#ifdef ADJOINT
      IF (.not.allocated(ad_PARfrac)) THEN
        allocate ( ad_PARfrac(Ngrids) )
      END IF
#endif
      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
      END IF
      IF (.not.allocated(PhyMRD)) THEN
        allocate ( PhyMRD(Ngrids) )
      END IF
      IF (.not.allocated(PhyMRN)) THEN
        allocate ( PhyMRN(Ngrids) )
      END IF
      IF (.not.allocated(Vm_NO3)) THEN
        allocate ( Vm_NO3(Ngrids) )
      END IF
      IF (.not.allocated(wDet)) THEN
        allocate ( wDet(Ngrids) )
      END IF
#ifdef TANGENT
      IF (.not.allocated(tl_wDet)) THEN
        allocate ( tl_wDet(Ngrids) )
      END IF
#endif
#ifdef ADJOINT
      IF (.not.allocated(ad_wDet)) THEN
        allocate ( ad_wDet(Ngrids) )
      END IF
#endif
      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
      END IF
#ifdef TANGENT
      IF (.not.allocated(tl_wPhy)) THEN
        allocate ( tl_wPhy(Ngrids) )
      END IF
#endif
#ifdef ADJOINT
      IF (.not.allocated(ad_wPhy)) THEN
        allocate ( ad_wPhy(Ngrids) )
      END IF
#endif
      IF (.not.allocated(ZooEED)) THEN
        allocate ( ZooEED(Ngrids) )
      END IF
      IF (.not.allocated(ZooEEN)) THEN
        allocate ( ZooEEN(Ngrids) )
      END IF
      IF (.not.allocated(ZooGR)) THEN
        allocate ( ZooGR(Ngrids) )
      END IF
      IF (.not.allocated(ZooMRD)) THEN
        allocate ( ZooMRD(Ngrids) )
      END IF
      IF (.not.allocated(ZooMRN)) THEN
        allocate ( ZooMRN(Ngrids) )
      END IF
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
