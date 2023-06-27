      MODULE mod_biology
!
!git $Id$
!svn $Id: npzd_Powell_mod.h 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Powell et al. (2006) ecosystem model:                !
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
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

#ifdef ANA_BIOLOGY
      real(r8), allocatable :: BioIni(:,:)
#endif
      real(r8), allocatable :: AttPhy(:)       ! m2/mmole
      real(r8), allocatable :: AttSW(:)        ! 1/m
      real(r8), allocatable :: DetRR(:)        ! 1/day
      real(r8), allocatable :: K_NO3(:)        ! mmol/m3
      real(r8), allocatable :: Ivlev(:)        ! nondimensional
      real(r8), allocatable :: PARfrac(:)      ! nondimensional
#ifdef TANGENT
      real(r8), allocatable :: tl_PARfrac(:)
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_PARfrac(:)
#endif
      real(r8), allocatable :: PhyIS(:)        ! m2/W
      real(r8), allocatable :: PhyMRD(:)       ! 1/day
      real(r8), allocatable :: PhyMRN(:)       ! 1/day
      real(r8), allocatable :: Vm_NO3(:)       ! 1/day
      real(r8), allocatable :: wDet(:)         ! m/day
#ifdef TANGENT
      real(r8), allocatable :: tl_wDet(:)
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_wDet(:)
#endif
      real(r8), allocatable :: wPhy(:)         ! m/day
#ifdef TANGENT
      real(r8), allocatable :: tl_wPhy(:)
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_wPhy(:)
#endif
      real(r8), allocatable :: ZooEED(:)       ! nondimensional
      real(r8), allocatable :: ZooEEN(:)       ! nondimensional
      real(r8), allocatable :: ZooGR(:)        ! 1/day
      real(r8), allocatable :: ZooMRD(:)       ! 1/day
      real(r8), allocatable :: ZooMRN(:)       ! 1/day
!
      CONTAINS
!
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
      NBT=4
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(AttPhy)) THEN
        allocate ( AttPhy(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(DetRR)) THEN
        allocate ( DetRR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(K_NO3)) THEN
        allocate ( K_NO3(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Ivlev)) THEN
        allocate ( Ivlev(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef TANGENT
      IF (.not.allocated(tl_PARfrac)) THEN
        allocate ( tl_PARfrac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

#ifdef ADJOINT
      IF (.not.allocated(ad_PARfrac)) THEN
        allocate ( ad_PARfrac(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyMRD)) THEN
        allocate ( PhyMRD(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(PhyMRN)) THEN
        allocate ( PhyMRN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Vm_NO3)) THEN
        allocate ( Vm_NO3(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(wDet)) THEN
        allocate ( wDet(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef TANGENT
      IF (.not.allocated(tl_wDet)) THEN
        allocate ( tl_wDet(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

#ifdef ADJOINT
      IF (.not.allocated(ad_wDet)) THEN
        allocate ( ad_wDet(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#ifdef TANGENT
      IF (.not.allocated(tl_wPhy)) THEN
        allocate ( tl_wPhy(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

#ifdef ADJOINT
      IF (.not.allocated(ad_wPhy)) THEN
        allocate ( ad_wPhy(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

      IF (.not.allocated(ZooEED)) THEN
        allocate ( ZooEED(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooEEN)) THEN
        allocate ( ZooEEN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooGR)) THEN
        allocate ( ZooGR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooMRD)) THEN
        allocate ( ZooMRD(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(ZooMRN)) THEN
        allocate ( ZooMRN(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
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
!
      RETURN
      END SUBROUTINE initialize_biology

      END MODULE mod_biology
