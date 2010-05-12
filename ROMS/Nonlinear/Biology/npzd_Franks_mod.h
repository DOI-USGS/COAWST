!
!svn $Id: npzd_Franks_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
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
      integer, dimension(Ngrids) :: BioIter

#ifdef ANA_BIOLOGY
      real(r8), allocatable :: BioIni(:,:)
#endif
      real(r8), dimension(Ngrids) :: DetRR           ! 1/day
      real(r8), dimension(Ngrids) :: K_ext           ! 1/m
      real(r8), dimension(Ngrids) :: K_NO3           ! 1/(mmol/m3)
      real(r8), dimension(Ngrids) :: K_Phy           ! mmol/m3
      real(r8), dimension(Ngrids) :: PhyMR           ! 1/day
      real(r8), dimension(Ngrids) :: Vm_NO3          ! 1/day
      real(r8), dimension(Ngrids) :: wDet            ! m/day
#ifdef TANGENT
      real(r8), dimension(Ngrids) :: tl_wDet
#endif
#ifdef ADJOINT
      real(r8), dimension(Ngrids) :: ad_wDet
#endif
      real(r8), dimension(Ngrids) :: ZooGR           ! 1/day
      real(r8), dimension(Ngrids) :: ZooMR           ! 1/day
      real(r8), dimension(Ngrids) :: ZooMD           ! 1/day
      real(r8), dimension(Ngrids) :: ZooGA           ! nondimensional
      real(r8), dimension(Ngrids) :: ZooEC           ! nondimensional

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
      NBT=4
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

      RETURN
      END SUBROUTINE initialize_biology
