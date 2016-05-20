!
!svn $Id: sediment_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for sediment model:                                      !
!  =============================                                       !
!                                                                      !
!   Csed            Sediment concentration (kg/m3), used during        !
!                     analytical initialization.                       !
!   Erate           Surface erosion rate (kg/m2/s).                    !
!   Sd50            Median sediment grain diameter (m).                !
!   Srho            Sediment grain density (kg/m3).                    !
!   SedIter         Maximum number of iterations.                      !
!   Wsed            Particle settling velocity (m/s).                  !
!   poros           Porosity (non-dimensional: 0.0-1.0):               !
!                     Vwater/(Vwater+Vsed).                            !
!   tau_ce          Kinematic critical shear for erosion (m2/s2).      !
!   tau_cd          Kinematic critical shear for deposition (m2/s2).   !
!                                                                      !
!   bedload_coeff   Bedload rate coefficient (nondimensional).         !
!   minlayer_thick  Minimum thickness for 2nd layer (m).               !
!   newlayer_thick  New layer deposit thickness criteria (m).          !
!   morph_fac       Morphological scale factor (nondimensional).       !
!                                                                      !
!  BED properties indices:                                             !
!  ======================                                              !
!                                                                      !
!   MBEDP           Number of bed properties (array dimension).        !
!   idBmas(:)       Sediment mass index.                               !
!   idSbed(:)       IO indices for bed properties variables.           !
!   idfrac(:)       Sediment class fraction (non-dimensional).         !
!   ithck           Sediment layer thickness (m).                      !
!   iaged           Sediment layer age (s).                            !
!   iporo           Sediment layer porosity (non-dimensional).         !
!   idiff           Sediment layer bio-diffusivity (m2/s).             !
!   ibtcr           Sediment critical stress for erosion (Pa).         !
!                                                                      !
!  BOTTOM properties indices:                                          !
!  =========================                                           !
!                                                                      !
!   MBOTP           Number of bottom properties (array dimension).     !
!   idBott(:)       IO indices for bottom properties variables.        !
!   isd50           Median sediment grain diameter (m).                !
!   idens           Median sediment grain density (kg/m3).             !
!   iwsed           Mean settling velocity (m/s).                      !
!   itauc           Mean critical erosion stress (m2/s2).              !
!   irlen           Sediment ripple length (m).                        !
!   irhgt           Sediment ripple height (m).                        !
!   ibwav           Bed wave excursion amplitude (m).                  !
!   izdef           Default bottom roughness (m).                      !
!   izapp           Apparent bottom roughness (m).                     !
!   izNik           Nikuradse bottom roughness (m).                    !
!   izbio           Biological bottom roughness (m).                   !
!   izbfm           Bed form bottom roughness (m).                     !
!   izbld           Bed load bottom roughness (m).                     !
!   izwbl           Bottom roughness used wave BBL (m).                !
!   iactv           Active layer thickness for erosive potential (m).  !
!   ishgt           Sediment saltation height (m).                     !
!   imaxD           Maximum inundation depth.                          !
!   isgrH           Seagrass height.                                   !
!   isgrD           Seagrass shoot density.                            !
!   idnet           Erosion or deposition.                             !
!   idoff           Offset for calculation of dmix erodibility         !
!                     profile (m).                                     !
!   idslp           Slope  for calculation of dmix or erodibility      !
!                     profile.                                         !
!   idtim           Time scale for restoring erodibility profile (s).  !
!   idbmx           Bed biodifusivity maximum.                         !
!   idbmm           Bed biodifusivity minimum.                         !
!   idbzs           Bed biodifusivity zs.                              !
!   idbzm           Bed biodifusivity zm.                              !
!   idbzp           Bed biodifusivity phi.                             !
!   idprp           Cohesive behavior.                                 !
!                                                                      !
!   nTbiom          Number of hours for depth integration              !
!=======================================================================
!
      USE mod_param
!
      implicit none

      integer, allocatable  :: idSbed(:)  ! bed    properties IDs
      integer, allocatable  :: idBott(:)  ! bottom properties IDs
!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!
      integer, allocatable :: idsed(:)    ! Cohesive and non-cohesive
      integer, allocatable :: idmud(:)    ! Cohesive sediment
      integer, allocatable :: isand(:)    ! Non-cohesive sediment
!
!-----------------------------------------------------------------------
!  Set bed property variables 
!-----------------------------------------------------------------------
!
      integer :: MBEDP                     ! Number of bed properties  
      integer :: ithck, iaged, iporo, idiff 
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      integer :: ibtcr 
#endif
!
!-----------------------------------------------------------------------
!  Set bottom property variables 
!-----------------------------------------------------------------------
!
      integer :: MBOTP                     ! Number of bottom properties
      integer :: isd50, idens, iwsed, itauc
      integer :: irlen, irhgt, ibwav, izdef
      integer :: izapp, izNik, izbio, izbfm
      integer :: izbld, izwbl, iactv, ishgt
      integer :: imaxD, idnet  
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      integer :: idoff, idslp, idtim, idbmx
      integer :: idbmm, idbzs, idbzm, idbzp 
#endif 
#if defined MIXED_BED
      integer :: idprp 
#endif
!
!  Sediment metadata indices vectors.
!
      integer, allocatable :: idBmas(:)    ! class mass indices
      integer, allocatable :: idfrac(:)    ! class fraction indices
      integer, allocatable :: idUbld(:)    ! bed load u-points
      integer, allocatable :: idVbld(:)    ! bed load v-points
!
!-----------------------------------------------------------------------
!  Input sediment parameters.
!-----------------------------------------------------------------------
!
      real(r8), allocatable :: newlayer_thick(:)   ! deposit thickness criteria
      real(r8), allocatable :: minlayer_thick(:)   ! 2nd layer thickness criteria
      real(r8), allocatable :: bedload_coeff(:)    ! bedload rate coefficient

      real(r8), allocatable :: Csed(:,:)       ! initial concentration
      real(r8), allocatable :: Erate(:,:)      ! erosion rate
      real(r8), allocatable :: Sd50(:,:)       ! mediam grain diameter
      real(r8), allocatable :: Srho(:,:)       ! grain density
      real(r8), allocatable :: Wsed(:,:)       ! settling velocity
      real(r8), allocatable :: poros(:,:)      ! porosity
      real(r8), allocatable :: tau_ce(:,:)     ! shear for erosion
      real(r8), allocatable :: tau_cd(:,:)     ! shear for deposition
      real(r8), allocatable :: morph_fac(:,:)  ! morphological factor

#if defined COHESIVE_BED || defined MIXED_BED
      real(r8), allocatable :: tcr_min(:)      ! minimum shear for erosion
      real(r8), allocatable :: tcr_max(:)      ! maximum shear for erosion
      real(r8), allocatable :: tcr_slp(:)      ! Tau_crit profile slope
      real(r8), allocatable :: tcr_off(:)      ! Tau_crit profile offset
      real(r8), allocatable :: tcr_tim(:)      ! Tau_crit consolidation rate
#endif

#if defined MIXED_BED
      real(r8), allocatable :: transC(:)       ! cohesive transition
      real(r8), allocatable :: transN(:)       ! noncohesive transition
#endif
#if defined SED_BIODIFF
      real(r8), allocatable :: Dbmx(:)         ! Dbmax  Maximum biodiffusivity
      real(r8), allocatable :: Dbmm(:)         ! Dbmin  Minimum biodiffusivity
      real(r8), allocatable :: Dbzs(:)         ! Dbzs   Depth of maximum biodiff
      real(r8), allocatable :: Dbzm(:)         ! Dbzm   Depth end exp biodiff
      real(r8), allocatable :: Dbzp(:)         ! Dbzp   Depth of minimum biodiff
#endif
#if defined SED_FLOCS && defined SED_DEFLOC
      real(r8), allocatable :: mud_frac_eq(:,:) ! Equilibrium fractional class distribution
      real(r8), allocatable :: t_dfloc(:)       ! Time scale of bed deflocculation
#endif

      CONTAINS

      SUBROUTINE initialize_sediment
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the sediment model.   !
!  It allocates and assigns sediment tracers indices.                  !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
      integer :: counter1, counter2 
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Set bed properties indices.
!-----------------------------------------------------------------------
!
      counter1 = 1           ! Initializing counter 
      ithck    = counter1    ! layer thickness
      counter1 = counter1+1
      iaged    = counter1    ! layer age 
      counter1 = counter1+1 
      iporo    = counter1    ! layer porosity 
      counter1 = counter1+1
      idiff    = counter1    ! layer bio-diffusivity 
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      counter1 = counter1+1 
      ibtcr    = counter1    ! layer critical stress 
#endif 
!
!-----------------------------------------------------------------------
!  Set bottom properties indices.
!-----------------------------------------------------------------------
! 
      counter2 = 1           ! Initializing counter 
      isd50    = counter2    ! Median sediment grain diameter (m).
      counter2 = counter2+1
      idens    = counter2    ! Median sediment grain density (kg/m3).
      counter2 = counter2+1
      iwsed    = counter2    ! Mean settling velocity (m/s).  
      counter2 = counter2+1
      itauc    = counter2    ! Mean critical erosion stress (m2/s2).
      counter2 = counter2+1
      irlen    = counter2    ! Sediment ripple length (m).
      counter2 = counter2+1
      irhgt    = counter2    ! Sediment ripple height (m).
      counter2 = counter2+1
      ibwav    = counter2    ! Bed wave excursion amplitude (m).
      counter2 = counter2+1
      izdef    = counter2    ! Default bottom roughness (m).
      counter2 = counter2+1
      izapp    = counter2    ! Apparent bottom roughness (m).
      counter2 = counter2+1
      izNik    = counter2    ! Nikuradse bottom roughness (m).
      counter2 = counter2+1  
      izbio    = counter2    ! Biological bottom roughness (m).
      counter2 = counter2+1  
      izbfm    = counter2    ! Bed form bottom roughness (m).
      counter2 = counter2+1
      izbld    = counter2    ! Bed load bottom roughness (m).
      counter2 = counter2+1
      izwbl    = counter2    ! Bottom roughness used wave BBL (m).
      counter2 = counter2+1  
      iactv    = counter2    ! Active layer thickness for erosive potential (m).
      counter2 = counter2+1  
      ishgt    = counter2    ! Sediment saltation height (m).
      counter2 = counter2+1   
      imaxD    = counter2    ! Maximum inundation depth.
      counter2 = counter2+1
      idnet    = counter2    ! Erosion/deposition
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      counter2 = counter2+1
      idoff    = counter2    ! Offset for calculation of dmix erodibility profile (m). 
      counter2 = counter2+1                                 
      idslp    = counter2    ! Slope  for calculation of dmix or erodibility profile.  
      counter2 = counter2+1                                       
      idtim    = counter2    ! Time scale for restoring erodibility profile (s).
      counter2 = counter2+1
      idbmx    = counter2    ! Bed biodifusivity maximum.
      counter2 = counter2+1    
      idbmm    = counter2    ! Bed biodifusivity minimum.
      counter2 = counter2+1
      idbzs    = counter2    ! Bed biodifusivity zs.
      counter2 = counter2+1     
      idbzm    = counter2    ! Bed biodifusivity zm.
      counter2 = counter2+1
      idbzp    = counter2    ! Bed biodifusivity phi.
#endif 
#if defined MIXED_BED
      counter2 = counter2+1
      idprp    = counter2    ! Cohesive behavior. 
#endif 
!
!  Allocate bed & bottom properties 
!
      MBEDP   = counter1
      IF (.not.allocated(idSbed)) THEN
        allocate ( idSbed(MBEDP) )
      END IF
!
      MBOTP   = counter2 
      IF (.not.allocated(idBott)) THEN
        allocate ( idBott(MBOTP) )
      END IF
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate various nested grid depended parameters
!
      IF (.not.allocated(newlayer_thick)) THEN
        allocate ( newlayer_thick(Ngrids) )
        newlayer_thick = IniVal
      END IF
      IF (.not.allocated(minlayer_thick)) THEN
        allocate ( minlayer_thick(Ngrids) )
        minlayer_thick = IniVal
      END IF
      IF (.not.allocated(bedload_coeff)) THEN
        allocate ( bedload_coeff(Ngrids) )
        bedload_coeff = IniVal
      END IF

#if defined COHESIVE_BED || defined MIXED_BED
      IF (.not.allocated(tcr_min)) THEN
        allocate ( tcr_min(Ngrids) )
        tcr_min = IniVal
      END IF
      IF (.not.allocated(tcr_max)) THEN
        allocate ( tcr_max(Ngrids) )
        tcr_max = IniVal
      END IF
      IF (.not.allocated(tcr_slp)) THEN
        allocate ( tcr_slp(Ngrids) )
        tcr_slp = IniVal
      END IF
      IF (.not.allocated(tcr_off)) THEN
        allocate ( tcr_off(Ngrids) )
        tcr_off = IniVal
      END IF
      IF (.not.allocated(tcr_tim)) THEN
        allocate ( tcr_tim(Ngrids) )
        tcr_tim = IniVal
      END IF
#endif

#if defined SED_FLOCS && defined SED_DEFLOC
      IF (.not.allocated(t_dfloc)) THEN
        allocate ( t_dfloc(Ngrids) )
        t_dfloc = IniVal
      END IF
#endif

#if defined MIXED_BED
      IF (.not.allocated(transC)) THEN
        allocate ( transC(Ngrids) )
        transC = IniVal
      END IF
      IF (.not.allocated(transN)) THEN
        allocate ( transN(Ngrids) )
        transN = IniVal
      END IF
#endif
!
!  Allocate sediment tracers indices vectors.
!
      IF (.not.allocated(idsed)) THEN
        allocate ( idsed(MAX(1,NST)) )
      END IF
      IF (.not.allocated(idmud)) THEN
        allocate ( idmud(MAX(1,NCS)) )
      END IF
      IF (.not.allocated(isand)) THEN
        allocate ( isand(MAX(1,NNS)) )
      END IF
      IF (.not.allocated(idBmas)) THEN
        allocate ( idBmas(NST) )
      END IF
      IF (.not.allocated(idfrac)) THEN
        allocate ( idfrac(NST) )
      END IF
      IF (.not.allocated(idUbld)) THEN
        allocate ( idUbld(NST) )
      END IF
      IF (.not.allocated(idVbld)) THEN
        allocate ( idVbld(NST) )
      END IF
!
!  Set cohesive and noncohesive suspended sediment tracers
!  identification indices.
!
      ic=NAT+NPT
      DO i=1,NCS
        ic=ic+1
        idmud(i)=ic
        idsed(i)=idmud(i)
      END DO
      DO i=1,NNS
        ic=ic+1
        isand(i)=ic
        idsed(NCS+i)=isand(i)
      END DO

      RETURN
      END SUBROUTINE initialize_sediment
