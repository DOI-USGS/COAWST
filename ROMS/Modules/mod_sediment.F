#include "cppdefs.h"
      MODULE mod_sediment
#if defined SEDIMENT || defined BBL_MODEL
!
!svn $Id: mod_sediment.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for sediment model:                                      !
!                                                                      !
!   Csed     Sediment concentration (kg/m3), used during analytical    !
!              initialization.                                         !
!   Erate    Surface erosion rate (kg/m2/s).                           !
!   Sd50     Median sediment grain diameter (m).                       !
!   Srho     Sediment grain density (kg/m3).                           !
!   SedIter  Maximum number of iterations.                             !
!   Wsed     Particle settling velocity (m/s).                         !
!   poros    Porosity (non-dimensional: 0.0-1.0):                      !
!              Vwater/(Vwater+Vsed).                                   !
!   tau_ce   Kinematic critical shear for erosion (m2/s2).             !
!   tau_cd   Kinematic critical shear for deposition (m2/s2).          !
!   morph_fac  Morphological scale factor (nondimensional).             !
!   newlayer_thick    New layer deposit thickness criteria (m).        !
!   minlayer_thick    Minimum thickness for 2nd layer (m)              !
!   bedload_coeff     Bedload rate coefficient (nondimensional).       !
!                                                                      !
!  BED properties indices:                                             !
!                                                                      !
!   MBEDP    Number of bed properties (array dimension).               !
!   idBmas   Sediment mass index.                                      !
!   idSbed   IO indices for bed properties variables.                  !
!   idfrac   sediment class fraction (non-dimensional).                !
!   ithck    Sediment layer thickness (m).                             !
!   iaged    Sediment layer age (s).                                   !
!   iporo    Sediment layer porosity (non-dimensional).                !
!   idiff    Sediment layer bio-diffusivity (m2/s).                    !
!   ibtcr    Sediment critical stress for erosion (Pa)
!                                                                      !
!  BOTTOM properties indices:                                          !
!                                                                      !
!   MBOTP    Number of bottom properties (array dimension).            !
!   idBott   IO indices for bottom properties variables.               !
!   isd50    Median sediment grain diameter (m).                       !
!   idens    Median sediment grain density (kg/m3).                    !
!   iwsed    Mean settling velocity (m/s).                             !
!   itauc    Mean critical erosion stress (m2/s2).                     !
!   irlen    Sediment ripple length (m).                               !
!   irhgt    Sediment ripple height (m).                               !
!   ibwav    Bed wave excursion amplitude (m).                         !
!   izdef    Default bottom roughness (m).                             !
!   izapp    Apparent bottom roughness (m).                            !
!   izNik    Nikuradse bottom roughness (m).                           !
!   izbio    Biological bottom roughness (m).                          !
!   izbfm    Bed form bottom roughness (m).                            !
!   izbld    Bed load bottom roughness (m).                            !
!   izwbl    Bottom roughness used wave BBL (m).                       ! 
!   iactv    Active layer thickness for erosive potential (m).         !
!   ishgt    Sediment saltation height (m).                            !
!   idflux   erosion flux                                              !
!   idnet    erosion or deposition                                     !
!   idoff    Offset for calculation of dmix erodibility profile (m)    !
!   idslp    Slope for calculation of dmix or erodibility profile      !
!   idtim    Time scale for restoring erodibility profile (s)          !
!   idbmx    bed biodifusivity max                                     !
!   idbmm    bed biodifusivity m                                       !
!   idbzs    bed biodifusivity zs                                      !
!   idbzm    bed biodifusivity zm                                      !
!   idbzp    bed biodifusivity phi                                     !
!   idprp    cohesive behavior                                         !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
        integer, parameter :: MBEDP = 5    ! layer tau crit
# else
        integer, parameter :: MBEDP = 4    ! Bed Properties:
# endif
        integer, parameter :: ithck = 1    ! layer thickness
        integer, parameter :: iaged = 2    ! layer age
        integer, parameter :: iporo = 3    ! layer porosity
        integer, parameter :: idiff = 4    ! layer bio-diffusivity
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
        integer, parameter :: ibtcr = 5    ! layer tau crit
# endif
# if defined MIXED_BED
        integer, parameter :: MBOTP = 27   ! Bottom Properties:
# elif defined COHESIVE_BED || defined SED_BIODIFF
        integer, parameter :: MBOTP = 26
# else
        integer, parameter :: MBOTP = 18   ! Bottom Properties:
# endif
        integer, parameter :: isd50 = 1    ! mean grain diameter
        integer, parameter :: idens = 2    ! mean grain density
        integer, parameter :: iwsed = 3    ! mean settle velocity
        integer, parameter :: itauc = 4    ! critical erosion stress
        integer, parameter :: irlen = 5    ! ripple length
        integer, parameter :: irhgt = 6    ! ripple height
        integer, parameter :: ibwav = 7    ! wave excursion amplitude
        integer, parameter :: izdef = 8    ! default bottom roughness
        integer, parameter :: izapp = 9    ! apparent bottom roughness
        integer, parameter :: izNik = 10   ! Nikuradse bottom roughness
        integer, parameter :: izbio = 11   ! biological bottom roughness
        integer, parameter :: izbfm = 12   ! bed form bottom roughness
        integer, parameter :: izbld = 13   ! bed load bottom roughness
        integer, parameter :: izwbl = 14   ! wave bottom roughness
        integer, parameter :: iactv = 15   ! active layer thickness
        integer, parameter :: ishgt = 16   ! saltation height
        integer, parameter :: idflux= 17   ! erosion flux
        integer, parameter :: idnet = 18   ! erosion or deposition
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
        integer, parameter :: idoff = 19   ! tau critical offset
        integer, parameter :: idslp = 20   ! tau critical slope
        integer, parameter :: idtim = 21   ! erodibility time scale
        integer, parameter :: idbmx = 22   ! diffusivity db_max
        integer, parameter :: idbmm = 23   ! diffusivity db_m
        integer, parameter :: idbzs = 24   ! diffusivity db_zs
        integer, parameter :: idbzm = 25   ! diffusivity db_zm
        integer, parameter :: idbzp = 26   ! diffusivity db_zphi
# endif
# if defined MIXED_BED
        integer, parameter :: idprp = 27   ! cohesive behavior
# endif

        integer  :: idBott(MBOTP)          ! bottom properties IDs
        integer  :: idSbed(MBEDP)          ! bed properties IDs

        integer, allocatable :: idBmas(:)  ! class mass indices
        integer, allocatable :: idfrac(:)  ! class fraction indices
        integer, allocatable :: idUbld(:)  ! bed load u-points
        integer, allocatable :: idVbld(:)  ! bed load v-points

        real(r8) :: newlayer_thick(Ngrids) ! deposit thickness criteria
        real(r8) :: minlayer_thick(Ngrids) ! 2nd layer thickness criteria
        real(r8) :: bedload_coeff(Ngrids)  ! bedload rate coefficient

        real(r8), allocatable :: Csed(:,:)       ! initial concentration
        real(r8), allocatable :: Erate(:,:)      ! erosion rate
        real(r8), allocatable :: Sd50(:,:)       ! mediam grain diameter
        real(r8), allocatable :: Srho(:,:)       ! grain density
        real(r8), allocatable :: Wsed(:,:)       ! settling velocity
        real(r8), allocatable :: poros(:,:)      ! porosity
        real(r8), allocatable :: tau_ce(:,:)     ! shear for erosion
        real(r8), allocatable :: tau_cd(:,:)     ! shear for deposition
        real(r8), allocatable :: morph_fac(:,:)  ! morphological factor

# if defined COHESIVE_BED || defined MIXED_BED
        real(r8) :: tcr_min(Ngrids) ! minimum shear for erosion
        real(r8) :: tcr_max(Ngrids) ! maximum shear for erosion
        real(r8) :: tcr_slp(Ngrids) ! Tau_crit profile slope
        real(r8) :: tcr_off(Ngrids) ! Tau_crit profile offset
        real(r8) :: tcr_tim(Ngrids) ! Tau_crit consolidation rate
# endif
# if defined MIXED_BED
        real(r8) :: transC(Ngrids)  ! cohesive transition
        real(r8) :: transN(Ngrids)  ! non-cohesive transition
# endif


#endif
 
      END MODULE mod_sediment
