!
!svn $Id: vegetation_mod.h 429 2019-06-10 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================= John C. Warner =====!   
!================================================= Neil K. Ganju ======!   
!================================================= Alexis Beudin ======!   
!================================================= Tarandeep S. Kalra =!
!=======================================================================
!                                                                      !
!  Vegetation Model Kernel Input Variables:                            !
!  =======================================                             !
!  NVEG          Number of vegetation types                            !
!  NVEGP         Number of vegetation array properties                 !
!  CD_VEG        Drag coefficient from each veg type                   ! 
!  E_VEG         Youngs modulus from each veg type                     !
!  VEG_MASSDEN   Mass density from each veg type                       !
!  VEGHMIXCOEF   Viscosity coefficient from vegetation boundary        ! 
!                                                                      !
!  Marsh Wave Induced Thrust Model Input Variables:                    !
!  ===============================================                     !
!  KFAC_MARSH    Marsh sediment erodibility coefficient                !
!  SCARP_HGHT    Absolute change in scarp hght to convert marsh to     !
!                open water cell (only to be used for high res. model) !
!                                                                      !
!  Marsh Vertical Growth Model Input Variables:                        !
!  ===============================================                     !
!  PAR_FAC1         Marsh parabolic curve growth parameter 1           !
!  PAR_FAC2         Marsh parabolic curve growth parameter 2           !
!TDAYS_MARSH_GROWTH Growing number of days for marsh                   !
! MARSH_BULK_DENS   Bulk density for marsh organic sediment            !
!  NUGP             Fraction of below ground biomass                   !
!  BMAX             Peak biomass                                       !
!  CHIREF           Fraction of recalcitrant Carbon                    !
!  ALPHA_PDENS      Growth parameter 1 for marsh veg. density          !
!  BETA_PDENS       Growth parameter 2 for marsh veg. density          !
!  ALPHA_PHGHT      Growth parameter 1 for marsh veg. height           !
!  BETA_PHGHT       Growth parameter 2 for marsh veg. height           !
!  ALPHA_PDIAM      Growth parameter 1 for marsh veg. diameter         !
!  BETA_PDIAM       Growth parameter 2 for marsh veg. diameter         !
!                                                                      !
!  Plant Property indices:                                             !
!  ======================                                              !
!  pdens         Density                                               !
!  phght         Height                                                !
!  pdiam         Diameter                                              !
!  pthck         Thickness                                             !
!                                                                      !
!  Plant Property indices:                                             !
!  ======================                                              !
!  idvprp        Indices storing plant properties                      ! 
!                                                                      !
!  Plant Property Output IDs:                                          !
!  ==========================                                          !
!  ipdens         Id to output plant density                           !
!  iphght         Id to output plant height                            !
!  ipdiam         Id to output plant diameter                          !
!  ipthck         Id to output plant thickness                         !
!  ipupbm         Id to output above ground biomass                    !
!  ipdwbm         Id to output below ground biomass                    !
!  idWdvg         Id to output wave dissipation from vegetation        !
!  idCdvg         Id to output spectral Cd from waves vegetation       !
!                                                                      !
!  Marsh wave induced erosion Output:                                  !
!  ==========================                                          !
!  idTims        Store masking marsh from marsh cells                  ! 
!  idTtot        Total thrust from all directions due to waves         !
!  idTmfo        Marsh sediment flux from marsh cells                  ! 
!  idTmmr        Amount of marsh retreat from all directions           !
!                                                                      !
!  Marsh vertical growth model:                                        !
!  ==============================                                      !
!  idTmtr        Mean tidal range  (MHHW-MLLW)                         !  
!  idTmhw        Mean high high water (MHWW)                           !
!  idTmbp        Below ground biomass for marsh growth                 !
!  idTmvg        Rate of marsh vertical growth                         ! 
!  idTmvt        Amount of marsh vertical growth                       ! 
!======================================================================!
!                                                                      !
      USE mod_param
      USE mod_sediment
!
      implicit none
!
#if defined VEG_DRAG || defined VEG_BIOMASS  
      integer :: NVEG, NVEGP
      integer :: counter
      integer :: phght, pdens, pdiam, pthck
      integer :: ipdens, iphght, ipdiam, ipthck
#endif
!
#ifdef VEG_STREAMING 
      integer :: idWdvg, idCdvg
#endif 
#if defined VEG_DRAG || defined VEG_BIOMASS  
      integer, allocatable :: idvprp(:)
#endif
!    
#ifdef MARSH_DYNAMICS
      integer ::  idTims
# if defined MARSH_WAVE_THRUST 
      integer ::  idTtot
# endif
# if defined MARSH_SED_EROSION 
      integer, allocatable ::  idTmfo(:)
# endif 
# if defined MARSH_RETREAT
      integer ::  idTmmr
# endif 
# if defined MARSH_TIDAL_RANGE_CALC
      integer ::  NTIMES_MARSH
      integer ::  idTmtr
# endif
# if defined MARSH_VERT_GROWTH
      integer :: idTmhw, idTmlw  
      integer :: idTmbp
      integer :: idTmvg
      integer :: idTmvt
# endif 
#endif 
!
#if defined VEG_DRAG || defined VEG_BIOMASS  
      real(r8), allocatable :: E_VEG(:,:)
      real(r8), allocatable :: CD_VEG(:,:)
      real(r8), allocatable :: VEG_MASSDENS(:,:)
      real(r8), allocatable :: VEGHMIXCOEF(:,:)
#endif
! 
#ifdef MARSH_DYNAMICS
# if defined MARSH_SED_EROSION
      real(r8), allocatable :: KFAC_MARSH(:)
#  if defined MARSH_RETREAT
      real(r8), allocatable :: SCARP_HGHT(:)
#  endif 
# endif 
# if defined MARSH_VERT_GROWTH 
      real(r8), allocatable :: PAR_FAC1(:), PAR_FAC2(:)
      real(r8), allocatable :: TDAYS_MARSH_GROWTH(:)
!      real(r8), allocatable :: MARSH_BULK_DENS(:)
      real(r8), allocatable :: NUGP(:)
      real(r8), allocatable :: BMAX(:)
      real(r8), allocatable :: CHIREF(:)
#  if defined MARSH_BIOMASS_VEG
      real(r8), allocatable :: ALPHA_PDENS(:), BETA_PDENS(:)
      real(r8), allocatable :: ALPHA_PHGHT(:), BETA_PHGHT(:)
      real(r8), allocatable :: ALPHA_PDIAM(:), BETA_PDIAM(:)
#  endif 
# endif 
#endif 
!
      CONTAINS 
! 
      SUBROUTINE initialize_vegetation
!
      USE mod_param
      USE mod_sediment
!
      implicit none 
!
!     Setup property indices 
!
#if defined VEG_DRAG || defined VEG_BIOMASS  
       counter = 1
       pdens   = counter 
       counter = counter+1 
       phght   = counter
       counter = counter+1 
       pdiam   = counter
       counter = counter+1 
       pthck   = counter
#endif 
!#ifdef VEG_BIOMASS 
!       counter = counter+1 
!       pabbm   = counter
!       counter = counter+1 
!       pbgbm   = counter 
!#endif 
#if defined VEG_DRAG || defined VEG_BIOMASS  
       NVEGP = counter
       IF (.not.allocated(idvprp)) THEN
         allocate ( idvprp(NVEGP) )
       END IF
#endif 
#ifdef MARSH_DYNAMICS
# if defined MARSH_SED_EROSION 
       IF (.not.allocated(idTmfo)) THEN
         allocate ( idTmfo(NST) )
       END IF
# endif 
#endif 
      END SUBROUTINE initialize_vegetation
