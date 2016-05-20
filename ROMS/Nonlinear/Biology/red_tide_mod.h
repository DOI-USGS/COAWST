!
!svn $Id: red_tide_mod.h 791 2016-05-05 22:39:42Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Powell et al. (2006) ecosystem model with iron       !
!  limitation:                                                         !
!                                                                      !
!  AttS         Mean light diffuse attenuation coefficient in the      !
!                 sediment (1/m)                                       !
!  AttW         Mean light diffuse attenuation coefficient in the      !
!                 water column (1/m)                                   !
!  DIN_Cdepth   Constant Dissolved Inorganic Nutrient concentration    !
!                 below of growth crical depth (millimole/m3)          !
!  Dg           Mean depth (m) of sediment over which cysts are able   !
!                 to germinate and contribute to the bloom             !
!  E_dark       Light level for germination under dark conditions      !
!                 (watt/m2)                                            !
!  E_light      Light level for germination under "light" conditions   !
!                 (Watt/m2)                                            !
!  Gmax         Maximum growth rate at optimal temperature and salinty !
!                 (1/day)                                              !
!  G_eff        Growth efficiency (m2/Watts/day)                       !
!  G_r          Maintanenance respiration rate (1/day)                 !
!  Kn           Half-saturation constant for nutrient limited growth   !
!                 (millimoles/m3)                                      !
!  Mor          Spatially and temporarily averaged mortality rate      !
!                 (1/day)                                              !
!  Tmin_growth  Coldest temperature limit used to compute temperature- !
!                 dependent growth factor from cubic polynomial fit    !
!                 based on available data (Celsius)                    !
!  srad_Cdepth  Averaged solar shortwave radiation used to compute     !
!                 critical depth in the growth function (Watts/m2)     !
!  wDino        Dinoflagellate (Alexandrium Fundyense) vertical        !
!                 swimming rate (positive; m/day)                      !
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
      integer :: idAsrf                 ! Averaged shortwave radiation
      integer :: idCyst                 ! Bottom cyst concentration
      integer :: idODIN                 ! Dissolved Inorganic Nutrient
      integer :: iDino                  ! Dinoflagellate concentration
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)


      real(r8), allocatable :: AttS(:)         ! 1/m
      real(r8), allocatable :: AttW(:)         ! 1/m
      real(r8), allocatable :: DIN_Cdepth(:)   ! millimoles/m3
      real(r8), allocatable :: Dg(:)           ! m
      real(r8), allocatable :: E_dark(:)       ! Watts/m2
      real(r8), allocatable :: E_light(:)      ! Watts/m2
      real(r8), allocatable :: Gmax(:)         ! 1/day
      real(r8), allocatable :: G_eff(:)        ! m2/Watts/day
      real(r8), allocatable :: G_r(:)          ! 1/day
      real(r8), allocatable :: Kn(:)           ! millimoles/m3
      real(r8), allocatable :: Mor(:)          ! 1/day
      real(r8), allocatable :: Tmin_growth(:)  ! Celsius
      real(r8), allocatable :: srad_Cdepth(:)  ! Watts/m2
      real(r8), allocatable :: wDino(:)        ! m/day
!
!  Mid-day of each month (YearDay = 0 for Jan 1, 00:00:00).
!
      real(r8), dimension(12) :: Month_MidDay =                         &
     &          (/ 15.5_r8,  45.0_r8,  74.5_r8, 105.0_r8,               &
     &            135.5_r8, 166.0_r8, 196.5_r8, 227.5_r8,               &
     &            258.0_r8, 288.5_r8, 319.0_r8, 349.5_r8 /)
!
! Monthly median germination potential.
!
      real(r8), dimension(12) :: GP =                                   &
     &          (/ 21.90_r8, 11.25_r8, 78.0_r8, 85.0_r8,                &
     &             96.8_r8,  93.0_r8,  60.0_r8, 50.0_r8,                &
     &             10.0_r8,  11.5_r8,  17.0_r8, 34.5_r8 /)
!
! Normalized montly mean germination potential.
!
      real(r8), dimension(12) :: GPN

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

      real(r8) :: GPmax
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
      NBT=1
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AttS)) THEN
        allocate ( AttS(Ngrids) )
      END IF
      IF (.not.allocated(AttW)) THEN
        allocate ( AttW(Ngrids) )
      END IF
      IF (.not.allocated(DIN_Cdepth)) THEN
        allocate ( DIN_Cdepth(Ngrids) )
      END IF
      IF (.not.allocated(Dg)) THEN
        allocate ( Dg(Ngrids) )
      END IF
      IF (.not.allocated(E_dark)) THEN
        allocate ( E_dark(Ngrids) )
      END IF
      IF (.not.allocated(E_light)) THEN
        allocate ( E_light(Ngrids) )
      END IF
      IF (.not.allocated(Gmax)) THEN
        allocate ( Gmax(Ngrids) )
      END IF
      IF (.not.allocated(G_eff)) THEN
        allocate ( G_eff(Ngrids) )
      END IF
      IF (.not.allocated(G_r)) THEN
        allocate ( G_r(Ngrids) )
      END IF
      IF (.not.allocated(Kn)) THEN
        allocate ( Kn(Ngrids) )
      END IF
      IF (.not.allocated(Mor)) THEN
        allocate ( Mor(Ngrids) )
      END IF
      IF (.not.allocated(Tmin_growth)) THEN
        allocate ( Tmin_growth(Ngrids) )
      END IF
      IF (.not.allocated(srad_Cdepth)) THEN
        allocate ( srad_Cdepth(Ngrids) )
      END IF
      IF (.not.allocated(wDino)) THEN
        allocate ( wDino(Ngrids) )
      END IF
#ifdef TANGENT
      IF (.not.allocated(tl_wDino)) THEN
        allocate ( tl_wDino(Ngrids) )
      END IF
#endif
#ifdef ADJOINT
      IF (.not.allocated(ad_wDino)) THEN
        allocate ( ad_wDino(Ngrids) )
      END IF
#endif
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
      iDino=ic+1
!
!  Compute normalized montly germination poterntial.
!
      GPmax=MAXVAL(GP)
      DO i=1,12
        GPN(i)=GP(i)/GPmax
      END DO

      RETURN
      END SUBROUTINE initialize_biology
