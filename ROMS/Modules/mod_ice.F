#include "cppdefs.h"
      MODULE mod_ice
#if defined ICE_MODEL || defined NCEP_FLUXES
        USE mod_kinds

        implicit none

        TYPE T_ICE

          real(r8), pointer :: ui(:,:,:)
          real(r8), pointer :: vi(:,:,:)
          real(r8), pointer :: uie(:,:,:)
          real(r8), pointer :: vie(:,:,:)
          real(r8), pointer :: ai(:,:,:)
          real(r8), pointer :: hi(:,:,:)
          real(r8), pointer :: hsn(:,:,:)
# ifdef MELT_PONDS
          real(r8), pointer :: apond(:,:,:)
          real(r8), pointer :: hpond(:,:,:)
# endif
          real(r8), pointer :: ageice(:,:,:)
          real(r8), pointer :: io_mflux(:,:)
          real(r8), pointer :: wfr(:,:)
          real(r8), pointer :: wai(:,:)
          real(r8), pointer :: wao(:,:)
          real(r8), pointer :: wio(:,:)
          real(r8), pointer :: wro(:,:)
          real(r8), pointer :: wdiv(:,:)
          real(r8), pointer :: tis(:,:)
          real(r8), pointer :: ti(:,:,:)
          real(r8), pointer :: enthalpi(:,:,:)
          real(r8), pointer :: hage(:,:,:)
          real(r8), pointer :: utau_iw(:,:)
          real(r8), pointer :: chu_iw(:,:)
          real(r8), pointer :: spd_iw(:,:)
          real(r8), pointer :: coef_ice_heat(:,:)
          real(r8), pointer :: rhs_ice_heat(:,:)
          real(r8), pointer :: pice(:,:)
#ifdef ICE_STRENGTH_QUAD
          real(r8), pointer :: pstar_grid(:,:)
#endif
          real(r8), pointer :: zetai(:,:)
          real(r8), pointer :: eta(:,:)
          real(r8), pointer :: tauaiu(:,:)
          real(r8), pointer :: tauaiv(:,:)
#ifdef ICE_LANDFAST
          real(r8), pointer :: h_loadu(:,:)
          real(r8), pointer :: h_loadv(:,:)
#endif
          real(r8), pointer :: uwater(:,:)
          real(r8), pointer :: vwater(:,:)
          real(r8), pointer :: sealev(:,:)
          real(r8), pointer :: s0mk(:,:)
          real(r8), pointer :: t0mk(:,:)
          real(r8), pointer :: sig11(:,:,:)
          real(r8), pointer :: sig22(:,:,:)
          real(r8), pointer :: sig12(:,:,:)

# if defined ICE_BIO  
          real(r8), pointer ::IcePhL(:,:,:)
          real(r8), pointer ::IceNO3(:,:,:)
          real(r8), pointer ::IceNH4(:,:,:)
          integer, pointer ::IceLog(:,:,:)
# endif

        END TYPE T_ICE

        TYPE (T_ICE), allocatable :: ICE(:)

      CONTAINS

      SUBROUTINE allocate_ice (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002-2014 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( ICE(Ngrids) )
!
      allocate ( ICE(ng) % ui(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % vi(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % uie(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % vie(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % ai(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hi(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hsn(LBi:UBi,LBj:UBj,2) )
# ifdef MELT_PONDS
      allocate ( ICE(ng) % apond(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hpond(LBi:UBi,LBj:UBj,2) )
# endif
      allocate ( ICE(ng) % ageice(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % io_mflux(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wfr(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wai(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wao(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wio(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wro(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wdiv(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tis(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % ti(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % enthalpi(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hage(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % utau_iw(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % chu_iw(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % spd_iw(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % coef_ice_heat(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % rhs_ice_heat(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % pice(LBi:UBi,LBj:UBj) )
#ifdef ICE_STRENGTH_QUAD
      allocate ( ICE(ng) % pstar_grid(LBi:UBi,LBj:UBj) )
#endif
      allocate ( ICE(ng) % zetai(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % eta(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tauaiu(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tauaiv(LBi:UBi,LBj:UBj) )
#ifdef ICE_LANDFAST
      allocate ( ICE(ng) % h_loadu(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % h_loadv(LBi:UBi,LBj:UBj) )
#endif
      allocate ( ICE(ng) % uwater(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % vwater(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % sealev(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % s0mk(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % t0mk(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % sig11(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % sig22(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % sig12(LBi:UBi,LBj:UBj,2) )

# if defined ICE_BIO  
      allocate (ICE(ng) % IcePhL(LBi:UBi,LBj:UBj,2) )
      allocate (ICE(ng) % IceNO3(LBi:UBi,LBj:UBj,2) )
      allocate (ICE(ng) % IceNH4(LBi:UBi,LBj:UBj,2) )
      allocate (ICE(ng) % IceLog(LBi:UBi,LBj:UBj,2) )
# endif

      RETURN
      END SUBROUTINE allocate_ice

      SUBROUTINE initialize_ice (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ICE(ng) % ui(i,j,1) = IniVal
          ICE(ng) % ui(i,j,2) = IniVal
          ICE(ng) % vi(i,j,1) = IniVal
          ICE(ng) % vi(i,j,2) = IniVal
          ICE(ng) % uie(i,j,1) = IniVal
          ICE(ng) % uie(i,j,2) = IniVal
          ICE(ng) % vie(i,j,1) = IniVal
          ICE(ng) % vie(i,j,2) = IniVal
          ICE(ng) % ai(i,j,1) = IniVal
          ICE(ng) % ai(i,j,2) = IniVal
# ifdef SOGLOBEC
          ICE(ng) % hi(i,j,1) = 0.5_r8
          ICE(ng) % hi(i,j,2) = 0.5_r8
# else
          ICE(ng) % hi(i,j,1) = IniVal
          ICE(ng) % hi(i,j,2) = IniVal
# endif
          ICE(ng) % hsn(i,j,1) = IniVal
          ICE(ng) % hsn(i,j,2) = IniVal
# ifdef MELT_PONDS
          ICE(ng) % apond(i,j,1) = IniVal
          ICE(ng) % apond(i,j,2) = IniVal
          ICE(ng) % hpond(i,j,1) = IniVal
          ICE(ng) % hpond(i,j,2) = IniVal
# endif
          ICE(ng) % ageice(i,j,1) = IniVal
          ICE(ng) % ageice(i,j,2) = IniVal
          ICE(ng) % io_mflux(i,j) = IniVal
          ICE(ng) % wfr(i,j) = IniVal
          ICE(ng) % wai(i,j) = IniVal
          ICE(ng) % wao(i,j) = IniVal
          ICE(ng) % wio(i,j) = IniVal
          ICE(ng) % wro(i,j) = IniVal
          ICE(ng) % wdiv(i,j) = IniVal
          ICE(ng) % tis(i,j) = IniVal
          ICE(ng) % ti(i,j,1) = IniVal
          ICE(ng) % ti(i,j,2) = IniVal
          ICE(ng) % enthalpi(i,j,1) = IniVal
          ICE(ng) % enthalpi(i,j,2) = IniVal
          ICE(ng) % hage(i,j,1) = IniVal
          ICE(ng) % hage(i,j,2) = IniVal
          ICE(ng) % utau_iw(i,j) = IniVal
          ICE(ng) % chu_iw(i,j) = IniVal
          ICE(ng) % spd_iw(i,j) = IniVal
          ICE(ng) % coef_ice_heat(i,j) = IniVal
          ICE(ng) % rhs_ice_heat(i,j) = IniVal
          ICE(ng) % pice(i,j) = IniVal
#ifdef ICE_STRENGTH_QUAD
          ICE(ng) % pstar_grid(i,j) = IniVal
#endif
          ICE(ng) % zetai(i,j) = IniVal
          ICE(ng) % eta(i,j) = IniVal
          ICE(ng) % tauaiu(i,j) = IniVal
          ICE(ng) % tauaiv(i,j) = IniVal
#ifdef ICE_LANDFAST
          ICE(ng) % h_loadu(i,j) = IniVal
          ICE(ng) % h_loadv(i,j) = IniVal
#endif
          ICE(ng) % uwater(i,j) = IniVal
          ICE(ng) % vwater(i,j) = IniVal
          ICE(ng) % sealev(i,j) = IniVal
          ICE(ng) % s0mk(i,j) = IniVal
          ICE(ng) % t0mk(i,j) = IniVal
          ICE(ng) % sig11(i,j,1) = IniVal
          ICE(ng) % sig11(i,j,2) = IniVal
          ICE(ng) % sig22(i,j,1) = IniVal
          ICE(ng) % sig22(i,j,2) = IniVal
          ICE(ng) % sig12(i,j,1) = IniVal
          ICE(ng) % sig12(i,j,2) = IniVal

# if defined ICE_BIO  
          ICE(ng) % IcePhL(i,j,1) = IniVal
          ICE(ng) % IcePhL(i,j,2) = IniVal
          ICE(ng) % IceNO3(i,j,1) = IniVal
          ICE(ng) % IceNO3(i,j,2) = IniVal
          ICE(ng) % IceNH4(i,j,1) = IniVal
          ICE(ng) % IceNH4(i,j,2) = IniVal
          ICE(ng) % IceLog(i,j,1) = -1
          ICE(ng) % IceLog(i,j,2) = -1
# endif
        END DO
      END DO

      RETURN
      END SUBROUTINE initialize_ice
#elif defined CICE_MODEL
        USE mod_kinds

        implicit none

        TYPE T_ICE
          real(r8), pointer :: wfr(:,:)
          real(r8), pointer :: ai(:,:)
          real(r8), pointer :: hi(:,:)
          real(r8), pointer :: tis(:,:)
        END TYPE T_ICE

        TYPE (T_ICE), allocatable :: ICE(:)

      CONTAINS

      SUBROUTINE allocate_ice (ng, LBi, UBi, LBj, UBj)
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( ICE(Ngrids) )
!
      allocate ( ICE(ng) % wfr(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % ai(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % hi(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tis(LBi:UBi,LBj:UBj) )

      RETURN
      END SUBROUTINE allocate_ice

      SUBROUTINE initialize_ice (ng, tile)
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ICE(ng) % wfr(i,j) = IniVal
          ICE(ng) % ai(i,j) = IniVal
          ICE(ng) % hi(i,j) = IniVal
          ICE(ng) % tis(i,j) = IniVal
        END DO
      END DO
!
      RETURN
      END SUBROUTINE initialize_ice
#endif
      END MODULE mod_ice
