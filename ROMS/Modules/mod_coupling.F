#include "cppdefs.h"
      MODULE mod_coupling
#ifdef SOLVE3D
!
!svn $Id: mod_coupling.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  DU_avg1   Time averaged U-flux for 2D equations (m3/s).             !
!  DU_avg2   Time averaged U-flux for 3D equations coupling (m3/s).    !
!  DV_avg1   Time averaged V-flux for 2D equations (m3/s).             !
!  DV_avg2   Time averaged V-flux for 3D equations coupling (m3/s).    !
!  Zt_avg1   Free-surface averaged over all short time-steps (m).      !
!  rhoA      Normalized vertical averaged density.                     !
!  rhoS      Normalized vertical averaged density perturbation.        !
!  rufrc     Right-hand-side forcing term for 2D U-momentum (m4/s2)    !
!  rvfrc     Right-hand-side forcing term for 2D V-momentum (m4/s2)    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_COUPLING
!
!  Nonlinear model state.
!
          real(r8), pointer :: DU_avg1(:,:)
          real(r8), pointer :: DU_avg2(:,:)
          real(r8), pointer :: DV_avg1(:,:)
          real(r8), pointer :: DV_avg2(:,:)
          real(r8), pointer :: Zt_avg1(:,:)
          real(r8), pointer :: rufrc(:,:)
          real(r8), pointer :: rvfrc(:,:)
# ifdef VAR_RHO_2D
          real(r8), pointer :: rhoA(:,:)
          real(r8), pointer :: rhoS(:,:)
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_DU_avg1(:,:)
          real(r8), pointer :: tl_DU_avg2(:,:)
          real(r8), pointer :: tl_DV_avg1(:,:)
          real(r8), pointer :: tl_DV_avg2(:,:)
          real(r8), pointer :: tl_Zt_avg1(:,:)
          real(r8), pointer :: tl_rufrc(:,:)
          real(r8), pointer :: tl_rvfrc(:,:)
#  ifdef VAR_RHO_2D
          real(r8), pointer :: tl_rhoA(:,:)
          real(r8), pointer :: tl_rhoS(:,:)
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_DU_avg1(:,:)
          real(r8), pointer :: ad_DU_avg2(:,:)
          real(r8), pointer :: ad_DV_avg1(:,:)
          real(r8), pointer :: ad_DV_avg2(:,:)
          real(r8), pointer :: ad_Zt_avg1(:,:)
          real(r8), pointer :: ad_rufrc(:,:)
          real(r8), pointer :: ad_rvfrc(:,:)
#  ifdef VAR_RHO_2D
          real(r8), pointer :: ad_rhoA(:,:)
          real(r8), pointer :: ad_rhoS(:,:)
#  endif
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
          real(r8), pointer :: DU_avg1G(:,:,:)
          real(r8), pointer :: DU_avg2G(:,:,:)
          real(r8), pointer :: DV_avg1G(:,:,:)
          real(r8), pointer :: DV_avg2G(:,:,:)
          real(r8), pointer :: rufrcG(:,:,:)
          real(r8), pointer :: rvfrcG(:,:,:)
# endif

        END TYPE T_COUPLING

        TYPE (T_COUPLING), allocatable :: COUPLING(:)

        CONTAINS

      SUBROUTINE allocate_coupling (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!   grids.                                                             !
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
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( COUPLING(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( COUPLING(ng) % DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rvfrc(LBi:UBi,LBj:UBj) )

# ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rhoS(LBi:UBi,LBj:UBj) )
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( COUPLING(ng) % tl_DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_rvfrc(LBi:UBi,LBj:UBj) )

#  ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % tl_rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_rhoS(LBi:UBi,LBj:UBj) )
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( COUPLING(ng) % ad_DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_rvfrc(LBi:UBi,LBj:UBj) )

#  ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % ad_rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_rhoS(LBi:UBi,LBj:UBj) )
#  endif
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      allocate ( COUPLING(ng) % DU_avg1G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % DU_avg2G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % DV_avg1G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % DV_avg2G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % rufrcG(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % rvfrcG(LBi:UBi,LBj:UBj,2) )
# endif

      RETURN
      END SUBROUTINE allocate_coupling

      SUBROUTINE initialize_coupling (ng, tile, model)
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
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# else
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
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            COUPLING(ng) % DU_avg1(i,j) = IniVal
            COUPLING(ng) % DU_avg2(i,j) = IniVal

            COUPLING(ng) % DV_avg1(i,j) = IniVal
            COUPLING(ng) % DV_avg2(i,j) = IniVal

            COUPLING(ng) % Zt_avg1(i,j) = IniVal

            COUPLING(ng) % rufrc(i,j) = IniVal
            COUPLING(ng) % rvfrc(i,j) = IniVal

# ifdef VAR_RHO_2D
            COUPLING(ng) % rhoA(i,j) = IniVal
            COUPLING(ng) % rhoS(i,j) = IniVal
# endif
          END DO
        END DO
      END IF

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            COUPLING(ng) % tl_DU_avg1(i,j) = IniVal
            COUPLING(ng) % tl_DU_avg2(i,j) = IniVal

            COUPLING(ng) % tl_DV_avg1(i,j) = IniVal
            COUPLING(ng) % tl_DV_avg2(i,j) = IniVal

            COUPLING(ng) % tl_Zt_avg1(i,j) = IniVal

            COUPLING(ng) % tl_rufrc(i,j) = IniVal
            COUPLING(ng) % tl_rvfrc(i,j) = IniVal

#  ifdef VAR_RHO_2D
            COUPLING(ng) % tl_rhoA(i,j) = IniVal
            COUPLING(ng) % tl_rhoS(i,j) = IniVal
#  endif
          END DO
        END DO
      END IF
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            COUPLING(ng) % ad_DU_avg1(i,j) = IniVal
            COUPLING(ng) % ad_DU_avg2(i,j) = IniVal

            COUPLING(ng) % ad_DV_avg1(i,j) = IniVal
            COUPLING(ng) % ad_DV_avg2(i,j) = IniVal

            COUPLING(ng) % ad_Zt_avg1(i,j) = IniVal

            COUPLING(ng) % ad_rufrc(i,j) = IniVal
            COUPLING(ng) % ad_rvfrc(i,j) = IniVal

#  ifdef VAR_RHO_2D
            COUPLING(ng) % ad_rhoA(i,j) = IniVal
            COUPLING(ng) % ad_rhoS(i,j) = IniVal
#  endif
          END DO
        END DO
      END IF
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      IF (model.eq.0) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            COUPLING(ng) % DU_avg1G(i,j,1) = IniVal
            COUPLING(ng) % DU_avg1G(i,j,2) = IniVal
            COUPLING(ng) % DU_avg2G(i,j,1) = IniVal
            COUPLING(ng) % DU_avg2G(i,j,2) = IniVal

            COUPLING(ng) % DV_avg1G(i,j,1) = IniVal
            COUPLING(ng) % DV_avg1G(i,j,2) = IniVal
            COUPLING(ng) % DV_avg2G(i,j,1) = IniVal
            COUPLING(ng) % DV_avg2G(i,j,2) = IniVal

            COUPLING(ng) % rufrcG(i,j,1) = IniVal
            COUPLING(ng) % rufrcG(i,j,2) = IniVal
            COUPLING(ng) % rvfrcG(i,j,1) = IniVal
            COUPLING(ng) % rvfrcG(i,j,2) = IniVal
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE initialize_coupling
#endif
      END MODULE mod_coupling
