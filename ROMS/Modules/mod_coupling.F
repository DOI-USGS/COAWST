#include "cppdefs.h"
      MODULE mod_coupling
#ifdef SOLVE3D
!
!git $Id$
!svn $Id: mod_coupling.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
!
        implicit none
!
        PUBLIC :: allocate_coupling
        PUBLIC :: deallocate_coupling
        PUBLIC :: initialize_coupling
!
!-----------------------------------------------------------------------
!  Define T_COUPLING structure.
!-----------------------------------------------------------------------
!
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
!
        TYPE (T_COUPLING), allocatable :: COUPLING(:)
!
      CONTAINS
!
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( COUPLING(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Nonlinear model state.
!
      allocate ( COUPLING(ng) % DU_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % DU_avg2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % DV_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % DV_avg2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % Zt_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % rufrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % rvfrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % rhoA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % rhoS(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( COUPLING(ng) % tl_DU_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_DU_avg2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_DV_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_DV_avg2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_Zt_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_rufrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_rvfrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % tl_rhoA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % tl_rhoS(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( COUPLING(ng) % ad_DU_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_DU_avg2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_DV_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_DV_avg2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_Zt_avg1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_rufrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_rvfrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % ad_rhoA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( COUPLING(ng) % ad_rhoS(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      allocate ( COUPLING(ng) % DU_avg1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( COUPLING(ng) % DU_avg2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( COUPLING(ng) % DV_avg1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( COUPLING(ng) % DV_avg2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( COUPLING(ng) % rufrcG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( COUPLING(ng) % rvfrcG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
!
      RETURN
      END SUBROUTINE allocate_coupling
!
      SUBROUTINE deallocate_coupling (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : Ngrids
# ifdef SUBOBJECT_DEALLOCATION
      USE destroy_mod, ONLY : destroy
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_coupling"

# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF (.not.destroy(ng, COUPLING(ng)%DU_avg1, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%DU_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%DU_avg2, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%DU_avg2')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%DV_avg1, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%DV_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%DV_avg2, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%DV_avg2')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%Zt_avg1, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%Zt_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%rufrc, MyFile,                  &
     &                 __LINE__, 'COUPLING(ng)%rufrc')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%rvfrc, MyFile,                  &
     &                 __LINE__, 'COUPLING(ng)%rvfrc')) RETURN

#  ifdef VAR_RHO_2D
      IF (.not.destroy(ng, COUPLING(ng)%rhoA, MyFile,                   &
     &                 __LINE__, 'COUPLING(ng)%rhoA')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%rhoS, MyFile,                   &
     &                 __LINE__, 'COUPLING(ng)%rhoS')) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF (.not.destroy(ng, COUPLING(ng)%tl_DU_avg1, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%tl_DU_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_DU_avg2, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%tl_DU_avg2')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_DV_avg1, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%tl_DV_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_DV_avg2, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%tl_DV_avg2')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_Zt_avg1, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%tl_Zt_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_rufrc, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%tl_rufrc')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_rvfrc, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%tl_rvfrc')) RETURN

#   ifdef VAR_RHO_2D
      IF (.not.destroy(ng, COUPLING(ng)%tl_rhoA, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%tl_rhoA')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%tl_rhoS, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%tl_rhoS')) RETURN
#   endif
#  endif

#  ifdef ADJOINT
!
!  Adjoint model state.
!
      IF (.not.destroy(ng, COUPLING(ng)%ad_DU_avg1, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%ad_DU_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_DU_avg2, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%ad_DU_avg2')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_DV_avg1, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%ad_DV_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_DV_avg2, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%ad_DV_avg2')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_Zt_avg1, MyFile,             &
     &                 __LINE__, 'COUPLING(ng)%ad_Zt_avg1')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_rufrc, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%ad_rufrc')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_rvfrc, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%ad_rvfrc')) RETURN

#   ifdef VAR_RHO_2D
      IF (.not.destroy(ng, COUPLING(ng)%ad_rhoA, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%ad_rhoA')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%ad_rhoS, MyFile,                &
     &                 __LINE__, 'COUPLING(ng)%ad_rhoS')) RETURN
#   endif
#  endif

#  if defined FORWARD_READ && \
     (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      IF (.not.destroy(ng, COUPLING(ng)%DU_avg1G, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%DU_avg1G')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%DU_avg2G, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%DU_avg2G')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%DV_avg1G, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%DV_avg1G')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%DV_avg2G, MyFile,               &
     &                 __LINE__, 'COUPLING(ng)%DV_avg2G')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%rufrcG, MyFile,                 &
     &                 __LINE__, 'COUPLING(ng)%rufrcG')) RETURN

      IF (.not.destroy(ng, COUPLING(ng)%rvfrcG, MyFile,                 &
     &                 __LINE__, 'COUPLING(ng)%rvfrcG')) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type COUPLING structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(COUPLING)) deallocate ( COUPLING )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_coupling
!
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
!
      RETURN
      END SUBROUTINE initialize_coupling
#endif
      END MODULE mod_coupling
