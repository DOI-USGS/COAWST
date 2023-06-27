#include "cppdefs.h"
      MODULE mod_bbl
#ifdef BBL_MODEL
!
!git $Id$
!svn $Id: mod_bbl.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Ubot         Wind-induced, bed wave orbital U-velocity (m/s) at     !
!                 RHO-points.                                          !
!  Ur           Bottom U-momentum above bed (m/s) at RHO-points.       !
!  Vbot         Wind-induced, bed wave orbital V-velocity (m/s) at     !
!                 RHO-points.                                          !
!  Vr           Bottom V-momentum above bed (m/s) at RHO-points.       !
!  bustrc       Kinematic bottom stress (m2/s2) due currents in the    !
!                 XI-direction at RHO-points.                          !
!  bustrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the XI-direction at horizontal RHO-points.     !
!  bustrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the XI-direction at RHO-points.      !
!  bvstrc       Kinematic bottom stress (m2/s2) due currents in the    !
!                 ETA-direction at RHO-points.                         !
!  bvstrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the ETA-direction at horizontal RHO-points.    !
!  bvstrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the ETA-direction RHO-points.        !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none
!
        PUBLIC :: allocate_bbl
        PUBLIC :: deallocate_bbl
        PUBLIC :: initialize_bbl
!
!-----------------------------------------------------------------------
!  Define T_BBL structure.
!-----------------------------------------------------------------------
!
        TYPE T_BBL

          integer,  pointer :: Iconv(:,:)

          real(r8), pointer :: Ubot(:,:)
          real(r8), pointer :: Ur(:,:)
          real(r8), pointer :: Vbot(:,:)
          real(r8), pointer :: Vr(:,:)
          real(r8), pointer :: bustrc(:,:)
          real(r8), pointer :: bvstrc(:,:)
          real(r8), pointer :: bustrw(:,:)
          real(r8), pointer :: bvstrw(:,:)
          real(r8), pointer :: bustrcwmax(:,:)
          real(r8), pointer :: bvstrcwmax(:,:)

        END TYPE T_BBL
!
        TYPE (T_BBL), allocatable :: BBL(:)
!
      CONTAINS
!
      SUBROUTINE allocate_bbl (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
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
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BBL(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Bottom boundary layer arrays.
!
      allocate ( BBL(ng) % Iconv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % Ubot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % Ur(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % Vbot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % Vr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % bustrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % bvstrc(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % bustrw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % bvstrw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % bustrcwmax(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( BBL(ng) % bvstrcwmax(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      RETURN
      END SUBROUTINE allocate_bbl
!
      SUBROUTINE deallocate_bbl (ng)
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
     &  __FILE__//", deallocate_bbl"

# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_BBL structure
!  separately.
!-----------------------------------------------------------------------
!
      IF (.not.destroy(ng, BBL(ng)%Iconv, MyFile,                       &
     &                 __LINE__, 'BBL(ng)%Iconv')) RETURN

      IF (.not.destroy(ng, BBL(ng)%Ubot, MyFile,                        &
     &                 __LINE__, 'BBL(ng)%Ubot')) RETURN

      IF (.not.destroy(ng, BBL(ng)%Ur, MyFile,                          &
     &                 __LINE__, 'BBL(ng)%Ur')) RETURN

      IF (.not.destroy(ng, BBL(ng)%Vbot, MyFile,                        &
     &                 __LINE__, 'BBL(ng)%Vbot')) RETURN

      IF (.not.destroy(ng, BBL(ng)%Vr, MyFile,                          &
     &                 __LINE__, 'BBL(ng)%Vr')) RETURN

      IF (.not.destroy(ng, BBL(ng)%bustrc, MyFile,                      &
     &                 __LINE__, 'BBL(ng)%bustrc')) RETURN

      IF (.not.destroy(ng, BBL(ng)%bvstrc, MyFile,                      &
     &                 __LINE__, 'BBL(ng)%bvstrc')) RETURN

      IF (.not.destroy(ng, BBL(ng)%bustrw, MyFile,                      &
     &                 __LINE__, 'BBL(ng)%bustrw')) RETURN

      IF (.not.destroy(ng, BBL(ng)%bvstrw, MyFile,                      &
     &                 __LINE__, 'BBL(ng)%bvstrw')) RETURN

      IF (.not.destroy(ng, BBL(ng)%bustrcwmax, MyFile,                  &
     &                 __LINE__, 'BBL(ng)%bustrcwmax')) RETURN

      IF (.not.destroy(ng, BBL(ng)%bvstrcwmax, MyFile,                  &
     &                 __LINE__, 'BBL(ng)%bvstrcwmax')) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type BBL structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(BBL)) deallocate ( BBL )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_bbl
!
      SUBROUTINE initialize_bbl (ng, tile)
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
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.

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
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          BBL(ng) % Iconv(i,j) = 0

          BBL(ng) % Ubot(i,j) = IniVal
          BBL(ng) % Ur(i,j) = IniVal

          BBL(ng) % Vbot(i,j) = IniVal
          BBL(ng) % Vr(i,j) = IniVal

          BBL(ng) % bustrc(i,j) = IniVal
          BBL(ng) % bvstrc(i,j) = IniVal

          BBL(ng) % bustrw(i,j) = IniVal
          BBL(ng) % bvstrw(i,j) = IniVal

          BBL(ng) % bustrcwmax(i,j) = IniVal
          BBL(ng) % bvstrcwmax(i,j) = IniVal
        END DO
      END DO
!
      RETURN
      END SUBROUTINE initialize_bbl
#endif
      END MODULE mod_bbl
