#include "cppdefs.h"
      MODULE mod_extract

#ifdef GRID_EXTRACT
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Grid extration module for I/O purposes.                             !
!                                                                      !
!  Hz         Thicknesses (m) of vertical RHO-points.                  !
!  CosAngler  Cosine of curvilinear angle, angler.                     !
!  SinAngler  Sine of curvilinear angle, angler.                       !
!  angler     Angle (radians) between XI-axis and true EAST at         !
!               RHO-points.                                            !
!  dmde       ETA-derivative of inverse metric factor pm,              !
!               d(1/pm)/d(ETA).                                        !
!  dndx       XI-derivative  of inverse metric factor pn,              !
!               d(1/pn)/d(XI).                                         !
!  f          Coriolis parameter (1/s).                                !
!  h          Bottom depth (m) at RHO-points.                          !
!  latp       Latitude (degrees_north) at PSI-points.                  !
!  latr       Latitude (degrees_north) at RHO-points.                  !
!  latu       Latitude (degrees_north) at U-points.                    !
!  latv       Latitude (degrees_north) at V-points.                    !
!  lonp       Longitude (degrees_east) at PSI-points.                  !
!  lonr       Longitude (degrees_east) at RHO-points.                  !
!  lonu       Longitude (degrees_east) at U-points.                    !
!  lonv       Longitude (degrees_east) at V-points.                    !
!  pm         Coordinate transformation metric "m" (1/meters)          !
!               associated with the differential distances in XI.      !
!  pn         Coordinate transformation metric "n" (1/meters)          !
!               associated with the differential distances in ETA.     !
# ifdef MASKING
!  pmask      Slipperiness time-independent mask at PSI-points:        !
!               (0=Land, 1=Sea, 2=no-slip).                            !
!  rmask      Time-independent mask at RHO-points (0=Land, 1=Sea).     !
!  umask      Time-independent mask at U-points (0=Land, 1=Sea).       !
!  vmask      Time-independent mask at V-points (0=Land, 1=Sea).       !
# endif
!  xp         XI-coordinates (m) at PSI-points.                        !
!  xr         XI-coordinates (m) at RHO-points.                        !
!  xu         XI-coordinates (m) at U-points.                          !
!  xv         XI-coordinates (m) at V-points.                          !
!  yp         ETA-coordinates (m) at PSI-points.                       !
!  yr         ETA-coordinates (m) at RHO-points.                       !
!  yu         ETA-coordinates (m) at U-points.                         !
!  yv         ETA-coordinates (m) at V-points.                         !
!  z_r        Actual depths (m) at horizontal RHO-points and           !
!               vertical RHO-points.                                   !
!                                                                      !
! Packed global arrays needed for field extraction by interpolation:   !
!
!  Gx_psi     Global X-coordinate at PSI-points.                       !
!  Gx_rho     Global X-coordinate at RHO-points.                       !
!  Gx_u       Global X-coordinate at U-points.                         !
!  Gx_v       Global X-coordinate at V-points.                         !
!  Gy_psi     Global Y-coordinate at PSI-points.                       !
!  Gy_rho     Global Y-coordinate at RHO-points.                       !
!  Gy_u       Global Y-coordinate at U-points.                         !
!  Gy_v       Global Y-coordinate at V-points.                         !

!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none
!
        PUBLIC :: allocate_extract
        PUBLIC :: deallocate_extract
        PUBLIC :: initialize_extract
!
!-----------------------------------------------------------------------
!  Define T_EXTRACT structure.
!-----------------------------------------------------------------------
!
        TYPE T_EXTRACT
!
!  Parameters.
!
          real(r8), allocatable :: Hmin(:)
          real(r8), allocatable :: Hmax(:)
          real(r8), allocatable :: LonMin(:)
          real(r8), allocatable :: LonMax(:)
          real(r8), allocatable :: LatMin(:)
          real(r8), allocatable :: LatMax(:)
!
!  Nonlinear model state.
!
          real(r8), pointer :: angler(:,:)
          real(r8), pointer :: CosAngler(:,:)
          real(r8), pointer :: SinAngler(:,:)

# if defined CURVGRID && defined UV_ADV
          real(r8), pointer :: dmde(:,:)
          real(r8), pointer :: dndx(:,:)
# endif
          real(r8), pointer :: f(:,:)
          real(r8), pointer :: h(:,:)
          real(r8), pointer :: latp(:,:)
          real(r8), pointer :: latr(:,:)
          real(r8), pointer :: latu(:,:)
          real(r8), pointer :: latv(:,:)
          real(r8), pointer :: lonp(:,:)
          real(r8), pointer :: lonr(:,:)
          real(r8), pointer :: lonu(:,:)
          real(r8), pointer :: lonv(:,:)
          real(r8), pointer :: MyLon(:,:)
          real(r8), pointer :: omn(:,:)
          real(r8), pointer :: om_p(:,:)
          real(r8), pointer :: om_r(:,:)
          real(r8), pointer :: om_u(:,:)
          real(r8), pointer :: om_v(:,:)
          real(r8), pointer :: on_p(:,:)
          real(r8), pointer :: on_r(:,:)
          real(r8), pointer :: on_u(:,:)
          real(r8), pointer :: on_v(:,:)
          real(r8), pointer :: pm(:,:)
          real(r8), pointer :: pn(:,:)
          real(r8), pointer :: xp(:,:)
          real(r8), pointer :: xr(:,:)
          real(r8), pointer :: xu(:,:)
          real(r8), pointer :: xv(:,:)
          real(r8), pointer :: yp(:,:)
          real(r8), pointer :: yr(:,:)
          real(r8), pointer :: yu(:,:)
          real(r8), pointer :: yv(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: Hz(:,:,:)
          real(r8), pointer :: z_r(:,:,:)
          real(r8), pointer :: z_v(:,:,:)
          real(r8), pointer :: z_w(:,:,:)
# endif
# ifdef MASKING
          real(r8), pointer :: pmask(:,:)
          real(r8), pointer :: rmask(:,:)
          real(r8), pointer :: umask(:,:)
          real(r8), pointer :: vmask(:,:)
# endif
!
! Set global geometry 1D arrays, packed in column-major order,
! fractional I- and J-coordinates of the extracted data with
! respect the donor grid ti facilitate bilinear or bicubic
! interpolations. The 3D fields are interpolated level by level,
! so only horizontal fractional coordinates are needed for each
! staggered location.
!
          real(r8), pointer :: Gx_psi(:)
          real(r8), pointer :: Gx_rho(:)
          real(r8), pointer :: Gx_u(:)
          real(r8), pointer :: Gx_v(:)

          real(r8), pointer :: Gy_psi(:)
          real(r8), pointer :: Gy_rho(:)
          real(r8), pointer :: Gy_u(:)
          real(r8), pointer :: Gy_v(:)

# ifdef MASKING
          real(r8), pointer :: Gmask_psi(:)
          real(r8), pointer :: Gmask_rho(:)
          real(r8), pointer :: Gmask_u(:)
          real(r8), pointer :: Gmask_v(:)
# endif

          real(r8), pointer :: Iout_psi(:)
          real(r8), pointer :: Iout_rho(:)
          real(r8), pointer :: Iout_u(:)
          real(r8), pointer :: Iout_v(:)

          real(r8), pointer :: Jout_psi(:)
          real(r8), pointer :: Jout_rho(:)
          real(r8), pointer :: Jout_u(:)
          real(r8), pointer :: Jout_v(:)

        END TYPE T_EXTRACT
!
        TYPE (T_EXTRACT), allocatable :: EXTRACT(:)
!
      CONTAINS
!
      SUBROUTINE allocate_extract (ng, Extract_Flag, LBi, UBi, LBj, UBj)
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
!  Local variable declarations.
!
      integer, intent(in) :: ng, Extract_Flag
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer  :: my_size
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( EXTRACT(Ngrids) )
!
!  Parameters.
!
      IF (.not.allocated(EXTRACT(ng) % Hmin))                           &
     &         allocate ( EXTRACT(ng) % Hmin(Ngrids) )

      IF (.not.allocated(EXTRACT(ng) % Hmax))                           &
     &         allocate ( EXTRACT(ng) % Hmax(Ngrids) )

      IF (.not.allocated(EXTRACT(ng) % LonMin))                         &
     &         allocate ( EXTRACT(ng) % LonMin(Ngrids) )

      IF (.not.allocated(EXTRACT(ng) % LonMax))                         &
     &         allocate ( EXTRACT(ng) % LonMax(Ngrids) )

      IF (.not.allocated(EXTRACT(ng) % LatMin))                         &
     &         allocate ( EXTRACT(ng) % LatMin(Ngrids) )

      IF (.not.allocated(EXTRACT(ng) % LatMax))                         &
     &         allocate ( EXTRACT(ng) % LatMax(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Nonlinear model state.
!
      allocate ( EXTRACT(ng) % angler(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % CosAngler(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % SinAngler(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# if defined CURVGRID && defined UV_ADV
      allocate ( EXTRACT(ng) % dmde(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % dndx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

      allocate ( EXTRACT(ng) % f(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % h(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % latp(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % latr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % latu(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % latv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % lonp(LBi:UBi,LBj:UBj))
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % lonr(LBi:UBi,LBj:UBj))
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % lonu(LBi:UBi,LBj:UBj))
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % lonv(LBi:UBi,LBj:UBj))
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % Mylon(LBi:UBi,LBj:UBj))
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % pm(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % pn(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % xp(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % xr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % xu(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % xv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % yp(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % yr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % yu(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % yv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# ifdef SOLVE3D
      allocate ( EXTRACT(ng) % Hz(LBi:UBi,LBj:UBj,N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

      allocate ( EXTRACT(ng) % z_r(LBi:UBi,LBj:UBj,N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

      allocate ( EXTRACT(ng) % z_v(LBi:UBi,LBj:UBj,N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

      allocate ( EXTRACT(ng) % z_w(LBi:UBi,LBj:UBj,0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
# endif

# ifdef MASKING
      allocate ( EXTRACT(ng) % pmask(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % rmask(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % umask(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( EXTRACT(ng) % vmask(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif
!
!  Allocate global arrays for interpolation
!
      IF (Extract_Flag.ge.1) THEN
        my_size=(xtr_IOBOUNDS(ng)%IUB_psi-xtr_IOBOUNDS(ng)%ILB_psi+1)*  &
     &          (xtr_IOBOUNDS(ng)%JUB_psi-xtr_IOBOUNDS(ng)%JLB_psi+1)

        allocate ( EXTRACT(ng) % Gx_psi(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Gy_psi(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

# ifdef MASKING
        allocate ( EXTRACT(ng) % Gmask_psi(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
# endif

        allocate ( EXTRACT(ng) % Iout_psi(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Jout_psi(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
!
        my_size=(xtr_IOBOUNDS(ng)%IUB_rho-xtr_IOBOUNDS(ng)%ILB_rho+1)*  &
     &          (xtr_IOBOUNDS(ng)%JUB_rho-xtr_IOBOUNDS(ng)%JLB_rho+1)

        allocate ( EXTRACT(ng) % Gx_rho(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Gy_rho(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

# ifdef MASKING
        allocate ( EXTRACT(ng) % Gmask_rho(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
# endif

        allocate ( EXTRACT(ng) % Iout_rho(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Jout_rho(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
!
        my_size=(xtr_IOBOUNDS(ng)%IUB_u-xtr_IOBOUNDS(ng)%ILB_u+1)*      &
     &          (xtr_IOBOUNDS(ng)%JUB_u-xtr_IOBOUNDS(ng)%JLB_u+1)

        allocate ( EXTRACT(ng) % Gx_u(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Gy_u(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

# ifdef MASKING
        allocate ( EXTRACT(ng) % Gmask_u(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
# endif

        allocate ( EXTRACT(ng) % Iout_u(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Jout_u(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
!
        my_size=(xtr_IOBOUNDS(ng)%IUB_v-xtr_IOBOUNDS(ng)%ILB_v+1)*      &
     &          (xtr_IOBOUNDS(ng)%JUB_v-xtr_IOBOUNDS(ng)%JLB_v+1)

        allocate ( EXTRACT(ng) % Gx_v(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Gy_v(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

# ifdef MASKING
        allocate ( EXTRACT(ng) % Gmask_v(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
# endif

        allocate ( EXTRACT(ng) % Iout_v(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)

        allocate ( EXTRACT(ng) % Jout_v(my_size) )
        Dmem(ng)=Dmem(ng)+REAL(my_size,r8)
      END IF
!
      RETURN
      END SUBROUTINE allocate_extract
!
      SUBROUTINE deallocate_extract (ng)
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
     &  __FILE__//", deallocate_extract"

# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_GRID structure
!  separately.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF (.not.destroy(ng, EXTRACT(ng)%angler, MyFile,                  &
     &                 __LINE__, 'EXTRACT(ng)%angler')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%CosAngler, MyFile,               &
     &                 __LINE__, 'EXTRACT(ng)%CosAngler')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%SinAngler, MyFile,               &
     &                 __LINE__, 'EXTRACT(ng)%SinAngler')) RETURN

#  if defined CURVGRID && defined UV_ADV
      IF (.not.destroy(ng, EXTRACT(ng)%dmde, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%dmde')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%dndx, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%dndx')) RETURN
#  endif

      IF (.not.destroy(ng, EXTRACT(ng)%f, MyFile,                       &
     &                 __LINE__, 'EXTRACT(ng)%f')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%h, MyFile,                       &
     &                 __LINE__, 'EXTRACT(ng)%h')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%latp, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%latp')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%latr, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%latr')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%latu, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%latu')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%latv, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%latv')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%lonp, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%lonp')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%lonr, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%lonr')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%lonu, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%lonu')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%lonv, MyFile,                    &
     &                 __LINE__, 'EXTRACT(ng)%lonv')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%Mylon, MyFile,                   &
     &                 __LINE__, 'EXTRACT(ng)%Mylon')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%pm, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%pm')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%pn, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%pn')) RETURN


      IF (.not.destroy(ng, EXTRACT(ng)%xp, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%xp')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%xr, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%xr')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%xu, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%xu')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%xv, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%xv')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%yp, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%yp')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%yr, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%yr')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%yu, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%yu')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%yv, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%yv')) RETURN

#  ifdef SOLVE3D
      IF (.not.destroy(ng, EXTRACT(ng)%Hz, MyFile,                      &
     &                 __LINE__, 'EXTRACT(ng)%Hz')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%z_r, MyFile,                     &
     &                 __LINE__, 'EXTRACT(ng)%z_r')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%z_v, MyFile,                     &
     &                 __LINE__, 'GRID(ng)%z_v')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%z_w, MyFile,                     &
     &                 __LINE__, 'GRID(ng)%z_w')) RETURN
#  endif

#  ifdef MASKING
      IF (.not.destroy(ng, EXTRACT(ng)%pmask, MyFile,                   &
     &                 __LINE__, 'EXTRACT(ng)%pmask')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%rmask, MyFile,                   &
     &                 __LINE__, 'EXTRACT(ng)%rmask')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%umask, MyFile,                   &
     &                 __LINE__, 'EXTRACT(ng)%umask')) RETURN

      IF (.not.destroy(ng, EXTRACT(ng)%vmask, MyFile,                   &
     &                 __LINE__, 'EXTRACT(ng)%vmask')) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type GRID structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(EXTRACT)) deallocate ( EXTRACT )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_extract
!
      SUBROUTINE initialize_extract (ng, tile, model)
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
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
# ifdef SOLVE3D
      integer :: k
# endif

      real(r8), parameter :: IniVal = 0.0_r8
      real(r8) :: IniMetricVal

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=xtr_BOUNDS(ng)%LBi(tile)
      Imax=xtr_BOUNDS(ng)%UBi(tile)
      Jmin=xtr_BOUNDS(ng)%LBj(tile)
      Jmax=xtr_BOUNDS(ng)%UBj(tile)
# else
      IF (xtr_DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=xtr_BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=xtr_BOUNDS(ng)%Istr(tile)
      END IF
      IF (xtr_DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=xtr_BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=xtr_BOUNDS(ng)%Iend(tile)
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=xtr_BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=xtr_BOUNDS(ng)%Jstr(tile)
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=xtr_BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=xtr_BOUNDS(ng)%Jend(tile)
      END IF
# endif
!
!  Set initialization value that it is special in nexting to just
!  load contact points that have not been initialized from the
!  regular physical grid. This is done to make sure that all these
!  important metric values have been set-up correctly.
!
# ifdef NESTING
      IniMetricVal=spval                   ! very large value
# else
      IniMetricVal=IniVal
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
            EXTRACT(ng) % angler(i,j) = IniMetricVal
            EXTRACT(ng) % CosAngler(i,j) = IniVal
            EXTRACT(ng) % SinAngler(i,j) = IniVal

# if defined CURVGRID && defined UV_ADV
            EXTRACT(ng) % dmde(i,j) = IniMetricVal
            EXTRACT(ng) % dndx(i,j) = IniMetricVal
# endif
            EXTRACT(ng) % f(i,j) = IniMetricVal
            EXTRACT(ng) % h(i,j) = IniMetricVal

            EXTRACT(ng) % latp(i,j) = IniVal
            EXTRACT(ng) % latr(i,j) = IniMetricVal
            EXTRACT(ng) % latu(i,j) = IniMetricVal
            EXTRACT(ng) % latv(i,j) = IniMetricVal
            EXTRACT(ng) % lonp(i,j) = IniVal
            EXTRACT(ng) % lonr(i,j) = IniMetricVal
            EXTRACT(ng) % lonu(i,j) = IniMetricVal
            EXTRACT(ng) % lonv(i,j) = IniMetricVal
            EXTRACT(ng) % MyLon(i,j) = IniMetricVal

            EXTRACT(ng) % pm(i,j) = IniMetricVal
            EXTRACT(ng) % pn(i,j) = IniMetricVal

            EXTRACT(ng) % xp(i,j) = IniVal
            EXTRACT(ng) % xr(i,j) = IniMetricVal
            EXTRACT(ng) % xu(i,j) = IniMetricVal
            EXTRACT(ng) % xv(i,j) = IniMetricVal
            EXTRACT(ng) % yp(i,j) = IniVal
            EXTRACT(ng) % yr(i,j) = IniMetricVal
            EXTRACT(ng) % yu(i,j) = IniMetricVal
            EXTRACT(ng) % yv(i,j) = IniMetricVal

# ifdef MASKING
            EXTRACT(ng) % pmask(i,j) = IniVal
            EXTRACT(ng) % rmask(i,j) = IniMetricVal
            EXTRACT(ng) % umask(i,j) = IniMetricVal
            EXTRACT(ng) % vmask(i,j) = IniMetricVal
# endif
          END DO

# ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              EXTRACT(ng) % Hz(i,j,k) = IniVal
              EXTRACT(ng) % z_r(i,j,k) = IniVal
              EXTRACT(ng) % z_v(i,j,k) = IniVal
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              EXTRACT(ng) % z_w(i,j,k) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
!
      RETURN
      END SUBROUTINE initialize_extract
#endif
      END MODULE mod_extract
