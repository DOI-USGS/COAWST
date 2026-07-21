#include "cppdefs.h"
      MODULE pack_field_mod
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!=======================================================================
!                                                                      !
!  These routines pack output variables data into a 1D array. If       !
!  distributed memory, it gathers tiled data from all processors       !
!  and serializes.                                                     !
!                                                                      !
!  Then, if applicable, they extract fields at different geometry by   !
!  decimation or interpolation.                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod,   ONLY : mp_collect,                          &
     &                             mp_gather2d,                         &
# ifdef GRID_EXTRACT
     &                             mp_gather2d_xtr,                     &
# endif
     &                             mp_gather3d
#endif
#ifdef GRID_EXTRACT
# ifdef ADJUST_BOUNDARY
      USE extract_field_mod, ONLY: extract_boundary
# endif
      USE extract_field_mod, ONLY: extract_field
#endif
!
      implicit none
!
#ifdef ADJUST_BOUNDARY
      PUBLIC :: pack_boundary2d
      PUBLIC :: pack_boundary3d
#endif
      PUBLIC :: pack_field2d
      PUBLIC :: pack_field3d
      PUBLIC :: pack_field4d
!
      CONTAINS

#ifdef ADJUST_BOUNDARY
!
      SUBROUTINE pack_boundary2d (ng, model, tile,                      &
     &                            gtype, ncvname, tindex,               &
     &                            Extract_Flag,                         &
     &                            LBij, UBij, Nrec,                     &
     &                            Bscl, Bdat,                           &
     &                            start, total, Npts, Bwrk)
!
!=======================================================================
!                                                                      !
!  Packs 2D boundary data into 1D array to be written into output      !
!  NetCDF file elsewhere.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction flag interpolation/decimation (integer)  !
!     LBij         Donor IJ-dimension Lower bound (integer)            !
!     UBij         Donor IJ-dimension Upper bound (integer)            !
!     Nrec         Number of boundary records (integer)                !
!     Bscl         Factor to scale boundary data before writing (real) !
!     Bdat         2D boundary data to process (real)                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Bwrk (integer)        !
!     Bwrk         Packed 2D boundary data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)   :: ng, model, tile
      integer, intent(in)   :: gtype, tindex
      integer, intent(in)   :: Extract_Flag
      integer, intent(in)   :: LBij, UBij, Nrec
      integer, intent(in)   :: Npts
      integer, intent(out)  :: start(:), total(:)
!
      real(dp), intent(in)  :: Bscl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: Bdat(LBij:,:,:)
# else
      real(r8), intent(in)  :: Bdat(LBij:LBij,4,Nrec)
# endif
      real(r8), intent(out) :: Bwrk(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      logical, dimension(4) :: bounded
!
      integer :: bc, i, ib, ic, ir, j, rc
      integer :: IorJ, Imin, Imax, Jmin, Jmax
!
      real(r8), parameter :: Bspv = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pack_boundary2d"
!
!------------------------------------------------------------------------
!  Pack 2D boundary data.
!------------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification.
!
# ifdef DISTRIBUTE
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (p2dvar, p3dvar)               ! PSI points field
!       ---------------------
!
          Imin=BOUNDS(ng)%Istr (tile)
          Imax=BOUNDS(ng)%Iend (tile)
          Jmin=BOUNDS(ng)%Jstr (tile)
          Jmax=BOUNDS(ng)%Jend (tile)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          Imin=BOUNDS(ng)%Istr (tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%Jstr (tile)
          Jmax=BOUNDS(ng)%JendR(tile)
!
!       ---------------------
        CASE DEFAULT                        ! RHO points field
!       ---------------------
!
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
      END SELECT
# else
      Imin=LBij
      Imax=UBij
      Jmin=LBij
      Jmax=UBij
# endif
!
      IorJ=IOBOUNDS(ng)%IorJ
!
!  Set switch to process boundary data by their associated tiles.
!
      bounded(iwest )=DOMAIN(ng)%Western_Edge (tile)
      bounded(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      bounded(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      bounded(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!-----------------------------------------------------------------------
!  Pack and scale output data.
!-----------------------------------------------------------------------
!
      Bwrk=Bspv
!
      DO ir=1,Nrec
        rc=(ir-1)*IorJ*4
        DO ib=1,4
          IF (bounded(ib)) THEN
            bc=(ib-1)*IorJ+rc
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              DO j=Jmin,Jmax
                ic=1+(j-LBij)+bc
                Bwrk(ic)=Bdat(j,ib,ir)*Bscl
# ifdef POSITIVE_ZERO
                IF (ABS(Bwrk(ic)).eq.0.0_r8) THEN
                  Bwrk(ic)=0.0_r8                 ! impose positive zero
                END IF
# endif
              END DO
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              DO i=Imin,Imax
                ic=1+(i-LBij)+bc
                Bwrk(ic)=Bdat(i,ib,ir)*Bscl
# ifdef POSITIVE_ZERO
                IF (ABS(Bwrk(ic)).eq.0.0_r8) THEN
                  Bwrk(ic)=0.0_r8                 ! impose positive zero
                END IF
# endif
              END DO
            END IF
          END IF
        END DO
      END DO

# ifdef DISTRIBUTE
!
!  If distributed-memory set-up, collect data from all spawned
!  processes.
!
      CALL mp_collect (ng, model, Npts, Bspv, Bwrk)
# endif
!
!-----------------------------------------------------------------------
!  If there is no extracting data, set the start and total vectors
!  needed by the NetCDF library for writing.
!-----------------------------------------------------------------------
!
      IF (Extract_Flag.le.0) THEN
        start(1)=1
        total(1)=IorJ
        start(2)=1
        total(2)=4
        start(3)=1
        total(3)=Nrec
        start(4)=tindex
        total(4)=1

# ifdef GRID_EXTRACT
!
!-----------------------------------------------------------------------
!  Otherwise, extract field by decimation or interpolation.
!-----------------------------------------------------------------------
!
      ELSE IF (Extract_Flag.ge.1) THEN
        CALL extract_boundary (ng, model, tile,                         &
     &                         gtype, ncvname, tindex,                  &
     &                         Extract_flag,                            &
     &                         Imin, Imax, Jmin, Jmax, Nrec,            &
     &                         start, total,                            &
     &                         Npts, Bwrk)
# endif
      END IF
!
      RETURN
      END SUBROUTINE pack_boundary2d
!
      SUBROUTINE pack_boundary3d (ng, model, tile,                      &
     &                            gtype, ncvname, tindex,               &
     &                            Extract_Flag,                         &
     &                            LBij, UBij, LBk, UBk, Nrec,           &
     &                            Bscl, Bdat,                           &
     &                            start, total, Npts, Bwrk)
!
!=======================================================================
!                                                                      !
!  Packs 3D boundary data into 1D array to be written into output      !
!  NetCDF file elsewhere.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction flag interpolation/decimation (integer)  !
!     LBij         Donor IJ-dimension Lower bound (integer)            !
!     UBij         Donor IJ-dimension Upper bound (integer)            !
!     LBk          Donor K-dimension lower bound (integer)             !
!     UBk          Donor K-dimension upper bound (integer)             !
!     Nrec         Number of boundary records (integer)                !
!     Bscl         Factor to scale boundary data before writing (real) !
!     Bdat         3D boundary data to process (real)                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Bwrk (integer)        !
!     Bwrk         Packed 3D boundary data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)   :: ng, model, tile
      integer, intent(in)   :: gtype, tindex
      integer, intent(in)   :: Extract_Flag
      integer, intent(in)   :: LBij, UBij, LBk, UBk, Nrec
      integer, intent(in)   :: Npts
      integer, intent(out)  :: start(:), total(:)
!
      real(dp), intent(in)  :: Bscl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: Bdat(LBij:,LBk:,:,:)
# else
      real(r8), intent(in)  :: Bdat(LBij:LBij,LBk:UBk,4,Nrec)
# endif
      real(r8), intent(out) :: Bwrk(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      logical, dimension(4) :: bounded
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: IorJ, IJKlen, Klen
      integer :: bc, i, ib, ic, ir, j, k, kc, rc
!
      real(r8), parameter :: Bspv = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pack_boundary3d"
!
!------------------------------------------------------------------------
!  Pack 3D boundary data.
!------------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification.
!
# ifdef DISTRIBUTE
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (p2dvar, p3dvar)               ! PSI points field
!       ---------------------
!
          Imin=BOUNDS(ng)%Istr (tile)
          Imax=BOUNDS(ng)%Iend (tile)
          Jmin=BOUNDS(ng)%Jstr (tile)
          Jmax=BOUNDS(ng)%Jend (tile)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          Imin=BOUNDS(ng)%Istr (tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%Jstr (tile)
          Jmax=BOUNDS(ng)%JendR(tile)
!
!       ---------------------
        CASE DEFAULT                        ! RHO points field
!       ---------------------
!
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
      END SELECT
# else
      Imin=LBij
      Imax=UBij
      Jmin=LBij
      Jmax=UBij
# endif
!
      IorJ=IOBOUNDS(ng)%IorJ
      Klen=UBk-LBk+1
      IJKlen=IorJ*Klen
!
!  Set switch to process boundary data by their associated tiles.
!
      bounded(iwest )=DOMAIN(ng)%Western_Edge (tile)
      bounded(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      bounded(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      bounded(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!-----------------------------------------------------------------------
!  Pack and scale output data.
!-----------------------------------------------------------------------
!
      Bwrk=Bspv
!
      DO ir=1,Nrec
        rc=(ir-1)*IJKlen*4
        DO ib=1,4
          IF (bounded(ib)) THEN
            bc=(ib-1)*IJKlen+rc
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              DO k=LBk,UBk
                kc=(k-LBk)*IorJ+bc
                DO j=Jmin,Jmax
                  ic=1+(j-LBij)+kc
                  Bwrk(ic)=Bdat(j,k,ib,ir)*Bscl
# ifdef POSITIVE_ZERO
                  IF (ABS(Bwrk(ic)).eq.0.0_r8) THEN
                    Bwrk(ic)=0.0_r8               ! impose positive zero
                  END IF
# endif
                END DO
              END DO
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              DO k=LBk,UBk
                kc=(k-LBk)*IorJ+bc
                DO i=Imin,Imax
                  ic=1+(i-LBij)+kc
                  Bwrk(ic)=Bdat(i,k,ib,ir)*Bscl
# ifdef POSITIVE_ZERO
                  IF (ABS(Bwrk(ic)).eq.0.0_r8) THEN
                    Bwrk(ic)=0.0_r8               ! impose positive zero
                  END IF
# endif
                END DO
              END DO
            END IF
          END IF
        END DO
      END DO

# ifdef DISTRIBUTE
!
!  If distributed-memory set-up, collect data from all spawned
!  processes.
!
      CALL mp_collect (ng, model, Npts, Bspv, Bwrk)
# endif
!
!-----------------------------------------------------------------------
!  If there is no extracting data, set the start and total vectors
!  needed by the NetCDF library for writing.
!-----------------------------------------------------------------------
!
      IF (Extract_Flag.le.0) THEN
        start(1)=1
        total(1)=IorJ
        start(2)=1
        total(2)=Klen
        start(3)=1
        total(3)=4
        start(4)=1
        total(4)=Nrec
        start(5)=tindex
        total(5)=1

# ifdef GRID_EXTRACT
!
!-----------------------------------------------------------------------
!  Otherwise, extract field by decimation or interpolation.
!-----------------------------------------------------------------------
!
      ELSE IF (Extract_Flag.ge.1) THEN
        CALL extract_boundary (ng, model, tile,                         &
     &                         gtype, ncvname, tindex,                  &
     &                         Extract_flag,                            &
     &                         Imin, Imax, Jmin, Jmax, LBk, UBk, Nrec,  &
     &                         start, total,                            &
     &                         Npts, Bwrk)
# endif
      END IF
!
      RETURN
      END SUBROUTINE pack_boundary3d
#endif
!
      SUBROUTINE pack_field2d (ng, model, tile,                         &
     &                         gtype, ifield, tindex,                   &
     &                         LandFill, Extract_Flag,                  &
     &                         LBi, UBi, LBj, UBj,                      &
#ifdef MASKING
     &                         Amask,                                   &
#endif
     &                         Ascl, Adat,                              &
     &                         start, total, Npts, Awrk)
!
!=======================================================================
!                                                                      !
!  Packs 2D field data into 1D array to be written into output NetCDF  !
!  file elsewhere.                                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     LandFill     Switch to set fill value in land areas (logical)    !                                         !
!     Extract_Flag Extraction flag interpolation/decimation (integer)  !
!     LBi          Donor field I-dimension Lower bound (integer)       !
!     UBi          Donor field I-dimension Upper bound (integer)       !
!     LBj          Donor field J-dimension Lower bound (integer)       !
!     UBj          Donor field J-dimension Upper bound (integer)       !
#ifdef MASKING
!     Amask        Land/Sea mask, if any (real 2D array)               !
#endif
!     Ascl         Factor to scale field after reading (real)          !
!     Adat         2D field to process (real)                          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Awrk.                 !
!     Awrk         Packed 2D field data (real)                         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in)  :: LandFill
!
      integer, intent(in)  :: ng, model, tile
      integer, intent(in)  :: gtype, ifield, tindex
      integer, intent(in)  :: Extract_Flag
      integer, intent(in)  :: LBi, UBi, LBj, UBj
      integer, intent(out) :: Npts
      integer, intent(out) :: start(:), total(:)
!
      real(dp), intent(in) :: Ascl
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Adat(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Adat(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
      integer :: i, ic, j
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Isize, Jsize, IJsize, MyType
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pack_field2d"

!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
#ifdef GRID_EXTRACT
!
!  Here, the ELSE statements (Extract_Flag < 0) indicate the processing
!  of data already on the extracted grid, like geometry variables
!  requiring no extraction since they are read during initialization.
#endif
!
      MyType=gtype
!
      SELECT CASE (ABS(MyType))
!
!       ---------------------
        CASE (p2dvar, p3dvar)               ! PSI points field
!       ---------------------
!
          IF (Extract_Flag.ge.0) THEN
            Imin=IOBOUNDS(ng)%ILB_psi
            Imax=IOBOUNDS(ng)%IUB_psi
            Jmin=IOBOUNDS(ng)%JLB_psi
            Jmax=IOBOUNDS(ng)%JUB_psi
#ifdef GRID_EXTRACT
          ELSE                              ! Already on extract grid
            Imin=xtr_IOBOUNDS(ng)%ILB_psi
            Imax=xtr_IOBOUNDS(ng)%IUB_psi
            Jmin=xtr_IOBOUNDS(ng)%JLB_psi
            Jmax=xtr_IOBOUNDS(ng)%JUB_psi
#endif
          END IF
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          IF (Extract_Flag.ge.0) THEN
            Imin=IOBOUNDS(ng)%ILB_rho
            Imax=IOBOUNDS(ng)%IUB_rho
            Jmin=IOBOUNDS(ng)%JLB_rho
            Jmax=IOBOUNDS(ng)%JUB_rho
#ifdef GRID_EXTRACT
          ELSE                              ! Already on extract grid
            Imin=xtr_IOBOUNDS(ng)%ILB_rho
            Imax=xtr_IOBOUNDS(ng)%IUB_rho
            Jmin=xtr_IOBOUNDS(ng)%JLB_rho
            Jmax=xtr_IOBOUNDS(ng)%JUB_rho
#endif
          END IF
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          IF (Extract_Flag.ge.0) THEN
            Imin=IOBOUNDS(ng)%ILB_u
            Imax=IOBOUNDS(ng)%IUB_u
            Jmin=IOBOUNDS(ng)%JLB_u
            Jmax=IOBOUNDS(ng)%JUB_u
#ifdef GRID_EXTRACT
          ELSE                              ! Already on extract grid
            Imin=xtr_IOBOUNDS(ng)%ILB_u
            Imax=xtr_IOBOUNDS(ng)%IUB_u
            Jmin=xtr_IOBOUNDS(ng)%JLB_u
            Jmax=xtr_IOBOUNDS(ng)%JUB_u
#endif
          END IF
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          IF (Extract_Flag.ge.0) THEN
            Imin=IOBOUNDS(ng)%ILB_v
            Imax=IOBOUNDS(ng)%IUB_v
            Jmin=IOBOUNDS(ng)%JLB_v
            Jmax=IOBOUNDS(ng)%JUB_v
#ifdef GRID_EXTRACT
          ELSE                              ! Already on extract grid
            Imin=xtr_IOBOUNDS(ng)%ILB_v
            Imax=xtr_IOBOUNDS(ng)%IUB_v
            Jmin=xtr_IOBOUNDS(ng)%JLB_v
            Jmax=xtr_IOBOUNDS(ng)%JUB_v
#endif
          END IF
!
!       ---------------------
        CASE DEFAULT                        ! RHO points field
!       ---------------------
!
          IF (Extract_Flag.ge.0) THEN
            Imin=IOBOUNDS(ng)%ILB_rho
            Imax=IOBOUNDS(ng)%IUB_rho
            Jmin=IOBOUNDS(ng)%JLB_rho
            Jmax=IOBOUNDS(ng)%JUB_rho
#ifdef GRID_EXTRACT
          ELSE                              ! Already on extract grid
            Imin=xtr_IOBOUNDS(ng)%ILB_rho
            Imax=xtr_IOBOUNDS(ng)%IUB_rho
            Jmin=xtr_IOBOUNDS(ng)%JLB_rho
            Jmax=xtr_IOBOUNDS(ng)%JUB_rho
#endif
          END IF
      END SELECT
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      IJsize=Isize*Jsize

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order.
# ifdef MASKING
#  ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
      IF (Extract_Flag.ge.0) THEN
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,                &
     &                    tindex, gtype, Ascl,                          &
# ifdef MASKING
     &                    Amask,                                        &
# endif
     &                    Adat, Npts, Awrk, LandFill)
# ifdef GRID_EXTRACT
      ELSE                                  ! Already on extract grid
        CALL mp_gather2d_xtr (ng, model, LBi, UBi, LBj, UBj,            &
     &                        tindex, gtype, Ascl,                      &
#  ifdef MASKING
     &                        Amask,                                    &
#  endif
     &                        Adat, Npts, Awrk, LandFill)
# endif
      END IF

#else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order.
# ifdef MASKING
#  ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
        ic=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ic=ic+1
            Awrk(ic)=Adat(i,j)*Ascl
# ifdef MASKING
            IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk(ic)=spval
            END IF
# endif
          END DO
        END DO
        Npts=ic
# ifdef MASKING
      ELSE
        ic=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            IF (Amask(i,j).gt.0.0_r8) THEN
              ic=ic+1
              Awrk(ic)=Adat(i,j)*Ascl
            END IF
          END DO
        END DO
        Npts=ic
# endif
      END IF
#endif
!
!-----------------------------------------------------------------------
!  If there is no extracting data, set the start and total vectors
!  needed by the NetCDF library for writing.
!-----------------------------------------------------------------------
!
      IF (Extract_Flag.le.0) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Isize
          start(2)=1
          total(2)=Jsize
          start(3)=tindex
          total(3)=1
#ifdef MASKING
        ELSE
          start(1)=1
          total(1)=Npts
          start(2)=tindex
          total(2)=1
#endif
        END IF

#ifdef GRID_EXTRACT
!
!-----------------------------------------------------------------------
!  Otherwise, extract field by decimation or interpolation.
!-----------------------------------------------------------------------
!
      ELSE IF (Extract_Flag.ge.1) THEN
        CALL extract_field (ng, model, tile,                            &
     &                      gtype, ifield, tindex,                      &
     &                      Extract_Flag,                               &
     &                      Imin, Imax, Jmin, Jmax,                     &
     &                      start, total,                               &
     &                      Npts, Awrk)
#endif
      END IF
!
      RETURN
      END SUBROUTINE pack_field2d
!
      SUBROUTINE pack_field3d (ng, model, tile,                         &
     &                         gtype, ifield, tindex,                   &
     &                         LandFill, Extract_Flag,                  &
     &                         LBi, UBi, LBj, UBj, LBk, UBk,            &
#ifdef MASKING
     &                         Amask,                                   &
#endif
     &                         Ascl, Adat,                              &
     &                         start, total, Npts, Awrk)
!
!=======================================================================
!                                                                      !
!  Packs 3D field data into 1D array to be written into output NetCDF  !
!  file elsewhere.                                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     LandFill     Switch to set fill value in land areas (logical)    !
!     Extract_Flag Extraction flag interpolation/decimation (integer)  !
!     LBi          Field I-dimension Lower bound (integer)             !
!     UBi          Field I-dimension Upper bound (integer)             !
!     LBj          Field J-dimension Lower bound (integer)             !
!     UBj          Field J-dimension Upper bound (integer)             !
!     LBk          Field K-dimension Lower bound (integer)             !
!     UBk          Field K-dimension Upper bound (integer)             !
#ifdef MASKING
!     Amask        Land/Sea mask, if any (real 2D array)               !
#endif
!     Ascl         Factor to scale field after reading (real)          !
!     Adat         3D field to process (real)                          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Awrk.                 !
!     Awrk         Packed 3D field data (real)                         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in)  :: LandFill
!
      integer, intent(in)  :: ng, model, tile
      integer, intent(in)  :: gtype, ifield, tindex
      integer, intent(in)  :: Extract_Flag
      integer, intent(in)  :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(out) :: Npts
      integer, intent(out) :: start(:), total(:)
!
      real(dp), intent(in) :: Ascl
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Adat(LBi:,LBj:,LBk:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Adat(LBi:UBi,LBj:UBj,LBk:UBk)
#endif
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
      integer :: i, j, k, ic
      integer :: Imin, Imax, Jmin, Jmax, Koff
      integer :: Isize, Jsize, Ksize, IJsize, MyType
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pack_field3d"
!
!------------------------------------------------------------------------
!  Pack 3D field data.
!------------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      MyType=gtype
!
      SELECT CASE (ABS(MyType))
!
!       ---------------------
        CASE (p3dvar)                       ! PSI points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
!
!       ---------------------
        CASE (r3dvar)                       ! RHO points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
!
!       ---------------------
        CASE (u3dvar)                       ! U points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
!
!       ---------------------
        CASE (v3dvar)                       ! V points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_v
          Imax=IOBOUNDS(ng)%IUB_v
          Jmin=IOBOUNDS(ng)%JLB_v
          Jmax=IOBOUNDS(ng)%JUB_v
!
!       ---------------------
        CASE DEFAULT                        ! RHO points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
      END SELECT
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      Ksize=UBk-LBk+1
      IJsize=Isize*Jsize
!
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order.
# ifdef MASKING
#  ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
      CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                  tindex, gtype, Ascl,                            &
# ifdef MASKING
     &                  Amask,                                          &
# endif
     &                  Adat, Npts, Awrk, LandFill)
#else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order.
# ifdef MASKING
#  ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              Awrk(ic)=Adat(i,j,k)*Ascl
# ifdef MASKING
              IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
                Awrk(ic)=spval
              END IF
# endif
            END DO
          END DO
        END DO
        Npts=ic
# ifdef MASKING
      ELSE
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                ic=ic+1
                Awrk(ic)=Adat(i,j,k)*Ascl
              END IF
            END DO
          END DO
        END DO
        Npts=ic
# endif
      END IF
#endif
!
!-----------------------------------------------------------------------
!  If there is no extracting data, set the start and total vectors
!  needed by the NetCDF library for writing.
!-----------------------------------------------------------------------
!
      IF (Extract_Flag.le.0) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Isize
          start(2)=1
          total(2)=Jsize
          start(3)=1
          total(3)=Ksize
          start(4)=tindex
          total(4)=1
#ifdef MASKING
        ELSE
          start(1)=1
          total(1)=Npts
          start(2)=tindex
          total(2)=1
#endif
        END IF

#ifdef GRID_EXTRACT
!
!-----------------------------------------------------------------------
!  Otherwise, extract field by decimation or interpolation.
!-----------------------------------------------------------------------
!
      ELSE IF (Extract_Flag.ge.1) THEN
        CALL extract_field (ng, model, tile,                            &
     &                      gtype, ifield, tindex,                      &
     &                      Extract_Flag,                               &
     &                      Imin, Imax, Jmin, Jmax, LBk, UBk,           &
     &                      start, total,                               &
     &                      Npts, Awrk)
#endif
      END IF

#ifdef POSITIVE_ZERO
!
!-----------------------------------------------------------------------
!  Impose positive zero. Recall that IEEE allows both negative and
!  positive zeros, which affects bit-by-bit comparison betwee different
!  tile partitions.
!-----------------------------------------------------------------------
!
      DO ic=1,Npts
        IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
          Awrk(ic)=0.0_r8
        END IF
      END DO
#endif
!
      RETURN
      END SUBROUTINE pack_field3d
!
      SUBROUTINE pack_field4d (ng, model, tile,                         &
     &                         gtype, ifield, tindex,                   &
     &                         LandFill, Extract_Flag,                  &
     &                         LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt,  &
     &                         fourth,                                  &
#ifdef MASKING
     &                         Amask,                                   &
#endif
     &                         Ascl, Adat,                              &
     &                         start, total, Npts, Awrk)
!
!=======================================================================
!                                                                      !
!  Packs 4D field data into 1D array to be written into output NetCDF  !
!  file elsewhere. The field data is processed by 3D slices to reduce  !
!  memory requirements.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     LandFill     Switch to set fill value in land areas (logical)    !
!     Extract_Flag Extraction flag interpolation/decimation (integer)  !
!     LBi          Field I-dimension Lower bound (integer)             !
!     UBi          Field I-dimension Upper bound (integer)             !
!     LBj          Field J-dimension Lower bound (integer)             !
!     UBj          Field J-dimension Upper bound (integer)             !
!     LBk          Field K-dimension Lower bound (integer)             !
!     UBk          Field K-dimension Upper bound (integer)             !
!     LBt          Field fourth-dimension Lower bound (integer)        !
!     UBt          Field fourth-dimension Upper bound (integer)        !
!     fourth       Fourth dimension index to process (integer)         !
#ifdef MASKING
!     Amask        Land/Sea mask, if any (real 2D array)               !
#endif
!     Ascl         Factor to scale field after reading (real)          !
!     Adat         4D field to process (real)                          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Awrk.                 !
!     Awrk         Packed 4D field data (real)                         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in)  :: LandFill
!
      integer, intent(in)  :: ng, model, tile
      integer, intent(in)  :: gtype, ifield, tindex
      integer, intent(in)  :: Extract_Flag
      integer, intent(in)  :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt
      integer, intent(in)  :: fourth
      integer, intent(out) :: Npts
      integer, intent(out) :: start(:), total(:)
!
      real(dp), intent(in) :: Ascl
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Adat(LBi:,LBj:,LBk:,LBt:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Adat(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
#endif
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
      integer :: i, j, k, ic
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, Koff, Loff
      integer :: Isize, Jsize, Ksize, IJsize, MyType
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pack_field4d"
!
!------------------------------------------------------------------------
!  Pack 4D field data by 3D slices.
!------------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      MyType=gtype
!
      SELECT CASE (ABS(MyType))
!
!       ---------------------
        CASE (p2dvar, p3dvar)               ! PSI points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_v
          Imax=IOBOUNDS(ng)%IUB_v
          Jmin=IOBOUNDS(ng)%JLB_v
          Jmax=IOBOUNDS(ng)%JUB_v
!
!       ---------------------
        CASE DEFAULT                        ! RHO points field
!       ---------------------
!
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
      END SELECT
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      Ksize=UBk-LBk+1
      IJsize=Isize*Jsize
!
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
!
      IF (LBt.eq.0) THEN
        Loff=1
      ELSE
        Loff=0
      END IF

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order.
# ifdef MASKING
#  ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
!  Process the data by 3D slices.
!
      CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                  tindex, gtype, Ascl,                            &
# ifdef MASKING
     &                  Amask,                                          &
# endif
     &                  Adat(:,:,:,fourth), Npts, Awrk, LandFill)
#else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order.
# ifdef MASKING
#  ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
!  The data is processed in 3D slices.
!
      IF (gtype.gt.0) THEN
        ic=0
        Npts=IJsize*Ksize
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              Awrk(ic)=Adat(i,j,k,fourth)*Ascl
# ifdef MASKING
              IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
                Awrk(ic)=spval
              END IF
# endif
            END DO
          END DO
        END DO
# ifdef MASKING
      ELSE
        Npts=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                Npts=Npts+1
                Awrk(Npts)=Adat(i,j,k,fourth)*Ascl
              END IF
            END DO
          END DO
        END DO
# endif
      END IF
#endif
!
!-----------------------------------------------------------------------
!  If there is no extracting data, set the start and total vectors
!  needed by the NetCDF library for writing.
!-----------------------------------------------------------------------
!
      IF (Extract_Flag.le.0) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Isize
          start(2)=1
          total(2)=Jsize
          start(3)=1
          total(3)=Ksize
          start(4)=fourth+Loff
          total(4)=1
          start(5)=tindex
          total(5)=1
#ifdef MASKING
        ELSE
          start(1)=1+(fourth+Loff-1)*Npts
          total(1)=Npts
          start(2)=tindex
          total(2)=1
#endif
        END IF

#ifdef GRID_EXTRACT
!
!-----------------------------------------------------------------------
!  Otherwise, extract field by decimation or interpolation.
!-----------------------------------------------------------------------
!
      ELSE IF (Extract_Flag.ge.1) THEN
        CALL extract_field (ng, model, tile,                            &
     &                      gtype, ifield, tindex,                      &
     &                      Extract_Flag,                               &
     &                      Imin, Imax, Jmin, Jmax, LBk, UBk,           &
     &                      start, total,                               &
     &                      Npts, Awrk)
#endif
      END IF

#ifdef POSITIVE_ZERO
!
!-----------------------------------------------------------------------
!  Impose positive zero. Recall that IEEE allows both negative and
!  positive zeros, which affects bit-by-bit comparison betwee different
!  tile partitions.
!-----------------------------------------------------------------------
!
      DO ic=1,Npts
        IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
          Awrk(ic)=0.0_r8
        END IF
      END DO
#endif
!
      RETURN
      END SUBROUTINE pack_field4d
!
      END MODULE pack_field_mod
