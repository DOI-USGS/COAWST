#include "cppdefs.h"
      MODULE extract_field_mod

#ifdef GRID_EXTRACT
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!=======================================================================
!                                                                      !
!  These routines extract output fields at coaser resolution using     !
!  decimation or horizontal interpolation. Both extraction strategies  !
!  required a grid NetCDF (XTRNAME) provided at execution.             !
!                                                                      !
!  If Extract_Flag > 1:   Decimation                                   !
!                                                                      !
!  It decimates the field solution at the prescribed integer factor,   !
!  Extract_Flag . For example, if Extract_Flag=2 (recommended), the    !
!  output fieldsare written at every other point, resulting in coarser !
!  data resolution. This strategy is advantageous in mixed resolution, !
!  split 4D-Var applications where the outer loop background (prior)   !
!  trajectory may be computed at a higher resolution than in the inner !
!  loop minimization to accelerate the calculations. For decimation to !
!  work, the number of parent grid RHO-points (0: Lm+1, 0:Mm+1) must   !
!  be multiples of dec_fator. That is,                                 !
!                                                                      !
!     MOD(Lm+1, Extract_Flag) = 0                                      !
!     MOD(Mm+1, Extract_Flag) = 0                                      !
!                                                                      !
!  If Extract_Flag = 1:   Interpolation                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_extract
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod,       ONLY : mp_gather2d, mp_gather2d_xtr
# endif
      USE roms_interpolate_mod, ONLY : cinterp2d, linterp2d, hindices
      USE strings_mod,          ONLY : FoundError
!
      implicit none
!
!  Interfaces for sam name routine overloading.
!
# ifdef ADJUST_BOUNDARY
      INTERFACE extract_boundary
        MODULE PROCEDURE extract_boundary2d
        MODULE PROCEDURE extract_boundary3d
      END INTERFACE extract_boundary
!
# endif
      INTERFACE extract_field
        MODULE PROCEDURE extract_field2d
        MODULE PROCEDURE extract_field3d
        MODULE PROCEDURE extract_field4d
      END INTERFACE extract_field
!
      PUBLIC  :: interp_coords
!
# ifdef ADJUST_BOUNDARY
      PRIVATE :: average_boundary2d
      PRIVATE :: average_boundary3d
# endif
      PRIVATE :: average_field2d
      PRIVATE :: average_field3d
# ifdef ADJUST_BOUNDARY
      PRIVATE :: decimate_boundary2d
      PRIVATE :: decimate_boundary3d
# endif
      PRIVATE :: decimate_field2d
      PRIVATE :: decimate_field3d
      PRIVATE :: decimate_field4d
      PRIVATE :: interp_field2d
      PRIVATE :: interp_field3d
      PRIVATE :: interp_field4d
      PRIVATE :: regrid_field2d
      PRIVATE :: regrid_field3d
!
!  Module parameters.
!
      logical, parameter :: ParallelOutput = .FALSE.
!
      integer, parameter :: Bilinear = 0    ! bilinear interpolation
      integer, parameter :: Bicubic  = 1    ! bicubic  interpolation
!
      CONTAINS
!
# ifdef ADJUST_BOUNDARY
!
      SUBROUTINE extract_boundary2d (ng, model, tile,                   &
     &                               gtype, ncvname, tindex,            &
     &                               Extract_flag,                      &
     &                               Imin, Imax, Jmin, Jmax,            &
     &                               Nrec,                              &
     &                               start, total,                      &
     &                               Npts, Bdat)
!
!=======================================================================
!                                                                      !
!  Extracts/samples boundary data from donor 2D array.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor boundary starting data I-index (integer)      !
!     Imax         Donor boundary ending   data I-index (integer)      !
!     Jmin         Donor boundary starting data J-index (integer)      !
!     Jmax         Donor boundary ending   data J-index (integer)      !
!     Nrec         Number of boundary records (integer)                !
!     Npts         Length of packed boundary data (integer)            !
!     Bdat         Packed global donor 2D boundary data (real, 1D)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Bdat         Extracted 2D boundary data (real, 1D)               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, tindex, Extract_flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)    :: Nrec
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Bdat(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_boundary2d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate or interpolate 2D boundary data.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        CALL decimate_boundary2d (ng, model, tile,                      &
     &                            gtype, ncvname, tindex, Extract_Flag, &
     &                            Imin, Imax, Jmin, Jmax,               &
     &                            Nrec,                                 &
     &                            start, total,                         &
     &                            Npts, Bdat)
      ELSE IF (Extract_Flag.eq.1) THEN
      END IF
!
      RETURN
      END SUBROUTINE extract_boundary2d
!
      SUBROUTINE extract_boundary3d (ng, model, tile,                   &
     &                               gtype, ncvname, tindex,            &
     &                               Extract_flag,                      &
     &                               Imin, Imax, Jmin, Jmax, LBk, UBk,  &
     &                               Nrec,                              &
     &                               start, total,                      &
     &                               Npts, Bdat)
!
!=======================================================================
!                                                                      !
!  Extracts/samples boundary data from donor 2D array.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor boundary starting data I-index (integer)      !
!     Imax         Donor boundary ending   data I-index (integer)      !
!     Jmin         Donor boundary starting data J-index (integer)      !
!     Jmax         Donor boundary ending   data J-index (integer)      !
!     LBk          Donor boundary K-dimension lower bound (integer)    !
!     UBk          Donor boundary K-dimension upper bound (integer)    !
!     Nrec         Number of boundary records (integer)                !
!     Npts         Length of packed boundary data (integer)            !
!     Bdat         Packed global donor 2D boundary data (real, 1D)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Bdat         Extracted 3D boundary data (real, 1D)               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, tindex, Extract_flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, LBk, UBk
      integer,  intent(in)    :: Nrec
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Bdat(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      integer :: ij
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_boundary3d"
!
!------------------------------------------------------------------------
!  Decimate or interpolate 2D boundary data.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        CALL decimate_boundary3d (ng, model, tile,                      &
     &                            gtype, ncvname, tindex, Extract_Flag, &
     &                            Imin, Imax, Jmin, Jmax, LBk, UBk,     &
     &                            Nrec,                                 &
     &                            start, total,                         &
     &                            Npts, Bdat)
      ELSE IF (Extract_Flag.eq.1) THEN
      END IF
!
      RETURN
      END SUBROUTINE extract_boundary3d
# endif
!
      SUBROUTINE extract_field2d (ng, model, tile,                      &
     &                            gtype, ifield, tindex, Extract_Flag,  &
     &                            Imin, Imax, Jmin, Jmax,               &
     &                            start, total,                         &
     &                            Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It decimates or interpolates data from donor 2D fields, packed as   !
!  a 1D global array, to specified extract grid geometry.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Npts         Length of packed field data (integer)
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_field2d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate or interpolate 2D field.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        CALL decimate_field2d (ng, model, tile,                         &
     &                         gtype, ifield, tindex, Extract_Flag,     &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         start, total,                            &
     &                         Npts, Fdat)
      ELSE IF (Extract_Flag.eq.1) THEN
        IF (ParallelOutput) THEN
          CALL interp_field2d (ng, model, tile,                         &
     &                         gtype, ifield, tindex,                   &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         start, total,                            &
     &                         Npts, Fdat)
        ELSE
          CALL interp_field2d_global (ng, model, tile,                  &
     &                                gtype, ifield, tindex,            &
     &                                Imin, Imax, Jmin, Jmax,           &
     &                                start, total,                     &
     &                                Npts, Fdat)
        END IF
      END IF
!
      RETURN
      END SUBROUTINE extract_field2d
!
      SUBROUTINE extract_field3d (ng, model, tile,                      &
     &                            gtype, ifield, tindex, Extract_Flag,  &
     &                            Imin, Imax, Jmin, Jmax, Kmin, Kmax,   &
     &                            start, total,                         &
     &                            Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It decimates or interpolates data from donor 3D fields, packed as   !
!  a 1D global array, to specified extract grid geometry.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_field3d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate or interpolate 3D field.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        CALL decimate_field3d (ng, model, tile,                         &
     &                         gtype, ifield, tindex, Extract_Flag,     &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         Kmin, Kmax,                              &
     &                         start, total,                            &
     &                         Npts, Fdat)
      ELSE IF (Extract_Flag.eq.1) THEN
        IF (ParallelOutput) THEN
          CALL interp_field3d (ng, model, tile,                         &
     &                         gtype, ifield, tindex,                   &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         Kmin, Kmax,                              &
     &                         start, total,                            &
     &                         Npts, Fdat)
        ELSE
          CALL interp_field3d_global (ng, model, tile,                  &
     &                                gtype, ifield, tindex,            &
     &                                Imin, Imax, Jmin, Jmax,           &
     &                                Kmin, Kmax,                       &
     &                                start, total,                     &
     &                                Npts, Fdat)
        END IF
      END IF
!
      RETURN
      END SUBROUTINE extract_field3d
!
      SUBROUTINE extract_field4d (ng, model, tile,                      &
     &                            gtype, ifield, tindex, Extract_Flag,  &
     &                            Imin, Imax, Jmin, Jmax, Kmin, Kmax,   &
     &                            Fourth, Loff,                         &
     &                            start, total,                         &
     &                            Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It decimates or interpolates data from donor 4D fields, packed as   !
!  a 1D global array, to specified extract grid geometry. The field    !
!  data is processed by 3D slices to reduce memory requirements in     !
!  the calling routine.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     fourth       Donor Fouth dimension index to process (integer)    !
!     Loff         Fourth dimension couter offset (integer)            !
!     Fdat         Packed global donor 3D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Fdat (integer)        !
!     Fdat         Extracted 3D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)    :: fourth, Loff
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_field4d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate or interpolate 4D field.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        CALL decimate_field4d (ng, model, tile,                         &
     &                         gtype, ifield, tindex, Extract_Flag,     &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         Kmin, Kmax, Fourth, Loff,                &
     &                         start, total,                            &
     &                         Npts, Fdat)
      ELSE IF (Extract_Flag.eq.1) THEN
        IF (ParallelOutput) THEN
          CALL interp_field4d (ng, model, tile,                         &
     &                         gtype, ifield, tindex,                   &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         Kmin, Kmax, Fourth, Loff,                &
     &                         start, total,                            &
     &                         Npts, Fdat)
        ELSE
          CALL interp_field4d_global (ng, model, tile,                  &
     &                                gtype, ifield, tindex,            &
     &                                Imin, Imax, Jmin, Jmax,           &
     &                                Kmin, Kmax, Fourth, Loff,         &
     &                                start, total,                     &
     &                                Npts, Fdat)
        END IF
      END IF
!
      RETURN
      END SUBROUTINE extract_field4d

# ifdef ADJUST_BOUNDARY
!
      SUBROUTINE average_boundary2d (ng, model, tile,                   &
     &                               gtype, ncvname, Extract_Flag,      &
     &                               Imin, Imax, Jmin, Jmax, Nrec,      &
     &                               Npts, Bdat, Bavg)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 2D boundary packed as 1D global array. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor boundary starting data I-index (integer)      !
!     Imax         Donor boundary ending   data I-index (integer)      !
!     Jmin         Donor boundary starting data J-index (integer)      !
!     Jmax         Donor boundary ending   data J-index (integer)      !
!     Nrec         Number of boundary records (integer)                !
!     Npts         Length of packed boundary data (integer)            !
!     Bdat         Packed global donor 2D boundary data (real; 1D)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Favg         Averaged 2D boundary data (real 1D array)           !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile
      integer,  intent(in)  :: gtype, Extract_Flag
      integer,  intent(in)  :: Imin, Imax, Jmin, Jmax, Nrec
      integer,  intent(in)  :: Npts
!
      real(r8), intent(in)  :: Bdat(:)
      real(r8), intent(out) :: Bavg(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      logical, dimension(4) :: bounded
!
      integer  :: bc, i, ib, ij, ir, j, rc
      integer  :: IorJ
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", average_boundary2d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Average 2D U- and V-boundary data at the appropiate location.
!------------------------------------------------------------------------
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
!  The logic below works well for Extract_Flag=2.
!
      IF (Extract_Flag.ge.2) THEN
!
        SELECT CASE (gtype)
!
!         ---------------------
          CASE (u2dvar, u3dvar)             ! U points field:
!         ---------------------               Average data at RHO-points
!                                             for southern and northern
            Bavg=Bdat                       ! boundaries only
            DO ir=1,Nrec
              rc=(ir-1)*IorJ*4
              DO ib=1,4
                bc=(ib-1)*IorJ+rc
                IF (bounded(ib).and.                                    &
     &              ((ib.eq.isouth).or.(ib.eq.inorth))) THEN
                  DO i=Imin,Imax-1
                    ij=i+bc
                    Bavg(ij)=0.5_r8*(Bdat(ij)+Bdat(ij+1))
                  END DO
                END IF
              END DO
            END DO
!
!         ---------------------
          CASE (v2dvar, v3dvar)             ! V points field:
!         ---------------------               Average data at RHO-points
!                                             for western and eastern
            Bavg=Bdat                       ! boundaries only
            DO ir=1,Nrec
              rc=(ir-1)*IorJ*4
              DO ib=1,4
                bc=(ib-1)*IorJ+rc
                IF (bounded(ib).and.                                    &
     &              ((ib.eq.iwest).or.(ib.eq.ieast))) THEN
                  DO j=Jmin,Jmax-1
                    ij=j+bc
                    Bavg(ij)=0.5_r8*(Bdat(ij)+Bdat(ij+1))
                  END DO
                END IF
              END DO
            END DO
!
        END SELECT
      END IF
!
      RETURN
      END SUBROUTINE average_boundary2d
!
      SUBROUTINE average_boundary3d (ng, model, tile,                   &
     &                               gtype, ncvname, Extract_Flag,      &
     &                               Imin, Imax, Jmin, Jmax, LBk, UBk,  &
     &                               Nrec,                              &
     &                               Npts, Bdat, Bavg)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 3D boundary packed as 1D global array. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor boundary starting data I-index (integer)      !
!     Imax         Donor boundary ending   data I-index (integer)      !
!     Jmin         Donor boundary starting data J-index (integer)      !
!     Jmax         Donor boundary ending   data J-index (integer)      !
!     LBk          Donor boundary K-dimension lower bound (integer)    !
!     UBk          Donor boundary K-dimension upper bound (integer)    !
!     Nrec         Number of boundary records (integer)                !
!     Npts         Length of packed boundary data (integer)            !
!     Bdat         Packed global donor 2D boundary data (real; 1D)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Favg         Averaged 3D boundary data (real 1D array)           !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile
      integer,  intent(in)  :: gtype, Extract_Flag
      integer,  intent(in)  :: Imin, Imax, Jmin, Jmax, LBk, UBk, Nrec
      integer,  intent(in)  :: Npts
!
      real(r8), intent(in)  :: Bdat(:)
      real(r8), intent(out) :: Bavg(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      logical, dimension(4) :: bounded
!
      integer :: bc, i, ib, ij, ir, j, k, kc, rc
      integer :: IJKlen, IorJ, Klen
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", average_boundary3d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Average 3D U- and V-boundary data at the appropiate location.
!------------------------------------------------------------------------
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
!  The logic below works well for Extract_Flag=2.
!
      IF (Extract_Flag.ge.2) THEN
        Klen=UBk-LBk+1
        IJKlen=IorJ*Klen
!
        SELECT CASE (gtype)
!
!         ---------------------
          CASE (u3dvar)                     ! U points field:
!         ---------------------               Average data at RHO-points
!                                             for southern and northern
            Bavg=Bdat                       ! boundaries only
            DO ir=1,Nrec
              rc=(ir-1)*IorJ*4
              DO ib=1,4
                bc=(ib-1)*IJKlen+rc
                IF (bounded(ib).and.                                    &
     &              ((ib.eq.isouth).or.(ib.eq.inorth))) THEN
                  DO k=LBk,UBk
                    kc=(k-LBk)*IorJ+bc
                    DO i=Imin,Imax-1
                      ij=i+kc
                      Bavg(ij)=0.5_r8*(Bdat(ij)+Bdat(ij+1))
                    END DO
                  END DO
                END IF
              END DO
            END DO
!
!         ---------------------
          CASE (v2dvar, v3dvar)             ! V points field:
!         ---------------------               Average data at RHO-points
!                                             for western and eastern
            Bavg=Bdat                       ! boundaries only
            DO ir=1,Nrec
              rc=(ir-1)*IorJ*4
              DO ib=1,4
                bc=(ib-1)*IJKlen+rc
                IF (bounded(ib).and.                                    &
     &              ((ib.eq.iwest).or.(ib.eq.ieast))) THEN
                  DO k=LBk,UBk
                    kc=(k-LBk)*IorJ+bc
                    DO j=Jmin,Jmax-1
                      ij=j+kc
                      Bavg(ij)=0.5_r8*(Bdat(ij)+Bdat(ij+1))
                    END DO
                  END DO
                END IF
              END DO
            END DO
!
        END SELECT
      END IF
!
      RETURN
      END SUBROUTINE average_boundary3d
#endif
!
      SUBROUTINE average_field2d (ng, model, tile,                      &
     &                            gtype, ifield, Extract_Flag,          &
     &                            Imin, Imax, Jmin, Jmax,               &
     &                            Npts, Fdat, Favg)
!
!=======================================================================
!                                                                      !
!  It averages staggered 2D U- and V-fields at the appropiate location !
!  to facilitate decimation.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Favg         Averaged 2D field data (real 1D array)              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile
      integer,  intent(in)  :: gtype, ifield, Extract_Flag
      integer,  intent(in)  :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)  :: Npts
!
      real(r8), intent(in)  :: Fdat(:)
      real(r8), intent(out) :: Favg(:)
!
!  Local variable declarations.
!
      integer  :: i, j, ij
      integer  :: Ioff, Joff, Isize, Jsize
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", average_field2d"
!
!------------------------------------------------------------------------
!  Average 2D U- and V-field data at the appropiate location.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        Isize=Imax-Imin+1
        Jsize=Jmax-Jmin+1
!
        SELECT CASE (gtype)
!
!         ---------------------
          CASE (u2dvar, u3dvar)             ! U points field
!         ---------------------
!
            Favg=Fdat
            Ioff=0
            Joff=1
            DO j=Jmin,Jmax                  ! Average data at RHO points
              DO i=Imin,Imax-1
                ij=i+Ioff+(j-1+Joff)*Isize
                Favg(ij)=0.5_r8*(Fdat(ij)+Fdat(ij+1))
              END DO
            END DO
!
!         ---------------------
          CASE (v2dvar, v3dvar)             ! V points field
!         ---------------------
!
            Favg=Fdat
            Ioff=1
            Joff=0
            DO j=Jmin,Jmax-1                ! Average data at RHO points
              DO i=Imin,Imax
                ij=i+Ioff+(j-1+Joff)*Isize
                Favg(ij)=0.5_r8*(Fdat(ij)+Fdat(ij+Isize))
              END DO
            END DO
!
        END SELECT
      END IF
!
      RETURN
      END SUBROUTINE average_field2d
!
      SUBROUTINE average_field3d (ng, model, tile,                      &
     &                            gtype, ifield, Extract_Flag,          &
     &                            Imin, Imax, Jmin, Jmax, Kmin, Kmax,   &
     &                            Npts, Fdat, Favg)
!
!=======================================================================
!                                                                      !
!  It averages staggered 3D U- and V-fields at the appropiate location !
!  to facilitate decimation.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Favg         Averaged 2D field data (real 1D array)              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile
      integer,  intent(in)  :: gtype, ifield, Extract_Flag
      integer,  intent(in)  :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)  :: Npts
!
      real(r8), intent(in)  :: Fdat(:)
      real(r8), intent(out) :: Favg(:)
!
!  Local variable declarations.
!
      integer  :: i, j, k, ij, ijk
      integer  :: Ioff, Joff, Koff
      integer  :: Isize, Jsize, Ksize, IJsize
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", average_field3d"
!
!------------------------------------------------------------------------
!  Average 3D U- and V-field data at the appropiate location.
!------------------------------------------------------------------------
!
      IF (Extract_Flag.ge.2) THEN
        Isize=Imax-Imin+1
        Jsize=Jmax-Jmin+1
        Ksize=Kmax-Kmin+1
        IJsize=Isize*Jsize
!
        SELECT CASE (ABS(gtype))
!
!         -------------
          CASE (u3dvar)                     ! U points field
!         -------------
!
            Favg=Fdat
            Ioff=0
            Joff=1
            Koff=0
            DO k=Kmin,Kmax                  ! Average data at RHO points
              DO j=Jmin,Jmax
                DO i=Imin,Imax-1
                  ij=i+Ioff+(j-1+Joff)*Isize
                  ijk=ij+(k-1+Koff)*IJsize
                  Favg(ijk)=0.5_r8*(Fdat(ijk)+Fdat(ijk+1))
                END DO
              END DO
            END DO
!
!         -------------
          CASE (v3dvar)                     ! V points field
!         -------------
!
            Favg=Fdat
            Ioff=1
            Joff=0
            Koff=0
            DO k=Kmin,Kmax                  ! Average data at RHO points
              DO j=Jmin,Jmax-1
                DO i=Imin,Imax
                  ij=i+Ioff+(j-1+Joff)*Isize
                  ijk=ij+(k-1+Koff)*IJsize
                  Favg(ijk)=0.5_r8*(Fdat(ijk)+Fdat(ijk+Isize))
                END DO
              END DO
            END DO
!
        END SELECT
      END IF
!
      RETURN
      END SUBROUTINE average_field3d

# ifdef ADJUST_BOUNDARY
!
      SUBROUTINE decimate_boundary2d (ng, model, tile,                  &
     &                                gtype, ncvname, tindex,           &
     &                                Extract_Flag,                     &
     &                                Imin, Imax, Jmin, Jmax,           &
     &                                Nrec,                             &
     &                                start, total,                     &
     &                                Npts, Bdat)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 2D boundary packed as 1D global array. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor boundary starting data I-index (integer)      !
!     Imax         Donor boundary ending   data I-index (integer)      !
!     Jmin         Donor boundary starting data J-index (integer)      !
!     Jmax         Donor boundary ending   data J-index (integer)      !
!     Nrec         Number of boundary records (integer)                !
!     Npts         Length of packed boundary data (integer)            !
!     Bdat         Packed global donor 2D boundary data (real; 1D)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Bdat         Decimated 2D boundary data (real; 1D)               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)    :: Nrec, Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Bdat(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      integer  :: bc, ib, ic, ij, ir, ifactor, rc
      integer  :: Idim, IJdim, IorJ, Mpts
!
      real(r8) :: Bwrk(SIZE(Bdat))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", decimate_boundary2d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate input packed 2D boundary by requested factor.
!------------------------------------------------------------------------
!
      ifactor=ABS(Extract_Flag)
      IorJ=IOBOUNDS(ng)%IorJ
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Bwrk=Bdat
          Bdat=0.0_r8
          ic=0
          DO ir=1,Nrec
            rc=(ir-1)*IorJ*4
            DO ib=1,4
              bc=(ib-1)*IorJ+rc
              DO ij=1,IorJ,ifactor
                ic=ic+1
                Bdat(ic)=Bwrk(ij)
              END DO
              IF ((ir.eq.1).and.(ib.eq.1)) IJdim=ic
            END DO
          END DO
          Mpts=ic
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          CALL average_boundary2d (ng, model, tile,                     &
     &                             gtype, ncvname, Extract_Flag,        &
     &                             Imin, Imax, Jmin, Jmax, Nrec,        &
     &                             Npts, Bdat, Bwrk)
!
          Bdat=0.0_r8
          ic=0
          DO ir=1,Nrec
            rc=(ir-1)*IorJ*4
            DO ib=1,4
              bc=(ib-1)*IorJ+rc
              DO ij=1,IorJ,ifactor
                ic=ic+1
                Bdat(ic)=Bwrk(ij)
              END DO
              IF ((ir.eq.1).and.(ib.eq.1)) IJdim=ic
            END DO
          END DO
          Mpts=ic
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          CALL average_boundary2d (ng, model, tile,                     &
     &                             gtype, ncvname, Extract_Flag,        &
     &                             Imin, Imax, Jmin, Jmax, Nrec,        &
     &                             Npts, Bdat, Bwrk)
!
          Bdat=0.0_r8
          ic=0
          DO ir=1,Nrec
            rc=(ir-1)*IorJ*4
            DO ib=1,4
              bc=(ib-1)*IorJ+rc
              DO ij=1,IorJ,ifactor
                ic=ic+1
                Bdat(ic)=Bwrk(ij)
              END DO
              IF ((ir.eq.1).and.(ib.eq.1)) IJdim=ic
            END DO
          END DO
          Mpts=ic
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for decimation:',          &
     &                        TRIM(ncvname)
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      start(1)=1
      total(1)=IJdim
      start(2)=1
      total(2)=4
      start(3)=1
      total(3)=Nrec
      start(4)=tindex
      total(4)=1
!
  10  FORMAT (' DECIMATE_BOUNDARY2D - Staggered variable, gtype = ',    &
     &        i0,/,20x,a,2x,a)
!
      RETURN
      END SUBROUTINE decimate_boundary2d
!
      SUBROUTINE decimate_boundary3d (ng, model, tile,                  &
     &                                gtype, ncvname, tindex,           &
     &                                Extract_Flag,                     &
     &                                Imin, Imax, Jmin, Jmax, LBk, UBk, &
     &                                Nrec,                             &
     &                                start, total,                     &
     &                                Npts, Bdat)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 3D boundary packed as 1D global array. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ncvname      NetCDF variable name (string)                       !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor boundary starting data I-index (integer)      !
!     Imax         Donor boundary ending   data I-index (integer)      !
!     Jmin         Donor boundary starting data J-index (integer)      !
!     Jmax         Donor boundary ending   data J-index (integer)      !
!     LBk          Donor boundary K-dimension lower bound (integer)    !
!     UBk          Donor boundary K-dimension upper bound (integer)    !
!     Nrec         Number of boundary records (integer)                !
!     Npts         Length of packed boundary data (integer)            !
!     Bdat         Packed global donor 3D boundary data (real; 1D)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Bdat         Decimated 3D boundary data (real; 1D)               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, LBk, UBk
      integer,  intent(in)    :: Nrec, Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Bdat(:)
!
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations.
!
      integer  :: bc, ib, ic, ij, ir, ifactor, k, rc
      integer  :: Idim, IJdim, IJKlen, IorJ, Klen, Mpts
!
      real(r8) :: Bwrk(SIZE(Bdat))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", decimate_boundary3d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate input packed 3D boundary by requested factor.
!------------------------------------------------------------------------
!
      ifactor=ABS(Extract_Flag)
      IorJ=IOBOUNDS(ng)%IorJ
      Klen=UBk-LBk+1
      IJKlen=IorJ*Klen
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Bwrk=Bdat
          Bdat=0.0_r8
          ic=0
          DO ir=1,Nrec
            rc=(ir-1)*IJKlen*4
            DO ib=1,4
              bc=(ib-1)*IJKlen+rc
              DO k=LBk,UBk
                DO ij=1,IorJ,ifactor
                  ic=ic+1
                  Bwrk(ic)=Bdat(ij)
                END DO
                IF ((ir.eq.1).and.(ib.eq.1).and.(k.eq.LBk)) IJdim=ic
              END DO
            END DO
          END DO
          Mpts=ic
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          CALL average_boundary3d (ng, model, tile,                     &
     &                             gtype, ncvname, Extract_Flag,        &
     &                             Imin, Imax, Jmin, Jmax, LBk, UBk,    &
     &                             Nrec,                                &
     &                             Npts, Bdat, Bwrk)
!
          Bdat=0.0_r8
          ic=0
          DO ir=1,Nrec
            rc=(ir-1)*IJKlen*4
            DO ib=1,4
              bc=(ib-1)*IJKlen+rc
              DO k=LBk,UBk
                DO ij=1,IorJ,ifactor
                  ic=ic+1
                  Bwrk(ic)=Bdat(ij)
                END DO
                IF ((ir.eq.1).and.(ib.eq.1).and.(k.eq.LBk)) IJdim=ic
              END DO
            END DO
          END DO
          Mpts=ic
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          CALL average_boundary3d (ng, model, tile,                     &
     &                             gtype, ncvname, Extract_Flag,        &
     &                             Imin, Imax, Jmin, Jmax, LBk, UBk,    &
     &                             Nrec,                                &
     &                             Npts, Bdat, Bwrk)
!
          Bdat=0.0_r8
          ic=0
          DO ir=1,Nrec
            rc=(ir-1)*IJKlen*4
            DO ib=1,4
              bc=(ib-1)*IJKlen+rc
              DO k=LBk,UBk
                DO ij=1,IorJ,ifactor
                  ic=ic+1
                  Bwrk(ic)=Bdat(ij)
                END DO
                IF ((ir.eq.1).and.(ib.eq.1).and.(k.eq.LBk)) IJdim=ic
              END DO
            END DO
          END DO
          Mpts=ic
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for decimation:',          &
     &                        TRIM(ncvname)
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      start(1)=1
      total(1)=IJdim
      start(2)=1
      total(2)=Klen
      start(3)=1
      total(3)=4
      start(4)=1
      total(4)=Nrec
      start(5)=tindex
      total(5)=1
!
  10  FORMAT (' DECIMATE_BOUNDARY3D - Staggered variable, gtype = ',    &
     &        i0,/,20x,a,2x,a)
!
      RETURN
      END SUBROUTINE decimate_boundary3d
# endif
!
      SUBROUTINE decimate_field2d (ng, model, tile,                     &
     &                             gtype, ifield, tindex, Extract_Flag, &
     &                             Imin, Imax, Jmin, Jmax,              &
     &                             start, total,                        &
     &                             Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 2D fiels packed as 1D global array.    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Decimated 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      integer  :: i, j, ij, ic, jc, ifactor
      integer  :: Idim, Jdim, Ioff, Joff, Isize, Jsize, Mpts
!
      real(r8) :: Fwrk(SIZE(Fdat))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", decimate_field2d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate input packed 2D field by requested factor.
!------------------------------------------------------------------------
!
      ifactor=ABS(Extract_Flag)
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Fwrk=Fdat
          Fdat=0.0_r8
          Ioff=1
          Joff=1
          ic=0
          jc=0
          DO j=Jmin,Jmax,ifactor            ! decimate
            jc=jc+1
            DO i=Imin,Imax,ifactor
              ij=i+Ioff+(j-1+Joff)*Isize
              ic=ic+1
              Fdat(ic)=Fwrk(ij)
            END DO
            IF (j.eq.Jmin) Idim=ic
          END DO
          Jdim=jc
          Mpts=ic
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          CALL average_field2d (ng, model, tile,                        &
     &                          gtype, ifield, Extract_Flag,            &
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          Npts, Fdat, Fwrk)
!
          Fdat=0.0_r8
          Ioff=0
          Joff=1
          ic=0
          jc=0
          DO j=Jmin,Jmax,ifactor            ! decimate
            jc=jc+1
            DO i=Imin,Imax,ifactor
              ij=i+Ioff+(j-1+Joff)*Isize
              ic=ic+1
# ifdef MASKING
              Fdat(ic)=Fwrk(ij)*EXTRACT(ng)%Gmask_u(ic)
# else
              Fdat(ic)=Fwrk(ij)
# endif
            END DO
            IF (j.eq.Jmin) Idim=ic
          END DO
          Jdim=jc
          Mpts=ic
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          CALL average_field2d (ng, model, tile,                        &
     &                          gtype, ifield, Extract_Flag,            &
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          Npts, Fdat, Fwrk)
!
          Fdat=0.0_r8
          Ioff=1
          Joff=0
          ic=0
          jc=0
          DO j=Jmin,Jmax,ifactor            ! decimate
            jc=jc+1
            DO i=Imin,Imax,ifactor
              ij=i+Ioff+(j-1+Joff)*Isize
              ic=ic+1
# ifdef MASKING
              Fdat(ic)=Fwrk(ij)*EXTRACT(ng)%Gmask_v(ic)
# else
              Fdat(ic)=Fwrk(ij)
# endif
            END DO
            IF (j.eq.Jmin) Idim=ic
          END DO
          Jdim=jc
          Mpts=ic
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for decimation:',          &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Idim
        start(2)=1
        total(2)=Jdim
        start(3)=tindex
        total(3)=1
# ifdef MASKING
      ELSE
        start(1)=1
        total(1)=Mpts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' DECIMATE_FIELD2D - Staggered variable, gtype = ', i0,   &
     &        /,20x,a,2x,a)
!
      RETURN
      END SUBROUTINE decimate_field2d
!
      SUBROUTINE decimate_field3d (ng, model, tile,                     &
     &                             gtype, ifield, tindex, Extract_Flag, &
     &                             Imin, Imax, Jmin, Jmax, Kmin, Kmax,  &
     &                             start, total,                        &
     &                             Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 2D fiels packed as 1D global array.    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      integer  :: i, j, k, ij, ijk, ic, jc, ifactor, mc
      integer  :: Idim, Jdim, Ioff, Joff, Koff
      integer  :: Isize, IJsize, Jsize, Ksize, Mpts
!
      real(r8) :: Fwrk(SIZE(Fdat))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", decimate_field3d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate input packed 3D field by requested factor.
!------------------------------------------------------------------------
!
      ifactor=ABS(Extract_Flag)
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      Ksize=Kmax-Kmin+1
      IJsize=Isize*Jsize
!
      SELECT CASE (gtype)
!
!       -------------
        CASE (r3dvar)                       ! RHO points field
!       -------------
!
          Fwrk=Fdat
          Fdat=0.0_r8
          Ioff=1
          Joff=1
          Koff=0
          ic=0
          jc=0
          DO k=Kmin,Kmax
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
                Fdat(ic)=Fwrk(ijk)
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE (u3dvar)                       ! U points field
!       -------------
!
          CALL average_field3d (ng, model, tile,                        &
     &                          gtype, ifield, Extract_Flag,            &
     &                          Imin, Imax, Jmin, Jmax, Kmin, Kmax,     &
     &                          Npts, Fdat, Fwrk)
!
          Fdat=0.0_r8
          Ioff=0
          Joff=1
          Koff=0
          ic=0
          jc=0
          DO k=Kmin,Kmax
# ifdef MASKING
            mc=0
# endif
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
# ifdef MASKING
                mc=mc+1
                Fdat(ic)=Fwrk(ijk)*EXTRACT(ng)%Gmask_u(mc)
# else
                Fdat(ic)=Fwrk(ijk)
# endif
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE (v3dvar)                       ! V points field
!       -------------
!
          CALL average_field3d (ng, model, tile,                        &
     &                          gtype, ifield, Extract_Flag,            &
     &                          Imin, Imax, Jmin, Jmax, Kmin, Kmax,     &
     &                          Npts, Fdat, Fwrk)
!
          Fdat=0.0_r8
          Ioff=1
          Joff=0
          Koff=0
          ic=0
          jc=0
          DO k=Kmin,Kmax
# ifdef MASKING
            mc=0
# endif
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
# ifdef MASKING
                mc=mc+1
                Fdat(ic)=Fwrk(ijk)*EXTRACT(ng)%Gmask_v(mc)
# else
                Fdat(ic)=Fwrk(ijk)
# endif
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE (w3dvar)                       ! W points field
!       -------------
!
          Fwrk=Fdat
          Fdat=0.0_r8
          Ioff=1
          Joff=1
          Koff=1
          ic=0
          jc=0
          DO k=Kmin,Kmax
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
                Fdat(ic)=Fwrk(ijk)
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE DEFAULT
!       -------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for decimation:',          &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Idim
        start(2)=1
        total(2)=Jdim
        start(3)=1
        total(3)=Ksize
        start(4)=tindex
        total(4)=1
# ifdef MASKING
      ELSE
        start(1)=1
        total(1)=Mpts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' DECIMATE_FIELD3D - Staggered variable, gtype = ', i0,   &
     &        /,20x,a,2x,a)
!
      RETURN
      END SUBROUTINE decimate_field3d
!
      SUBROUTINE decimate_field4d (ng, model, tile,                     &
     &                             gtype, ifield, tindex, Extract_Flag, &
     &                             Imin, Imax, Jmin, Jmax, Kmin, Kmax,  &
     &                             Fourth, Loff,                        &
     &                             start, total,                        &
     &                             Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It decimates data from donor 4D field packed as 1D global array.    !
!  The data is processed by 3D slices to reduce memory requirements    !
!  in the calling routine.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     fourth       Donor Fouth dimension index to process (integer)    !
!     Loff         Fourth dimension couter offset (integer)            !
!     Fdat         Packed global donor 3D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 3D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex, Extract_Flag
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)    :: fourth, Loff
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      integer  :: i, j, k, ij, ijk, ic, jc, ifactor, mc
      integer  :: Idim, Jdim, Ioff, Joff, Koff
      integer  :: Isize, IJsize, Jsize, Ksize, Mpts
!
      real(r8) :: Fwrk(SIZE(Fdat))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", decimate_field4d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Decimate input packed 4D field by requested factor.
!------------------------------------------------------------------------
!
      ifactor=ABS(Extract_Flag)
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      Ksize=Kmax-Kmin+1
      IJsize=Isize*Jsize
!
      SELECT CASE (gtype)
!
!       -------------
        CASE (r3dvar)                       ! RHO points field
!       -------------
!
          Fwrk=Fdat
          Fdat=0.0_r8
          Ioff=1
          Joff=1
          Koff=0
          ic=0
          jc=0
          DO k=Kmin,Kmax
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
                Fdat(ic)=Fwrk(ijk)
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE (u3dvar)                       ! U points field
!       -------------
!
          CALL average_field3d (ng, model, tile,                        &
     &                          gtype, ifield, Extract_Flag,            &
     &                          Imin, Imax, Jmin, Jmax, Kmin, Kmax,     &
     &                          Npts, Fdat, Fwrk)
!
          Fdat=0.0_r8
          Ioff=0
          Joff=1
          Koff=0
          ic=0
          jc=0
          DO k=Kmin,Kmax
# ifdef MASKING
            mc=0
# endif
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
# ifdef MASKING
                mc=mc+1
                Fdat(ic)=Fwrk(ijk)*EXTRACT(ng)%Gmask_u(mc)
# else
                Fdat(ic)=Fwrk(ijk)
# endif
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE (v3dvar)                       ! V points field
!       -------------
!
          CALL average_field3d (ng, model, tile,                        &
     &                          gtype, ifield, Extract_Flag,            &
     &                          Imin, Imax, Jmin, Jmax, Kmin, Kmax,     &
     &                          Npts, Fdat, Fwrk)
!
          Fdat=0.0_r8
          Ioff=1
          Joff=0
          Koff=0
          ic=0
          jc=0
          DO k=Kmin,Kmax
# ifdef MASKING
            mc=0
# endif
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
# ifdef MASKING
                mc=mc+1
                Fdat(ic)=Fwrk(ijk)*EXTRACT(ng)%Gmask_v(mc)
# else
                Fdat(ic)=Fwrk(ijk)
# endif
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE (w3dvar)                       ! W points field
!       -------------
!
          Fwrk=Fdat
          Fdat=0.0_r8
          Ioff=1
          Joff=1
          Koff=1
          ic=0
          jc=0
          DO k=Kmin,Kmax
            DO j=Jmin,Jmax,ifactor
              DO i=Imin,Imax,ifactor
                ij=i+Ioff+(j-1+Joff)*Isize
                ijk=ij+(k-1+Koff)*IJsize
                ic=ic+1
                Fdat(ic)=Fwrk(ijk)
              END DO
              IF ((j.eq.Jmin).and.(k.eq.Kmin)) Idim=ic
              IF (k.eq.Kmin) jc=jc+1
            END DO
          END DO
          Jdim=jc
          Mpts=ic
!
!       -------------
        CASE DEFAULT
!       -------------
!
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)), gtype
  10      FORMAT (' DECIMATE_FIELD4D - Staggered variable, gtype = ',   &
     &            i0,/,20x,'not supported for field: ',a)
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Idim
        start(2)=1
        total(2)=Jdim
        start(3)=1
        total(3)=Ksize
        start(4)=fourth+Loff
        total(4)=1
        start(5)=tindex
        total(5)=1
# ifdef MASKING
      ELSE
        start(1)=1+(fourth+Loff-1)*Mpts
        total(1)=Mpts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
      RETURN
      END SUBROUTINE decimate_field4d
!
      SUBROUTINE interp_coords (ng, tile, model, gtype,                 &
     &                          Ainp,                                   &
# ifdef MASKING
     &                          Minp,                                   &
# endif
     &                          Xinp, Yinp,                             &
# ifdef MASKING
     &                          Mout,                                   &
# endif
     &                          Xout, Yout,                             &
     &                          Iout, Jout)
!
!=======================================================================
!                                                                      !
!  Suppose fields are extracted with the interpolation method. In that !
!  case, this routine computes the fractional I- and J-coordinates     !
!  (Iout, Jout) of the extracted data in terms of the donor grid       !
!  application to facilitate the interpolation. Ipos and Jpos are the  !
!  donor interpolant left-bottom cell corner coordinates containing    !
!  each extracted field grid point.                                    !
!                                                                      !
!  The 3D fields are interpolated level by level, so only horizontal   !
!  fractional coordinates are needed for each staggered location.      !
!                                                                      !
!  Since the interpolation is done during the writing of the output    !
!  extraction fields, the strategy here is to set the geometry arrays  !
!  as 1D global arrays to facilitate the interpolation between the     !
!  donor grid-tiled data in parallel applications. The interpolation   !
!  stencil is trivial in fractional global coordinates for bilinear    !
!  or bicubic methodologies.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ainp         Donor grid global curvilinear angle (1D real)       !
# ifdef MASKING
!     Minp         Donor grid global land/sea mask (1D real)           !
# endif
!     Xinp         Donor grid global X-coordinate  (1D real)           !
!     Yinp         Donor grid global Y-coordinate  (1D real)           !
# ifdef MASKING
!     Mout         Interpolated data global land/sea mask (1D real)    !
# endif
!     Xout         Interpolated data global X-coordinates (1D real)    !
!     Yout         Interpolated data global Y-coordinates (1D real)    !
!     Iout         Global fractional I-coordinate  (1D real)           !
!     Jout         Global fractional J-coordinate  (1D real)           !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, tile, model, gtype
!
      real(r8), intent(out) :: Ainp(:)
# ifdef MASKING
      real(r8), intent(out) :: Minp(:)
# endif
      real(r8), intent(out) :: Xinp(:)
      real(r8), intent(out) :: Yinp(:)
# ifdef MASKING
      real(r8), intent(out) :: Mout(:)
# endif
      real(r8), intent(out) :: Xout(:)
      real(r8), intent(out) :: Yout(:)
      real(r8), intent(out) :: Iout(:)
      real(r8), intent(out) :: Jout(:)
!
!  Local variable declarations.
!
      integer :: LBi_inp, UBi_inp, LBj_inp, UBj_inp
      integer :: LBi_out, UBi_out, LBj_out, UBj_out
      integer :: Is_inp, Ie_inp, Js_inp, Je_inp
      integer :: Is_out, Ie_out, Js_out, Je_out
      integer :: Isize, Jsize, Msize, Nsize
      integer :: Cgrid, Npts, ghost, i, ic, j
!
      real(dp) :: scale
!
      real(r8), allocatable :: angle(:,:)
# ifdef MASKING
      real(r8), pointer :: mask_inp(:,:), mask_out(:,:)
# endif
      real(r8), pointer :: Xi(:,:), Yi(:,:)
      real(r8), pointer :: Xo(:,:), Yo(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", interp_coord"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute fractional coordinates (Ipos, Jpos).
!-----------------------------------------------------------------------
!
      scale=1.0_dp
!
!  Get input donor and output extract grids arrays tiled bounds.
!
      LBi_inp=BOUNDS(ng)%LBi(tile)
      UBi_inp=BOUNDS(ng)%UBi(tile)
      LBj_inp=BOUNDS(ng)%LBj(tile)
      UBj_inp=BOUNDS(ng)%UBj(tile)
!
      LBi_out=xtr_BOUNDS(ng)%LBi(tile)
      UBi_out=xtr_BOUNDS(ng)%UBi(tile)
      LBj_out=xtr_BOUNDS(ng)%LBj(tile)
      UBj_out=xtr_BOUNDS(ng)%UBj(tile)
!
!  Get input donor and output extract grids I/O bounds and coordinates
!  according to staggered locations.
!
      SELECT CASE (gtype)
!
!       -----------------------------
        CASE (p2dvar, p3dvar)                    ! PSI-points
!       -----------------------------
!
          Cgrid=1
!
          Is_inp=IOBOUNDS(ng)%ILB_psi
          Ie_inp=IOBOUNDS(ng)%IUB_psi
          Js_inp=IOBOUNDS(ng)%JLB_psi
          Je_inp=IOBOUNDS(ng)%JUB_psi
!
          Is_out=xtr_IOBOUNDS(ng)%ILB_psi
          Ie_out=xtr_IOBOUNDS(ng)%IUB_psi
          Js_out=xtr_IOBOUNDS(ng)%JLB_psi
          Je_out=xtr_IOBOUNDS(ng)%JUB_psi
!
          IF (spherical) THEN
            Xi => GRID(ng)%lonp
            Yi => GRID(ng)%latp
            Xo => EXTRACT(ng)%lonp
            Yo => EXTRACT(ng)%latp
          ELSE
            Xi => GRID(ng)%xp
            Yi => GRID(ng)%yp
            Xo => EXTRACT(ng)%xp
            Yo => EXTRACT(ng)%yp
          END IF
# ifdef MASKING
          mask_inp => GRID(ng)%pmask
          mask_out => EXTRACT(ng)%pmask
# endif
!
!       -----------------------------
        CASE (r2dvar, r3dvar, w3dvar)            ! RHO-points
!       -----------------------------
!
          Cgrid=2
!
          Is_inp=IOBOUNDS(ng)%ILB_rho
          Ie_inp=IOBOUNDS(ng)%IUB_rho
          Js_inp=IOBOUNDS(ng)%JLB_rho
          Je_inp=IOBOUNDS(ng)%JUB_rho
!
          Is_out=xtr_IOBOUNDS(ng)%ILB_rho
          Ie_out=xtr_IOBOUNDS(ng)%IUB_rho
          Js_out=xtr_IOBOUNDS(ng)%JLB_rho
          Je_out=xtr_IOBOUNDS(ng)%JUB_rho
!
          IF (spherical) THEN
            Xi => GRID(ng)%lonr
            Yi => GRID(ng)%latr
            Xo => EXTRACT(ng)%lonr
            Yo => EXTRACT(ng)%latr
          ELSE
            Xi => GRID(ng)%xr
            Yi => GRID(ng)%yr
            Xo => EXTRACT(ng)%xr
            Yo => EXTRACT(ng)%yr
          END IF
# ifdef MASKING
          mask_inp => GRID(ng)%rmask
          mask_out => EXTRACT(ng)%rmask
# endif
!
!       -----------------------------
        CASE (u2dvar, u3dvar)                    ! U-points
!       -----------------------------
!
          Cgrid=3
!
          Is_inp=IOBOUNDS(ng)%ILB_u
          Ie_inp=IOBOUNDS(ng)%IUB_u
          Js_inp=IOBOUNDS(ng)%JLB_u
          Je_inp=IOBOUNDS(ng)%JUB_u
!
          Is_out=xtr_IOBOUNDS(ng)%ILB_u
          Ie_out=xtr_IOBOUNDS(ng)%IUB_u
          Js_out=xtr_IOBOUNDS(ng)%JLB_u
          Je_out=xtr_IOBOUNDS(ng)%JUB_u
!
          IF (spherical) THEN
            Xi => GRID(ng)%lonu
            Yi => GRID(ng)%latu
            Xo => EXTRACT(ng)%lonu
            Yo => EXTRACT(ng)%latu
          ELSE
            Xi => GRID(ng)%xu
            Yi => GRID(ng)%yu
            Xo => EXTRACT(ng)%xu
            Yo => EXTRACT(ng)%yu
          END IF
# ifdef MASKING
          mask_inp => GRID(ng)%umask
          mask_out => EXTRACT(ng)%umask

# endif
!
!       -----------------------------
        CASE (v2dvar, v3dvar)                    ! V-points
!       -----------------------------
!
          Cgrid=4
!
          Is_inp=IOBOUNDS(ng)%ILB_v
          Ie_inp=IOBOUNDS(ng)%IUB_v
          Js_inp=IOBOUNDS(ng)%JLB_v
          Je_inp=IOBOUNDS(ng)%JUB_v
!
          Is_out=xtr_IOBOUNDS(ng)%ILB_v
          Ie_out=xtr_IOBOUNDS(ng)%IUB_v
          Js_out=xtr_IOBOUNDS(ng)%JLB_v
          Je_out=xtr_IOBOUNDS(ng)%JUB_v
!
          IF (spherical) THEN
            Xi => GRID(ng)%lonv
            Yi => GRID(ng)%latv
            Xo => EXTRACT(ng)%lonv
            Yo => EXTRACT(ng)%latv
          ELSE
            Xi => GRID(ng)%xv
            Yi => GRID(ng)%yv
            Xo => EXTRACT(ng)%xv
            Yo => EXTRACT(ng)%yv
          END IF
# ifdef MASKING
          mask_inp => GRID(ng)%vmask
          mask_out => EXTRACT(ng)%vmask
# endif
!
      END SELECT
!
      Isize=Ie_inp-Is_inp+1
      Jsize=Je_inp-Js_inp+1
      Msize=Ie_out-Is_out+1
      Nsize=Je_out-Js_out+1
!
!  Set donor grid curvilinear rotation angle.
!
      IF (.not.allocated(angle)) THEN
        allocate ( angle(LBi_inp:UBi_inp,LBj_inp:UBj_inp) )
      END IF
!
      IF (Cgrid.eq.1) THEN                                ! PSI-points
        DO j=JstrP,JendT
          DO i=IstrP,IendT
            angle(i,j)=0.25_r8*(GRID(ng)%angler(i-1,j-1)+               &
     &                          GRID(ng)%angler(i-1,j  )+               &
     &                          GRID(ng)%angler(i  ,j-1)+               &
     &                          GRID(ng)%angler(i  ,j  ))
          END DO
        END DO
      ELSE IF (Cgrid.eq.2) THEN                           ! RHO-points
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            angle(i,j)=GRID(ng)%angler(i,j)
          END DO
        END DO
      ELSE IF (Cgrid.eq.3) THEN                           ! U-points
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            angle(i,j)=0.5_r8*(GRID(ng)%angler(i-1,j)+                  &
     &                         GRID(ng)%angler(i  ,j))
          END DO
        END DO
      ELSE IF (Cgrid.eq.4) THEN                           ! V-points
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            angle(i,j)=0.5_r8*(GRID(ng)%angler(i,j-1)+                  &
     &                         GRID(ng)%angler(i,j  ))
          END DO
        END DO
      END IF

# ifdef DISTRIBUTE
!
!  If distributed-memory, collect donor grid location data from all
!  spawned nodes and store it into a global scratch 1D array, packed
!  in column-major order to facilitate interpolation.

#  ifdef MASKING
!
      CALL mp_gather2d (ng, model,                                      &
     &                  LBi_inp, UBi_inp, LBj_inp, UBj_inp,             &
     &                  0, gtype, scale,                                &
     &                  mask_inp,                                       &
     &                  mask_inp, Npts, Minp, .FALSE.)
#  endif
!
      CALL mp_gather2d (ng, model,                                      &
     &                  LBi_inp, UBi_inp, LBj_inp, UBj_inp,             &
     &                  0, gtype, scale,                                &
#  ifdef MASKING
     &                  mask_inp,                                       &
#  endif
     &                  angle, Npts, Ainp, .FALSE.)
!
      CALL mp_gather2d (ng, model,                                      &
     &                  LBi_inp, UBi_inp, LBj_inp, UBj_inp,             &
     &                  0, gtype, scale,                                &
#  ifdef MASKING
     &                  mask_inp,                                       &
#  endif
     &                  Xi, Npts, Xinp, .FALSE.)
!
      CALL mp_gather2d (ng, model,                                      &
     &                  LBi_inp, UBi_inp, LBj_inp, UBj_inp,             &
     &                  0, gtype, scale,                                &
#  ifdef MASKING
     &                  mask_inp,                                       &
#  endif
     &                  Yi, Npts, Yinp, .FALSE.)

!
      CALL mp_gather2d (ng, model,                                      &
     &                  LBi_inp, UBi_inp, LBj_inp, UBj_inp,             &
     &                  0, gtype, scale,                                &
#  ifdef MASKING
     &                  mask_inp,                                       &
#  endif
     &                  Xi, Npts, Xinp, .FALSE.)
!
      CALL mp_gather2d (ng, model,                                      &
     &                  LBi_inp, UBi_inp, LBj_inp, UBj_inp,             &
     &                  0, gtype, scale,                                &
#  ifdef MASKING
     &                  mask_inp,                                       &
#  endif
     &                  Yi, Npts, Yinp, .FALSE.)

#  ifdef MASKING
!
      CALL mp_gather2d_xtr (ng, model,                                  &
     &                      LBi_out, UBi_out, LBj_out, UBj_out,         &
     &                      0, gtype, scale,                            &
     &                      mask_out,                                   &
     &                      mask_out, Npts, Mout, .FALSE.)
#  endif
!
      CALL mp_gather2d_xtr (ng, model,                                  &
     &                      LBi_out, UBi_out, LBj_out, UBj_out,         &
     &                      0, gtype, scale,                            &
#  ifdef MASKING
     &                      mask_out,                                   &
#  endif
     &                      Xo, Npts, Xout, .FALSE.)
!
      CALL mp_gather2d_xtr (ng, model,                                  &
     &                      LBi_out, UBi_out, LBj_out, UBj_out,         &
     &                      0, gtype, scale,                            &
#  ifdef MASKING
     &                      mask_out,                                   &
#  endif
     &                      Yo, Npts, Yout, .FALSE.)

# else

!
!  If serial or shared-memory applications, pack donor grid data into a
!  global 1D array in column-major order.
!
      Ainp=PACK(angle, .TRUE.)
#  ifdef MASKING
      Minp=PACK(mask_inp, .TRUE.)
#  endif
      Xinp=PACK(Xi, .TRUE.)
      Yinp=PACK(Yi, .TRUE.)
!
#  ifdef MASKING
      Mout=PACK(mask_out, .TRUE.)
#  endif
      Xout=PACK(Xo, .TRUE.)
      Yout=PACK(Yo, .TRUE.)
# endif

# ifdef DISTRIBUTE
!
!  Find the donor grid cell's fractional indices (Iout, Jout) containing
!  the extracted field points. Thus, computing the interpolation weights
!  is trivial. Here, each mpi-process calculates the global 1D-packed
!  Iout and Jout indices once to avoid collecting communications when
!  writing the extracted data.
!
      CALL hindices (ng,                                                &
     &               1, Isize, 1, Jsize,                                &
     &               1, Isize, 1, Jsize,                                &
     &               Ainp, Xinp, Yinp,                                  &
     &               1, Msize, 1, Nsize,                                &
     &               1, Msize, 1, Nsize,                                &
     &               Xout, Yout,                                        &
     &               Iout, Jout,                                        &
     &               0.0_r8, .FALSE.)
# else
!
!  Find the fractional indices (Iout,Jout) of the donor grid cells
!  containing extracted field points.
!
      CALL hindices (ng,                                                &
     &               1, Isize, 1, Jsize,                                &
     &               1, Isize, 1, Jsize,                                &
     &               Ainp, Xinp, Yinp,                                  &
     &               LBi_out, UBi_out, LBj_out, UBj_out,                &
     &               Is_out, Ie_out, Js_out, Je_out,                    &
     &               Xout, Yout,                                        &
     &               Iout, Jout,                                        &
     &               0.0_r8, .FALSE.)
# endif
!
!  Deallocate.
!
      IF (allocated(angle)) deallocate ( angle )
!
      RETURN
      END SUBROUTINE interp_coords
!
      SUBROUTINE interp_field2d (ng, model, tile,                       &
     &                           gtype, ifield, tindex,                 &
     &                           Imin, Imax, Jmin, Jmax,                &
     &                           start, total,                          &
     &                           Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It interpolates data from donor 2D fields, packed as a 1D global,   !
!  array, to specified extract grid geometry.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      logical  :: LandFill
!
      integer  :: ghost, i, j, ij, ic, jc, ifactor, method
      integer  :: Cgrid, Ilen, Jlen, Ioff, Joff, Isize, Jsize
      integer  :: Istr, Iend, Jstr, Jend
      integer  :: LBi, UBi, LBj, UBj
!
      real(r8) :: Fmin, Fmax

      real(r8) :: Fwrk(Npts)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_field2d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Interpolate
!------------------------------------------------------------------------
!
      LandFill=.TRUE.
      ghost=0
      method=Bilinear
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
!  At this point, the donor data is packed as global 1D arrays to
!  facilitate interpolating to the extracted grid tiled data.
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Cgrid=2
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize,                    &
     &                           GRID(ng) % Gx_rho,                     &
     &                           GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                           GRID(ng) % Gmask_rho,                  &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_rho,                &
     &                           EXTRACT(ng) % Jout_rho,                &
     &                           EXTRACT(ng) % lonr,                    &
     &                           EXTRACT(ng) % latr,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize,                    &
     &                           GRID(ng) % Gx_rho,                     &
     &                           GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                           GRID(ng) % Gmask_rho,                  &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_rho,                &
     &                           EXTRACT(ng) % Jout_rho,                &
     &                           EXTRACT(ng) % xr,                      &
     &                           EXTRACT(ng) % Yr,                      &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          Cgrid=3
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize,                    &
     &                           GRID(ng) % Gx_u,                       &
     &                           GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_u,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_u,                  &
     &                           EXTRACT(ng) % Jout_u,                  &
     &                           EXTRACT(ng) % lonu,                    &
     &                           EXTRACT(ng) % latu,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize,                    &
     &                           GRID(ng) % Gx_u,                       &
     &                           GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_u,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_u,                  &
     &                           EXTRACT(ng) % Jout_u,                  &
     &                           EXTRACT(ng) % xu,                      &
     &                           EXTRACT(ng) % yu,                      &
# ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          Cgrid=4
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize,                    &
     &                           GRID(ng) % Gx_v,                       &
     &                           GRID(ng) % Gy_v,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_v,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_v,                  &
     &                           EXTRACT(ng) % Jout_v,                  &
     &                           EXTRACT(ng) % lonv,                    &
     &                           EXTRACT(ng) % latv,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize,                    &
     &                           GRID(ng) % Gx_v,                       &
     &                           GRID(ng) % Gy_v,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_v,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_v,                  &
     &                           EXTRACT(ng) % Jout_v,                  &
     &                           EXTRACT(ng) % lonv,                    &
     &                           EXTRACT(ng) % latv,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for interpolation:',       &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=Jend-Jstr+1
        start(3)=tindex
        total(3)=1
# ifdef MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' INTERP_FIELD2D - Staggered variable, gtype = ', i0,     &
     &        /,18x,a,2x,a)
!
      RETURN
      END SUBROUTINE interp_field2d
!
      SUBROUTINE interp_field3d (ng, model, tile,                       &
     &                           gtype, ifield, tindex,                 &
     &                           Imin, Imax, Jmin, Jmax, Kmin, Kmax,    &
     &                           start, total,                          &
     &                           Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It interpolates data from donor 3D fields, packed as a 1D global,   !
!  array, to specified extract grid geometry.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Extract_Flag Extraction/decimation flag (integer)                !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Fwrk.                 !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile
      integer,  intent(in)  :: gtype, ifield, tindex
      integer,  intent(in)  :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)  :: Npts
      integer,  intent(out) :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      logical  :: LandFill
!
      integer  :: ghost, i, j, k, ij, ijk, ic, jc, ifactor, method
      integer  :: Ilen, Jlen, Ioff, Joff, Koff
      integer  :: Cgrid, Isize, IJsize, Jsize, Ksize
      integer  :: Istr, Iend, Jstr, Jend
      integer  :: LBi, UBi, LBj, UBj
!
      real(r8) :: Fmin, Fmax

      real(r8) :: Fwrk(Npts)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", interp_field3d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Interpolate
!------------------------------------------------------------------------
!
      LandFill=.TRUE.
      ghost=0
      method=Bilinear
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      Ksize=Kmax-Kmin+1
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r3dvar, w3dvar)               ! RHO and W points field
!       ---------------------
!
          Cgrid=2
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_rho,                     &
     &                           GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                           GRID(ng) % Gmask_rho,                  &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_rho,                &
     &                           EXTRACT(ng) % Jout_rho,                &
     &                           EXTRACT(ng) % lonr,                    &
     &                           EXTRACT(ng) % latr,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_rho,                     &
     &                           GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                           GRID(ng) % Gmask_rho,                  &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_rho,                &
     &                           EXTRACT(ng) % Jout_rho,                &
     &                           EXTRACT(ng) % xr,                      &
     &                           EXTRACT(ng) % Yr,                      &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ---------------------
        CASE (u3dvar)                       ! U points field
!       ---------------------
!
          Cgrid=3
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_u,                       &
     &                           GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_u,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_u,                  &
     &                           EXTRACT(ng) % Jout_u,                  &
     &                           EXTRACT(ng) % lonu,                    &
     &                           EXTRACT(ng) % latu,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_u,                       &
     &                           GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_u,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_u,                  &
     &                           EXTRACT(ng) % Jout_u,                  &
     &                           EXTRACT(ng) % xu,                      &
     &                           EXTRACT(ng) % yu,                      &
# ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ---------------------
        CASE (v3dvar)                       ! V points field
!       ---------------------
!
          Cgrid=4
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_v,                       &
     &                           GRID(ng) % Gy_v,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_v,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_v,                  &
     &                           EXTRACT(ng) % Jout_v,                  &
     &                           EXTRACT(ng) % lonv,                    &
     &                           EXTRACT(ng) % latv,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_v,                       &
     &                           GRID(ng) % Gy_v,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_v,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_v,                  &
     &                           EXTRACT(ng) % Jout_v,                  &
     &                           EXTRACT(ng) % lonv,                    &
     &                           EXTRACT(ng) % latv,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for interpolation:',       &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=Jend-Jstr+1
        start(3)=1
        total(3)=Ksize
        start(4)=tindex
        total(4)=1
# ifdef MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' INTERP_FIELD3D - Staggered variable, gtype = ', i0,     &
     &        /,18x,a,2x,a)
!
      RETURN
      END SUBROUTINE interp_field3d
!
      SUBROUTINE interp_field4d (ng, model, tile,                       &
     &                           gtype, ifield, tindex,                 &
     &                           Imin, Imax, Jmin, Jmax, Kmin, Kmax,    &
     &                           Fourth, Loff,                          &
     &                           start, total,                          &
     &                           Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It interpolates data from donor 2D fields, packed as a 1D global,   !
!  array, to specified extract grid geometry. The data is processed    !
!  by 3D slices to reduce memory requirements in the calling routine.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     fourth       Donor Fouth dimension index to process (integer)    !
!     Loff         Fourth dimension couter offset (integer)            !
!     Fdat         Packed global donor 3D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 3D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)    :: fourth, Loff
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      logical  :: LandFill
!
      integer  :: ghost, i, j, k, ij, ijk, ic, jc, ifactor, method
      integer  :: Cgrid, Ilen, Jlen, Ioff, Joff, Koff
      integer  :: Isize, IJsize, Jsize, Ksize
      integer  :: Istr, Iend, Jstr, Jend
      integer  :: LBi, UBi, LBj, UBj
!
      real(r8) :: Fmin, Fmax

      real(r8) :: Fwrk(Npts)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", interp_field4d"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Interpolate.
!------------------------------------------------------------------------
!
      LandFill=.TRUE.
      ghost=0
      method=Bilinear
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r3dvar, w3dvar)               ! RHO and W points field
!       ---------------------
!
          Cgrid=2
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_rho,                     &
     &                           GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                           GRID(ng) % Gmask_rho,                  &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_rho,                &
     &                           EXTRACT(ng) % Jout_rho,                &
     &                           EXTRACT(ng) % lonr,                    &
     &                           EXTRACT(ng) % latr,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_rho,                     &
     &                           GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                           GRID(ng) % Gmask_rho,                  &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_rho,                &
     &                           EXTRACT(ng) % Jout_rho,                &
     &                           EXTRACT(ng) % xr,                      &
     &                           EXTRACT(ng) % Yr,                      &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ---------------------
        CASE (u3dvar)                       ! U points field
!       ---------------------
!
          Cgrid=3
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_u,                       &
     &                           GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_u,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_u,                  &
     &                           EXTRACT(ng) % Jout_u,                  &
     &                           EXTRACT(ng) % lonu,                    &
     &                           EXTRACT(ng) % latu,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_u,                       &
     &                           GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_u,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_u,                  &
     &                           EXTRACT(ng) % Jout_u,                  &
     &                           EXTRACT(ng) % xu,                      &
     &                           EXTRACT(ng) % yu,                      &
# ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ---------------------
        CASE (v3dvar)                       ! V points field
!       ---------------------
!
          Cgrid=4
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          IF (spherical) THEN
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_v,                       &
     &                           GRID(ng) % Gy_v,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_v,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_v,                  &
     &                           EXTRACT(ng) % Jout_v,                  &
     &                           EXTRACT(ng) % lonv,                    &
     &                           EXTRACT(ng) % latv,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          ELSE
            CALL regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           1, Isize, 1, Jsize, Kmin, Kmax,        &
     &                           GRID(ng) % Gx_v,                       &
     &                           GRID(ng) % Gy_v,                       &
# ifdef MASKING
     &                           GRID(ng) % Gmask_v,                    &
# endif
     &                           Fwrk,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           EXTRACT(ng) % Iout_v,                  &
     &                           EXTRACT(ng) % Jout_v,                  &
     &                           EXTRACT(ng) % lonv,                    &
     &                           EXTRACT(ng) % latv,                    &
# ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
# endif
     &                           Npts, Fdat, Fmin, Fmax)
          END IF
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for interpolation:',       &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=Jend-Jstr+1
        start(3)=1
        total(3)=Ksize
        start(4)=fourth+Loff
        total(4)=1
        start(5)=tindex
        total(5)=1
# ifdef MASKING
      ELSE
        start(1)=1+(fourth+Loff-1)*Npts
        total(1)=Npts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' INTERP_FIELD4D - Staggered variable, gtype = ', i0,     &
     &        /,18x,a,2x,a)
!
      RETURN
      END SUBROUTINE interp_field4d
!
      SUBROUTINE interp_field2d_global (ng, model, tile,                &
     &                                  gtype, ifield, tindex,          &
     &                                  Imin, Imax, Jmin, Jmax,         &
     &                                  start, total,                   &
     &                                  Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It interpolates data from donor 2D fields, packed as a 1D global,   !
!  array, to specified extract grid geometry.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      logical  :: LandFill
!
      integer  :: ghost, i, j, ij, ic, jc, ifactor, method
      integer  :: Cgrid, Ilen, Jlen, Ioff, Joff, Isize, Jsize
      integer  :: Istr, Iend, Jstr, Jend
      integer  :: LBi, UBi, LBj, UBj
!
      real(r8) :: Fmin, Fmax

      real(r8) :: Fwrk(Npts)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", extract_field2d_global"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Interpolate
!------------------------------------------------------------------------
!
      LandFill=.TRUE.
      ghost=0
      method=Bilinear
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
!  At this point, the donor data is packed as global 1D arrays to
!  facilitate interpolating to the extracted grid tiled data.
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r2dvar, r3dvar)               ! RHO points field
!       ---------------------
!
          Cgrid=2
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                         method, LandFill,                      &
     &                         Imin, Imax, Jmin, Jmax,                &
     &                         GRID(ng) % Gx_rho,                     &
     &                         GRID(ng) % Gy_rho,                     &
# ifdef MASKING
     &                         GRID(ng) % Gmask_rho,                  &
# endif
     &                         Fwrk,                                  &
     &                         Istr, Iend, Jstr, Jend,                &
     &                         Istr, Iend, Jstr, Jend,                &
     &                         EXTRACT(ng) % Iout_rho,                &
     &                         EXTRACT(ng) % Jout_rho,                &
     &                         EXTRACT(ng) % Gx_rho,                  &
     &                         EXTRACT(ng) % Gy_rho,                  &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_rho,               &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ---------------------
        CASE (u2dvar, u3dvar)               ! U points field
!       ---------------------
!
          Cgrid=3
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                         method, LandFill,                      &
     &                         Imin, Imax, Jmin, Jmax,                &
     &                         GRID(ng) % Gx_u,                       &
     &                         GRID(ng) % Gy_u,                       &
# ifdef MASKING
     &                         GRID(ng) % Gmask_u,                    &
# endif
     &                         Fwrk,                                  &
     &                         Istr, Iend, Jstr, Jend,                &
     &                         Istr, Iend, Jstr, Jend,                &
     &                         EXTRACT(ng) % Iout_u,                  &
     &                         EXTRACT(ng) % Jout_u,                  &
     &                         EXTRACT(ng) % Gx_u,                    &
     &                         EXTRACT(ng) % Gx_u,                    &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_u,                 &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ---------------------
        CASE (v2dvar, v3dvar)               ! V points field
!       ---------------------
!
          Cgrid=4
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field2d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         Imin, Imax, Jmin, Jmax,                  &
     &                         GRID(ng) % Gx_v,                         &
     &                         GRID(ng) % Gy_v,                         &
# ifdef MASKING
     &                         GRID(ng) % Gmask_v,                      &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_v,                    &
     &                         EXTRACT(ng) % Jout_v,                    &
     &                         EXTRACT(ng) % Gx_v,                      &
     &                         EXTRACT(ng) % Gy_v,                      &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_v,                   &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for interpolation:',       &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=Jend-Jstr+1
        start(3)=tindex
        total(3)=1
# ifdef MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' INTERP_FIELD2D_GLOBAL - Staggered variable, gtype = ',  &
     &        i0,/,18x,a,2x,a)
!
      RETURN
      END SUBROUTINE interp_field2d_global
!
      SUBROUTINE interp_field3d_global (ng, model, tile,                &
     &                                  gtype, ifield, tindex,          &
     &                                  Imin, Imax, Jmin, Jmax,         &
     &                                  Kmin, Kmax,                     &
     &                                  start, total,                   &
     &                                  Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It interpolates data from donor 3D fields, packed as a 1D global,   !
!  array, to specified extract grid geometry.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     Fdat         Packed global donor 2D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Npts         Number of points processed in Fwrk.                 !
!     Fdat         Extracted 2D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile
      integer,  intent(in)  :: gtype, ifield, tindex
      integer,  intent(in)  :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)  :: Npts
      integer,  intent(out) :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      logical  :: LandFill
!
      integer  :: ghost, i, j, k, ij, ijk, ic, jc, ifactor, method
      integer  :: Ilen, Jlen, Ioff, Joff, Koff
      integer  :: Cgrid, Isize, IJsize, Jsize, Ksize
      integer  :: Istr, Iend, Jstr, Jend
      integer  :: LBi, UBi, LBj, UBj
!
      real(r8) :: Fmin, Fmax

      real(r8) :: Fwrk(Npts)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", interp_field3d_global"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Interpolate
!------------------------------------------------------------------------
!
      LandFill=.TRUE.
      ghost=0
      method=Bilinear
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
      Ksize=Kmax-Kmin+1
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r3dvar, w3dvar)               ! RHO and W points field
!       ---------------------
!
          Cgrid=2
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field3d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         Imin, Imax, Jmin, Jmax, Kmin, Kmax,      &
     &                         GRID(ng) % Gx_rho,                       &
     &                         GRID(ng) % Gy_rho,                       &
# ifdef MASKING
     &                         GRID(ng) % Gmask_rho,                    &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_rho,                  &
     &                         EXTRACT(ng) % Jout_rho,                  &
     &                         EXTRACT(ng) % Gx_rho,                    &
     &                         EXTRACT(ng) % Gy_rho,                    &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_rho,                 &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ---------------------
        CASE (u3dvar)                       ! U points field
!       ---------------------
!
          Cgrid=3
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field3d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         Imin, Imax, Jmin, Jmax, Kmin, Kmax,      &
     &                         GRID(ng) % Gx_u,                         &
     &                         GRID(ng) % Gy_u,                         &
# ifdef MASKING
     &                         GRID(ng) % Gmask_u,                      &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_u,                    &
     &                         EXTRACT(ng) % Jout_u,                    &
     &                         EXTRACT(ng) % Gx_u,                      &
     &                         EXTRACT(ng) % Gy_u,                      &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_u,                   &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ---------------------
        CASE (v3dvar)                       ! V points field
!       ---------------------
!
          Cgrid=4
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field3d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         Imin, Imax, Jmin, Jmax, Kmin, Kmax,      &
     &                         GRID(ng) % Gx_v,                         &
     &                         GRID(ng) % Gy_v,                         &
# ifdef MASKING
     &                         GRID(ng) % Gmask_v,                      &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_v,                    &
     &                         EXTRACT(ng) % Jout_v,                    &
     &                         EXTRACT(ng) % Gx_v,                      &
     &                         EXTRACT(ng) % Gy_v,                      &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_v,                   &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for interpolation:',       &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=Jend-Jstr+1
        start(3)=1
        total(3)=Ksize
        start(4)=tindex
        total(4)=1
# ifdef MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' INTERP_FIELD3D_GLOBAL - Staggered variable, gtype = ',  &
     &        i0,/,18x,a,2x,a)
!
      RETURN
      END SUBROUTINE interp_field3d_global
!
      SUBROUTINE interp_field4d_global (ng, model, tile,                &
     &                                  gtype, ifield, tindex,          &
     &                                  Imin, Imax, Jmin, Jmax,         &
     &                                  Kmin, Kmax, Fourth, Loff,       &
     &                                  start, total,                   &
     &                                  Npts, Fdat)
!
!=======================================================================
!                                                                      !
!  It interpolates data from donor 2D fields, packed as a 1D global,   !
!  array, to specified extract grid geometry. The data is processed    !
!  by 3D slices to reduce memory requirements in the calling routine.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     tile         Domain partition (integer)                          !
!     gtype        Staggered C-grid type (integer)                     !
!     ifield       Field metadata index (integer)                      !
!     tindex       Time record index to process (integer)              !
!     Imin         Donor field starting data I-index (integer)         !
!     Imax         Donor field ending   data I-index (integer)         !
!     Jmin         Donor field starting data J-index (integer)         !
!     Jmax         Donor field ending   data J-index (integer)         !
!     Kmin         Donor field starting data K-index (integer)         !
!     Kmax         Donor field ending   data K-index (integer)         !
!     fourth       Donor Fouth dimension index to process (integer)    !
!     Loff         Fourth dimension couter offset (integer)            !
!     Fdat         Packed global donor 3D field data (real 1D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     start        Start index where the first of the data values will !
!                    be written along each dimension (integer)         !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     Fdat         Extracted 3D field data (real)                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in)    :: ng, model, tile
      integer,  intent(in)    :: gtype, ifield, tindex
      integer,  intent(in)    :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer,  intent(in)    :: fourth, Loff
      integer,  intent(in)    :: Npts
      integer,  intent(out)   :: start(:), total(:)
!
      real(r8), intent(inout) :: Fdat(:)
!
!  Local variable declarations.
!
      logical  :: LandFill
!
      integer  :: ghost, i, j, k, ij, ijk, ic, jc, ifactor, method
      integer  :: Cgrid, Ilen, Jlen, Ioff, Joff, Koff
      integer  :: Isize, IJsize, Jsize, Ksize
      integer  :: Istr, Iend, Jstr, Jend
      integer  :: LBi, UBi, LBj, UBj
!
      real(r8) :: Fmin, Fmax

      real(r8) :: Fwrk(Npts)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", interp_field4d_global"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Interpolate.
!------------------------------------------------------------------------
!
      LandFill=.TRUE.
      ghost=0
      method=Bilinear
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
!
      LBi=xtr_BOUNDS(ng)%LBi(tile)
      UBi=xtr_BOUNDS(ng)%UBi(tile)
      LBj=xtr_BOUNDS(ng)%LBj(tile)
      UBj=xtr_BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (gtype)
!
!       ---------------------
        CASE (r3dvar, w3dvar)               ! RHO and W points field
!       ---------------------
!
          Cgrid=2
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field3d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         Imin, Imax, Jmin, Jmax, Kmin, Kmax,      &
     &                         GRID(ng) % Gx_rho,                       &
     &                         GRID(ng) % Gy_rho,                       &
# ifdef MASKING
     &                         GRID(ng) % Gmask_rho,                    &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_rho,                  &
     &                         EXTRACT(ng) % Jout_rho,                  &
     &                         EXTRACT(ng) % Gx_rho,                    &
     &                         EXTRACT(ng) % Gy_rho,                    &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_rho,                 &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ---------------------
        CASE (u3dvar)                       ! U points field
!       ---------------------
!
          Cgrid=3
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field3d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         Imin, Imax, Jmin, Jmax, Kmin, Kmax,      &
     &                         GRID(ng) % Gx_u,                         &
     &                         GRID(ng) % Gy_u,                         &
# ifdef MASKING
     &                         GRID(ng) % Gmask_u,                      &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_u,                    &
     &                         EXTRACT(ng) % Jout_u,                    &
     &                         EXTRACT(ng) % Gx_u,                      &
     &                         EXTRACT(ng) % Gy_u,                      &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_u,                   &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ---------------------
        CASE (v3dvar)                       ! V points field
!       ---------------------
!
          Cgrid=4
          Istr=xtr_BOUNDS(ng)%Imin(Cgrid,ghost,tile)
          Iend=xtr_BOUNDS(ng)%Imax(Cgrid,ghost,tile)
          Jstr=xtr_BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
          Jend=xtr_BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
!
          Fwrk=Fdat(1:Npts)
!
          CALL regrid_field3d (ng, model, tile, gtype, ifield,          &
     &                         method, LandFill,                        &
     &                         1, Isize, 1, Jsize, Kmin, Kmax,          &
     &                         GRID(ng) % Gx_v,                         &
     &                         GRID(ng) % Gy_v,                         &
# ifdef MASKING
     &                         GRID(ng) % Gmask_v,                      &
# endif
     &                         Fwrk,                                    &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         EXTRACT(ng) % Iout_v,                    &
     &                         EXTRACT(ng) % Jout_v,                    &
     &                         EXTRACT(ng) % Gx_v,                      &
     &                         EXTRACT(ng) % Gy_v,                      &
# ifdef MASKING
     &                         EXTRACT(ng) % Gmask_v,                   &
# endif
     &                         Npts, Fdat, Fmin, Fmax)
!
!       ------------
        CASE DEFAULT
!       ------------
!
          IF (Master) THEN
            WRITE (stdout,10) gtype,                                    &
     &                        'not supported for interpolation:',       &
     &                        TRIM(Vname(1,ifield))
          END IF
          exit_flag=3
!
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set start and total vectors needed to write into output NetCDF file.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=Jend-Jstr+1
        start(3)=1
        total(3)=Ksize
        start(4)=fourth+Loff
        total(4)=1
        start(5)=tindex
        total(5)=1
# ifdef MASKING
      ELSE
        start(1)=1+(fourth+Loff-1)*Npts
        total(1)=Npts
        start(2)=tindex
        total(2)=1
# endif
      END IF
!
  10  FORMAT (' INTERP_FIELD4D_GLOBAL - Staggered variable, gtype = ',  &
     &        i0,/,18x,a,2x,a)
!
      RETURN
      END SUBROUTINE interp_field4d_global
!
      SUBROUTINE regrid_field2d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           LBx, UBx, LBy, UBy,                    &
     &                           Xinp, Yinp,                            &
# ifdef MASKING
     &                           maskInp,                               &
# endif
     &                           Finp,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Iout, Jout,                            &
     &                           Xout, Yout,                            &
# ifdef MASKING
     &                           maskOut,                               &
# endif
     &                           Npts, Fout, Fmin, Fmax)
!
!=======================================================================
!                                                                      !
!  It extracts 2D fields from donor grid by interpolation.             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number (integer)                         !
!     model       Calling model identifier (integer)                   !
!     tile        Domain partition (integer)                           !
!     gtype       Staggered C-grid type (integer)                      !
!     ifield      Field metadata index (integer)                       !
!     method      Interpolation method (integer)                       !
!     LandFill    Switch to set fill value on land points (logical)    !
!     LBx         Donor data I-dimension lower bound (integer)         !
!     UBx         Donor data I-dimension upper bound (integer)         !
!     LBy         Donor data J-dimension lower bound (integer)         !
!     UBy         Donor data J-dimension upper bound (integer)         !
!     Xinp        Donor data X-locations (real array)                  !
!     Yinp        Donor data Y-locations (real array)                  !
# ifdef MASKING
!     maskInp     Donor data land/sea masking (real array)             !
# endif
!     Finp        Donor 2D field to interpolate from (real array)      !
!     Istr        Interpolsted field Starting tile I-index (integer)   !
!     Iend        Interpolated field ending   tile I-index (integer)   !
!     Jstr        Interpolated field starting tile J-index (integer)   !
!     Jend        Interpolated field ending   tile J-index (integer)   !
!     LBi         Interpolated field I-dimension Lower bound (integer) !
!     UBi         Interpolated field I-dimension Upper bound (integer) !
!     LBj         Interpolated field J-dimension Lower bound (integer) !
!     UBj         Interpolated field J-dimension Upper bound (integer) !
!     Iout        I-fractional Xinp grid cell containing Xout (real)   !
!     Jout        J-fractional Yinp grid cell containing Yout (real)   !
!     Xout        Interpolated field X-locations (real array)          !
!     Yout        Interpolated field Y-locations (real array)          !
# ifdef MASKING
!     maskOut     Interpolated field land/sea masking (real array)     !
# endif
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts         Number of points processed in Fout (integer)        !
!     Fout         Interpolated 2D field data (real array)             !
!     Fmin         Interpolated field minimum value (real)             !
!     Fmax         Interpolated field maximum value (real)             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,  intent(in)   :: LandFill
!
      integer,  intent(in)   :: ng, model, tile, gtype, ifield, method
      integer,  intent(in)   :: LBx, UBx, LBy, UBy
      integer,  intent(in)   :: LBi, UBi, LBj, UBj
      integer,  intent(in)   :: Istr, Iend, Jstr, Jend
      integer,  intent(in)   :: Npts
!
      real(r8), intent(in)   :: Xinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in)   :: Yinp(LBx:UBx,LBy:UBy)
# ifdef MASKING
      real(r8), intent(in)   :: maskInp(LBx:UBx,LBy:UBy)
# endif
      real(r8), intent(in)   :: Finp(LBx:UBx,LBy:UBy)
      real(r8), intent(in)   :: Iout(LBi:UBi,LBj:UBj)
      real(r8), intent(in)   :: Jout(LBi:UBi,LBj:UBj)
      real(r8), intent(in)   :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in)   :: Yout(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in)   :: maskOut(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out)  :: Fout(LBi:UBi,LBj:UBj)
      real(r8), intent(out)  :: Fmin, Fmax
!
!  Local variable declarations.
!
      integer :: i, ic, j
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", regrid_field2d"
!
!-----------------------------------------------------------------------
!  Interpolate 2D field.
!-----------------------------------------------------------------------
!
      SELECT CASE (ABS(method))
        CASE (Bilinear)
          CALL linterp2d (ng,                                           &
     &                    LBx, UBx, LBy, UBy,                           &
     &                    Xinp, Yinp, Finp,                             &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Istr, Iend, Jstr, Jend,                       &
     &                    Iout, Jout,                                   &
     &                    Xout, Yout,                                   &
     &                    Fout)
        CASE (Bicubic)
          CALL cinterp2d (ng,                                           &
     &                    LBx, UBx, LBy, UBy,                           &
     &                    Xinp, Yinp, Finp,                             &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Istr, Iend, Jstr, Jend,                       &
     &                    Iout, Jout,                                   &
     &                    Xout, Yout,                                   &
     &                    Fout)
        CASE DEFAULT
          IF (Master) WRITE(stdout,10) method
          exit_flag=5
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  If requested, set fill value in land areas.
!-----------------------------------------------------------------------
!
      Fmin=spval
      Fmax=-spval
      ic=0
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ic=ic+1
# ifdef MASKING
          IF ((maskOut(i,j).lt.1.0_r8).and.LandFill) THEN
            Fout(i,j)=spval
          ELSE
            Fmin=MIN(Fmin,Fout(i,j))
            Fmax=MAX(Fmax,Fout(i,j))
          END IF
# else
          Fmin=MIN(Fmin,Fout(i,j))
          Fmax=MAX(Fmax,Fout(i,j))
# endif
        END DO
      END DO
!
  10  FORMAT (' REGRID_FIEL2D - Illegal interpolation method =', i0)
!
      RETURN
      END SUBROUTINE regrid_field2d
!
      SUBROUTINE regrid_field3d (ng, model, tile, gtype, ifield,        &
     &                           method, LandFill,                      &
     &                           LBx, UBx, LBy, UBy, Kmin, Kmax,        &
     &                           Xinp, Yinp,                            &
# ifdef MASKING
     &                           maskInp,                               &
# endif
     &                           Finp,                                  &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Iout, Jout,                            &
     &                           Xout, Yout,                            &
# ifdef MASKING
     &                           maskOut,                               &
# endif
     &                           Npts, Fout, Fmin, Fmax)
!
!=======================================================================
!                                                                      !
!  It extracts 3D fields from donor grid by interpolation level-by-    !
!  level.                                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number (integer)                         !
!     model       Calling model identifier (integer)                   !
!     tile        Domain partition (integer)                           !
!     gtype       Staggered C-grid type (integer)                      !
!     ifield      Field metadata index (integer)                       !
!     method      Interpolation method (integer)                       !
!     LandFill    Switch to set fill value on land points (logical)    !
!     LBx         Donor data I-dimension lower bound (integer)         !
!     UBx         Donor data I-dimension upper bound (integer)         !
!     LBy         Donor data J-dimension lower bound (integer)         !
!     UBy         Donor data J-dimension upper bound (integer)         !
!     Kmin        Donor data K-dimension lower bound (integer)         !
!     Kmax        Donor data K-dimension upper bound (integer)         !
!     Xinp        Donor data X-locations (real array)                  !
!     Yinp        Donor data Y-locations (real array)                  !
# ifdef MASKING
!     maskInp     Donor data land/sea masking (real array)             !
# endif
!     Finp        Donor 3D field to interpolate from (real array)      !
!     Istr        Interpolsted field Starting tile I-index (integer)   !
!     Iend        Interpolated field ending   tile I-index (integer)   !
!     Jstr        Interpolated field starting tile J-index (integer)   !
!     Jend        Interpolated field ending   tile J-index (integer)   !
!     LBi         Interpolated field I-dimension Lower bound (integer) !
!     UBi         Interpolated field I-dimension Upper bound (integer) !
!     LBj         Interpolated field J-dimension Lower bound (integer) !
!     UBj         Interpolated field J-dimension Upper bound (integer) !
!     Iout        I-fractional Xinp grid cell containing Xout (real)   !
!     Jout        J-fractional Yinp grid cell containing Yout (real)   !
!     Xout        Interpolated field X-locations (real array)          !
!     Yout        Interpolated field Y-locations (real array)          !
# ifdef MASKING
!     maskOut     Interpolated field land/sea masking (real array)     !
# endif
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts         Number of points processed in Fout (integer)        !
!     Fout         Interpolated 3D field data (real array)             !
!     Fmin         Interpolated field minimum value (real)             !
!     Fmax         Interpolated field maximum value (real)             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,  intent(in)  :: landFill
!
      integer,  intent(in)  :: ng, model, tile, gtype, ifield, method
      integer,  intent(in)  :: LBx, UBx, LBy, UBy
      integer,  intent(in)  :: Kmin, Kmax
      integer,  intent(in)  :: Istr, Iend, Jstr, Jend
      integer,  intent(in)  :: LBi, UBi, LBj, UBj
      integer,  intent(in)  :: Npts
!
      real(r8), intent(in)  :: Xinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in)  :: Yinp(LBx:UBx,LBy:UBy)
# ifdef MASKING
      real(r8), intent(in)  :: maskInp(LBx:UBx,LBy:UBy)
# endif
      real(r8), intent(in)  :: Finp(LBx:UBx,LBy:UBy,Kmin:Kmax)
      real(r8), intent(in)  :: Iout(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: Jout(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: Yout(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in)  :: maskOut(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj,Kmin:Kmax)
      real(r8), intent(out) :: Fmin, Fmax
!
!  Local variable declarations.
!
      integer :: i, ic, j, k
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", regrid_field3d"
!
!-----------------------------------------------------------------------
!  Interpolate 2D field.
!-----------------------------------------------------------------------
!
      SELECT CASE (ABS(method))
        CASE (Bilinear)
          DO k=Kmin,Kmax
            CALL linterp2d (ng,                                         &
     &                      LBx, UBx, LBy, UBy,                         &
     &                      Xinp, Yinp, Finp(:,:,k),                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      Iout, Jout,                                 &
     &                      Xout, Yout,                                 &
     &                      Fout(LBi:,LBj:,k))
          END DO
        CASE (Bicubic)
          DO k=Kmin,Kmax
            CALL cinterp2d (ng,                                         &
     &                      LBx, UBx, LBy, UBy,                         &
     &                      Xinp, Yinp, Finp(:,:,k),                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      Iout, Jout,                                 &
     &                      Xout, Yout,                                 &
     &                      Fout(LBi:,LBj:,k))
          END DO
        CASE DEFAULT
          IF (Master) WRITE(stdout,10) method
          exit_flag=5
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  If requested, set fill value in land areas.
!-----------------------------------------------------------------------
!
      Fmin=spval
      Fmax=-spval
      ic=0
      DO k=Kmin,Kmax
        DO j=Jstr,Jend
          DO i=Istr,Iend
            ic=ic+1
# ifdef MASKING
            IF ((maskOut(i,j).lt.1.0_r8).and.LandFill) THEN
              Fout(i,j,k)=spval
            ELSE
              Fmin=MIN(Fmin,Fout(i,j,k))
              Fmax=MAX(Fmax,Fout(i,j,k))
            END IF
# else
            Fmin=MIN(Fmin,Fout(i,j,k))
            Fmax=MAX(Fmax,Fout(i,j,k))
# endif
          END DO
        END DO
      END DO
!
  10  FORMAT (' REGRID_FIEL3D - Illegal interpolation method =', i0)
!
      RETURN
      END SUBROUTINE regrid_field3d
!
#endif
      END MODULE extract_field_mod
