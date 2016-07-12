#include "cppdefs.h"
      MODULE nf_fwrite2d_mod
!
!svn $Id: nf_fwrite2d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 2D array into an  !
!  output NetCDF file.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvarid      NetCDF variable ID.                                 !
!     tindex       NetCDF time record index to write.                  !
!     gtype        Grid type. If negative, only write water points.    !
!     LBi          I-dimension Lower bound.                            !
!     UBi          I-dimension Upper bound.                            !
!     LBj          J-dimension Lower bound.                            !
!     UBj          J-dimension Upper bound.                            !
!     Amask        land/Sea mask, if any (real).                       !
!     Ascl         Factor to scale field before writing (real).        !
!     A            Field to write out (real).                          !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional).                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite  Error flag (integer).                                 !
!                                                                      !
#ifdef POSITIVE_ZERO
!  Starting F95 zero values can be signed (-0 or +0) following the     !
!  IEEE 754 floating point standard.  This may produce different       !
!  output data in serial and parallel applications. Since comparing    !
!  serial and parallel output is essential for tracking parallel       !
!  partition  bugs, "positive zero" is enforced.                       !
!                                                                      !
#endif
!=======================================================================
!
      implicit none

      CONTAINS

#if defined PARALLEL_OUT && defined DISTRIBUTE
!
!***********************************************************************
      FUNCTION nf_fwrite2d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, Ascl,                   &
# ifdef MASKING
     &                      Amask,                                      &
# endif
     &                      A, SetFillVal)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# if defined WRITE_WATER && defined MASKING
!
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj

      real(r8), intent(in) :: Ascl

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
# endif
      integer :: i, ic, j, jc, Npts
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ioff, Joff
      integer :: Istr, Iend
      integer :: Ilen, Isize, Jlen, Jsize, IJlen
      integer :: MyType, status

      integer, dimension(3) :: start, total

      integer :: nf_fwrite2d

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), allocatable :: Awrk(:)
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification.
!
      MyType=gtype

      SELECT CASE (ABS(MyType))
        CASE (p2dvar, p3dvar)
          Imin=BOUNDS(ng)%Istr (MyRank)
          Imax=BOUNDS(ng)%IendR(MyRank)
          Jmin=BOUNDS(ng)%Jstr (MyRank)
          Jmax=BOUNDS(ng)%JendR(MyRank)
          Isize=IOBOUNDS(ng)%xi_psi
          Jsize=IOBOUNDS(ng)%eta_psi
        CASE (r2dvar, r3dvar)
          Imin=BOUNDS(ng)%IstrR(MyRank)
          Imax=BOUNDS(ng)%IendR(MyRank)
          Jmin=BOUNDS(ng)%JstrR(MyRank)
          Jmax=BOUNDS(ng)%JendR(MyRank)
          Isize=IOBOUNDS(ng)%xi_rho
          Jsize=IOBOUNDS(ng)%eta_rho
        CASE (u2dvar, u3dvar)
          Imin=BOUNDS(ng)%Istr (MyRank)
          Imax=BOUNDS(ng)%IendR(MyRank)
          Jmin=BOUNDS(ng)%JstrR(MyRank)
          Jmax=BOUNDS(ng)%JendR(MyRank)
          Isize=IOBOUNDS(ng)%xi_u
          Jsize=IOBOUNDS(ng)%eta_u
        CASE (v2dvar, v3dvar)
          Imin=BOUNDS(ng)%IstrR(MyRank)
          Imax=BOUNDS(ng)%IendR(MyRank)
          Jmin=BOUNDS(ng)%Jstr (MyRank)
          Jmax=BOUNDS(ng)%JendR(MyRank)
          Isize=IOBOUNDS(ng)%xi_v
          Jsize=IOBOUNDS(ng)%eta_v
        CASE DEFAULT
          Imin=BOUNDS(ng)%IstrR(MyRank)
          Imax=BOUNDS(ng)%IendR(MyRank)
          Jmin=BOUNDS(ng)%JstrR(MyRank)
          Jmax=BOUNDS(ng)%JendR(MyRank)
          Isize=IOBOUNDS(ng)%xi_rho
          Jsize=IOBOUNDS(ng)%eta_rho
      END SELECT

      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IJlen=Ilen*Jlen

# ifdef MASKING
!
!  Set switch to replace land areas with fill value, spval.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Parallel I/O: Pack tile data into 1D array in column-major order.
# ifdef MASKING
!                Overwrite masked points with special value.
# endif
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
!
!  Set offsets due the NetCDF dimensions. Recall that some output
!  variables not always start at one.
!
        SELECT CASE (ABS(MyType))
          CASE (p2dvar, p3dvar)
            Ioff=0
            Joff=0
          CASE (r2dvar, r3dvar)
            Ioff=1
            Joff=1
          CASE (u2dvar, u3dvar)
            Ioff=0
            Joff=1
          CASE (v2dvar, v3dvar)
            Ioff=1
            Joff=0
          CASE DEFAULT
            Ioff=1
            Joff=1
        END SELECT
!
!  Allocate and initialize scratch work array.
!
        Npts=IJlen
        IF (.not.allocated(Awrk)) THEN
          allocate ( Awrk(Npts) )
          Awrk=IniVal
        END IF
!
!  Pack and scale tile data.
!
        ic=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ic=ic+1
            Awrk(ic)=A(i,j)*Ascl
#ifdef POSITIVE_ZERO
            IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
              Awrk(ic)=0.0_r8                     ! impose positive zero
            END IF
#endif
# ifdef MASKING
            IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk(ic)=spval
            END IF
# endif
          END DO
        END DO
!
!  Write out data: all parallel nodes write their own packed tile data.
!
        start(1)=Imin+Ioff
        total(1)=Ilen
        start(2)=Jmin+Joff
        total(2)=Jlen
        start(3)=tindex
        total(3)=1

        status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
        nf_fwrite2d=status
      END IF

# if defined WRITE_WATER && defined MASKING
!
!-----------------------------------------------------------------------
!  Parallel I/O: Remove land points and pack tile data into 1D array.
!-----------------------------------------------------------------------
!
      IF (gtype.lt.0) THEN
!
!  Set offsets due array packing into 1D array in column-major order.
!
        SELECT CASE (ABS(MyType))
          CASE (p2dvar, p3dvar)
            Ioff=0
            Joff=1
          CASE (r2dvar, r3dvar)
            Ioff=1
            Joff=0
          CASE (u2dvar, u3dvar)
            Ioff=0
            Joff=0
          CASE (v2dvar, v3dvar)
            Ioff=1
            Joff=1
          CASE DEFAULT
            Ioff=1
            Joff=0
        END SELECT
!
!  Allocate and initialize scratch work array.
!
        Npts=Isize*Jsize
        IF (.not.allocated(Awrk)) THEN
          allocate ( Awrk(Npts) )
          Awrk=IniVal
        END IF
!
!  Scale and gather data from all spawned nodes. Store data into a 1D
!  global array, packed in column-major order. Flag land point with
!  special value.
!
        DO j=Jmin,Jmax
          jc=(j-Joff)*Isize
          DO i=Imin,Imax
            ic=i+Ioff+jc
            Awrk(ic)=A(i,j)*Ascl
#ifdef POSITIVE_ZERO
            IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
              Awrk(ic)=0.0_r8                     ! impose positive zero
            END IF
#endif
            IF (Amask(i,j).eq.0.0_r8) THEN
              Awrk(ic)=spval
            END IF
          END DO
        END DO
!
!  Global reduction of work array.
!
        CALL mp_collect (ng, model, Npts, IniVal, Awrk)
!
!  Remove land points.
!
        ic=0
        DO i=1,Npts
          IF (Awrk(i).lt.spval) THEN
            ic=ic+1
            Awrk(ic)=Awrk(i)
          END IF
        END DO
        Npts=ic
!
!  Write out data: all parallel nodes write a section of the packed
!                  data.
!
        CALL tile_bounds_1d (ng, MyRank, Npts, Istr, Iend)

        start(1)=Istr
        total(1)=Iend-Istr+1
        start(2)=tindex
        total(2)=1

        status=nf90_put_var(ncid, ncvarid, Awrk(Istr:), start, total)
        nf_fwrite2d=status
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Deallocate scratch work array.
!-----------------------------------------------------------------------
!
      IF (allocated(Awrk)) THEN
        deallocate (Awrk)
      END IF

      RETURN
      END FUNCTION nf_fwrite2d

#else

!
!***********************************************************************
      FUNCTION nf_fwrite2d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, Ascl,                   &
# ifdef MASKING
     &                      Amask,                                      &
# endif
     &                      A, SetFillVal)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj

      real(r8), intent(in) :: Ascl

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
# endif
      integer :: i, j, ic, Npts
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen, MyType, status

      integer, dimension(3) :: start, total

      integer :: nf_fwrite2d

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      MyType=gtype

      SELECT CASE (ABS(MyType))
        CASE (p2dvar, p3dvar)
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
        CASE (r2dvar, r3dvar)
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
        CASE (u2dvar, u3dvar)
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
        CASE (v2dvar, v3dvar)
          Imin=IOBOUNDS(ng)%ILB_v
          Imax=IOBOUNDS(ng)%IUB_v
          Jmin=IOBOUNDS(ng)%JLB_v
          Jmax=IOBOUNDS(ng)%JUB_v
        CASE DEFAULT
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
      END SELECT

      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      IJlen=Ilen*Jlen

# ifdef MASKING
!
!  Set switch to replace land areas with fill value, spval.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
# endif
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Awrk=0.0_r8

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order.
#  ifdef MASKING
#   ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#   else
!  Overwrite masked points with special value.
#   endif
#  endif
!-----------------------------------------------------------------------
!
      CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,                  &
     &                  tindex, gtype, Ascl,                            &
#  ifdef MASKING
     &                  Amask,                                          &
#  endif
     &                  A, Npts, Awrk, SetFillVal)
# else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order.
#  ifdef MASKING
#   ifdef WRITE_WATER
!  Remove land points and pack water points into 1D-array.
#   else
!  Overwrite masked points with special value.
#   endif
#  endif
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
        ic=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ic=ic+1
            Awrk(ic)=A(i,j)*Ascl
#  ifdef MASKING
            IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk(ic)=spval
            END IF
#  endif
          END DO
        END DO
        Npts=IJlen
#  ifdef MASKING
      ELSE
        Npts=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            IF (Amask(i,j).gt.0.0_r8) THEN
              Npts=Npts+1
              Awrk(Npts)=A(i,j)*Ascl
            END IF
          END DO
        END DO
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      nf_fwrite2d=nf90_noerr
      IF (OutThread) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Ilen
          start(2)=1
          total(2)=Jlen
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
        status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
      END IF
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast IO error flag to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcasti (ng, model, status)
# endif
      nf_fwrite2d=status

      RETURN
      END FUNCTION nf_fwrite2d
#endif
      END MODULE nf_fwrite2d_mod
