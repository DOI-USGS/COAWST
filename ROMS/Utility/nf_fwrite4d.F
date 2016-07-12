#include "cppdefs.h"
      MODULE nf_fwrite4d_mod
!
!svn $Id: nf_fwrite4d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 4D array into an  !
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
!     LBk          K-dimension Lower bound.                            !
!     UBk          K-dimension Upper bound.                            !
!     LBt          Time-dimension Lower bound.                         !
!     UBt          Time-dimension Upoer bound.                         !
!     Amask        land/Sea mask, if any (real).                       !
!     Ascl         Factor to scale field before writing (real).        !
!     A            Field to write out (real).                          !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional).                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite4d  Error flag (integer).                               !
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
      FUNCTION nf_fwrite4d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt,     &
     &                      Ascl,                                       &
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

# ifdef MASKING
!
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt

      real(r8), intent(in) :: Ascl

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:,LBt:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
# endif
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
# endif
      integer :: i, ic, j, jc, k, kc, l, lc, Npts
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax
      integer :: Ioff, Joff, Koff, Loff
      integer :: Istr, Iend
      integer :: Ilen, Isize, Jlen, Jsize, IJlen, IJKlen, Klen, Llen
      integer :: MyType, status

      integer, dimension(5) :: start, total

      integer :: nf_fwrite4d

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
      Klen=UBk-LBk+1
      Llen=UBt-LBt+1

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

        IF (LBk.eq.0) THEN
          Koff=1
        ELSE
          Koff=0
        END IF

        IF (LBt.eq.0) THEN
          Loff=1
        ELSE
          Loff=0
        END IF
!
!  Allocate and initialize scratch work array.
!
        IJlen=Ilen*Jlen
        IJKlen=IJlen*Klen
        Npts=IJKlen*Llen

        IF (.not.allocated(Awrk)) THEN
          allocate ( Awrk(Npts) )
          Awrk=IniVal
        END IF
!
!  Pack and scale tile data.
!
        ic=0
        DO l=LBt,UBt
          DO k=LBk,UBk
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ic=ic+1
                Awrk(ic)=A(i,j,k,l)*Ascl
# ifdef POSITIVE_ZERO
                IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
                  Awrk(ic)=0.0_r8                 ! impose positive zero
                END IF
# endif
# ifdef MASKING
                IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
                  Awrk(ic)=spval
                END IF
# endif
              END DO
            END DO
          END DO
        END DO
!
!  Write out data: all parallel nodes write their own packed tile data.
!
        start(1)=Imin+Ioff
        total(1)=Ilen
        start(2)=Jmin+Joff
        total(2)=Jlen
        start(3)=LBk+Koff
        total(3)=Klen
        start(4)=LBt+Loff
        total(4)=Llen
        start(5)=tindex
        total(5)=1

        status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
        nf_fwrite4d=status
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

        IF (LBk.eq.0) THEN
          Koff=0
        ELSE
          Koff=1
        END IF

        IF (LBt.eq.0) THEN
          Loff=0
        ELSE
          Loff=1
        END IF
!
!  Allocate and initialize scratch work array.
!
        IJlen=Isize*Jsize
        IJKlen=IJlen*Klen
        Npts=IJKlen*Llen

        IF (.not.allocated(Awrk)) THEN
          allocate ( Awrk(Npts) )
          Awrk=IniVal
        END IF
!
!  Scale and gather data from all spawned nodes. Store data into a 1D
!  global array, packed in column-major order. Flag land point with
!  special value.
!
        DO l=LBt,UBt
          lc=(l-Loff)*IJKlen
          DO k=LBk,UBk
            kc=(k-Koff)*IJlen+lc
            DO j=Jmin,Jmax
              jc=(j-Joff)*Isize+kc
              DO i=Imin,Imax
                ic=i+Ioff+jc
                Awrk(ic)=A(i,j,k,l)*Ascl
#  ifdef POSITIVE_ZERO
                IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
                  Awrk(ic)=0.0_r8                 ! impose positive zero
                END IF
#  endif
                IF (Amask(i,j).eq.0.0_r8) THEN
                  Awrk(ic)=spval
                END IF
              END DO
            END DO
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
        nf_fwrite4d=status
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
      END FUNCTION nf_fwrite4d

#else

!
!***********************************************************************
      FUNCTION nf_fwrite4d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt,     &
     &                      Ascl,                                       &
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
      USE distribute_mod, ONLY : mp_bcasti, mp_gather3d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt

      real(r8), intent(in) :: Ascl

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:,LBt:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
# endif
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
# endif
      integer :: i, j, k, ic, fourth, Npts
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, Koff, Loff
      integer :: Ilen, Jlen, Klen, IJlen, MyType, status

      integer, dimension(5) :: start, total

      integer :: nf_fwrite4d

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Awrk
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
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen

      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF

      IF (LBt.eq.0) THEN
        Loff=1
      ELSE
        Loff=0
      END IF

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
!  Process data as 3D slices.
!
      DO fourth=LBt,UBt
        CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                    tindex, gtype, Ascl,                          &
#  ifdef MASKING
     &                    Amask,                                        &
#  endif
     &                    A(:,:,:,fourth), Npts, Awrk, SetFillVal)
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
        nf_fwrite4d=nf90_noerr
        IF (OutThread) THEN
          IF (gtype.gt.0) THEN
            start(1)=1
            total(1)=Ilen
            start(2)=1
            total(2)=Jlen
            start(3)=1
            total(3)=Klen
            start(4)=fourth+Loff
            total(4)=1
            start(5)=tindex
            total(5)=1
#  ifdef MASKING
          ELSE
            start(1)=1+(fourth+Loff-1)*Npts
            total(1)=Npts
            start(2)=tindex
            total(2)=1
#  endif
          END IF
#  ifdef POSITIVE_ZERO
          DO ic=1,Npts
            IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
              Awrk(ic)=0.0_r8                     ! impose positive zero
            END IF
          END DO
#  endif
          status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
          nf_fwrite4d=status
        END IF
      END DO
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
!  Process data as 3D slices.
!
      DO fourth=LBt,UBt
        IF (gtype.gt.0) THEN
          ic=0
          Npts=IJlen*Klen
          DO k=LBk,UBk
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ic=ic+1
                Awrk(ic)=A(i,j,k,fourth)*Ascl
#  ifdef MASKING
                IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
                  Awrk(ic)=spval
                END IF
#  endif
              END DO
            END DO
          END DO
#  ifdef MASKING
        ELSE
          Npts=0
          DO k=LBk,UBk
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                IF (Amask(i,j).gt.0.0_r8) THEN
                  Npts=Npts+1
                  Awrk(Npts)=A(i,j,k,fourth)*Ascl
                END IF
              END DO
            END DO
          END DO
#  endif
        END IF
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
        nf_fwrite4d=nf90_noerr
        IF (OutThread) THEN
          IF (gtype.gt.0) THEN
            start(1)=1
            total(1)=Ilen
            start(2)=1
            total(2)=Jlen
            start(3)=1
            total(3)=Klen
            start(4)=fourth+Loff
            total(4)=1
            start(5)=tindex
            total(5)=1
#  ifdef MASKING
          ELSE
            start(1)=1+(fourth+Loff-1)*Npts
            total(1)=Npts
            start(2)=tindex
            total(2)=1
#  endif
          END IF
#  ifdef POSITIVE_ZERO
          DO ic=1,Npts
            IF (ABS(Awrk(ic)).eq.0.0_r8) THEN
              Awrk(ic)=0.0_r8                     ! impose positive zero
            END IF
          END DO
#  endif
          status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
        END IF
      END DO
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast IO error flag to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcasti (ng, model, status)
# endif
      nf_fwrite4d=status

      RETURN
      END FUNCTION nf_fwrite4d
#endif
      END MODULE nf_fwrite4d_mod
