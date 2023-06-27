#include "cppdefs.h"
      MODULE nf_fwrite2d_mod
!
!git $Id$
!svn $Id: nf_fwrite2d.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes out a generic floating point 2D array into       !
!  an output file using either the standard NetCDF library or the      !
!  Parallel-IO (PIO) library.                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number                                  !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
#if defined PIO_LIB && defined DISTRIBUTE
!or   pioFile      PIO file descriptor structure, TYPE(file_desc_t)    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
#endif
!     ncvarid      NetCDF variable ID                                  !
#if defined PIO_LIB && defined DISTRIBUTE
!or   pioVar       PIO variable descriptor structure, TYPE(My_VarDesc) !
!                    pioVar%vd     variable descriptor TYPE(Var_Desc_t)!
!                    pioVar%dkind  variable data kind                  !
!                    pioVar%gtype  variable C-gridtype                 !
#endif
!     tindex       NetCDF time record index to write                   !
!     gtype        Grid type. If negative, only write water points     !
#if defined PIO_LIB && defined DISTRIBUTE
!or   pioDesc      IO data decomposition descriptor, TYPE(IO_desc_t)   !
#endif
!     LBi          I-dimension Lower bound                             !
!     UBi          I-dimension Upper bound                             !
!     LBj          J-dimension Lower bound                             !
!     UBj          J-dimension Upper bound                             !
!     Amask        land/Sea mask, if any (real)                        !
!     Ascl         Factor to scale field before writing (real)         !
!     Adat         Field to write out (real)                           !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (OPTIONAL)                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status       Error flag (integer)                                !
!     MinValue     Minimum value (real, OPTIONAL)                      !
!     MaxValue     Maximum value (real, OPTIONAL)                      !
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
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
      INTERFACE nf_fwrite2d
        MODULE PROCEDURE nf90_fwrite2d
#if defined PIO_LIB && defined DISTRIBUTE
        MODULE PROCEDURE pio_fwrite2d
#endif
      END INTERFACE nf_fwrite2d
!
      CONTAINS

#if defined PARALLEL_IO && defined DISTRIBUTE
!
!***********************************************************************
      FUNCTION nf90_fwrite2d (ng, model, ncid, ncvarid, tindex, gtype,  &
     &                        LBi, UBi, LBj, UBj, Ascl,                 &
# ifdef MASKING
     &                        Amask,                                    &
# endif
     &                        Adat, SetFillVal,                         &
     &                        MinValue, MaxValue) RESULT (status)
!***********************************************************************
!
      USE mod_netcdf

# if defined WRITE_WATER && defined MASKING
!
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(dp), intent(in) :: Ascl
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Adat(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Adat(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out), optional :: MinValue
      real(r8), intent(out), optional :: MaxValue
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
!
# endif
      integer :: i, ic, j, jc, Npts
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ioff, Joff
      integer :: Istr, Iend
      integer :: Ilen, Isize, Jlen, Jsize, IJlen
      integer :: MyType
      integer :: status

      integer, dimension(3) :: start, total
!
      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), allocatable :: Awrk(:)
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
!
!  Set first and last grid point according to staggered C-grid
!  classification.
!
      MyType=gtype
!
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
!
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
            Awrk(ic)=Adat(i,j)*Ascl
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
            Awrk(ic)=Adat(i,j)*Ascl
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
      END FUNCTION nf90_fwrite2d

#else

!
!***********************************************************************
      FUNCTION nf90_fwrite2d (ng, model, ncid, ncvarid, tindex, gtype,  &
     &                        LBi, UBi, LBj, UBj, Ascl,                 &
# ifdef MASKING
     &                        Amask,                                    &
# endif
     &                        Adat, SetFillVal,                         &
     &                        MinValue, MaxValue) RESULT (status)
!***********************************************************************
!
      USE mod_netcdf

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(dp), intent(in) :: Ascl
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Adat(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Adat(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out), optional :: MinValue
      real(r8), intent(out), optional :: MaxValue
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
!
# endif
      integer :: i, j, ic, Npts
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen, MyType
      integer :: status

      integer, dimension(3) :: start, total
!
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      MyType=gtype
!
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
     &                  Adat, Npts, Awrk, SetFillVal)
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
            Awrk(ic)=Adat(i,j)*Ascl
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
              Awrk(Npts)=Adat(i,j)*Ascl
            END IF
          END DO
        END DO
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  If applicable, compute output field minimum and maximum values.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MinValue)) THEN
        IF (OutThread) THEN
          MinValue=spval
          MaxValue=-spval
          DO i=1,Npts
            IF (ABS(Awrk(i)).lt.spval) THEN
              MinValue=MIN(MinValue,Awrk(i))
              MaxValue=MAX(MaxValue,Awrk(i))
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
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
!
      RETURN
      END FUNCTION nf90_fwrite2d
#endif

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      FUNCTION pio_fwrite2d (ng, model, pioFile, pioVar, tindex,        &
     &                       pioDesc,                                   &
     &                       LBi, UBi, LBj, UBj, Ascl,                  &
# ifdef MASKING
     &                       Amask,                                     &
# endif
     &                       Adat, SetFillVal,                          &
     &                       MinValue, MaxValue) RESULT (status)
!***********************************************************************
!
      USE mod_pio_netcdf
!
      USE distribute_mod, ONLY : mp_reduce
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, tindex
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(dp), intent(in) :: Ascl
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Adat(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Adat(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out), optional :: MinValue
      real(r8), intent(out), optional :: MaxValue
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (IO_Desc_t),   intent(inout) :: pioDesc
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations.
!
      logical :: LandFill, Lminmax

      logical,  pointer :: Lwater(:,:)
!
      integer :: i, j
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Cgrid, dkind, ghost, gtype
      integer :: status

      integer, dimension(3) :: start, total
!
      real(r8), dimension(2) :: rbuffer

      real(r4), pointer :: Awrk4(:,:)
      real(r8), pointer :: Awrk8(:,:)
!
      character (len= 3), dimension(2) :: op_handle
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
      status=PIO_noerr
!
!  Set first and last tile computational grid point according to the
!  staggered C-grid location. Ghost points are not included.
!
      ghost=0
      dkind=pioVar%dkind
      gtype=pioVar%gtype
!
      SELECT CASE (gtype)
        CASE (p2dvar, p3dvar)
          Cgrid=1                                         ! PSI-points
        CASE (r2dvar, r3dvar)
          Cgrid=2                                         ! RHO-points
        CASE (u2dvar, u3dvar)
          Cgrid=3                                         ! U-points
        CASE (v2dvar, v3dvar)
          Cgrid=4                                         ! V-points
        CASE DEFAULT
          Cgrid=2                                         ! RHO-points
      END SELECT
!
      Imin=BOUNDS(ng)%Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng)%Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng)%Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng)%Jmax(Cgrid,ghost,MyRank)
!
!  Set switch to compute minimum and maximum values.
!
      IF (PRESENT(MinValue)) THEN
        Lminmax=.TRUE.
        IF (.not.associated(Lwater)) THEN
           allocate ( Lwater(LBi:UBi,LBj:UBj) )
           Lwater=.TRUE.
        END IF
      ELSE
        Lminmax=.FALSE.
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
!-----------------------------------------------------------------------
!  Write out data into NetCDF file.
!-----------------------------------------------------------------------
!
!  Allocate, initialize and load data into local array used for
!  writing. Overwrite masked points with special value.
!
      IF (dkind.eq.PIO_double) THEN                 ! double precision
        IF (.not.associated(Awrk8)) THEN
          allocate ( Awrk8(LBi:UBi,LBj:UBj) )
          Awrk8=0.0_r8
        END IF
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk8(i,j)=Adat(i,j)*Ascl
# ifdef MASKING
            IF((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk8(i,j)=spval
              IF (Lminmax) Lwater(i,j)=.FALSE.
            END iF
# endif
          END DO
        END DO
        IF (Lminmax) THEN
          rbuffer(1)=MINVAL(Awrk8, MASK=Lwater)
          rbuffer(2)=MAXVAL(Awrk8, MASK=Lwater)
        END IF
      ELSE                                          ! single precision
        IF (.not.associated(Awrk4)) THEN
          allocate ( Awrk4(LBi:UBi,LBj:UBj) )
          Awrk4=0.0_r4
        END IF
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk4(i,j)=REAL(Adat(i,j)*Ascl, r4)

# ifdef MASKING
            IF((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk4(i,j)=REAL(spval, r4)
              IF (Lminmax) Lwater(i,j)=.FALSE.
            END iF
# endif
          END DO
        END DO
        IF (Lminmax) THEN
          rbuffer(1)=REAL(MINVAL(Awrk4, MASK=Lwater),r8)
          rbuffer(2)=REAL(MAXVAL(Awrk4, MASK=Lwater),r8)
        END IF
      END IF
!
!  Set unlimited time dimension record to write, if any.
!
      IF (tindex.gt.0) THEN
        CALL PIO_setframe (pioFile,                                     &
     &                     pioVar%vd,                                   &
     &                     INT(tindex, kind=PIO_OFFSET_KIND))
      END IF
!
!  Write out data into NetCDF.
!
      IF (dkind.eq.PIO_double) THEN                 ! double precision
        CALL PIO_write_darray (pioFile,                                 &
     &                         pioVar%vd,                               &
     &                         pioDesc,                                 &
     &                         Awrk8(Imin:Imax,Jmin:Jmax),              &
     &                         status)
      ELSE                                          ! single precision
        CALL PIO_write_darray (pioFile,                                 &
     &                         pioVar%vd,                               &
     &                         piodesc,                                 &
     &                         Awrk4(Imin:Imax,Jmin:Jmax),              &
     &                         status)
      END IF
!
!-----------------------------------------------------------------------
!  If applicable, compute global minimum and maximum values.
!-----------------------------------------------------------------------
!
      IF (Lminmax) THEN
        op_handle(1)='MIN'
        op_handle(2)='MAX'
        CALL mp_reduce (ng, model, 2, rbuffer, op_handle)
        MinValue=rbuffer(1)
        MaxValue=rbuffer(2)
        IF (associated(Lwater)) deallocate (Lwater)
      END IF
!
!  Deallocate local array.
!
      IF (dkind.eq.PIO_double) THEN
        IF (associated(Awrk8)) deallocate (Awrk8)
      ELSE
        IF (associated(Awrk4)) deallocate (Awrk4)
      END IF
!
      RETURN
      END FUNCTION pio_fwrite2d
#endif
!
      END MODULE nf_fwrite2d_mod
