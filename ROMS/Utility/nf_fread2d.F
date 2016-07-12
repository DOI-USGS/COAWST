#include "cppdefs.h"
      MODULE nf_fread2d_mod
!
!svn $Id: nf_fread2d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function reads in a generic floating point 2D array from an    !
!  input NetCDF file.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
!     ncvname    NetCDF variable name (string)                         !
!     ncvarid    NetCDF variable ID (integer)                          !
!     tindex     NetCDF time record index to read (integer)            !
!     gtype      C-grid type (integer)                                 !
!     Vsize      Variable dimensions in NetCDF file (integer 1D array) !
!     LBi        I-dimension Lower bound (integer)                     !
!     UBi        I-dimension Upper bound (integer)                     !
!     LBj        J-dimension Lower bound (integer)                     !
!     UBj        J-dimension Upper bound (integer)                     !
!     Ascl       Factor to scale field after reading (real).           !
!     Amask      Land/Sea mask, if any (real 2D array)                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Amin       Field minimum value (real)                            !
!     Amax       Field maximum value (real)                            !
!     A          Field to read in (real 2D array)                      !
!     nf_fread2d Error flag (integer)                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS

#if defined PARALLEL_IN && defined DISTRIBUTE
!
!***********************************************************************
      FUNCTION nf_fread2d (ng, model, ncname, ncid,                     &
     &                     ncvname, ncvarid,                            &
     &                     tindex, gtype, Vsize,                        &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Ascl, Amin, Amax,                            &
# ifdef MASKING
     &                     Amask,                                       &
# endif
     &                     A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_bcasti, mp_collect, mp_reduce
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Vsize(4)

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      logical :: interpolate

      logical, dimension(3) :: foundit

      integer :: i, ic, ij, j, jc, np, MyNpts, Npts
      integer :: Imin, Imax, Isize, Jmin, Jmax, Jsize, IJsize
      integer :: Istr, Iend, Jstr, Jend
      integer :: Ioff, Joff, IJoff
      integer :: Ilen, Itile, Jlen, Jtile, IJlen
      integer :: Cgrid, MyType, ghost, status, wtype

      integer, dimension(3) :: start, total

      integer :: nf_fread2d

      real(r8) :: Afactor, Aoffset, Aspval, cff

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), dimension(2) :: buffer
      real(r8), dimension(3) :: AttValue

      real(r8), allocatable :: Awrk(:,:)
# if defined MASKING && defined READ_WATER
      real(r8), allocatable :: A2d(:)
# endif
      real(r8), allocatable :: wrk(:)

      character (len= 3), dimension(2) :: op_handle
      character (len=12), dimension(3) :: AttName
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set the offsets for variables with starting
!  zero-index.  Recall the NetCDF does not support a zero-index.
!
!  Notice that (Imin,Jmin) and (Imax,Jmax) are the corner of the
!  computational tile. If ghost=0, ghost points are not processed.
!  They will be processed elsewhere by the appropriate call to any
!  of the routines in "mp_exchange.F".  If ghost=1, the ghost points
!  are read.
!
# ifdef NO_READ_GHOST
      ghost=0                        ! non-overlapping, no ghost points
# else
      IF (model.eq.iADM) THEN
        ghost=0                      ! non-overlapping, no ghost points
      ELSE
        ghost=1                      ! overlapping, read ghost points
      END IF
# endif

      MyType=gtype

      SELECT CASE (ABS(MyType))
        CASE (p2dvar, p3dvar)
          Cgrid=1
          Isize=IOBOUNDS(ng)%xi_psi
          Jsize=IOBOUNDS(ng)%eta_psi
        CASE (r2dvar, r3dvar)
          Cgrid=2
          Isize=IOBOUNDS(ng)%xi_rho
          Jsize=IOBOUNDS(ng)%eta_rho
        CASE (u2dvar, u3dvar)
          Cgrid=3
          Isize=IOBOUNDS(ng)%xi_u
          Jsize=IOBOUNDS(ng)%eta_u
        CASE (v2dvar, v3dvar)
          Cgrid=4
          Isize=IOBOUNDS(ng)%xi_v
          Jsize=IOBOUNDS(ng)%eta_v
        CASE DEFAULT
          Cgrid=2
          Isize=IOBOUNDS(ng)%xi_rho
          Jsize=IOBOUNDS(ng)%eta_rho
      END SELECT

      Imin=BOUNDS(ng)%Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng)%Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng)%Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng)%Jmax(Cgrid,ghost,MyRank)

      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
!
!  Determine if interpolating from coarse gridded data to model grid
!  is required.  This is only allowed for gridded 2D fields.  This is
!  convenient for atmospheric forcing datasets that are usually on
!  coarser grids. The user can provide coarser gridded data to avoid
!  very large input files.
!
      interpolate=.FALSE.
      IF (((Vsize(1).gt.0).and.(Vsize(1).ne.Isize)).or.                 &
     &    ((Vsize(2).gt.0).and.(Vsize(2).ne.Jsize))) THEN
        interpolate=.TRUE.
        Ilen=Vsize(1)
        Jlen=Vsize(2)
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, ncvarid, AttName,        &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) THEN
        nf_fread2d=ioerror
        RETURN
      END IF

      IF (.not.foundit(1)) THEN
        Afactor=1.0_r8
      ELSE
        Afactor=AttValue(1)
      END IF

      IF (.not.foundit(2)) THEN
        Aoffset=0.0_r8
      ELSE
        Aoffset=AttValue(2)
      END IF

      IF (.not.foundit(3)) THEN
        Aspval=spval_check
      ELSE
        Aspval=AttValue(3)
      END IF
!
!-----------------------------------------------------------------------
!  Parallel I/O: Read in tile data from requested field and scale it.
!                Processing both water and land points.
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

        Npts=Ilen*Jlen
!
!  Allocate scratch work arrays.
!
        IF (interpolate) THEN
          IF (.not.allocated(Awrk)) THEN
#ifdef GLOBAL_PERIODIC
            allocate ( Awrk(0:Ilen+1,Jlen+1) )
#else
            allocate ( Awrk(Ilen,Jlen) )
#endif
            Awrk=IniVal
          END IF
          IF (.not.allocated(wrk)) THEN
            allocate ( wrk(Npts) )
            wrk=0.0_r8
          END IF
        ELSE
          IF (.not.allocated(wrk)) THEN
            allocate ( wrk(TileSize(ng)) )
            wrk=0.0_r8
          END IF
        END IF
!
!  Read in data: all parallel nodes read their own tile data.
!
        IF (Interpolate) THEN
          CALL tile_bounds_2d (ng, MyRank, Ilen, Jlen, Itile, Jtile,    &
     &                         Istr, Iend, Jstr, Jend)
          start(1)=Istr
          total(1)=Iend-Istr+1
          start(2)=Jstr
          total(2)=Jend-Jstr+1
        ELSE
          start(1)=Imin+Ioff
          total(1)=Ilen
          start(2)=Jmin+Joff
          total(2)=Jlen
        END IF
        MyNpts=total(1)*total(2)
        start(3)=tindex
        total(3)=1

        status=nf90_get_var(ncid, ncvarid, wrk, start, total)
        nf_fread2d=status
!
!  Scale read data and process fill values, if any.  Compute minimum
!  and maximum values.
!
        IF (status.eq.nf90_noerr) THEN
          Amin=spval
          Amax=-spval
          DO i=1,MyNpts
            IF (ABS(wrk(i)).ge.ABS(Aspval)) THEN
              wrk(i)=0.0_r8                     ! masked with _FillValue
            ELSE
              wrk(i)=Ascl*(Afactor*wrk(i)+Aoffset)
              Amin=MIN(Amin,wrk(i))
              Amax=MAX(Amax,wrk(i))
            END IF
          END DO
!
!  Set minimum and maximum values: global reduction.
!
          buffer(1)=Amin
          op_handle(1)='MIN'
          buffer(2)=Amax
          op_handle(2)='MAX'
          CALL mp_reduce (ng, model, 2, buffer, op_handle)
          Amin=buffer(1)
          Amax=buffer(2)
!
!  Unpack read data. If not interpolating, the data is loaded into
!  model array.
!
          IF (interpolate) THEN
            ic=0
            DO j=Jstr,Jend
              DO i=Istr,Iend
                ic=ic+1
                Awrk(i,j)=wrk(ic)
              END DO
            END DO
#ifdef GLOBAL_PERIODIC
            DO j=Jstr,Jend
              Awrk(0,j) = Awrk(Iend,j)
              Awrk(Iend+1,j) = Awrk(1,j)
            END DO
              cff = 0
            DO i=Istr,Iend
              cff = cff + Awrk(i,Jend)
            END DO
              cff = cff/Jend
            DO i=0,Iend+1
              Awrk(i,Jend+1) = cff
            END DO
            Npts=(Ilen+2)*(Jlen+1)
#endif
            CALL mp_collect (ng, model, Npts, IniVal, Awrk)
          ELSE
            ic=0
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ic=ic+1
                A(i,j)=wrk(ic)
              END DO
            END DO
          END IF
        ELSE
          exit_flag=2
          ioerror=status
        END IF
      END IF

# if defined MASKING && defined READ_WATER
!
!-----------------------------------------------------------------------
!  Parallel I/O: Read in tile data from requested field and scale it.
!                Processing water points only. Interpolation is not
!                allowed.
!-----------------------------------------------------------------------
!
      IF (gtype.lt.0) THEN
!
!  Set number of points to process, grid type switch, and offsets due
!  array packing into 1D array in column-major order.
!
        SELECT CASE (ABS(MyType))
          CASE (p2dvar)
            Npts=IOBOUNDS(ng)%xy_psi
            wtype=p2dvar
            Ioff=0
            Joff=1
          CASE (r2dvar)
            Npts=IOBOUNDS(ng)%xy_rho
            wtype=r2dvar
            Ioff=1
            Joff=0
          CASE (u2dvar)
            Npts=IOBOUNDS(ng)%xy_u
            wtype=u2dvar
            Ioff=0
            Joff=0
          CASE (v2dvar)
            Npts=IOBOUNDS(ng)%xy_v
            wtype=v2dvar
            Ioff=1
            Joff=1
          CASE DEFAULT
            Npts=IOBOUNDS(ng)%xy_rho
            wtype=r2dvar
            Ioff=1
            Joff=0
        END SELECT
        IJsize=Isize*Jsize
!
!  Allocate scratch work arrays.
!
        IF (.not.allocated(A2d)) THEN
          allocate ( A2d(IJsize) )
          A2d=IniVal
        END IF
        IF (.not.allocated(wrk)) THEN
          allocate ( wrk(Npts) )
          wrk=IniVal
        END IF
!
!  Read in data: all parallel nodes read a segment of the 1D data.
!  Recall that water points are pack in the NetCDF file in a single
!  dimension.
!
        CALL tile_bounds_1d (ng, MyRank, Npts, Istr, Iend)

        start(1)=Istr
        total(1)=Iend-Istr+1
        start(2)=1
        total(2)=tindex

        status=nf90_get_var(ncid, ncvarid, wrk(Istr:), start, total)
        nf_fread2d=status
!
!  Global reduction of work array.  We need this because the packing
!  of the water point only affects the model tile partition.
!
        IF (status.eq.nf90_noerr) THEN
          CALL mp_collect (ng, model, Npts, IniVal, wrk)
!
!  Scale read data and process fill values, if any.  Compute minimum
!  and maximum values.
!
          Amin=spval
          Amax=-spval
          DO i=1,Npts
            IF (ABS(wrk(i)).ge.ABS(Aspval)) THEN
              wrk(i)=0.0_r8                     ! masked with _FillValue
            ELSE
              wrk(i)=Ascl*(Afactor*wrk(i)+Aoffset)
              Amin=MIN(Amin,wrk(i))
              Amax=MAX(Amax,wrk(i))
            END IF
          END DO
!
!  Unpack read data.  This is tricky in parallel I/O.  The cheapeast
!  thing to do is reconstruct a packed global array and then select
!  the appropriate values for the tile.
!
          DO np=1,Npts
            ij=SCALARS(ng)%IJwater(np,wtype)
            A2d(ij)=wrk(np)
          END DO

          DO j=Jmin,Jmax
            jc=(j-Joff)*Isize
            DO i=Imin,Imax
              ij=i+Ioff+jc
              A(i,j)=A2d(ij)
            END DO
          END DO
        ELSE
          exit_flag=2
          ioerror=status
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Parallel I/O: If interpolating from gridded data, read its associated
!                locations and interpolate.
!-----------------------------------------------------------------------
!
      IF (interpolate.and.(status.eq.nf90_noerr)) THEN
        SELECT CASE (ABS(MyType))
          CASE (p2dvar, p3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonp,                               &
     &                   GRID(ng) % latp,                               &
     &                   A)
          CASE (r2dvar, r3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonr,                               &
     &                   GRID(ng) % latr,                               &
     &                   A)
          CASE (u2dvar, u3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonu,                               &
     &                   GRID(ng) % latu,                               &
     &                   A)
          CASE (v2dvar, v3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonv,                               &
     &                   GRID(ng) % latv,                               &
     &                   A)

          CASE DEFAULT
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonr,                               &
     &                   GRID(ng) % latr,                               &
     &                   A)
        END SELECT
      END IF
!
!-----------------------------------------------------------------------
!  Deallocate scratch work arrays.
!-----------------------------------------------------------------------
!
      IF (interpolate) THEN
        IF (allocated(Awrk)) THEN
          deallocate ( Awrk )
        END IF
      END IF

# if defined MASKING && defined READ_WATER
      IF (allocated(A2d)) THEN
        deallocate (A2d)
      END IF
# endif

      IF (allocated(wrk)) THEN
        deallocate (wrk)
      END IF

      RETURN
      END FUNCTION nf_fread2d

#else

!
!***********************************************************************
      FUNCTION nf_fread2d (ng, model, ncname, ncid,                     &
     &                     ncvname, ncvarid,                            &
     &                     tindex, gtype, Vsize,                        &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Ascl, Amin, Amax,                            &
# ifdef MASKING
     &                     Amask,                                       &
# endif
     &                     A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_scatter2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Vsize(4)

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname

      real(r8), allocatable :: Awrk(:,:)
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      logical :: interpolate

      logical, dimension(3) :: foundit

      integer :: i, j, ic, Npts, NWpts, status, wtype
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen, MyType
# ifdef DISTRIBUTE
      integer :: Nghost
# endif
      integer, dimension(3) :: start, total

      integer :: nf_fread2d

      real(r8) :: Afactor, Aoffset, Aspval, cff

      real(r8), dimension(3) :: AttValue

      real(r8), allocatable :: wrk(:)

      character (len=12), dimension(3) :: AttName
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
        CASE (p2dvar)
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
        CASE (r2dvar)
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
        CASE (u2dvar)
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
        CASE (v2dvar)
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

# ifdef DISTRIBUTE
!
!  Set the number of tile ghost points, Nghost, to scatter in
!  distributed-memory applications. If Nghost=0, the ghost points
!  are not processed.  They will be processed elsewhere by the
!  appropriate call to any of the routines in "mp_exchange.F".
!
#  ifdef NO_READ_GHOST
      Nghost=0
#  else
      IF (model.eq.iADM) THEN
        Nghost=0
      ELSE
        Nghost=NghostPoints
      END IF
#  endif
# endif
!
!  Determine if interpolating from coarse gridded data to model grid
!  is required.  This is only allowed for gridded 2D fields.  This is
!  convenient for atmospheric forcing datasets that are usually on
!  coarser grids. The user can provide coarser gridded data to avoid
!  very large input files.
!
      interpolate=.FALSE.
      IF (((Vsize(1).gt.0).and.(Vsize(1).ne.Ilen)).or.                  &
     &    ((Vsize(2).gt.0).and.(Vsize(2).ne.Jlen))) THEN
        interpolate=.TRUE.
        Ilen=Vsize(1)
        Jlen=Vsize(2)
      END IF
      IJlen=Ilen*Jlen
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, ncvarid, AttName,        &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) THEN
        nf_fread2d=ioerror
        RETURN
      END IF

      IF (.not.foundit(1)) THEN
        Afactor=1.0_r8
      ELSE
        Afactor=AttValue(1)
      END IF

      IF (.not.foundit(2)) THEN
        Aoffset=0.0_r8
      ELSE
        Aoffset=AttValue(2)
      END IF

      IF (.not.foundit(3)) THEN
        Aspval=spval_check
      ELSE
        Aspval=AttValue(3)
      END IF

# if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      SELECT CASE (ABS(MyType))
        CASE (p2dvar)
          Npts=IOBOUNDS(ng)%xy_psi
          wtype=p2dvar
        CASE (r2dvar)
          Npts=IOBOUNDS(ng)%xy_rho
          wtype=r2dvar
        CASE (u2dvar)
          Npts=IOBOUNDS(ng)%xy_u
          wtype=u2dvar
        CASE (v2dvar)
          Npts=IOBOUNDS(ng)%xy_v
          wtype=v2dvar
        CASE DEFAULT
          Npts=IOBOUNDS(ng)%xy_rho
          wtype=r2dvar
      END SELECT
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
# endif
!
!  Set NetCDF dimension counters for processing requested field.
!
      IF (MyType.gt.0) THEN
        Npts=IJlen
        start(1)=1
        total(1)=Ilen
        start(2)=1
        total(2)=Jlen
        start(3)=tindex
        total(3)=1
# if defined READ_WATER && defined MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=1
        total(2)=tindex
# endif
      END IF
!
!  Allocate scratch work vector. The dimension of this vector is
!  unknown when interpolating input data to model grid. Notice
!  that the array length is increased by two because the minimum
!  and maximum values are appended in distributed-memory
!  communications.
!
      IF (.not.allocated(wrk)) THEN
        IF (interpolate) THEN
          IF (.not.allocated(Awrk)) THEN
#ifdef GLOBAL_PERIODIC
            allocate ( Awrk(0:Ilen+1,Jlen+1) )
#else
            allocate ( Awrk(Ilen,Jlen) )
#endif
            Awrk=0.0_r8
          END IF
          allocate ( wrk(Npts) )
        ELSE
          allocate ( wrk(Npts+2) )
        END IF
        wrk=0.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  Serial I/O: Read in requested field and scale it.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
      IF (InpThread) THEN
        status=nf90_get_var(ncid, ncvarid, wrk, start, total)
        IF (status.eq.nf90_noerr) THEN
          Amin=spval
          Amax=-spval
          DO i=1,Npts
            IF (ABS(wrk(i)).ge.ABS(Aspval)) THEN
              wrk(i)=0.0_r8                     ! masked with _FillValue
            ELSE
              wrk(i)=Ascl*(Afactor*wrk(i)+Aoffset)
              Amin=MIN(Amin,wrk(i))
              Amax=MAX(Amax,wrk(i))
            END IF
          END DO
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, status)
# endif
      IF (status.ne.nf90_noerr) THEN
        exit_flag=2
        ioerror=status
        nf_fread2d=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Serial I/O: If not interpolating, unpack read field.
!-----------------------------------------------------------------------
!
      IF (.not.interpolate) THEN
# ifdef DISTRIBUTE
        CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,               &
     &                     Nghost, MyType, Amin, Amax,                  &
#  if defined READ_WATER && defined MASKING
     &                     NWpts, SCALARS(ng)%IJwater(:,wtype),         &
#  endif
     &                     Npts, wrk, A)
# else
        IF (MyType.gt.0) THEN
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=wrk(ic)
            END DO
          END DO
#  if defined MASKING || defined READ_WATER
        ELSE
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                ic=ic+1
                A(i,j)=wrk(ic)
              ELSE
                A(i,j)=0.0_r8
              END IF
            END DO
          END DO
#  endif
        END IF
# endif
# ifdef DISTRIBUTE
      ELSE
        CALL mp_bcastf (ng, model, wrk)
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Serial I/O: If interpolating from gridded data, read its associated
!  locations and interpolate.
!-----------------------------------------------------------------------
!
      IF (interpolate) THEN
        ic=0
        DO j=1,Jlen
          DO i=1,Ilen
            ic=ic+1
            Awrk(i,j)=wrk(ic)
          END DO
        END DO
#ifdef GLOBAL_PERIODIC
        DO j=1,Jlen
          Awrk(0,j) = Awrk(Ilen,j)
          Awrk(Ilen+1,j) = Awrk(1,j)
        END DO
        cff = 0
        DO i=1,Ilen
          cff = cff + Awrk(i,Jlen)
        END DO
        cff = cff/Ilen
        DO i=0,Ilen+1
          Awrk(i,Jlen+1) = cff
        END DO
        Ilen = Ilen+2
        Jlen = Jlen+1
        Npts=Ilen*Jlen
#endif
!       CALL mp_collect (ng, model, Npts, 0.0_r8, Awrk)
        SELECT CASE (ABS(MyType))
          CASE (p2dvar, p3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonp,                               &
     &                   GRID(ng) % latp,                               &
     &                   A)
          CASE (r2dvar, r3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonr,                               &
     &                   GRID(ng) % latr,                               &
     &                   A)
          CASE (u2dvar, u3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonu,                               &
     &                   GRID(ng) % latu,                               &
     &                   A)
          CASE (v2dvar, v3dvar)
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonv,                               &
     &                   GRID(ng) % latv,                               &
     &                   A)

          CASE DEFAULT
            CALL regrid (ng, model, ncname, ncid, ncvname, ncvarid,     &
     &                   MyType, InterpFlag,                            &
     &                   Ilen, Jlen, Awrk, Amin, Amax,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   GRID(ng) % lonr,                               &
     &                   GRID(ng) % latr,                               &
     &                   A)
        END SELECT
      END IF
!
!-----------------------------------------------------------------------
!  Deallocate scratch work vector.
!-----------------------------------------------------------------------
!
      IF (allocated(wrk)) THEN
        deallocate (wrk)
      END IF
      IF (interpolate) THEN
        IF (allocated(Awrk)) THEN
          deallocate (Awrk)
        END IF
      END IF

      nf_fread2d=status

      RETURN
      END FUNCTION nf_fread2d
#endif
      END MODULE nf_fread2d_mod
