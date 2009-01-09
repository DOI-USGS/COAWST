#include "cppdefs.h"
#undef DEBUG
      FUNCTION nf_fread2d (ng, model, ncid, ncvarid, tindex, gtype,     &
     &                     Vsize, LBi, UBi, LBj, UBj, Ascl,             &
     &                     Amin, Amax,                                  &
#ifdef MASKING
     &                     Amask,                                       &
#endif
     &                     A)
!
!svn $Id: nf_fread2d.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function reads in a generic floating point 2D array from an    !
!  input NetCDF file.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ncid       NetCDF file ID.                                       !
!     ncvarid    NetCDF variable ID.                                   !
!     tindex     NetCDF time record index to read (integer).           !
!     gtype      Grid type.                                            !
!     Vsize      Variable dimensions in NetCDF file.                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     Ascl       Factor to scale field after reading (real).           !
!     Amask      Land/Sea mask, if any (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Amin       Field minimum value (real).                           !
!     Amax       Field maximum value (real).                           !
!     A          Field to read in (real).                              !
#ifdef DISTRIBUTE
!     gtype      Grid type.                                            !
#endif
!     nf_fread2d Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_scatter
# ifdef DEBUG
      USE distribute_mod, ONLY : mp_dump
# endif
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Vsize(4)

#ifndef DISTRIBUTE
      integer, intent(in) :: gtype
#else
      integer, intent(inout) :: gtype
#endif

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

#ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      logical :: interpolate
      integer :: i, j, ic, Npts, NWpts, status, wtype
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen, MyType, Nghost

      integer, dimension(3) :: start, total

      integer :: nf_fread2d

      real(r8) :: Aval

      real(r8), allocatable :: wrk(:)

#ifdef DEBUG
      character (len=20) :: Aname
#endif
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------

#ifdef DISTRIBUTE
!
!  In some parts of the code the IO is serial and the value of the
!  grid type is only known by the input thread so broadcast its value
!  to all nodes.
!
      CALL mp_bcasti (ng, model, gtype, 1)
#endif
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.  Notice that Nghost is set
!  to zero when processing an adjoint solution.  Therefore, the
!  ghost points will be not assigned in mp_scatter. This is the
!  correct adjoint solution reading.  The ghost points are then
!  assumed to be zero.
!
#ifdef COMPOSED_GRID
      Imin=LOWER_BOUND_I
      Imax=Lm(ng)+GHOST_POINTS
      Jmin=LOWER_BOUND_J
      Jmax=Mm(ng)+GHOST_POINTS
      MyType=gtype
      IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
        Imin=LOWER_BOUND_I+1
      END IF
      IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
        Jmin=LOWER_BOUND_J+1
      END IF
#elif defined REFINED_GRID
      IF (ng.eq.1) THEN
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        MyType=gtype
        IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
          Imin=1
        END IF
        IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
          Jmin=1
        END IF
      ELSE
        Imin=LOWER_BOUND_I
        Imax=Lm(ng)+GHOST_POINTS
        Jmin=LOWER_BOUND_J
        Jmax=Mm(ng)+GHOST_POINTS
        MyType=gtype
        IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
          Imin=LOWER_BOUND_I+1
        END IF
        IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
          Jmin=LOWER_BOUND_J+1
        END IF
      END IF
#else
      Imin=0
      Imax=Lm(ng)+1
      Jmin=0
      Jmax=Mm(ng)+1
      MyType=gtype
      IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.u2dvar)) THEN
        Imin=1
      END IF
      IF ((ABS(MyType).eq.p2dvar).or.(ABS(MyType).eq.v2dvar)) THEN
        Jmin=1
      END IF
#endif
      IF (model.eq.iADM) THEN
        Nghost=0
      ELSE  
        Nghost=GHOST_POINTS
      END IF
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
!
!  Determine if interpolating from coarse gridded data to model grid
!  is required.  This is only allowed for gridded 2D fields.  This is
!  convinient for atmospheric forcing datasets that are usually on
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

#ifdef REFINED_GRID
      Npts=IJlen
#endif

#if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      IF (ABS(MyType).eq.p2dvar) THEN
        Npts=Nxyp(ng)
        wtype=p2dvar
      ELSE IF (ABS(MyType).eq.u2dvar) THEN
        Npts=Nxyu(ng)
        wtype=u2dvar
      ELSE IF (ABS(MyType).eq.v2dvar) THEN
        Npts=Nxyv(ng)
        wtype=v2dvar
      ELSE
        Npts=Nxyr(ng)
        wtype=r2dvar
      END IF
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
#endif
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
#if defined READ_WATER && defined MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=1
        total(2)=tindex
#endif
      END IF
!
!  Allocate scratch work vector. The dimension of this vector is 
!  unknown when interpolating input data to model grid. Notice
!  that the array length is increased by two because the minimum
!  and maximum values are appended in distributed-memory
!  communications.
!
      IF (.not.allocated(wrk)) THEN
        allocate (wrk(Npts+2))
      END IF
!
!-----------------------------------------------------------------------
!  Read in requested field and scale it.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
      IF (InpThread) THEN
        status=nf90_get_var(ncid, ncvarid, wrk, start, total)
        IF (status.eq.nf90_noerr) THEN
          Amin=wrk(1)*Ascl
          Amax=wrk(1)*Ascl
          DO i=1,Npts
            wrk(i)=Ascl*wrk(i)
            Amin=MIN(Amin,wrk(i))
            Amax=MAX(Amax,wrk(i))
          END DO
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, status, 1)
#endif
      IF (status.ne.nf90_noerr) THEN
        nf_fread2d=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  If not interpolating, unpack read field.
!-----------------------------------------------------------------------
!
      IF (.not.interpolate) THEN
#ifdef DISTRIBUTE
        CALL mp_scatter (ng, model, LBi, UBi, LBj, UBj, 1, 1,           &
     &                   Nghost, MyType, Amin, Amax,                    &
# if defined READ_WATER && defined MASKING
     &                   NWpts, SCALARS(ng)%IJwater(1,wtype),           &
# endif
# ifdef REFINED_GRID
     &                   Npts, wrk, A)
# else
     &                   Npts+2, wrk, A)
# endif
#else
        IF (MyType.gt.0) THEN
          ic=0
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j)=wrk(ic)
            END DO
          END DO
# if defined MASKING || defined READ_WATER
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
# endif
        END IF
#endif
      END IF
!
!-----------------------------------------------------------------------
!  If interpolating from gridded data, read its associated locations
!  and interpolate.
!-----------------------------------------------------------------------
!
      IF (interpolate) THEN
        IF (ABS(MyType).eq.u2dvar) THEN
          CALL regrid (ng, model, ncid, ncvarid, MyType, InterpFlag,    &
     &                 Vsize(1), Vsize(2), wrk, Amin, Amax,             &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 GRID(ng) % lonu(LBi,LBj),                        &
     &                 GRID(ng) % latu(LBi,LBj),                        &
     &                 A(LBi,LBj))
        ELSE IF (ABS(MyType).eq.v2dvar) THEN
          CALL regrid (ng, model, ncid, ncvarid, MyType, InterpFlag,    &
     &                 Vsize(1), Vsize(2), wrk, Amin, Amax,             &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 GRID(ng) % lonv,                                 &
     &                 GRID(ng) % latv,                                 &
     &                 A(LBi,LBj))
        ELSE IF (ABS(MyType).eq.r2dvar) THEN
          CALL regrid (ng, model, ncid, ncvarid, MyType, InterpFlag,    &
     &                 Vsize(1), Vsize(2), wrk, Amin, Amax,             &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Imin, Imax, Jmin, Jmax,                          &
     &                 GRID(ng) % lonr(LBi,LBj),                        &
     &                 GRID(ng) % latr(LBi,LBj),                        &
     &                 A(LBi,LBj))
        END IF
      END IF
#if defined DISTRIBUTE && defined DEBUG
!
!-----------------------------------------------------------------------
!  If debugging, write distributed data into formatted files.
!-----------------------------------------------------------------------
!
      status=nf90_inquire_variable(ncid, ncvarid, Aname)
      CALL mp_dump (ng, MyRank, MyType, LBi, UBi, LBj, UBj, 1, 1, A,    &
     &              Aname)
#endif
!
!-----------------------------------------------------------------------
!  Deallocate scratch work vector.
!-----------------------------------------------------------------------
!
      IF (allocated(wrk)) THEN
        deallocate (wrk)
      END IF
      nf_fread2d=status

      RETURN
      END FUNCTION nf_fread2d
