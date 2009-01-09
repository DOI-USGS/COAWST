#include "cppdefs.h"
#undef DEBUG
      FUNCTION nf_fread4d (ng, model, ncid, ncvarid, tindex, gtype,     &
     &                     Vsize,                                       &
     &                     LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt,      &
     &                     Ascl, Amin, Amax,                            &
#ifdef MASKING
     &                     Amask,                                       &
#endif
     &                     A)
!
!svn $Id: nf_fread4d.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function reads in a generic floating point 4D array from an    !
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
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     LBt        Time-dimension Lower bound.                           !
!     UBt        Time-dimension Upper bound.                           !
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
!     nf_fread4d Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_scatter
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBt, UBt
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
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBt:UBt)
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, fourth, Npts, NWpts, status, wtype
      integer :: Imin, Imax, Jmin, Jmax, Loff
      integer :: Ilen, Jlen, Klen, IJlen, MyType, Nghost

      integer, dimension(5) :: start, total

      integer :: nf_fread4d

      real(r8) :: BIG = 1.0E+35_r8
      real(r8) :: Aval

#ifdef COMPOSED_GRID
      real(r8), dimension(2+(Lm(ng)+7)*(Mm(ng)+7)*(UBk-LBk+1)) :: wrk
#elif defined REFINED_GRID
      real(r8), allocatable :: wrk(:)
#else
      real(r8), dimension(2+(Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: wrk
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
      IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
        Imin=LOWER_BOUND_I+1
      END IF
      IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
        Jmin=LOWER_BOUND_J+1
      END IF
#elif defined REFINED_GRID
      IF (ng.eq.1) THEN
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        MyType=gtype
        IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
          Imin=1
        END IF
        IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
          Jmin=1
        END IF
      ELSE
        Imin=LOWER_BOUND_I
        Imax=Lm(ng)+GHOST_POINTS
        Jmin=LOWER_BOUND_J
        Jmax=Mm(ng)+GHOST_POINTS
        MyType=gtype
        IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
          Imin=LOWER_BOUND_I+1
        END IF
        IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
          Jmin=LOWER_BOUND_J+1
        END IF
      END IF
#else
      Imin=0
      Imax=Lm(ng)+1
      Jmin=0
      Jmax=Mm(ng)+1
      MyType=gtype
      IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.u3dvar)) THEN
        Imin=1
      END IF
      IF ((ABS(MyType).eq.p3dvar).or.(ABS(MyType).eq.v3dvar)) THEN
        Jmin=1
      END IF
#endif
      IF (LBt.eq.0) THEN
        Loff=1
      ELSE
        Loff=0
      END IF
      IF (model.eq.iADM) THEN
        Nghost=0
      ELSE  
        Nghost=GHOST_POINTS
      END IF
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen

#if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      IF (ABS(MyType).eq.p3dvar) THEN
        Npts=Nxyp(ng)
        wtype=p2dvar
      ELSE IF (ABS(MyType).eq.u3dvar) THEN
        Npts=Nxyu(ng)
        wtype=u2dvar
      ELSE IF (ABS(MyType).eq.v3dvar) THEN
        Npts=Nxyv(ng)
        wtype=v2dvar
      ELSE
        Npts=Nxyr(ng)
        wtype=r2dvar
      END IF
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
      Npts=Npts*Klen
#endif
!
!-----------------------------------------------------------------------
!  Read in requested field and scale it.
!-----------------------------------------------------------------------
#ifdef REFINED_GRID
      Npts=IJlen*Klen
      IF (.not.allocated(wrk)) THEN
        allocate (wrk(Npts+2))
      END IF
#endif
!
!  Proccess data as 3D slides.
!
      Amin=BIG
      Amax=-BIG
      DO fourth=LBt,UBt
        IF (MyType.gt.0) THEN
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
          Npts=IJlen*Klen
#if defined READ_WATER && defined MASKING
        ELSE
          start(1)=1+(fourth+Loff-1)*Npts
          total(1)=Npts
          start(2)=1
          total(2)=tindex
#endif
        END IF
        status=nf90_noerr
        IF (InpThread) THEN
          status=nf90_get_var(ncid, ncvarid, wrk, start, total)
          IF (status.eq.nf90_noerr) THEN
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
          nf_fread4d=status
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Unpack read field.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
        CALL mp_scatter (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                   Nghost, MyType, Amin, Amax,                    &
# if defined READ_WATER && defined MASKING
     &                   NWpts, SCALARS(ng)%IJwater(1,wtype),           &
# endif
     &                   Npts, wrk, A(LBi,LBj,LBk,fourth))
#else        
        IF (MyType.gt.0) THEN
          ic=0
          DO k=LBk,UBk
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ic=ic+1
                A(i,j,k,fourth)=wrk(ic)
              END DO
            END DO
          END DO
# if defined MASKING || defined READ_WATER
        ELSE
          ic=0
          DO k=LBk,UBk
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                IF (Amask(i,j).gt.0.0_r8) THEN
                  ic=ic+1
                  A(i,j,k,fourth)=wrk(ic)
                ELSE
                  A(i,j,k,fourth)=0.0_r8
                END IF
              END DO
            END DO
          END DO
# endif
        END IF
#endif
      END DO

#ifdef REFINED_GRID
      IF (allocated(wrk)) THEN
        deallocate (wrk)
      END IF
#endif

      nf_fread4d=status

      RETURN
      END FUNCTION nf_fread4d
