#include "cppdefs.h"
      MODULE regrid_mod
!
!git $Id$
!svn $Id: regrid.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates gridded data, Finp, to model locations    !
!  Xout and Yout.  Input data are allow to have longitudes ranging     !
!  from [-180 180] or [0 360].  Usually, the ROMS application have     !
!  longitudes in the range [-180 180].                                 !
!                                                                      !
!  The Imin, Imax, Jmin, and Jmax indices are the global values of     !
!  application domain in serial I/O. However, in parallel I/O they     !
!  are the tiled values.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
#if defined PIO_LIB && defined DISTRIBUTE
!or   pioFile    PIO file descriptor structure, TYPE(file_desc_t)      !
!                  pioFile%fh         file handler                     !
!                  pioFile%iosystem   IO system descriptor (struct)    !
#endif
!     ncvname    NetCDF variable name (string)                         !
!     ncvarid    NetCDF variable ID (integer)                          !
#if defined PIO_LIB && defined DISTRIBUTE
!or   pioDesc    PIO variable descriptor structure, TYPE(My_VarDesc)   !
!                  pioVar%vd     variable descriptor TYPE(Var_Desc_t)  !
!                  pioVar%dkind  variable data kind                    !
!                  pioVar%gtype  variable C-gridtype                   !
#endif
!     gtype      C-grid type (integer)                                 !
!     iflag      Interpolation flag (integer, 0: linear, 1: cubic)     !
!     Nx         X-dimension size for gridded data, Finp (integer)     !
!     Ny         Y-dimension size for gridded data, Finp (integer)     !
!     Finp       Gridded data to interpolate from (real)               !
!     Amin       Gridded data minimum value (integer)                  !
!     Amax       Gridded data maximum value (integer)                  !
!     LBi        Fout I-dimension Lower bound (integer)                !
!     UBi        Fout I-dimension Upper bound (integer)                !
!     LBj        Fout J-dimension Lower bound (integer)                !
!     UBj        Fout J-dimension Upper bound (integer)                !
!     Imin       Fout starting data I-index (integer)                  !
!     Imax       Fout ending   data I-index (integer)                  !
!     Jmin       Fout starting data J-index (integer)                  !
!     Jmax       Fout ending   data J-index (integer)                  !
# ifdef MASKING
!     MyMask     Fout Land/Sea mask (real)                             !
# endif
!     MyXout     Work X-locations (longitude) for regridding (real)    !
!     Xout       X-locations (longitude) to interpolate (real)         !
!     Yout       Y-locations (latitude)  to interpolate (real)         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated field (real)                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE roms_interpolate_mod
!
#ifdef DISTRIBUTE
      USE distribute_mod,    ONLY : mp_reduce
#endif
      USE get_varcoords_mod, ONLY : get_varcoords
#if defined DISTRIBUTE && defined REGRID_SHAPIRO
      USE mp_exchange_mod,   ONLY : mp_exchange2d
#endif
#ifdef REGRID_SHAPIRO
      USE shapiro_mod,       ONLY : shapiro2d_tile
#endif
      USE strings_mod,       ONLY : FoundError
!
      implicit none
!
!  Interface for same name routine overloading.
!
      PUBLIC :: regrid_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: regrid_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE regrid_nf90 (ng, model, ncname, ncid,                  &
     &                        ncvname, ncvarid, gtype, iflag,           &
     &                        Nx, Ny, Finp, Amin, Amax,                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Imin, Imax, Jmin, Jmax,                   &
#ifdef MASKING
     &                        MyMask,                                   &
#endif
     &                        MyXout, Xout, Yout, Fout)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, gtype, iflag
      integer, intent(in) :: Nx, Ny
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
!
      real(r8), intent(inout) :: Amin, Amax
!
      real(r8), intent(inout) :: Finp(Nx,Ny)

      real(r8), intent(inout) :: MyXout(LBi:UBi,LBj:UBj)

#ifdef MASKING
      real(r8), intent(in) :: MyMask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations
!
      logical :: EastLon, rectangular
!
      integer :: i, j
      integer :: Istr, Iend, Jstr, Jend
#ifdef REGRID_SHAPIRO
      integer :: IminS, ImaxS, JminS, JmaxS, tile
#endif
#ifdef DISTRIBUTE
      integer :: Cgrid, ghost
#endif
!
      real(r8), parameter :: IJspv = 0.0_r8

      real(r8) :: my_min, my_max, Xmin, Xmax, Ymin, Ymax
      real(r8) :: MyLonMin, MyLonMax

      real(r8), dimension(Nx,Ny) :: angle
      real(r8), dimension(Nx,Ny) :: Xinp
      real(r8), dimension(Nx,Ny) :: Yinp

      real(r8), dimension(LBi:UBi,LBj:UBj) :: Iout
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Jout

#ifdef DISTRIBUTE
      real(r8), dimension(2) :: rbuffer
!
      character (len=3), dimension(2) :: op_handle
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", regrid_nf90"
!
!-----------------------------------------------------------------------
!  Get input variable coordinates.
!-----------------------------------------------------------------------
!
      CALL get_varcoords (ng, model, ncname, ncid,                      &
     &                    ncvname, ncvarid, Nx, Ny,                     &
     &                    Xmin, Xmax, Xinp, Ymin, Ymax, Yinp,           &
     &                    rectangular)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set input gridded data rotation angle.
!
      DO i=1,Nx
        DO j=1,Ny
          angle(i,j)=0.0_r8
        END DO
      END DO
!
!  Initialize local fractional coordinates arrays to avoid
!  deframentation.
!
      Iout=0.0_r8
      Jout=0.0_r8
!
!  Copy longitude coordinate Xout to MyXout. If the longitude of the
!  data is from a global grid [0-360] or in degrees_east, convert Xout
!  to east longitudes (MyXout) to facilitate regridding. In such case,
!  positive multiples of 360 map to 360 and negative multiples of 360
!  map to zero using the MODULO intrinsic Fortran function.
!
      IF ((Xmin.ge.0.0_r8).and.(Xmax.gt.0.0_r8).and.                    &
     &    ((Xmax-Xmin).gt.315.0_r8)) THEN
        EastLon=.TRUE.
        MyLonMin=MODULO(LonMin(ng), 360.0_r8)
        IF ((MyLonMin.eq.0.0_r8).and.                                   &
     &          (LonMin(ng).gt.0.0_r8)) MyLonMin=360.0_r8
        MyLonMax=MODULO(LonMax(ng), 360.0_r8)
        IF ((MyLonMax.eq.0.0_r8).and.                                   &
     &          (LonMax(ng).gt.0.0_r8)) MyLonMax=360.0_r8
      ELSE
        EastLon=.FALSE.
        MyLonMin=LonMin(ng)
        MyLonMax=LonMax(ng)
      END IF
      IF (EastLon) THEN
        DO j=LBj,UBj
          DO i=LBi,UBi
            MyXout(i,j)=MODULO(Xout(i,j), 360.0_r8)   ! range [0 360]
            IF ((MyXout(i,j).eq.0.0_r8).and.                            &
     &          (Xout(i,j).gt.0.0_r8)) MyXout(i,j)=360.0_r8
          END DO
        END DO
      ELSE
        DO j=LBj,UBj
          DO i=LBi,UBi
            MyXout(i,j)=Xout(i,j)                     ! range [-180 180]
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Check if gridded data contains model grid.
!-----------------------------------------------------------------------
!
      IF ((MyLonMin  .lt.Xmin).or.                                      &
     &    (MyLonMax  .gt.Xmax).or.                                      &
     &    (LatMin(ng).lt.Ymin).or.                                      &
     &    (LatMax(ng).gt.Ymax)) THEN
        IF (Master) THEN
          WRITE (stdout,10) Xmin, Xmax, Ymin, Ymax,                     &
     &                      MyLonMin  , MyLonMax,                       &
     &                      LatMin(ng), LatMax(ng)
 10       FORMAT (/, ' REGRID - input gridded data does not contain',   &
     &               ' model grid:', /,                                 &
     &            /,10x,'Gridded:  LonMin = ',f9.4,' LonMax = ',f9.4,   &
     &            /,10x,'          LatMin = ',f9.4,' LatMax = ',f9.4,   &
     &            /,10x,'Model:    LonMin = ',f9.4,' LonMax = ',f9.4,   &
     &            /,10x,'          LatMin = ',f9.4,' LatMax = ',f9.4)
        END IF
        exit_flag=4
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Interpolate (bilinear or bicubic) to requested positions.
!-----------------------------------------------------------------------
!
!  Set tile starting and ending indices.
!
#ifdef DISTRIBUTE
      ghost=0                        ! non-overlapping, no ghost points

      SELECT CASE (ABS(gtype))
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT
          Cgrid=1
      END SELECT

      Istr=BOUNDS(ng)%Imin(Cgrid,ghost,MyRank)
      Iend=BOUNDS(ng)%Imax(Cgrid,ghost,MyRank)
      Jstr=BOUNDS(ng)%Jmin(Cgrid,ghost,MyRank)
      Jend=BOUNDS(ng)%Jmax(Cgrid,ghost,MyRank)
#else
      Istr=Imin
      Iend=Imax
      Jstr=Jmin
      Jend=Jmax
#endif
!
!  Find fractional indices (Iout,Jout) of the grid cells in Finp
!  containing positions to intepolate.
!
      CALL hindices (ng, 1, Nx, 1, Ny, 1, Nx, 1, Ny,                    &
     &               angle, Xinp, Yinp,                                 &
     &               LBi, UBi, LBj, UBj,                                &
     &               Istr, Iend, Jstr, Jend,                            &
     &               MyXout, Yout,                                      &
     &               Iout, Jout,                                        &
     &               IJspv, rectangular)

      IF (iflag.eq.linear) THEN
        CALL linterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Finp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, MyXout, Yout,                       &
     &                  Fout, my_min, my_max)
      ELSE IF (iflag.eq.cubic) THEN
        CALL cinterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Finp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, MyXout, Yout,                       &
     &                  Fout, my_min, my_max)
      END IF

#ifdef REGRID_SHAPIRO
!
!  Apply Shapiro filter to interpolated date. Often, input data is
!  coarse and field have small scale noise. So, it is desirable to
!  the 2-delta noise.
!
# ifdef DISTRIBUTE
      tile=MyRank
      CALL mp_exchange2d (ng, tile, 1, 1,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Fout)
# else
      tile=-1
# endif
# ifdef NESTING
      IminS=BOUNDS(ng)%Istr(tile)-4
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-4
      JmaxS=BOUNDS(ng)%Jend(tile)+3
# else
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
# endif
!
      CALL shapiro2d_tile (ng, tile, 1,                                 &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# ifdef MASKING
     &                     MyMask,                                      &
# endif
     &                     Fout)
#endif
!
!  Compute global interpolated field minimum and maximum values.
!  Notice that gridded data values are overwritten.
!
#ifdef DISTRIBUTE
      rbuffer(1)=my_min
      rbuffer(2)=my_max
      op_handle(1)='MIN'
      op_handle(2)='MAX'
      CALL mp_reduce (ng, model, 2, rbuffer, op_handle)
      Amin=rbuffer(1)
      Amax=rbuffer(2)
#else
      Amin=my_min
      Amax=my_max
#endif
!
      RETURN
      END SUBROUTINE regrid_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE regrid_pio (ng, model, ncname, pioFile,                &
     &                       ncvname, pioVar, gtype, iflag,             &
     &                       Nx, Ny, Finp, Amin, Amax,                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Imin, Imax, Jmin, Jmax,                    &
# ifdef MASKING
     &                       MyMask,                                    &
# endif
     &                       MyXout, Xout, Yout, Fout)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype, iflag
      integer, intent(in) :: Nx, Ny
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
!
      real(r8), intent(inout) :: Amin, Amax
!
      real(r8), intent(inout) :: Finp(Nx,Ny)

      real(r8), intent(inout) :: MyXout(LBi:UBi,LBj:UBj)

# ifdef MASKING
      real(r8), intent(in) :: MyMask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations
!
      logical :: EastLon, rectangular
!
      integer :: i, j
      integer :: Istr, Iend, Jstr, Jend
# ifdef REGRID_SHAPIRO
      integer :: IminS, ImaxS, JminS, JmaxS, tile
# endif
# ifdef DISTRIBUTE
      integer :: Cgrid, ghost
# endif
!
      real(r8), parameter :: IJspv = 0.0_r8

      real(r8) :: my_min, my_max, Xmin, Xmax, Ymin, Ymax
      real(r8) :: MyLonMin, MyLonMax

      real(r8), dimension(Nx,Ny) :: angle
      real(r8), dimension(Nx,Ny) :: Xinp
      real(r8), dimension(Nx,Ny) :: Yinp

      real(r8), dimension(LBi:UBi,LBj:UBj) :: Iout
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Jout

# ifdef DISTRIBUTE
      real(r8), dimension(2) :: rbuffer
!
      character (len=3), dimension(2) :: op_handle
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", regrid_pio"
!
!-----------------------------------------------------------------------
!  Get input variable coordinates.
!-----------------------------------------------------------------------
!
      CALL get_varcoords (ng, model, ncname, pioFile,                   &
     &                    ncvname, pioVar, Nx, Ny,                      &
     &                    Xmin, Xmax, Xinp, Ymin, Ymax, Yinp,           &
     &                    rectangular)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set input gridded data rotation angle.
!
      DO i=1,Nx
        DO j=1,Ny
          angle(i,j)=0.0_r8
        END DO
      END DO
!
!  Initialize local fractional coordinates arrays to avoid
!  deframentation.
!
      Iout=0.0_r8
      Jout=0.0_r8
!
!  Copy longitude coordinate Xout to MyXout. If the longitude of the
!  data is from a global grid [0-360] or in degrees_east, convert Xout
!  to east longitudes (MyXout) to facilitate regridding. In such case,
!  positive multiples of 360 map to 360 and negative multiples of 360
!  map to zero using the MODULO intrinsic Fortran function.
!
      IF ((Xmin.ge.0.0_r8).and.(Xmax.gt.0.0_r8).and.                    &
     &    ((Xmax-Xmin).gt.315.0_r8)) THEN
        EastLon=.TRUE.
        MyLonMin=MODULO(LonMin(ng), 360.0_r8)
        IF ((MyLonMin.eq.0.0_r8).and.                                   &
     &          (LonMin(ng).gt.0.0_r8)) MyLonMin=360.0_r8
        MyLonMax=MODULO(LonMax(ng), 360.0_r8)
        IF ((MyLonMax.eq.0.0_r8).and.                                   &
     &          (LonMax(ng).gt.0.0_r8)) MyLonMax=360.0_r8
      ELSE
        EastLon=.FALSE.
        MyLonMin=LonMin(ng)
        MyLonMax=LonMax(ng)
      END IF
      IF (EastLon) THEN
        DO j=LBj,UBj
          DO i=LBi,UBi
            MyXout(i,j)=MODULO(Xout(i,j), 360.0_r8)   ! range [0 360]
            IF ((MyXout(i,j).eq.0.0_r8).and.                            &
     &          (Xout(i,j).gt.0.0_r8)) MyXout(i,j)=360.0_r8
          END DO
        END DO
      ELSE
        DO j=LBj,UBj
          DO i=LBi,UBi
            MyXout(i,j)=Xout(i,j)                     ! range [-180 180]
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Check if gridded data contains model grid.
!-----------------------------------------------------------------------
!
      IF ((MyLonMin  .lt.Xmin).or.                                      &
     &    (MyLonMax  .gt.Xmax).or.                                      &
     &    (LatMin(ng).lt.Ymin).or.                                      &
     &    (LatMax(ng).gt.Ymax)) THEN
        IF (Master) THEN
          WRITE (stdout,10) Xmin, Xmax, Ymin, Ymax,                     &
     &                      MyLonMin  , MyLonMax,                       &
     &                      LatMin(ng), LatMax(ng)
 10       FORMAT (/, ' REGRID - input gridded data does not contain',   &
     &               ' model grid:', /,                                 &
     &            /,10x,'Gridded:  LonMin = ',f9.4,' LonMax = ',f9.4,   &
     &            /,10x,'          LatMin = ',f9.4,' LatMax = ',f9.4,   &
     &            /,10x,'Model:    LonMin = ',f9.4,' LonMax = ',f9.4,   &
     &            /,10x,'          LatMin = ',f9.4,' LatMax = ',f9.4)
        END IF
        exit_flag=4
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Interpolate (bilinear or bicubic) to requested positions.
!-----------------------------------------------------------------------
!
!  Set tile starting and ending indices.
!
# ifdef DISTRIBUTE
      ghost=0                        ! non-overlapping, no ghost points

      SELECT CASE (ABS(gtype))
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT
          Cgrid=1
      END SELECT

      Istr=BOUNDS(ng)%Imin(Cgrid,ghost,MyRank)
      Iend=BOUNDS(ng)%Imax(Cgrid,ghost,MyRank)
      Jstr=BOUNDS(ng)%Jmin(Cgrid,ghost,MyRank)
      Jend=BOUNDS(ng)%Jmax(Cgrid,ghost,MyRank)
# else
      Istr=Imin
      Iend=Imax
      Jstr=Jmin
      Jend=Jmax
# endif
!
!  Find fractional indices (Iout,Jout) of the grid cells in Finp
!  containing positions to intepolate.
!
      CALL hindices (ng, 1, Nx, 1, Ny, 1, Nx, 1, Ny,                    &
     &               angle, Xinp, Yinp,                                 &
     &               LBi, UBi, LBj, UBj,                                &
     &               Istr, Iend, Jstr, Jend,                            &
     &               MyXout, Yout,                                      &
     &               Iout, Jout,                                        &
     &               IJspv, rectangular)

      IF (iflag.eq.linear) THEN
        CALL linterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Finp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, MyXout, Yout,                       &
     &                  Fout, my_min, my_max)
      ELSE IF (iflag.eq.cubic) THEN
        CALL cinterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Finp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, MyXout, Yout,                       &
     &                  Fout, my_min, my_max)
      END IF

# ifdef REGRID_SHAPIRO
!
!  Apply Shapiro filter to interpolated date. Often, input data is
!  coarse and field have small scale noise. So, it is desirable to
!  the 2-delta noise.
!
#  ifdef DISTRIBUTE
      tile=MyRank
      CALL mp_exchange2d (ng, tile, 1, 1,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Fout)
#  else
      tile=-1
#  endif
#  ifdef NESTING
      IminS=BOUNDS(ng)%Istr(tile)-4
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-4
      JmaxS=BOUNDS(ng)%Jend(tile)+3
#  else
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
#  endif
!
      CALL shapiro2d_tile (ng, tile, 1,                                 &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     MyMask,                                      &
#  endif
     &                     Fout)
# endif
!
!  Compute global interpolated field minimum and maximum values.
!  Notice that gridded data values are overwritten.
!
# ifdef DISTRIBUTE
      rbuffer(1)=my_min
      rbuffer(2)=my_max
      op_handle(1)='MIN'
      op_handle(2)='MAX'
      CALL mp_reduce (ng, model, 2, rbuffer, op_handle)
      Amin=rbuffer(1)
      Amax=rbuffer(2)
# else
      Amin=my_min
      Amax=my_max
# endif
!
      RETURN
      END SUBROUTINE regrid_pio
#endif
      END MODULE regrid_mod
