#include "cppdefs.h"
      SUBROUTINE regrid (ng, model, ncid, ncvarid, gtype, iflag,        &
     &                   Nx, Ny, Ainp, Amin, Amax,                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   Xout, Yout, Aout)
!
!svn $Id: regrid.F 771 2008-09-30 02:30:00Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates gridded data, Ainp, to model locations    !
!  Xout and Yout.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ncid       NetCDF file ID.                                       !
!     ncvarid    NetCDF variable ID.                                   !
!     gtype      Grid type.                                            !
!     iflag      Interpolation flag (0: linear, 1: cubic).             !
!     Nx         X-dimension size for gridded data, Ainp.              !
!     Ny         Y-dimension size for gridded data, Ainp.              !
!     Ainp       Gridded data to interpolate from.                     !
!     Amin       Gridded data minimum value.                           !
!     Amax       Gridded data maximum value.                           !
!     LBi        I-dimension Lower bound of data to interpolate, Aout. !
!     UBi        I-dimension Upper bound of data to interpolate, Aout. !
!     LBj        J-dimension Lower bound of data to interpolate, Aout. !
!     UBj        J-dimension Upper bound of data to interpolate, Aout. !
!     Imin       Starting data I-index to interpolate, Aout.           !
!     Imax       Ending   data I-index to interpolate, Aout.           !
!     Jmin       Starting data J-index to interpolate, Aout.           !
!     Jmax       Ending   data J-index to interpolate, Aout.           !
!     Xout       X-locations to interpolate.                           !
!     Yout       Y-locations to interpolate.                           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aout       Interpolated field.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE interpolate_mod
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_reduce
#endif
!
      implicit none
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
      real(r8), intent(inout) :: Ainp(Nx,Ny)

      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Aout(LBi:UBi,LBj:UBj)
!
!  Local variable declarations
!
      logical :: rectangular
      integer :: i, j, Nghost
      integer :: Istr, Iend, Jstr, Jend
      integer :: Xminr, Xmaxr, Yminr, Ymaxr
#ifdef DISTRIBUTE
      integer :: Itile, Jtile, Npts, ic, mc, nc
#endif
      real(r8), parameter :: IJspv = 0.0_r8

      real(r8) :: my_min, my_max, Xmin, Xmax, Ymin, Ymax

      real(r8), dimension(Nx,Ny) :: angle
      real(r8), dimension(Nx,Ny) :: Xinp
      real(r8), dimension(Nx,Ny) :: Yinp

      real(r8), dimension(LBi:UBi,LBj:UBj) :: Iout
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Jout

      real(r8), dimension(7+3*Nx*Ny) :: Awrk

#ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
#endif

!
!-----------------------------------------------------------------------
!  Get input variable coordinates.
!-----------------------------------------------------------------------
!
      CALL get_varcoords (ng, model, ncid, ncvarid, Nx, Ny,             &
     &                    Xmin, Xmax, Xinp, Ymin, Ymax, Yinp,           &
     &                    rectangular)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set input gridded data rotation angle.
!
      DO i=1,Nx
        DO j=1,Ny
          angle(i,j)=0.0_r8
        END DO
      END DO

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory, broadcast gridded data to all nodes.
!-----------------------------------------------------------------------
!
!  Pack all gridded data so it can be broadcasted to all nodes in the
!  group.
!
      ic=0
      mc=Nx*Ny
      nc=2*mc
      IF (MyRank.eq.MyMaster) THEN
        DO j=1,Ny
          DO i=1,Nx
            ic=ic+1
            Awrk(ic   )=Xinp(i,j)
            Awrk(ic+mc)=Yinp(i,j)
            Awrk(ic+nc)=Ainp(i,j)
          END DO
        END DO
        ic=3*mc
        Awrk(ic+1)=Amin
        Awrk(ic+2)=Amax
        Awrk(ic+3)=Xmin
        Awrk(ic+4)=Xmax
        Awrk(ic+5)=Ymin
        Awrk(ic+6)=Ymax
        IF (rectangular) THEN
          Awrk(ic+7)=1.0_r8
        ELSE
          Awrk(ic+7)=0.0_r8
        END IF
      END IF
      Npts=3*mc+7
!
!  Broadcast gridded data to all nodes in the group, itself included.
!
      CALL mp_bcastf (ng, model, Awrk, Npts)
!
!  Unpack broadcasted data.
!
      IF (MyRank.ne.MyMaster) THEN
        ic=0
        DO j=1,Ny
          DO i=1,Nx
            ic=ic+1
            Xinp(i,j)=Awrk(ic   )
            Yinp(i,j)=Awrk(ic+mc)
            Ainp(i,j)=Awrk(ic+nc)
          END DO
        END DO
        ic=3*mc
        Amin=Awrk(ic+1)
        Amax=Awrk(ic+2)
        Xmin=Awrk(ic+3)
        Xmax=Awrk(ic+4)
        Ymin=Awrk(ic+5)
        Ymax=Awrk(ic+6)
        IF (Awrk(ic+7).gt.0.0_r8) THEN
          rectangular=.TRUE.
        ELSE
          rectangular=.FALSE.
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Check if gridded data contains model grid.
!-----------------------------------------------------------------------
!
      IF ((LonMin(ng).lt.Xmin).or.                                      &
     &    (LonMax(ng).gt.Xmax).or.                                      &
     &    (LatMin(ng).lt.Ymin).or.                                      &
     &    (LatMax(ng).gt.Ymax)) THEN
        IF (Master) THEN
          WRITE (stdout,10) Xmin, Xmax, Ymin, Ymax,                     &
     &                      LonMin(ng), LonMax(ng),                     &
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
#ifdef DISTRIBUTE

!  Set physical, non-overlaping (no ghost-points) ranges according to
!  tile rank.
!
      Nghost=0
      CALL get_bounds (ng, MyRank, gtype, Nghost, Itile, Jtile,         &
     &                 Istr, Iend, Jstr, Jend)
#else
      Istr=Imin
      Iend=Imax
      Jstr=Jmin
      Jend=Jmax
#endif
!
!  Find fractional indices (Iout,Jout) of the grid cells in Ainp
!  containing positions to intepolate.
!
      CALL hindices (ng, 1, Nx, 1, Ny, 1, Nx, 1, Ny,                    &
     &               angle, Xinp, Yinp,                                 &
     &               LBi, UBi, LBj, UBj,                                &
     &               Istr, Iend, Jstr, Jend,                            &
     &               Xout, Yout,                                        &
     &               Iout, Jout,                                        &
     &               IJspv, rectangular, Xminr, Xmaxr, Yminr, Ymaxr)

      IF (iflag.eq.linear) THEN
        CALL linterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Ainp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, Xout, Yout,                         &
     &                  Aout, my_min, my_max)
      ELSE IF (iflag.eq.cubic) THEN
        CALL cinterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Ainp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, Xout, Yout,                         &
     &                  Aout, my_min, my_max)
      END IF
!
!  Compute global interpolated field minimum and maximum values.
!  Notice that gridded data values are overwritten.
!
#ifdef DISTRIBUTE
      buffer(1)=my_min
      buffer(2)=my_max
      op_handle(1)='MIN'
      op_handle(2)='MAX'
      CALL mp_reduce (ng, model, 2, buffer, op_handle)
      Amin=buffer(1)
      Amax=buffer(2)
#else
      Amin=my_min
      Amax=my_max
#endif

      RETURN
      END SUBROUTINE regrid
