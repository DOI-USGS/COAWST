#include "cppdefs.h"

#undef CUBIC_MASKED    /* Bicubic interpolation with mask, needs work */

      MODULE roms_interpolate_mod
!
!git $Id$
!svn $Id: interpolate.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This module contains several all-purpose generic interpolations     !
!  routines:                                                           !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    cinterp2d     Bicubic  interpolation for any 2D field.            !
!    linterp2d     Bilinear interpolation for any 2D field.            !
!    hindices      Finds model grid cell for any datum.                !
!    try_range     Binary search of model grid cell for any datum.     !
!    inside        Closed polygon datum search.                        !
!                                                                      !
!  Overloading routines using "roms_interp_type" structure, S, that    !
!  can be used in the ROMS-JEDI interface.                             !
!                                                                      !
!    - Interpolate a ROMS 2D state field to observation locations      !
!      (ObsVetting is an optional output argument)                     !
!                                                                      !
!      roms_datum_interp (S, fld2d, datum(:), method)                  !
!                                                                      !
!    - Interpolate a ROMS 3D state field to observation locations      !
!      (ObsVetting is an optional output argument)                     !
!                                                                      !
!      roms_datum_interp (S, fld3d, zfld3d, zlocs(:), datum(:), method)!
!                                                                      !
!    - Interpolate a ROMS 3D state field to observation locations      !
!      level-by-level and returns datum(1:nlevs,1:nlocs)               !
!                                                                      !
!      roms_datum_interp (S, fld3d, zfld3d, datum(:,:), method)        !
!                                                                      !
!    - Interpolate a ROMS 2D source field to a new destination         !
!      location (2D remapping)                                         !
!                                                                      !
!      roms_horiz_interp (S, fld2d_src, fld2d_dst, method)             !
!                                                                      !
!    - Interpolate a ROMS 3D source field to a new destination         !
!      locations level-by-level (3D remapping)                         !
!                                                                      !
!      roms_horiz_interp (S, fld2d_src, fld2d_dst, method)             !
!                                                                      !
!  Managing routines for "roms_interp_type" structure, S:              !
!                                                                      !
!    - Creates a ROMS interpolation object                             !
!                                                                      !
!      roms_interp_create (ng, S)                                      !
!                                                                      !
!    - Deallocates pointer arrays in the interpolation structure       !
!                                                                      !
!      roms_interp_delete (S)                                          !
!                                                                      !
!    - Computes the horizontal fractional coordinates S%x_dst and      !
!      S%y_dst of the source cells containing the destination values   !
!                                                                      !
!      roms_interp_fractional (S)                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_assemble
#endif
      USE mod_iounits,    ONLY : stdout
      USE strings_mod,    ONLY : FoundError
!
      implicit none
!
!-----------------------------------------------------------------------
!  Set ROMS interpolation structure for Object Oriented Programming.
!-----------------------------------------------------------------------
!
      TYPE :: roms_interp_type

        logical :: rectangular = .FALSE.            ! plaid grid switch
!
        integer :: ng                               ! nested grid number
        integer :: model                            ! model identifier
!
        real(r8):: spval                            ! unbounded value
!
!  Source grid array declaration bounds, tile partition range,
!  longitude/latitude, curvilinear angle, and land/sea mask
!
        integer :: LBi_src, UBi_src, LBj_src, UBj_src
        integer :: Istr_src, Iend_src, Jstr_src, Jend_src
!
        real(r8) :: min_src, max_src
!
        real(r8), allocatable :: lon_src(:,:)
        real(r8), allocatable :: lat_src(:,:)

        real(r8), allocatable :: angle_src(:,:)
        real(r8), allocatable :: mask_src(:,:)
!
!  Destination grid array declaration bounds, tile partition range,
!  longitude/latitude, fractional coordinates, and land/sea mask.
!
        integer :: LBi_dst, UBi_dst, LBj_dst, UBj_dst
        integer :: Istr_dst, Iend_dst, Jstr_dst, Jend_dst
!
        real(r8) :: min_dst, max_dst
!
        real(r8), allocatable :: lon_dst(:,:)
        real(r8), allocatable :: lat_dst(:,:)

        real(r8), allocatable :: mask_dst(:,:)

        real(r8), allocatable :: x_dst(:,:)
        real(r8), allocatable :: y_dst(:,:)

      END TYPE roms_interp_type
!
!  Interpolation methods available
!
      integer, parameter :: BilinearMethod = 0    ! bilinear
      integer, parameter :: BicubicMethod = 1     ! bicubic
!
!-----------------------------------------------------------------------
!  Interface for same name routine overloading.
!-----------------------------------------------------------------------
!
     INTERFACE roms_datum_interp
       MODULE PROCEDURE roms_datum_interp_2d
       MODULE PROCEDURE roms_datum_interp_3d
       MODULE PROCEDURE roms_datum_column_interp
     END INTERFACE roms_datum_interp
!
     INTERFACE roms_horiz_interp
       MODULE PROCEDURE roms_horiz_interp_2d
       MODULE PROCEDURE roms_horiz_interp_3d
     END INTERFACE roms_horiz_interp
!
!-----------------------------------------------------------------------
!  Interpolation routines.
!-----------------------------------------------------------------------
!
      PUBLIC  :: linterp2d
      PUBLIC  :: cinterp2d
      PUBLIC  :: hindices
!
      PUBLIC  :: roms_interp_delete
      PUBLIC  :: roms_interp_fractional
!
      PRIVATE :: inside
      PRIVATE :: try_range
!
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE linterp2d (ng, LBx, UBx, LBy, UBy,                     &
     &                      Xinp, Yinp, Finp,                           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      Iout, Jout,                                 &
     &                      Xout, Yout,                                 &
     &                      Fout, MinVal, MaxVal)
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field, Finp, this routine linearly interpolate !
!  to locations (Xout,Yout).  To facilitate the  interpolation  within !
!  any irregularly gridded 2D field,  the fractional grid cell indices !
!  (Iout,Jout) with respect Finp are needed at input.  Notice that the !
!  routine "hindices" can be used to compute these indices.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     LBi        I-dimension Lower bound of data to interpolate, Fout. !
!     UBi        I-dimension Upper bound of data to interpolate, Fout. !
!     LBj        J-dimension Lower bound of data to interpolate, Fout. !
!     UBj        J-dimension Upper bound of data to interpolate, Fout. !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Iout       I-fractional Xinp grid cell containing Xout.          !
!     Jout       J-fractional Yinp grid cell containing Yout.          !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!     MinVal     Interpolated field minimum value (Optional).          !
!     MaxVal     Interpolated field maximum value (Optional).          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBx, UBx, LBy, UBy
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(r8), intent(in) :: Xinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Yinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Finp(LBx:UBx,LBy:UBy)

      real(r8), intent(in) :: Iout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Jout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out), optional :: MinVal
      real(r8), intent(out), optional :: MaxVal
!
!  Local variable declarations.
!
      integer :: i, i1, i2, j, j1, j2
!
      real(r8) :: p1, p2, q1, q2
      real(r8) :: Fmin, Fmax
!
!-----------------------------------------------------------------------
!  Linearly interpolate requested field
!-----------------------------------------------------------------------
!
      Fmin=1.0E+35_r8
      Fmax=-1.0E+35_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          i1=INT(Iout(i,j))
          i2=i1+1
          j1=INT(Jout(i,j))
          j2=j1+1
          IF (((LBx.le.i1).and.(i1.le.UBx)).and.                        &
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
            p2=REAL(i2-i1,r8)*(Iout(i,j)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Jout(i,j)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            Fout(i,j)=p1*q1*Finp(i1,j1)+                                &
     &                p2*q1*Finp(i2,j1)+                                &
     &                p2*q2*Finp(i2,j2)+                                &
     &                p1*q2*Finp(i1,j2)
            Fmin=MIN(Fmin,Fout(i,j))
            Fmax=MAX(Fmax,Fout(i,j))
          END IF
        END DO
      END DO
!
!  Return minimum and maximum values.
!
      IF (PRESENT(MinVal)) THEN
        MinVal=Fmin
      END IF
      IF (PRESENT(MaxVal)) THEN
        MaxVal=Fmax
      END IF
!
      RETURN
      END SUBROUTINE linterp2d
!
      SUBROUTINE cinterp2d (ng, LBx, UBx, LBy, UBy,                     &
     &                      Xinp, Yinp, Finp,                           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      Iout, Jout,                                 &
     &                      Xout, Yout,                                 &
     &                      Fout, MinVal, MaxVal)
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field,  Finp, at locations (Xinp,Yinp) this    !
!  routine performs bicubic interpolation at locations (Xout,Yout).    !
!  To facilitate the interpolation within any  irregularly  gridded    !
!  field, the fractional grid cell indices (Iout,Jout) with respect    !
!  Finp are needed at input. Notice that the routine "hindices" can    !
!  be used to compute these indices.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     LBi        I-dimension Lower bound of data to interpolate, Fout. !
!     UBi        I-dimension Upper bound of data to interpolate, Fout. !
!     LBj        J-dimension Lower bound of data to interpolate, Fout. !
!     UBj        J-dimension Upper bound of data to interpolate, Fout. !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Iout       I-fractional Xinp grid cell containing Xout.          !
!     Jout       J-fractional Yinp grid cell containing Yout.          !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!     MinVal     Interpolated field minimum value (OPTIONAL).          !
!     MaxVal     Interpolated field maximum value (OPTIONAL).          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBx, UBx, LBy, UBy
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(r8), intent(in) :: Xinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Yinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Finp(LBx:UBx,LBy:UBy)

      real(r8), intent(in) :: Iout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Jout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out), optional :: MinVal
      real(r8), intent(out), optional :: MaxVal
!
!  Local variable declarations.
!
      integer i, ic, iter, i1, i2, j, jc, j1, j2
!
      real(r8) :: a11, a12, a21, a22
      real(r8) :: e11, e12, e21, e22
      real(r8) :: cff, d1, d2, dfc, dx, dy, eta, xi, xy, yx
      real(r8) :: f0, fx, fxx, fxxx, fxxy, fxy, fxyy, fy, fyy, fyyy
!
      real(r8), parameter :: C01 = 1.0_r8/48.0_r8
      real(r8), parameter :: C02 = 1.0_r8/32.0_r8
      real(r8), parameter :: C03 = 0.0625_r8                  ! 1/16
      real(r8), parameter :: C04 = 1.0_r8/6.0_r8
      real(r8), parameter :: C05 = 0.25_r8
      real(r8), parameter :: C06 = 0.5_r8
      real(r8), parameter :: C07 = 0.3125_r8                  ! 5/16
      real(r8), parameter :: C08 = 0.625_r8                   ! 5/8
      real(r8), parameter :: C09 = 1.5_r8
      real(r8), parameter :: C10 = 13.0_r8/24.0_r8

      real(r8), parameter :: LIMTR = 3.0_r8
      real(r8), parameter :: spv = 0.0_r8            ! HGA need work

      real(r8) :: Fmin, Fmax
!
      real(r8), dimension(-1:2,-1:2) :: dfx, dfy, f
!
!-----------------------------------------------------------------------
!  Interpolates requested field locations (Xout,Yout).
!-----------------------------------------------------------------------
!
      Fmin=1.0E+35_r8
      Fmax=-1.0E+35_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          i1=INT(Iout(i,j))
          i2=i1+1
          j1=INT(Jout(i,j))
          j2=j1+1
          IF (((LBx.le.i1).and.(i1.le.UBx)).and.                        &
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
!
!  Determine local fractional coordinates (xi,eta) corresponding to
!  the target point (Xout,Yout) on the grid (Xinp,Yinp). Here, "xi"
!  and "eta" are defined, in such a way, that xi=eta=0 corresponds
!  to the middle of the cell (i1:i1+1,j1:j1+1), while xi=+/-1/2 and
!  eta=+/-1/2 (any combination +/- signs) corresponds to the four
!  corner points of the cell. Inside the cell it is assumed that
!  (Xout,Yout) are expressed via bi-linear functions of (xi,eta),
!  where term proportional to xi*eta does not vanish because
!  coordinate transformation may be at least weakly non-orthogonal
!  due to discretization errors. The associated non-linear system
!  is solved by iterative method of Newton.
!
            xy=Xinp(i2,j2)-Xinp(i1,j2)-Xinp(i2,j1)+Xinp(i1,j1)
            yx=Yinp(i2,j2)-Yinp(i1,j2)-Yinp(i2,j1)+Yinp(i1,j1)
            dx=Xout(i,j)-0.25_r8*(Xinp(i2,j2)+Xinp(i1,j2)+              &
     &                            Xinp(i2,j1)+Xinp(i1,j1))
            dy=Yout(i,j)-0.25_r8*(Yinp(i2,j2)+Yinp(i1,j2)+              &
     &                            Yinp(i2,j1)+Yinp(i1,j1))
!
!  The coordinate transformation matrix:
!
!           e11 e12
!           e21 e22
!
!  contains derivatives of (Xinp,Yinp) with respect to (xi,eta). Because
!  the coordinates may be non-orthogonal (at least due to discretization
!  errors), the nonlinear system
!
!           e11*xi+e12*eta+xy*xi*eta=dx
!           e21*xi+e22*eta+yx*xi*eta=dy
!
!  needs to be solved in order to retain symmetry.
!
            e11=0.5_r8*(Xinp(i2,j2)-Xinp(i1,j2)+Xinp(i2,j1)-Xinp(i1,j1))
            e12=0.5_r8*(Xinp(i2,j2)+Xinp(i1,j2)-Xinp(i2,j1)-Xinp(i1,j1))
            e21=0.5_r8*(Yinp(i2,j2)-Yinp(i1,j2)+Yinp(i2,j1)-Yinp(i1,j1))
            e22=0.5_r8*(Yinp(i2,j2)+Yinp(i1,j2)-Yinp(i2,j1)-Yinp(i1,j1))
!
            cff=1.0_r8/(e11*e22-e12*e21)
            xi=cff*(e22*dx-e12*dy)
            eta=cff*(e11*dy-e21*dx)
!
            DO iter=1,4
              d1=dx-e11*xi-e12*eta-xy*xi*eta
              d2=dy-e21*xi-e22*eta-yx*xi*eta
              a11=e11+xy*eta
              a12=e12+xy*xi
              a21=e21+yx*eta
              a22=e22+yx*xi
              cff=1.0_r8/(a11*a22-a12*a21)
              xi =xi +cff*(a22*d1-a12*d2)
              eta=eta+cff*(a11*d2-a21*d1)
            END DO

#ifndef CUBIC_MASKED
!
!  Genuinely two-dimensional, isotropic cubic interpolation scheme
!  using 12-point stencil.  In the code below the interpolated field,
!  Fout, is expanded into two-dimensional Taylor series of local
!  fractional coordinates "xi" and "eta", retaining all terms of
!  combined power up to third order (that is, xi, eta, xi^2, eta^2,
!  xi*eta, xi^3, eta^3, xi^2*eta, and xi*eta^2), with all
!  coefficients (i.e, derivatives) computed via           x  x
!  two-dimensional finite difference expressions          |  |
!  of "natural" order of accuracy: 4th-order for       x--x--x--x
!  the field itself and its first derivatives in          |  |
!  both directions; and 2nd-order for all higher-      x--x--x--x
!  order derivatives. The permissible range of            |  |
!  of coordinates is -1/2 < xi,eta < +1/2, which          x--x
!  covers the central cell on the stencil, while
!  xi=eta=0 corresponds to its center. This interpolation scheme has
!  the property that if xi,eta=+/-1/2 (any combination of +/- signs)
!  it reproduces exactly value of the function at the corresponding
!  corner of the central "working" cell. However, it does not pass
!  exactly through the  extreme points of the stencil, where either
!  xi=+/-3/2 or eta+/-3/2. And, unlike a split-directional scheme,
!  when interpolating along the line eta=+/-1/2 (similarly xi=+/-1/2),
!  it has non-zero contribution from points on the side from the line,
!  except if xi=-1/2; 0; +1/2 (similarly eta=-1/2; 0; +1/2).
!
            DO jc=-1,2
              DO ic=-1,2
                f(ic,jc)=Finp(MAX(1,MIN(UBx,i1+ic)),                    &
     &                        MAX(1,MIN(UBy,j1+jc)))
              END DO
            END DO

            f0=C07*(f(1,1)+f(1,0)+f(0,1)+f(0,0))-                       &
     &         C02*(f(2,0)+f(2,1)+f(1,2)+f(0,2)+                        &
     &              f(-1,1)+f(-1,0)+f(0,-1)+f(1,-1))

            fx=C08*(f(1,1)+f(1,0)-f(0,1)-f(0,0))-                       &
     &         C01*(f(2,1)+f(2,0)-f(-1,1)-f(-1,0))-                     &
     &         C03*(f(1,2)-f(0,2)+f(1,-1)-f(0,-1))

            fy=C08*(f(1,1)-f(1,0)+f(0,1)-f(0,0))-                       &
     &         C01*(f(1,2)+f(0,2)-f(1,-1)-f(0,-1))-                     &
     &         C03*(f(2,1)-f(2,0)+f(-1,1)-f(-1,0))

            fxy=f(1,1)-f(1,0)-f(0,1)+f(0,0)

            fxx=C05*(f(2,1)-f(1,1)-f(0,1)+f(-1,1)+                      &
     &               f(2,0)-f(1,0)-f(0,0)+f(-1,0))

            fyy=C05*(f(1,2)-f(1,1)-f(1,0)+f(1,-1)+                      &
     &               f(0,2)-f(0,1)-f(0,0)+f(0,-1))

            fxxx=C06*(f(2,1)+f(2,0)-f(-1,1)-f(-1,0))-                   &
     &           C09*(f(1,1)+f(1,0)-f(0,1)-f(0,0))

            fyyy=C06*(f(1,2)+f(0,2)-f(1,-1)-f(0,-1))-                   &
     &           C09*(f(1,1)-f(1,0)+f(0,1)-f(0,0))

            fxxy=C06*(f(2,1)-f(1,1)-f(0,1)+f(-1,1)-                     &
     &                f(2,0)+f(1,0)+f(0,0)-f(-1,0))

            fxyy=C06*(f(1,2)-f(1,1)-f(1,0)+f(1,-1)-                     &
     &                f(0,2)+f(0,1)+f(0,0)-f(0,-1))
#else
!
!  Algorithm below is equivalent to the one above, except that special
!  care is taken to avoid interpolation accross land. This is achieved
!  by shortening the stencil and reducing order of polynomial, if
!  extreme points of the stencil touch land. This is achieved by
!  expressing all f0,fx,fy,...,fxyy in terms of values of interpolated
!  field at the four corners of central cell (which already checked to
!  stay away from land), and eight one-sided differences dfx,dfy (see
!  below) in such a way that field values at the extreme points of the
!  12-point stencil do not participate directly into f0,fx,...,fxyy.
!  Should an extreme point of the stencil touch land, thus making it
!  impossible to compute the corresponding one-sided difference, this
!  difference is retracted toward the center of the stencil.
!
!  Optionally, a slope-limiting algorithm may be employed to prevent
!  spurious oscillations of the interpolant. This is a valuable property,
!  if dealing with rough data, however, as a side effect, it turns off
!  high-order interpolation in the vicinity of extrema.
!
!  The slope-limiting algorithm employed here checks that two consecutive
!  elementary differences, "dfx" and "dfc" have the same sign and differ
!  in magnitude by no more than factor of 3.
!
            f(0,0)=Finp(i1,j1)
            f(1,0)=Finp(i2,j1)
            f(0,1)=Finp(i1,j2)
            f(1,1)=Finp(i2,j2)
!
            dfc=f(1,1)-f(0,1)
            IF (i1+2.gt.UBx) THEN
              dfx(1,1)=dfc
            ELSE IF (Finp(i1+2,j2).eq.spv) THEN
              dfx(1,1)=dfc
            ELSE
              dfx(1,1)=Finp(i1+2,j2)-f(1,1)
# ifdef LIMTR
              IF ((dfx(1,1)*dfc).lt.0.0_r8) THEN
                dfx(1,1)=0.0_r8
              ELSE IF (ABS(dfx(1,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(1,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,0)-f(0,0)
            IF ((i1+2).gt.UBx) THEN
              dfx(1,0)=dfc
            ELSE IF (Finp(i1+2,j1).eq.spv) THEN
              dfx(1,0)=dfc
            ELSE
              dfx(1,0)=Finp(i1+2,j1)-f(1,0)
# ifdef LIMTR
              IF ((dfx(1,0)*dfc).lt.0.0_r8) THEN
                dfx(1,0)=0.0_r8
              ELSE IF (ABS(dfx(1,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(1,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,1)-f(0,1)
            IF (i1-1.lt.1) THEN
              dfx(0,1)=dfc
            ELSE IF (Finp(i1-1,j2).eq.spv) THEN
              dfx(0,1)=dfc
            ELSE
              dfx(0,1)=f(0,1)-Finp(i1-1,j2)
# ifdef LIMTR
              IF ((dfx(0,1)*dfc).lt.0.0_r8) THEN
                dfx(0,1)=0.0_r8
              ELSE IF (ABS(dfx(0,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(0,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,0)-f(0,0)
            IF (i1-1.lt.1) THEN
              dfx(0,0)=dfc
            ELSE IF (Finp(i1-1,j1).eq.spv) THEN
              dfx(0,0)=dfc
            ELSE
              dfx(0,0)=f(0,0)-Finp(i1-1,j1)
# ifdef LIMTR
              IF ((dfx(0,0)*dfc).lt.0.0_r8) THEN
                dfx(0,0)=0.0_r8
              ELSE IF (ABS(dfx(0,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(0,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,1)-f(1,0)
            IF (j1+2.gt.UBy) THEN
              dfy(1,1)=dfc
            ELSE IF (Finp(i2,j1+2).eq.spv) THEN
              dfy(1,1)=dfc
            ELSE
              dfy(1,1)=Finp(i2,j1+2)-f(1,1)
# ifdef LIMTR
              IF ((dfy(1,1)*dfc).lt.0.0_r8) THEN
                dfy(1,1)=0.0_r8
              ELSEIF (ABS(dfy(1,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(1,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(0,1)-f(0,0)
            IF (j1+2.gt.UBy) THEN
              dfy(0,1)=dfc
            ELSE IF (Finp(i1,j1+2).eq.spv) THEN
              dfy(0,1)=dfc
            ELSE
              dfy(0,1)=Finp(i1,j1+2)-f(0,1)
# ifdef LIMTR
              IF ((dfy(0,1)*dfc).lt.0.0_r8) THEN
                dfy(0,1)=0.0_r8
              ELSE IF (ABS(dfy(0,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(0,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,1)-f(1,0)
            IF (j1-1.lt.1) THEN
              dfy(1,0)=dfc
            ELSE IF (Finp(i2,j1-1).eq.spv) THEN
              dfy(1,0)=dfc
            ELSE
              dfy(1,0)=f(1,0)-Finp(i2,j1-1)
# ifdef LIMTR
              IF ((dfy(1,0)*dfc).lt.0.0_r8) THEN
                dfy(1,0)=0.0_r8
              ELSE IF (ABS(dfy(1,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(1,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(0,1)-f(0,0)
            IF (j1-1.lt.1) THEN
              dfy(0,0)=dfc
            ELSE IF (Finp(i1,j1-1).eq.spv) THEN
              dfy(0,0)=dfc
            ELSE
              dfy(0,0)=f(0,0)-Finp(i1,j1-1)
# ifdef LIMTR
              IF ((dfy(0,0)*dfc).lt.0.0_r8) THEN
                dfy(0,0)=0.0_r8
              ELSEIF (ABS(dfy(0,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(0,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            f0=C05*(f(1,1)+f(1,0)+f(0,1)+f(0,0))-                       &
     &         C02*(dfx(1,1)+dfx(1,0)-dfx(0,1)-dfx(0,0)+                &
     &              dfy(1,1)-dfy(1,0)+dfy(0,1)-dfy(0,0))

            fx=C10*(f(1,1)-f(0,1)+f(1,0)-f(0,0))-                       &
     &         C01*(dfx(1,1)+dfx(1,0)+dfx(0,1)+dfx(0,0))-               &
     &         C03*(dfy(1,1)-dfy(0,1)-dfy(1,0)+dfy(0,0))

            fy=C10*(f(1,1)-f(1,0)+f(0,1)-f(0,0))-                       &
     &         C01*(dfy(1,1)+dfy(0,1)+dfy(1,0)+dfy(0,0))-               &
     &         C03*(dfx(1,1)-dfx(1,0)-dfx(0,1)+dfx(0,0))

            fxy=f(1,1)-f(1,0)-f(0,1)+f(0,0)

            fxx=C05*(dfx(1,1)-dfx(0,1)+dfx(1,0)-dfx(0,0))

            fyy=C05*(dfy(1,1)-dfy(1,0)+dfy(0,1)-dfy(0,0))

            fxxx=C06*(dfx(1,1)+dfx(1,0)+dfx(0,1)+dfx(0,0))-             &
     &           f(1,1)+f(0,1)-f(1,0)+f(0,0)

            fyyy=C06*(dfy(1,1)+dfy(0,1)+dfy(1,0)+dfy(0,0))-             &
     &           f(1,1)+f(1,0)-f(0,1)+f(0,0)

            fxxy=C06*(dfx(1,1)-dfx(0,1)-dfx(1,0)+dfx(0,0))

            fxyy=C06*(dfy(1,1)-dfy(1,0)-dfy(0,1)+dfy(0,0))
#endif
            Fout(i,j)=f0+                                               &
     &                fx*xi+                                            &
     &                fy*eta+                                           &
     &                C06*fxx*xi*xi+                                    &
     &                fxy*xi*eta+                                       &
     &                C06*fyy*eta*eta+                                  &
     &                C04*fxxx*xi*xi*xi+                                &
     &                C06*fxxy*xi*xi*eta+                               &
     &                C04*fyyy*eta*eta*eta+                             &
     &                C06*fxyy*xi*eta*eta
            Fmin=MIN(Fmin,Fout(i,j))
            Fmax=MAX(Fmax,Fout(i,j))
          END IF
        END DO
      END DO
!
!  Return minimum and maximum values.
!
      IF (PRESENT(MinVal)) THEN
        MinVal=Fmin
      END IF
      IF (PRESENT(MaxVal)) THEN
        MaxVal=Fmax
      END IF
!
      RETURN
      END SUBROUTINE cinterp2d
!
      SUBROUTINE hindices (ng, LBi, UBi, LBj, UBj,                      &
     &                     Is, Ie, Js, Je,                              &
     &                     angler, Xgrd, Ygrd,                          &
     &                     LBm, UBm, LBn, UBn,                          &
     &                     Ms, Me, Ns, Ne,                              &
     &                     Xpos, Ypos, Ipos, Jpos,                      &
     &                     IJspv, rectangular)
!
!=======================================================================
!                                                                      !
!  Given any geographical locations Xpos and Ypos, this routine finds  !
!  the corresponding array cell indices (Ipos, Jpos) of gridded  data  !
!  Xgrd and Ygrd containing each requested location. This indices are  !
!  usually used for interpolation.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number.                                  !
!     LBi         I-dimension Lower bound of gridded data.             !
!     UBi         I-dimension Upper bound of gridded data.             !
!     LBj         J-dimension Lower bound of gridded data.             !
!     UBj         J-dimension Upper bound of gridded data.             !
!     Is          Starting gridded data I-index to search.             !
!     Ie          Ending   gridded data I-index to search.             !
!     Js          Starting gridded data J-index to search.             !
!     Je          Ending   gridded data J-index to search.             !
!     angler      Gridded data angle between X-axis and true EAST      !
!                   (radians).                                         !
!     Xgrd        Gridded data X-locations (usually, longitude).       !
!     Ygrd        Gridded data Y-locations (usually, latitude).        !
!     LBm         I-dimension Lower bound of requested locations.      !
!     UBm         I-dimension Upper bound of requested locations.      !
!     LBn         J-dimension Lower bound of requested locations.      !
!     UBn         J-dimension Upper bound of requested locations.      !
!     Ms          Starting requested locations I-index to search.      !
!     Me          Ending   requested locations I-index to search.      !
!     Ns          Starting requested locations J-index to search.      !
!     Ne          Ending   requested locations J-index to search.      !
!     Xpos        Requested X-locations to process (usually longitude).!
!     Ypos        Requested Y-locations to process (usually latitude). !
!     IJspv       Unbounded special value to assign.                   !
!     rectangular Logical switch indicating that gridded data has a    !
!                   plaid distribution.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ipos       Fractional I-cell index containing locations in data. !
!     Jpos       Fractional J-cell index containing locations in data. !
!                                                                      !
!  Calls:    Try_Range                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,  intent(in) :: rectangular
!
      integer,  intent(in) :: ng
      integer,  intent(in) :: LBi, UBi, LBj, UBj, Is, Ie, Js, Je
      integer,  intent(in) :: LBm, UBm, LBn, UBn, Ms, Me, Ns, Ne

      real(r8), intent(in) :: IJspv
!
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Ygrd(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: Xpos(LBm:UBm,LBn:UBn)
      real(r8), intent(in) :: Ypos(LBm:UBm,LBn:UBn)

      real(r8), intent(out) :: Ipos(LBm:UBm,LBn:UBn)
      real(r8), intent(out) :: Jpos(LBm:UBm,LBn:UBn)
!
!  Local variable declarations.
!
      logical :: found, foundi, foundj

      integer :: Imax, Imin, Jmax, Jmin, i, i0, j, j0, mp, np

      real(r8) :: aa2, ang, bb2, diag2, dx, dy, phi
      real(r8) :: xfac, xpp, yfac, ypp
!
!-----------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!-----------------------------------------------------------------------
!
      DO np=Ns,Ne
        DO mp=Ms,Me
          Ipos(mp,np)=IJspv
          Jpos(mp,np)=IJspv
!
!  The gridded data has a plaid distribution so the search is trivial.
!
          IF (rectangular) THEN
            foundi=.FALSE.
            I_LOOP : DO i=LBi,UBi-1
              IF ((Xgrd(i  ,1).le.Xpos(mp,np)).and.                     &
     &            (Xgrd(i+1,1).gt.Xpos(mp,np))) THEN
                Imin=i
                foundi=.TRUE.
                EXIT I_LOOP
              END IF
            END DO I_LOOP
            foundj=.FALSE.
            J_LOOP : DO j=LBj,UBj-1
              IF ((Ygrd(1,j  ).le.Ypos(mp,np)).and.                     &
     &            (Ygrd(1,j+1).gt.Ypos(mp,np))) THEN
                Jmin=j
                foundj=.TRUE.
                EXIT J_LOOP
              END IF
            END DO J_LOOP
            found=foundi.and.foundj
!
!  Check each position to find if it falls inside the whole domain.
!  Once it is stablished that it inside, find the exact cell to which
!  it belongs by successively dividing the domain by a half (binary
!  search).
!
          ELSE
            found=try_range(ng, LBi, UBi, LBj, UBj,                     &
     &                      Xgrd, Ygrd,                                 &
     &                      Is, Ie, Js, Je,                             &
     &                      Xpos(mp,np), Ypos(mp,np))
            IF (found) THEN
              Imin=Is
              Imax=Ie
              Jmin=Js
              Jmax=Je
              DO while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
                IF ((Imax-Imin).gt.1) THEN
                  i0=(Imin+Imax)/2
                  found=try_range(ng, LBi, UBi, LBj, UBj,               &
     &                            Xgrd, Ygrd,                           &
     &                            Imin, i0, Jmin, Jmax,                 &
     &                            Xpos(mp,np), Ypos(mp,np))
                  IF (found) THEN
                    Imax=i0
                  ELSE
                    Imin=i0
                  END IF
                END IF
                IF ((Jmax-Jmin).gt.1) THEN
                  j0=(Jmin+Jmax)/2
                  found=try_range(ng, LBi, UBi, LBj, UBj,               &
     &                            Xgrd, Ygrd,                           &
     &                            Imin, Imax, Jmin, j0,                 &
     &                            Xpos(mp,np), Ypos(mp,np))
                  IF (found) THEN
                    Jmax=j0
                  ELSE
                    Jmin=j0
                  END IF
                END IF
              END DO
              found=(Is.le.Imin).and.(Imin.le.Ie).and.                  &
     &              (Is.le.Imax).and.(Imax.le.Ie).and.                  &
     &              (Js.le.Jmin).and.(Jmin.le.Je).and.                  &
     &              (Js.le.Jmax).and.(Jmax.le.Je)
            END IF
          END IF
!
!  Knowing the correct cell, calculate the exact indices, accounting
!  for a possibly rotated grid.  If spherical, convert all positions
!  to meters first.
!
          IF (found) THEN
            IF (spherical) THEN
              yfac=Eradius*deg2rad
              xfac=yfac*COS(Ypos(mp,np)*deg2rad)
              xpp=(Xpos(mp,np)-Xgrd(Imin,Jmin))*xfac
              ypp=(Ypos(mp,np)-Ygrd(Imin,Jmin))*yfac
            ELSE
              xfac=1.0_r8
              yfac=1.0_r8
              xpp=Xpos(mp,np)-Xgrd(Imin,Jmin)
              ypp=Ypos(mp,np)-Ygrd(Imin,Jmin)
            END IF
!
!  Use Law of Cosines to get cell parallelogram "shear" angle.
!
            diag2=((Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+      &
     &            ((Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
            aa2=((Xgrd(Imin,Jmin)-Xgrd(Imin+1,Jmin))*xfac)**2+          &
     &          ((Ygrd(Imin,Jmin)-Ygrd(Imin+1,Jmin))*yfac)**2
            bb2=((Xgrd(Imin,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+          &
     &          ((Ygrd(Imin,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
            phi=ASIN((diag2-aa2-bb2)/(2.0_r8*SQRT(aa2*bb2)))
!
!  Transform float position into curvilinear coordinates. Assume the
!  cell is rectanglar, for now.
!
            ang=angler(Imin,Jmin)
            dx=xpp*COS(ang)+ypp*SIN(ang)
            dy=ypp*COS(ang)-xpp*SIN(ang)
!
!  Correct for parallelogram.
!
            dx=dx+dy*TAN(phi)
            dy=dy/COS(phi)
!
!  Scale with cell side lengths to translate into cell indices.
!
            dx=MIN(MAX(0.0_r8,dx/SQRT(aa2)),1.0_r8)
            dy=MIN(MAX(0.0_r8,dy/SQRT(bb2)),1.0_r8)
            Ipos(mp,np)=REAL(Imin,r8)+dx
            Jpos(mp,np)=REAL(Jmin,r8)+dy
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE hindices
!
      LOGICAL FUNCTION try_range (ng, LBi, UBi, LBj, UBj, Xgrd, Ygrd,   &
     &                            Imin, Imax, Jmin, Jmax, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given a grided domain with matrix coordinates Xgrd and Ygrd, this   !
!  function finds if the point (Xo,Yo)  is inside the box defined by   !
!  the requested corners (Imin,Jmin) and (Imax,Jmax). It will return   !
!  logical switch  try_range=.TRUE.  if (Xo,Yo) is inside, otherwise   !
!  it will return false.                                               !
!                                                                      !
!  Calls:   inside                                                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: ng, LBi, UBi, LBj, UBj
      integer,  intent(in) :: Imin, Imax, Jmin, Jmax

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Xgrd(LBi:,LBj:)
      real(r8), intent(in) :: Ygrd(LBi:,LBj:)
#else
      real(r8), intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Ygrd(LBi:UBi,LBj:UBj)
#endif

      real(r8), intent(in) :: Xo, Yo
!
!  Local variable declarations.
!
      integer ::  Nb, i, j, shft, ic

      real(r8), dimension(2*(Jmax-Jmin+Imax-Imin)+1) :: Xb, Yb
!
!-----------------------------------------------------------------------
!  Define closed polygon.
!-----------------------------------------------------------------------
!
!  Note that the last point (Xb(Nb),Yb(Nb)) does not repeat first
!  point (Xb(1),Yb(1)).  Instead, in function inside, it is implied
!  that the closing segment is (Xb(Nb),Yb(Nb))-->(Xb(1),Yb(1)). In
!  fact, function inside sets Xb(Nb+1)=Xb(1) and Yb(Nb+1)=Yb(1).
!
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      DO i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      END DO
      shft=1-Jmin+Imax-Imin
      DO j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      END DO
      shft=1+Jmax-Jmin+2*Imax-Imin
      DO i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      END DO
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      DO j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      END DO
!
!-----------------------------------------------------------------------
!  Check if point (Xo,Yo) is inside of the defined polygon.
!-----------------------------------------------------------------------
!
      try_range=inside(Nb, Xb, Yb, Xo, Yo)
!
      RETURN
      END FUNCTION try_range
!
      LOGICAL FUNCTION inside (Nb, Xb, Yb, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in   ) :: Nb

      real(r8), intent(in   ) :: Xo, Yo

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Xb(:), Yb(:)
#else
      real(r8), intent(inout) :: Xb(Nb+1), Yb(Nb+1)
#endif
!
!  Local variable declarations.
!
      logical :: ItIsInside
!
      integer, parameter :: Nstep =128

      integer :: crossings, i, inc, k, kk, nc

      integer, dimension(Nstep) :: Sindex
!
      real(r8) :: dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      ItIsInside=.FALSE.
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0_r8).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              ItIsInside=.TRUE.               ! terminate inner DO-loop
              EXIT
            ELSE IF (((dx1.eq.0.0_r8).and.(Yo.ge.Yb(k  ))).or.          &
     &              ((dx2.eq.0.0_r8).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0_r8).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0_r8)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
        IF (ItIsInside) EXIT                  ! terminate outer DO-loop
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
      IF (crossings.eq.0) THEN
        inside=.FALSE.
      ELSE
        inside=.TRUE.
      END IF
!
      RETURN
      END FUNCTION inside
!
      SUBROUTINE roms_datum_interp_2d (S, fld_src, datum, method,       &
     &                                 ObsVetting)
!
!=======================================================================
!                                                                      !
!  Interpolates a ROMS 2D field to observation locations.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S          Interpolation structure, TYPE(roms_interp_type),      !
!                  initialized elsewhere                               !
!     fld_src    2D source field data (real)                           !
!     method     Interpolation method flag (integer)                   !
!                  BilinearMethod => bilinear interpolation            !
!                  BicubicMethod  => bicubic  interpolation            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     datum      Interpolated values at datum locations (vector).      !
!     ObsVetting Observation screenning flag, 0: reject or 1: accept   !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in   ) :: method
!
      real(r8), intent(in   ) :: fld_src(:,:)
      real(r8), intent(out  ) :: datum(:)

      real(r8), intent(inout), optional :: ObsVetting(:)
!
      TYPE (roms_interp_type), intent(inout) :: S
!
!  Local variable declarations.
!
      integer :: ic, iobs, i1, i2, j1, j2, nlocs
!
      real(r8) :: p1, p2, q1, q2, wsum

      real(r8), dimension(4) :: Hmat
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", roms_datum_interp_2d"
!
!-----------------------------------------------------------------------
!  Interpolates a 2D ROMS state field to observation locations.
!-----------------------------------------------------------------------
!
      nlocs=SIZE(datum, DIM=1)
!
      DO iobs=1,nlocs
        i1=INT(S%x_dst(iobs,1))
        j1=INT(S%y_dst(iobs,1))
        IF (((S%Istr_src.le.i1).and.(i1.lt.S%Iend_src)).and.            &
     &      ((S%Jstr_src.le.j1).and.(j1.lt.S%Jend_src))) THEN
          i2=i1+1
          j2=j1+1
          p2=REAL(i2-i1,r8)*(S%x_dst(iobs,1)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(S%y_dst(iobs,1)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          Hmat(1)=p1*q1
          Hmat(2)=p2*q1
          Hmat(3)=p2*q2
          Hmat(4)=p1*q2
#ifdef MASKING
          Hmat(1)=Hmat(1)*S%mask_src(i1,j1)
          Hmat(2)=Hmat(2)*S%mask_src(i2,j1)
          Hmat(3)=Hmat(3)*S%mask_src(i2,j2)
          Hmat(4)=Hmat(4)*S%mask_src(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4
              Hmat(ic)=Hmat(ic)*wsum
            END DO
          END IF
#endif
          datum(iobs)=Hmat(1)*fld_src(i1,j1)+                          &
     &                Hmat(2)*fld_src(i2,j1)+                          &
     &                Hmat(3)*fld_src(i2,j2)+                          &
     &                Hmat(4)*fld_src(i1,j2)
!
!  Set observation screening flag, 0: reject or 1: accept.
!
          IF (PRESENT(ObsVetting)) THEN
#ifdef MASKING
            IF (wsum.gt.0.0_r8) ObsVetting(iobs)=1.0_r8
#else
            ObsVetting(iobs)=1.0_r8
#endif
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE roms_datum_interp_2d
!
      SUBROUTINE roms_datum_interp_3d (S, fld_src, z_src, z_locs,      &
                                       datum, method, obsVetting)
!
!=======================================================================
!                                                                      !
!  Interpolates a 3D ROMS field to observation horizontal and vertical !
!  locations.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S          Interpolation structure, TYPE(roms_interp_type),      !
!                  initialized elsewhere                               !
!     fld_src    3D source field data (real)                           !
!     z_src      3D source field depth (real; meter; negative)         !
!     z_locs     Datum depth locations (real vector). Negative values  !
!                  are depths in meters while positive values are in   !
!                  vertical fractional coordinates of 1 to N.          !
!     method     Interpolation method flag (integer)                   !
!                  BilinearMethod => bilinear interpolation            !
!                  BicubicMethod  => bicubic  interpolation            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     datum      Interpolated values at observation locations (vector) !
!     ObsVetting Observation screenning flag, 0: reject or 1: accept   !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent( in) :: method
!
      real(r8), intent( in) :: fld_src(:,:,:)
      real(r8), intent( in) :: z_src(:,:,:)
      real(r8), intent( in) :: z_locs(:)
      real(r8), intent(out) :: datum(:)

      real(r8), intent(inout), optional :: ObsVetting(:)
!
      TYPE (roms_interp_type), intent(inout) :: S
!
!  Local variable declarations.
!
      integer :: ic, iobs, i1, i2, j1, j2, k, k1, k2, nlevs, nlocs
!
      real(r8) :: Zbot, Ztop, dz, p1, p2, q1, q2, r1, r2
      real(r8) :: w11, w12, w21, w22, wsum

      real(r8), dimension(8) :: Hmat
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", roms_datum_interp_3d"
!
!-----------------------------------------------------------------------
!  Interpolate from 3D ROMS state field at observation locations.
!-----------------------------------------------------------------------
!
      nlevs=SIZE(fld_src, DIM=3)
      nlocs=SIZE(datum, DIM=1)
!
      DO iobs=1,nlocs
        i1=INT(S%x_dst(iobs,1))
        j1=INT(S%y_dst(iobs,1))
        IF (((S%Istr_src.le.i1).and.(i1.lt.S%Iend_src)).and.            &
     &      ((S%Jstr_src.le.j1).and.(j1.lt.S%Jend_src))) THEN
          i2=i1+1
          j2=j1+1
          p2=REAL(i2-i1,r8)*(S%x_dst(iobs,1)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(S%y_dst(iobs,1)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          w11=p1*q1
          w21=p2*q1
          w22=p2*q2
          w12=p1*q2
          IF (z_locs(iobs).gt.0.0_r8) THEN
            k1=MAX(1,INT(z_locs(iobs)))        ! Positions in fractional
            k2=MIN(INT(z_locs(iobs))+1,nlevs)  ! levels
            r2=REAL(k2-k1,r8)*(z_locs(iobs)-REAL(k1,r8))
            r1=1.0_r8-r2
          ELSE
            Ztop=z_src(i1,j1,nlevs)
            Zbot=z_src(i1,j1,1    )
            IF (z_locs(iobs).ge.Ztop) THEN
              r1=0.0_r8                        ! If shallower, ignore.
              r2=0.0_r8
              IF (PRESENT(ObsVetting)) ObsVetting(iobs)=0.0_r8
            ELSE IF (Zbot.ge.z_locs(iobs)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
              IF (PRESENT(ObsVetting)) ObsVetting(iobs)=0.0_r8
            ELSE
              DO k=nlevs,2,-1                  ! Otherwise, interpolate
                Ztop=z_src(i1,j1,k  )          ! to fractional level
                Zbot=z_src(i1,j1,k-1)
                IF ((Ztop.gt.z_locs(iobs)).and.                         &
                    (z_locs(iobs).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=z_src(i1,j1,k2)-z_src(i1,j1,k1)
              r2=(z_locs(iobs)-z_src(i1,j1,k1))/dz
              r1=1.0_r8-r2
            END IF
          END IF
          IF ((r1+r2).gt.0.0_r8) THEN
            Hmat(1)=w11*r1
            Hmat(2)=w21*r1
            Hmat(3)=w22*r1
            Hmat(4)=w12*r1
            Hmat(5)=w11*r2
            Hmat(6)=w21*r2
            Hmat(7)=w22*r2
            Hmat(8)=w12*r2
#ifdef MASKING
            Hmat(1)=Hmat(1)*S%mask_src(i1,j1)
            Hmat(2)=Hmat(2)*S%mask_src(i2,j1)
            Hmat(3)=Hmat(3)*S%mask_src(i2,j2)
            Hmat(4)=Hmat(4)*S%mask_src(i1,j2)
            Hmat(5)=Hmat(5)*S%mask_src(i1,j1)
            Hmat(6)=Hmat(6)*S%mask_src(i2,j1)
            Hmat(7)=Hmat(7)*S%mask_src(i2,j2)
            Hmat(8)=Hmat(8)*S%mask_src(i1,j2)
            wsum=0.0_r8
            DO ic=1,8
              wsum=wsum+Hmat(ic)
            END DO
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              DO ic=1,8
                Hmat(ic)=Hmat(ic)*wsum
              END DO
            END IF
#endif
            datum(iobs)=Hmat(1)*fld_src(i1,j1,k1)+                      &
     &                  Hmat(2)*fld_src(i2,j1,k1)+                      &
     &                  Hmat(3)*fld_src(i2,j2,k1)+                      &
     &                  Hmat(4)*fld_src(i1,j2,k1)+                      &
     &                  Hmat(5)*fld_src(i1,j1,k2)+                      &
     &                  Hmat(6)*fld_src(i2,j1,k2)+                      &
     &                  Hmat(7)*fld_src(i2,j2,k2)+                      &
     &                  Hmat(8)*fld_src(i1,j2,k2)
!
!  Set observation screening flag, 0: reject or 1: accept.
!
            IF (PRESENT(ObsVetting)) THEN
#ifdef MASKING
              IF (wsum.gt.0.0_r8) ObsVetting(iobs)=1.0_r8
#else
              ObsVetting(iobs)=1.0_r8
#endif
#ifndef ALLOW_BOTTOM_OBS
!
!  Reject observations that lie in the lower bottom grid cell (k=1) to
!  avoid clustering due shallowing of bathymetry during smoothing and
!  coarse level half-thickness (-h < z_locs < z_src(:,:,1)) in deep
!  water.
!
              IF ((z_locs(iobs).gt.0.0_r8).and.                         &
                  (z_locs(iobs).le.1.0_r8)) THEN
                ObsVetting(iobs)=0.0_r8
              END IF
#endif
            END IF
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE roms_datum_interp_3d
!
      SUBROUTINE roms_datum_column_interp (S, fld_src, datum, method)
!
!=======================================================================
!                                                                      !
!  Horizontally interpolates a 3D ROMS field to datum locations        !
!  level-by-level. The vertical interpolation will be carried out      !
!  elsewhere from the full column of values.  This is what is          !
!  required by the JEDI-UFO operator.                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S          Interpolation structure, TYPE(roms_interp_type),      !
!                  initialized elsewhere.                              !
!     fld_src    3D field source data (real).                          !
!     method     Interpolation method flag (integer)                   !
!                  BilinearMethod => bilinear interpolation            !
!                  BicubicMethod  => bicubic  interpolation            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     datum      Interpolated values at datum locations (matrix):      !
!                  datum(nlocs,nlevs)                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (roms_interp_type), intent(inout) :: S
!
      real(r8), intent( in) :: fld_src(:,:,:)
      real(r8), intent(out) :: datum(:,:)
!
      integer,  intent( in) :: method
!
!  Local variable declarations.
!
      integer :: ic, iobs, i1, i2, j1, j2, k, nlevs, nlocs
!
      real(r8) :: p1, p2, q1, q2, wsum

      real(r8), dimension(4) :: Hmat
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", roms_datum_interp_2d"
!
!-----------------------------------------------------------------------
!  Interpolate from 2D state field to datum locations.
!-----------------------------------------------------------------------
!
      nlevs=SIZE(fld_src, DIM=3)
      nlocs=SIZE(datum, DIM=1)
!
      DO iobs=1,nlocs
        i1=INT(S%x_dst(iobs,1))
        j1=INT(S%y_dst(iobs,1))
        IF (((S%Istr_src.le.i1).and.(i1.lt.S%Iend_src)).and.            &
     &      ((S%Jstr_src.le.j1).and.(j1.lt.S%Jend_src))) THEN
          i2=i1+1
          j2=j1+1
          p2=REAL(i2-i1,r8)*(S%x_dst(iobs,1)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(S%y_dst(iobs,1)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          Hmat(1)=p1*q1
          Hmat(2)=p2*q1
          Hmat(3)=p2*q2
          Hmat(4)=p1*q2
#ifdef MASKING
          Hmat(1)=Hmat(1)*S%mask_src(i1,j1)
          Hmat(2)=Hmat(2)*S%mask_src(i2,j1)
          Hmat(3)=Hmat(3)*S%mask_src(i2,j2)
          Hmat(4)=Hmat(4)*S%mask_src(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4
              Hmat(ic)=Hmat(ic)*wsum
            END DO
          END IF
#endif
          DO k=1,nlevs
            datum(iobs,k)=Hmat(1)*fld_src(i1,j1,k)+                    &
     &                    Hmat(2)*fld_src(i2,j1,k)+                    &
     &                    Hmat(3)*fld_src(i2,j2,k)+                    &
     &                    Hmat(4)*fld_src(i1,j2,k)
          END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE roms_datum_column_interp
!
      SUBROUTINE roms_horiz_interp_2d (S, fld_src, fld_dst, method)
!
!=======================================================================
!                                                                      !
!  Interpolate a ROMS 2D source field to a new destination locations.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S          Interpolation structure, TYPE(roms_interp_type),      !
!                  initialized elsewhere.                              !
!     fld_src    2D field source data (real).                          !
!     method     Interpolation method flag (integer)                   !
!                  BilinearMethod => bilinear interpolation            !
!                  BicubicMethod  => bicubic  interpolation            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     fld_dst    2D interpolated field data (real).                    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (roms_interp_type), intent(inout) :: S
!
      real(r8), intent(in ) :: fld_src(:,:)
      real(r8), intent(out) :: fld_dst(:,:)
!
      integer,  intent(in ) :: method
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", roms_horiz_interp_2d"
!
!-----------------------------------------------------------------------
!  Interpolate 2D field.
!-----------------------------------------------------------------------
!
      SELECT CASE (method)
        CASE (BilinearMethod)
          CALL linterp2d (S%ng,                                         &
     &                    S%LBi_src, S%UBi_src, S%LBj_src, S%UBj_src,   &
     &                    S%lon_src, S%lat_src, fld_src,                &
     &                    S%LBi_dst, S%UBi_dst, S%LBj_dst, S%UBj_dst,   &
     &                    S%Istr_dst, S%Iend_dst,                       &
     &                    S%Jstr_dst, S%Jend_dst,                       &
     &                    S%x_dst, S%y_dst, S%lon_dst, S%lat_dst,       &
     &                    fld_dst,                                      &
     &                    MinVal = S%min_dst,                           &
     &                    MaxVal = S%max_dst)
        CASE (BicubicMethod)
          CALL cinterp2d (S%ng,                                         &
     &                    S%LBi_src, S%UBi_src, S%LBj_src, S%UBj_src,   &
     &                    S%lon_src, S%lat_src, fld_src,                &
     &                    S%LBi_dst, S%UBi_dst, S%LBj_dst, S%UBj_dst,   &
     &                    S%Istr_dst, S%Iend_dst,                       &
     &                    S%Jstr_dst, S%Jend_dst,                       &
     &                    S%x_dst, S%y_dst, S%lon_dst, S%lat_dst,       &
     &                    fld_dst,                                      &
     &                    MinVal = S%min_dst,                           &
     &                    MaxVal = S%max_dst)
        CASE DEFAULT
          IF (Master) WRITE(stdout,10) method
          exit_flag=5
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' ROMS_HORIZ_INTERP_2D - Illegal interpolation method =', &
     &        i0)
!
      RETURN
      END SUBROUTINE roms_horiz_interp_2d
!
      SUBROUTINE roms_horiz_interp_3d (S, fld_src, fld_dst, method)
!
!=======================================================================
!                                                                      !
!  Interpolate a ROMS 3D source field to a new destination locations   !
!  level-by-level. There is no vertical interpolation.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S          Interpolation structure, TYPE(roms_interp_type),      !
!                  initialized elsewhere.                              !
!     fld_src    3D field source data (real).                          !
!     method     Interpolation method flag (integer)                   !
!                  BilinearMethod => bilinear interpolation            !
!                  BicubicMethod  => bicubic  interpolation            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     fld_dst    3D interpolated field data (real).                    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (roms_interp_type), intent(inout) :: S
!
      real(r8), intent( in) :: fld_src(:,:,:)
      real(r8), intent(out) :: fld_dst(:,:,:)
!
      integer,  intent( in) :: method
!
!  Local variable declarations.
!
      integer :: k, nk
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", roms_horiz_interp_3d"
!
!-----------------------------------------------------------------------
!  Interpolate 3D field level-by-level.
!-----------------------------------------------------------------------
!
      nk=SIZE(fld_src, DIM=3)!
!
      SELECT CASE (method)
        CASE (BilinearMethod)
          DO k=1,nk
            CALL linterp2d (S%ng,                                       &
     &                      S%LBi_src, S%UBi_src, S%LBj_src, S%UBj_src, &
     &                      S%lon_src, S%lat_src, fld_src(:,:,k),       &
     &                      S%LBi_dst, S%UBi_dst, S%LBj_dst, S%UBj_dst, &
     &                      S%Istr_dst, S%Iend_dst,                     &
     &                      S%Jstr_dst, S%Jend_dst,                     &
     &                      S%x_dst, S%y_dst, S%lon_dst, S%lat_dst,     &
     &                      fld_dst(:,:,k),                             &
     &                      MinVal = S%min_dst,                         &
     &                      MaxVal = S%max_dst)
          END DO
        CASE (BicubicMethod)
          DO k=1,nk
            CALL cinterp2d (S%ng,                                       &
     &                      S%LBi_src, S%UBi_src, S%LBj_src, S%UBj_src, &
     &                      S%lon_src, S%lat_src, fld_src(:,:,k),       &
     &                      S%LBi_dst, S%UBi_dst, S%LBj_dst, S%UBj_dst, &
     &                      S%Istr_dst, S%Iend_dst,                     &
     &                      S%Jstr_dst, S%Jend_dst,                     &
     &                      S%x_dst, S%y_dst, S%lon_dst, S%lat_dst,     &
     &                      fld_dst(:,:,k),                             &
     &                      MinVal = S%min_dst,                         &
     &                      MaxVal = S%max_dst)
          END DO
        CASE DEFAULT
          IF (Master) WRITE(stdout,10) method
          exit_flag=5
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' ROMS_HORIZ_INTERP_3D - Illegal interpolation method =', &
     &        i0)
!
      RETURN
      END SUBROUTINE roms_horiz_interp_3d
!
      SUBROUTINE roms_interp_create (ng, tile, gtype, S)
!
!=======================================================================
!                                                                      !
!  Creates ROMS intepolation Object.  The interpolation structure has  !
!  two components: source and destination. Here, the source component  !
!  is initialized. The destination component is initialized elsewhere  !
!  for flexibility.                                                    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype
!
      TYPE (roms_interp_type), intent(inout) :: S
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Creates source component of interpolation object.
!-----------------------------------------------------------------------
!
!  If applicable, destroy interpolation object.
!
      IF (allocated(S%lon_src)) THEN
        CALL roms_interp_delete (S)
      END IF
!
!  Allocate source component arrays.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      allocate ( S%lon_src(LBi:UBi,LBj:UBj) )
      allocate ( S%lat_src(LBi:UBi,LBj:UBj) )
      allocate ( S%angle_src(LBi:UBi,LBj:UBj) )
#ifdef MASKING
      allocate ( S%mask_src(LBi:UBi,LBj:UBj) )
#endif
!
!  Associate and assign values to source component.
!
      S%ng = ng
      S%model = iNLM
      S%spval = 0.0_r8
!
      SELECT CASE (gtype)
        CASE (r2dvar, r3dvar)
          S%lon_src   = GRID(ng)%lonr
          S%lat_src   = GRID(ng)%latr
          S%angle_src = GRID(ng)%angler
#ifdef MASKING
          S%mask_src  = GRID(ng)%rmask
#endif
          S%Istr_src  = BOUNDS(ng)%IstrR(tile)
          S%Iend_src  = BOUNDS(ng)%IendR(tile)
          S%Jstr_src  = BOUNDS(ng)%JstrR(tile)
          S%Jend_src  = BOUNDS(ng)%JendR(tile)
        CASE (u2dvar, u3dvar)
          S%lon_src   = GRID(ng)%lonu
          S%lat_src   = GRID(ng)%latu
          S%angle_src = GRID(ng)%angler
#ifdef MASKING
          S%mask_src  = GRID(ng)%umask
#endif
          S%Istr_src  = BOUNDS(ng)%Istr (tile)
          S%Iend_src  = BOUNDS(ng)%IendR(tile)
          S%Jstr_src  = BOUNDS(ng)%JstrR(tile)
          S%Jend_src  = BOUNDS(ng)%JendR(tile)
        CASE (v2dvar, v3dvar)
          S%lon_src   = GRID(ng)%lonv
          S%lat_src   = GRID(ng)%latv
          S%angle_src = GRID(ng)%angler
#ifdef MASKING
          S%mask_src  = GRID(ng)%vmask
#endif
          S%Istr_src  = BOUNDS(ng)%IstrR(tile)
          S%Iend_src  = BOUNDS(ng)%IendR(tile)
          S%Jstr_src  = BOUNDS(ng)%Jstr (tile)
          S%Jend_src  = BOUNDS(ng)%JendR(tile)
      END SELECT
!
!
      S%LBi_src = LBi
      S%UBi_src = UBi
      S%LBj_src = LBj
      S%UBj_src = UBj
!
      RETURN
      END SUBROUTINE roms_interp_create
!
      SUBROUTINE roms_interp_delete (S)
!
!=======================================================================
!                                                                      !
!  Deallocates pointer arrays in the intepolation structure, S.        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (roms_interp_type), intent(inout) :: S
!
!-----------------------------------------------------------------------
!  Deallocate interpolation structure.
!-----------------------------------------------------------------------
!
      IF (allocated(S%angle_src)) deallocate (S%angle_src)
      IF (allocated(S%lon_src))   deallocate (S%lon_src)
      IF (allocated(S%lat_src))   deallocate (S%lat_src)
      IF (allocated(S%mask_src))  deallocate (S%mask_src)

      IF (allocated(S%lon_dst))   deallocate (S%lon_dst)
      IF (allocated(S%lat_dst))   deallocate (S%lat_dst)
      IF (allocated(S%mask_dst))  deallocate (S%mask_dst)
      IF (allocated(S%x_dst))     deallocate (S%x_dst)
      IF (allocated(S%y_dst))     deallocate (S%y_dst)
!
      RETURN
      END SUBROUTINE roms_interp_delete
!
      SUBROUTINE roms_interp_fractional (S)
!
!=======================================================================
!                                                                      !
!  Computes the horizontal fractional coordinates (S%x_dst, S%y_dst)   !
!  of the source cells containing the destination values.              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      TYPE (roms_interp_type), intent(inout) :: S
!
!  Local variable declarations.
!
      integer :: Npts
!
!-----------------------------------------------------------------------
!  Determine the horizontal fractional coordinates: S%x_dst and S%y_dst.
!-----------------------------------------------------------------------
!
      CALL hindices (S%ng,                                              &
     &               S%LBi_src,   S%UBi_src,  S%LBj_src,  S%UBj_src,    &
     &               S%Istr_src , S%Iend_src, S%Jstr_src, S%Jend_src,   &
     &               S%angle_src, S%lon_src,  S%lat_src,                &
     &               S%LBi_dst,   S%UBi_dst,  S%LBj_dst,  S%UBj_dst,    &
     &               S%Istr_dst,  S%Iend_dst, S%Jstr_dst, S%Jend_dst,   &
     &               S%lon_dst,   S%lat_dst,  S%x_dst,    S%y_dst,      &
     &               S%spval,     S%rectangular)

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect fractional coordinates.
!-----------------------------------------------------------------------
!
      Npts=SIZE(S%x_dst)
      CALL mp_assemble (S%ng, S%model, Npts, S%spval, S%x_dst)
      CALL mp_assemble (S%ng, S%model, Npts, S%spval, S%y_dst)
#endif
!
      RETURN
      END SUBROUTINE roms_interp_fractional
!
      END MODULE roms_interpolate_mod
