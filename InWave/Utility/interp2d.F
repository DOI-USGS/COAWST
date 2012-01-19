#include "cppdefs.h"
!
!************************************************************************
      SUBROUTINE interp2d (LBx, UBx, LBy, UBy,                          &
     &                     Xinp, Yinp, Finp,                            &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Istr, Iend, Jstr, Jend,                      &
     &                     Xout, Yout,                                  &
     &                     Fout)
!***********************************************************************
!
#ifdef INWAVE_MODEL
!
!svn $Id: swan_reader.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010
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
!     Fmin       Interpolated field minimum value.                     !
!     Fmax       Interpolated field maximum value.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: LBx, UBx, LBy, UBy
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(r8), intent(in) :: Xinp(LBx:UBx)
      real(r8), intent(in) :: Yinp(LBy:UBy)
      real(r8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Xout(LBi:UBi)
      real(r8), intent(in) :: Yout(LBj:UBj)
      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, i1, i2, j, j1, j2
      integer:: Iout(LBi:UBi)
      integer:: Jout(LBj:UBj)
      real(r8):: p1, p2, q1, q2
      real(r8):: cff7, cff8, cff9
!
!-----------------------------------------------------------------------
!  Linearly interpolate requested field
!-----------------------------------------------------------------------
!
      DO i=Istr,Iend
       IF((Xout(i).lt.Xinp(LBx)))THEN
        Iout(i)=0
       ELSE
        IF((Xout(i).gt.Xinp(UBx-1)))THEN
         Iout(i)=0
        ELSE
         DO j=LBx,UBx-1
          IF((Xout(i).gt.Xinp(j)).and.(Xout(i).le.Xinp(j+1)))THEN
           Iout(i)=j
          ENDIF
         ENDDO
        ENDIF
       ENDIF
      ENDDO

      DO i=Jstr,Jend
       IF((Yout(i).lt.Yinp(LBy)))THEN
        Jout(i)=0
       ELSE
        IF((Yout(i).gt.Yinp(UBy-1)))THEN
         Jout(i)=0
        ELSE
         DO j=LBy,UBy-1
          IF((Yout(i).gt.Yinp(j)).and.(Yout(i).le.Yinp(j+1)))THEN
           Jout(i)=j
          ENDIF
         ENDDO
        ENDIF
       ENDIF
      ENDDO

      DO j=Jstr,Jend
       DO i=Istr,Iend
        IF(Iout(i).ne.0.and.Jout(j).ne.0)THEN
         i1=Iout(i)
         i2=i1+1
         j1=Jout(j)
         j2=j1+1
         p2=(Xinp(i2)-Xout(i))/(Xinp(i2)-Xinp(i1))
         q2=(Yinp(j2)-Yout(j))/(Yinp(j2)-Yinp(j1))
         p1=1.0_r8-p2
         q1=1.0_r8-q2
         Fout(i,j)=p1*q1*Finp(i1,j1)+                                   &
     &             p2*q1*Finp(i2,j1)+                                   &
     &             p2*q2*Finp(i2,j2)+                                   &
     &             p1*q2*Finp(i1,j2)
        ELSE
         Fout(i,j)=0.0_r8
        ENDIF
       END DO
      END DO
#endif
      RETURN
      END SUBROUTINE interp2d