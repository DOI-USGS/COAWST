#include "cppdefs.h"
      SUBROUTINE ludcmp (a, n, np, indx, d)
!
!svn $Id: ludcmp.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Given an N x N matrix A, with physical dimension NP, this routine   !
!  replaces it by the  LU  decomposition of a rowwise permutation of   !
!  itself.   A and N are input.   A is output, arranged according to   !
!  Crout algorithm;   INDX is an output vector which records the row   !
!  permutation effected by the partial pivoting;  D is output as  +1   !
!  or  -1  depending on  whether the  number of row interchanges was   !
!  even or odd,  respectively.  IER is output as 1 or 0 depending on   !
!  whether the matrix A is singular or not, respectively. It is used   !
!  in  combination with LUBKSB to solve linear equations or invert a   !
!  matrix.                                                             !
!                                                                      !
!  Reference:   Press, W.H, et al., 1989: Numerical Recipes, The Art   !
!               of Scientific Computing, pp 31-37.                     !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
!  Imported variable declarations.
!
      integer, intent(in) :: n, np

      integer, intent(out) :: indx(n)

      real(r8), intent(inout) :: a(np,np)
      real(r8), intent(out) :: d
!
!  Local variable declarations.
!
      integer i, ier, imax, j, k

      real(r8), parameter :: tiny = 1.0E-20_r8

      real(r8), dimension(n) :: vv

      real(r8) :: aamax, dum, MySum
!
!-----------------------------------------------------------------------
!  Replace matrix A by its LU decomposition.
!-----------------------------------------------------------------------
!
      ier=0
      d=1.0_r8
      DO i=1,n
        aamax=0.0_r8
        DO j=1,n
          IF (ABS(a(i,j)).gt.aamax) aamax=ABS(a(i,j))
        END DO
        IF (aamax.eq.0.0_r8) THEN
          ier=1
        ELSE
          vv(i)=1.0_r8/aamax
        END IF
      END DO
      IF (ier.eq.1) RETURN
      DO j=1,n
        IF (j.gt.1) THEN
          DO i=1,j-1
            MySum=a(i,j)
            IF (i.gt.1) THEN
              DO k=1,i-1
                MySum=MySum-a(i,k)*a(k,j)
              END DO
              a(i,j)=MySum
            END IF
          END DO
        END IF
        aamax=0.0_r8
        DO i=j,n
          MySum=a(i,j)
          IF (j.gt.1) THEN
            DO k=1,j-1
              MySum=MySum-a(i,k)*a(k,j)
            END DO
            a(i,j)=MySum
          END IF
          dum=vv(i)*ABS(MySum)
          IF (dum.ge.aamax) THEN
            imax=i
            aamax=dum
          END IF
        END DO
        IF (j.ne.imax) THEN
          DO k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          END DO
          d=-d
          vv(imax)=vv(j)
        END IF
        indx(j)=imax
        IF (j.ne.n) THEN
          IF (a(j,j).eq.0.0_r8) a(j,j)=tiny
          dum=1./a(j,j)
          DO i=j+1,n
            a(i,j)=a(i,j)*dum
          END DO
        END IF
      END DO
      IF (a(n,n).eq.0.0_r8) a(n,n)=tiny

      RETURN
      END SUBROUTINE ludcmp
