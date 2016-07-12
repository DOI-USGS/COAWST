#include "cppdefs.h"
      SUBROUTINE lubksb (a, n, np, indx, b)
!
!svn $Id: lubksb.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Solves the set of N linear equations  A X = B. Here  A is input,    !
!  not as the matrix A  but rather as its LU decomposition,  set by    !
!  routine LUDCMP. INDX is input as the permutation vector returned    !
!  by  LUDCMP.  B is input as  the  right-hand  side vector B,  and    !
!  returns with the solution vector X. A,N,NP,INDX are not modified    !
!  by this  routine and  can be left in place for  successive calls    !
!  with  different  right-hand sides  B.  This  routine  takes into    !
!  account  the  possiblility  that  B  will  begin with  many zero    !
!  elements, so is efficient for use in matrix inversion.              !
!                                                                      !
!  Reference:  Press, W.H, et al., 1989: Numerical Recipes, The Art    !
!              of Scientific Computing, pp 31-37.                      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
!  Imported variable declarations.
!
      integer, intent(in) :: n, np

      integer, intent(in) :: indx(n)

      real(r8), intent(in) :: a(np,np)
      real(r8), intent(inout) :: b(n)
!
!  Local variable declarations.
!
      integer :: i, ii, j, ll

      real(r8) :: MySum
!
!-----------------------------------------------------------------------
!  Solve set of linear equation by LU decomposition.
!-----------------------------------------------------------------------
!
      ii=0
      DO i=1,n
        ll=indx(i)
        MySum=b(ll)
        b(ll)=b(i)
        IF (ii.ne.0) THEN
          DO j=ii,i-1
            MySum=MySum-a(i,j)*b(j)
          END DO
        ELSE IF (MySum.ne.0.0_r8) THEN
          ii=i
        END IF
        b(i)=MySum
      END DO
      DO i=n,1,-1
        MySum=b(i)
        IF (i.lt.n) THEN
          DO j=i+1,n
            MySum=MySum-a(i,j)*b(j)
          END DO
        END IF
        b(i)=MySum/a(i,i)
      END DO

      RETURN
      END SUBROUTINE lubksb
