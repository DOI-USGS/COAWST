#include "cppdefs.h"

#if defined WEAK_CONSTRAINT && defined MINRES
        SUBROUTINE sqlq (innLoop, a, tau, y)
!
!svn $Id: sqlq.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs an LQ factorization of the square matrix A.   !
!  The algorithm is based on that used in the LAPACK routine DGELQF.   !
!                                                                      !
!  NOTE: A, TAU and Y are overwritten on exit.                         !
!                                                                      !
!  On exit, the elements on and below the diagonal of the array        !
!  contain the innLoop-by-innLoop lower trapezoidal matrix L (L is     !
!  lower triangular; the elements above the diagonal,                  !
!  with the array TAU, represent the orthogonal matrix Q as a          !
!  product of elementary reflectors.                                   !
!  The matrix Q is represented as a product of elementary reflectors   !
!                                                                      !
!     Q = H(innLoop) . . . H(2) H(1)                                   !
!                                                                      !
!  Each H(i) has the form                                              !
!                                                                      !
!     H(i) = I - tau * v * v'                                          !
!                                                                      !
!  where tau is a real scalar, and v is a real vector with             !
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:innLoop) is stored on exit in      !
!  A(i,i+1:innLoop), and tau in TAU(i).                                !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: innLoop

      real(r8), dimension(innLoop,innLoop), intent(inout) :: a
      real(r8), dimension(innLoop), intent(inout) :: tau, y
!
!  Local variable declarations.
!
      integer :: i, j, ii, jj

      real(r8) :: znorm, zbeta, zaii, ztemp
!
!-----------------------------------------------------------------------
!  LQ factorization of a square matrix.
!-----------------------------------------------------------------------
!
      DO i=1,innLoop
        tau(i)=0.0_r8
      END DO
!
      DO ii=1,innLoop-1
!
!  Generate elementary reflector H(ii) to annihilate A(ii,ii+1:innLoop).
!
        znorm=0.0_r8
        DO j=ii+1,innLoop
          znorm=znorm+a(ii,j)*a(ii,j)
        END DO
        znorm=SQRT(znorm)
        zbeta=SQRT(znorm*znorm+a(ii,ii)*a(ii,ii))
        zbeta=-SIGN(zbeta,a(ii,ii))
        tau(ii)=(zbeta-a(ii,ii))/zbeta
        DO j=ii+1,innLoop
          a(ii,j)=a(ii,j)/(a(ii,ii)-zbeta)
        END DO
        a(ii,ii)=zbeta
!
!  Apply H(ii) to A(ii+1:innLoop,ii:innLoop) from the right.
!
        zaii=a(ii,ii)
        a(ii,ii)=1.0_r8
        IF (tau(ii).ne.0.0_r8) THEN
          DO i=1,innLoop
            y(i)=0.0_r8
          END DO
          jj=ii-1
          DO j=1,innLoop-jj
            ztemp=a(ii,j+jj)
            DO i=1,innLoop-ii
              y(i)=y(i)+ztemp*a(ii+i,j+jj)
            END DO
          END DO
!
          DO j=1,innLoop-jj
            ztemp=-tau(ii)*a(ii,j+jj)
            DO i=1,innLoop-ii
              a(ii+i,j+jj)=a(ii+i,j+jj)+y(i)*ztemp
            END DO
          END DO
        END IF
        a(ii,ii)=zaii
      END DO
!
      RETURN
      END SUBROUTINE sqlq
#else
      SUBROUTINE sqlq
      RETURN
      END SUBROUTINE sqlq
#endif


