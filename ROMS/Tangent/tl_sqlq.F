#include "cppdefs.h"

#if (defined SENSITIVITY_4DVAR || defined TL_W4DPSAS || \
     defined TL_W4DVAR) && defined MINRES

        SUBROUTINE tl_sqlq(innLoop, a, tl_a, tau, tl_tau, y, tl_y)
!
!svn $Id: tl_sqlq.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs the tangent linearization of a LQ             !
!  factorization of the square matrix A.                               !
!                                                                      !
!  NOTE: A, tl_A, TAU, tl_TAU, Y and tl_Y are overwritten on exit.     !
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

      real(r8), dimension(innLoop,innLoop), intent(inout) :: a, tl_a
      real(r8), dimension(innLoop), intent(inout) :: tau, tl_tau
      real(r8), dimension(innLoop), intent(inout) :: y, tl_y
!
!  Local variable declarations.
!
      integer :: i, j, ii, jj

      real(r8) :: znorm, zbeta, zaii, ztemp
      real(r8) :: tl_znorm, tl_zbeta, tl_zaii, tl_ztemp
!
!-----------------------------------------------------------------------
!  Tangent linearization of a LQ factorization of a square matrix.
!-----------------------------------------------------------------------
!
      DO i=1,innLoop
        tau(i)=0.0_r8
        tl_tau(i)=0.0_r8
      END DO
!
      DO ii=1,innLoop-1
!
!  Generate elementary reflector H(ii) to annihilate A(ii,ii+1:innLoop).
!
        znorm=0.0_r8
        tl_znorm=0.0_r8
        DO j=ii+1,innLoop
          znorm=znorm+a(ii,j)*a(ii,j)
          tl_znorm=tl_znorm+2.0_r8*a(ii,j)*tl_a(ii,j)
        END DO
        znorm=SQRT(znorm)
        IF (znorm.NE.0.0_r8) THEN
          tl_znorm=0.5_r8*tl_znorm/znorm
        ELSE
          tl_znorm=0.0_r8
        END IF
        zbeta=SQRT(znorm*znorm+a(ii,ii)*a(ii,ii))
        IF (zbeta.NE.0.0_r8) THEN
          tl_zbeta=(tl_znorm*znorm+tl_a(ii,ii)*a(ii,ii))/zbeta
        ELSE
          tl_zbeta=0.0_r8
        END IF
!>      zbeta=-SIGN(zbeta,a(ii,ii))
!>
        tl_zbeta=-SIGN(1.0_r8,zbeta)*SIGN(1.0_r8,a(ii,ii))*tl_zbeta
        zbeta=-SIGN(zbeta,a(ii,ii))
        tau(ii)=(zbeta-a(ii,ii))/zbeta
        tl_tau(ii)=(tl_zbeta-tl_a(ii,ii))/zbeta-tl_zbeta*tau(ii)/zbeta
        DO j=ii+1,innLoop
          a(ii,j)=a(ii,j)/(a(ii,ii)-zbeta)
          tl_a(ii,j)=tl_a(ii,j)/(a(ii,ii)-zbeta)-                       &
     &               (tl_a(ii,ii)-tl_zbeta)*a(ii,j)/(a(ii,ii)-zbeta)
        END DO
        a(ii,ii)=zbeta
        tl_a(ii,ii)=tl_zbeta
!
!  Apply H(ii) to A(ii+1:innLoop,ii:innLoop) from the right.
!
        zaii=a(ii,ii)
        tl_zaii=tl_a(ii,ii)
        a(ii,ii)=1.0_r8
        tl_a(ii,ii)=0.0_r8
        IF (tau(ii).ne.0.0_r8) THEN
          DO i=1,innLoop
            y(i)=0.0_r8
            tl_y(i)=0.0_r8
          END DO
          jj=ii-1
          DO j=1,innLoop-jj
            ztemp=a(ii,j+jj)
            tl_ztemp=tl_a(ii,j+jj)
            DO i=1,innLoop-ii
              y(i)=y(i)+ztemp*a(ii+i,j+jj)
              tl_y(i)=tl_y(i)+tl_ztemp*a(ii+i,j+jj)+                    &
     &                        ztemp*tl_a(ii+i,j+jj)
            END DO
          END DO
!
          DO j=1,innLoop-jj
            ztemp=-tau(ii)*a(ii,j+jj)
            tl_ztemp=-tl_tau(ii)*a(ii,j+jj)-tau(ii)*tl_a(ii,j+jj)
            DO i=1,innLoop-ii
              a(ii+i,j+jj)=a(ii+i,j+jj)+y(i)*ztemp
              tl_a(ii+i,j+jj)=tl_a(ii+i,j+jj)+tl_y(i)*ztemp+            &
     &                                        y(i)*tl_ztemp
            END DO
          END DO
        END IF
        a(ii,ii)=zaii
        tl_a(ii,ii)=tl_zaii
      END DO
!
      RETURN
      END SUBROUTINE tl_sqlq
#else
      SUBROUTINE tl_sqlq
      RETURN
      END SUBROUTINE tl_sqlq
#endif
