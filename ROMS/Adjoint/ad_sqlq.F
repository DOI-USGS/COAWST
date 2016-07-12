#include "cppdefs.h"

#if (defined SENSITIVITY_4DVAR || defined TL_W4DPSAS || \
     defined TL_W4DVAR)        && defined MINRES

      SUBROUTINE ad_sqlq (innLoop, a, ad_a, tau, ad_tau, y, ad_y)
!
!svn $Id: ad_sqlq.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs the adjoint of a LQ factorization of the      !
!  square matrix A.                                                    !
!                                                                      !
!  NOTE: A, ad_A, TAU, ad_TAU, Y and ad_Y are overwritten on exit.     !
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

      real(r8), dimension(innLoop,innLoop), intent(inout) :: a, ad_a
      real(r8), dimension(innLoop), intent(inout) :: tau, ad_tau
      real(r8), dimension(innLoop), intent(inout) :: y, ad_y
!
!  Local variable declarations.
!
      integer :: i, j, ii, jj, m, kk, iflag

      real(r8) :: znorm, zbeta, zaii, ztemp, zbetas
      real(r8) :: adfac, ad_znorm, ad_zbeta, ad_zaii, ad_ztemp

      real(r8), dimension(innLoop,innLoop) :: as
      real(r8), dimension(innLoop) :: arow
!
!-----------------------------------------------------------------------
!  Adjoint of a LQ factorization of a square matrix.
!-----------------------------------------------------------------------
!
      ad_znorm=0.0_r8
      ad_zbeta=0.0_r8
      ad_zaii=0.0_r8
      ad_ztemp=0.0_r8
      DO i=1,innLoop
        ad_y(i)=0.0_r8
        arow(i)=0.0_r8
      END DO
!
!  Save a copy of a.
!
      DO j=1,innLoop
        DO i=1,innLoop
          as(i,j)=a(i,j)
        END DO
      END DO
!
      DO ii=innLoop-1,1,-1
!>      tl_a(ii,ii)=tl_zaii
!>
        ad_zaii=ad_zaii+ad_a(ii,ii)
        ad_a(ii,ii)=0.0_r8
!
!  Recompute basic state arrays.
!
        DO j=1,innLoop
          DO i=1,innLoop
            a(i,j)=as(i,j)
          END DO
        END DO
!
        iflag=1
        kk=ii
        m=innLoop
        call reclqbs (iflag, kk, m, a, tau, y)

        IF (tau(ii).ne.0.0_r8) THEN
          jj=ii-1
          DO j=1,innLoop-jj
            ztemp=-tau(ii)*a(ii,j+jj)
            DO i=1,innLoop-ii
!>            tl_a(ii+i,j+jj)=tl_a(ii+i,j+jj)+tl_y(i)*ztemp+            &
!>   &                                        y(i)*tl_ztemp
!>
              ad_ztemp=ad_ztemp+y(i)*ad_a(ii+i,j+jj)
              ad_y(i)=ad_y(i)+ztemp*ad_a(ii+i,j+jj)
            END DO
!>          tl_ztemp=-tl_tau(ii)*a(ii,j+jj)-tau(ii)*tl_a(ii,j+jj)
!>
            ad_tau(ii)=ad_tau(ii)-a(ii,j+jj)*ad_ztemp
            ad_a(ii,j+jj)=ad_a(ii,j+jj)-tau(ii)*ad_ztemp
            ad_ztemp=0.0_r8
          END DO
!
          DO j=1,innLoop-jj
            ztemp=a(ii,j+jj)
            DO i=1,innLoop-ii
!>            tl_y(i)=tl_y(i)+tl_ztemp*a(ii+i,j+jj)+                    &
!>   &                        ztemp*tl_a(ii+i,j+jj)
!>
              ad_ztemp=ad_ztemp+a(ii+i,j+jj)*ad_y(i)
              ad_a(ii+i,j+jj)=ad_a(ii+i,j+jj)+ztemp*ad_y(i)
            END DO
!>          tl_ztemp=tl_a(ii,j+jj)
!>
            ad_a(ii,j+jj)=ad_a(ii,j+jj)+ad_ztemp
            ad_ztemp=0.0_r8
          END DO
!
          DO i=1,innLoop
!>          tl_y(i)=0.0_r8
            ad_y(i)=0.0_r8
          END DO
        END IF
!>      tl_a(ii,ii)=0.0_r8
!>
        ad_a(ii,ii)=0.0_r8
!>      tl_zaii=tl_a(ii,ii)
!>
        ad_a(ii,ii)=ad_a(ii,ii)+ad_zaii
        ad_zaii=0.0_r8
!>      tl_a(ii,ii)=tl_zbeta
!>
        ad_zbeta=ad_zbeta+ad_a(ii,ii)
        ad_a(ii,ii)=0.0_r8
!
!  Recompute basic state arrays.
!
        DO j=1,innLoop
          DO i=1,innLoop
            a(i,j)=as(i,j)
          END DO
        END DO
!
        iflag=2
        kk=ii
        m=innLoop
        call reclqbs (iflag, kk, m, a, tau, y)
!
        znorm=0.0_r8
        DO j=ii+1,innLoop
          znorm=znorm+a(ii,j)*a(ii,j)
        END DO
        znorm=SQRT(znorm)
        zbeta=SQRT(znorm*znorm+a(ii,ii)*a(ii,ii))
        zbetas=zbeta
        zbeta=-SIGN(zbeta,a(ii,ii))
        tau(ii)=(zbeta-a(ii,ii))/zbeta
!
        DO j=ii+1,innLoop
          arow(j)=a(ii,j)
          a(ii,j)=a(ii,j)/(a(ii,ii)-zbeta)
!>        tl_a(ii,j)=tl_a(ii,j)/(a(ii,ii)-zbeta)-                       &
!>   &               (tl_a(ii,ii)-tl_zbeta)*a(ii,j)/(a(ii,ii)-zbeta)
!>
          adfac=ad_a(ii,j)/(a(ii,ii)-zbeta)
          ad_zbeta=ad_zbeta+a(ii,j)*adfac
          ad_a(ii,ii)=ad_a(ii,ii)-a(ii,j)*adfac
          ad_a(ii,j )=adfac
        END DO
!>      tl_tau(ii)=(tl_zbeta-tl_a(ii,ii))/zbeta-tl_zbeta*tau(ii)/zbeta
!>
        adfac=ad_tau(ii)/zbeta
        ad_zbeta=ad_zbeta+adfac-tau(ii)*adfac
        ad_a(ii,ii)=ad_a(ii,ii)-adfac
        ad_tau(ii)=0.0_r8
!>      tl_zbeta=-SIGN(1.0_r8,zbeta)*SIGN(1.0_r8,a(ii,ii))*tl_zbeta
!>
        ad_zbeta=-SIGN(1.0_r8,zbetas)*SIGN(1.0_r8,a(ii,ii))*ad_zbeta
        IF (zbetas.NE.0.0_r8) THEN
!>        tl_zbeta=(tl_znorm*znorm+tl_a(ii,ii)*a(ii,ii))/zbetas
!>
          adfac=ad_zbeta/zbetas
          ad_znorm=ad_znorm+znorm*adfac
          ad_a(ii,ii)=ad_a(ii,ii)+a(ii,ii)*adfac
          ad_zbeta=0.0_r8
        ELSE
!>        tl_zbeta=0.0_r8
!>
          ad_zbeta=0.0_r8
        END IF
        IF (znorm.NE.0.0_r8) THEN
!>        tl_znorm=0.5_r8*tl_znorm/znorm
!>
          ad_znorm=0.5_r8*ad_znorm/znorm
        ELSE
!>        tl_znorm=0.0_r8
!>
          ad_znorm=0.0_r8
        END IF
        DO j=ii+1,innLoop
!>        tl_znorm=tl_znorm+2.0_r8*a(ii,j)*tl_a(ii,j)
!>
          ad_a(ii,j)=ad_a(ii,j)+2.0_r8*arow(j)*ad_znorm
        END DO
!>      tl_znorm=0.0_r8
!>
        ad_znorm=0.0_r8
      END DO
      DO i=1,innLoop
!>      tl_tau(i)=0.0_r8
!>
        ad_tau(i)=0.0_r8
      END DO
!
      RETURN
      END SUBROUTINE ad_sqlq

      SUBROUTINE reclqbs (iflag, kk, innLoop, a, tau, y)
!
!***********************************************************************
!                                                                      !
!  This routine recomputes the intermediate arrays created during      !
!  the LQ factorization of A.                                          !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: innLoop, iflag, kk
      real(r8), dimension(innLoop,innLoop), intent(inout) :: a
      real(r8), dimension(innLoop), intent(inout) :: tau, y
!
!  Local variable declarations.
!
      integer :: i, j, ii, jj
      real(r8) :: znorm, zbeta, zaii, ztemp
!
!-----------------------------------------------------------------------
!  Recompute intermiediate arrays for the LQ factorization.
!-----------------------------------------------------------------------
!
      DO i=1,innLoop
        tau(i)=0.0_r8
      END DO
!
      DO ii=1,kk
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
!
        IF ((iflag.eq.2).and.(ii.eq.kk)) RETURN
!
        DO j=ii+1,innLoop
          a(ii,j)=a(ii,j)/(a(ii,ii)-zbeta)
        END DO
        a(ii,ii)=zbeta
!
!  Apply H(ii) to A(ii+1:innLoop,ii:innLoop) from the right.
!
        zaii = a(ii,ii)
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
          IF ((iflag.eq.1).and.(ii.eq.kk)) RETURN
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
      END SUBROUTINE reclqbs
#else
      SUBROUTINE ad_sqlq
      RETURN
      END SUBROUTINE ad_sqlq
#endif
