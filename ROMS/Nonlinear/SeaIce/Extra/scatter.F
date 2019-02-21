#include "cppdefs.h"
      SUBROUTINE scatter (ng, model, LBi, UBi, LBj, UBj, gtype,         &
     &           Acomb, A, ai, Laice)
!=======================================================================
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     LBi          I-dimension Lower bound.                            !
!     UBi          I-dimension Upper bound.                            !
!     LBj          J-dimension Lower bound.                            !
!     UBj          J-dimension Upper bound.                            !
!     Acomb        Combined field (real).                              !
!     ai           Tiled ice area field.                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A            One tile of field (real).                           !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_ncparam, ONLY : r2dvar, u2dvar, v2dvar

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_scatter2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, gtype
      logical, intent(in) :: Laice

      real(r8), intent(in) :: Acomb(0:Lm(ng)+1,0:Mm(ng)+1)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j, Npts, np
      integer :: NWpts
      real(r8), dimension(2 + (Lm(ng)+2)*(Mm(ng)+2)) :: wrk
      real(r8) :: wrk2(LBi:UBi,LBj:UBj)
      real(r8) :: Amin, Amax, ai_vel

#ifdef DISTRIBUTE
      Npts=(Lm(ng)+2)*(Mm(ng)+2)
      np=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          np=np+1
          wrk(np) = Acomb(i,j)
        END DO
      END DO
!
! Everything from CICE is on rho points, not going to work with READ_WATER
!
      CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,                 &
     &                   NghostPoints, r2dvar, Amin, Amax,              &
# if defined READ_WATER && defined MASKING
     &                   NWpts, SCALARS(ng)%IJwater(1,wtype),           &
# endif
     &                   Npts, wrk, wrk2)
      SELECT CASE (gtype)
        CASE (u2dvar)
          DO j=LBj,MIN(UBj,Mm(ng)+1)
            DO i=LBi+1,MIN(UBi,Lm(ng)+1)
              ai_vel = 0.5*(ai(i,j)+ai(i-1,j))
              A(i,j) = 0.5*ai_vel*(wrk2(i,j)+wrk2(i-1,j)) +             &
     &                 (1-ai_vel)*A(i,j)
            END DO
          END DO
        CASE (v2dvar)
          DO j=LBj+1,MIN(UBj,Mm(ng)+1)
            DO i=LBi,MIN(UBi,Lm(ng)+1)
              ai_vel = 0.5*(ai(i,j)+ai(i,j-1))
              A(i,j) = 0.5*ai_vel*(wrk2(i,j)+wrk2(i,j-1)) +             &
     &                 (1-ai_vel)*A(i,j)
            END DO
          END DO
        CASE DEFAULT                              ! RHO-points
          IF (Laice) THEN
            DO j=LBj,MIN(UBj,Mm(ng)+1)
              DO i=LBi,MIN(UBi,Lm(ng)+1)
                A(i,j) = wrk2(i,j)
              END DO
            END DO
          ELSE
            DO j=LBj,MIN(UBj,Mm(ng)+1)
              DO i=LBi,MIN(UBi,Lm(ng)+1)
                A(i,j) = ai(i,j)*wrk2(i,j) + (1-ai(i,j))*A(i,j)
              END DO
            END DO
          END IF
      END SELECT
#else
!
!-----------------------------------------------------------------------
! Just return array (not quite)
!-----------------------------------------------------------------------
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            A(i,j) = ai(i,j)*wrk2(i,j) + (1-ai(i,j))*Acomb(i,j)
          END DO
        END DO
#endif
!
      RETURN
      END SUBROUTINE scatter
