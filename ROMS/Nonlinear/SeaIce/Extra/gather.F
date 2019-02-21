#include "cppdefs.h"
      FUNCTION gather (ng, model, LBi, UBi, LBj, UBj, gtype, A, Acomb)
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
!     A            One tile of field (real).                           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Acomb        Combined field to return (real).                    !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
      USE mod_grid
#ifdef MASKING
      USE mod_ncparam
#endif

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastf, mp_gather2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, gtype

      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Acomb(0:Lm(ng)+1,0:Mm(ng)+1)
!
!  Local variable declarations.
!
      integer :: i, j, Npts, np, Imin, Jmin
      integer :: gather, MyType
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: AoutI

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order. If applicable, remove land points.
!-----------------------------------------------------------------------
!
! The tindex argument to mp_gather2d is used as a switch for setting
! land to spval. I know, it makes no sense at all.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Imin = 1
          Jmin = 1
        CASE (u2dvar, u3dvar)
          Imin = 1
          Jmin = 0
        CASE (v2dvar, v3dvar)
          Imin = 0
          Jmin = 1
        CASE (r2dvar, r3dvar)
          Imin = 0
          Jmin = 0
        CASE DEFAULT
          Imin = 0
          Jmin = 0
      END SELECT

# ifdef MASKING
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                    0, gtype, 1.0_r8,                             &
     &                    GRID(ng) % pmask,                             &
     &                    A, Npts, AoutI)
        CASE (u2dvar, u3dvar)
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                    0, gtype, 1.0_r8,                             &
     &                    GRID(ng) % umask,                             &
     &                    A, Npts, AoutI)
        CASE (v2dvar, v3dvar)
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                    0, gtype, 1.0_r8,                             &
     &                    GRID(ng) % vmask,                             &
     &                    A, Npts, AoutI)
        CASE (r2dvar, r3dvar)
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                    0, gtype, 1.0_r8,                             &
     &                    GRID(ng) % rmask,                             &
     &                    A, Npts, AoutI)
        CASE DEFAULT
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                    0, gtype, 1.0_r8,                             &
     &                    GRID(ng) % rmask,                             &
     &                    A, Npts, AoutI)
      END SELECT
# else
        CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj, gtype,         &
     &                    0, 1.0_r8,                                    &
     &                    A, Npts, AoutI)
# endif
      np=0
      DO j=Jmin,Mm(ng)+1
        DO i=Imin,Lm(ng)+1
          np=np+1
          Acomb(i,j) = AoutI(np)
        END DO
      END DO
      call mp_bcastf(ng, model, Acomb)
#else
!
!-----------------------------------------------------------------------
! Just return array
!-----------------------------------------------------------------------
!
#endif
      RETURN
      END FUNCTION gather
