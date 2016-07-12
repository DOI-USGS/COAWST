#include "cppdefs.h"
      MODULE shapiro_mod
!
!svn $Id: shapiro.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group         Kate Hedstrom   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package contains shapiro filter routines for order 2 and       !
!  reduced order at the boundary and mask edges.                       !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    shapirp2d_tile       Shapiro filter for 2D fields.                !
!    shapirp3d_tile       Shapiro filter for 3D fields.                !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE shapiro2d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
#ifdef MASKING
     &                           Amask,                                 &
#endif
     &                           A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk2

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Shapiro filter requested 2D field.
!-----------------------------------------------------------------------
!
!  This subroutine will apply a Shapiro filter of order 2 (defined
!  as twice the order in Shapiro (1970), with N even) to an array, A.
!  The order of the filter is reduced at the boundaries and at the
!  mask edges, if any.
!
!  Initialize filter in the Y-direction.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
#ifdef MASKING
          Awrk1(i,j)=0.25_r8*                                           &
     &               (A(i,j-1)*Amask(i,j-1)+                            &
     &                A(i,j+1)*Amask(i,j+1)-                            &
     &                2.0_r8*A(i,j)*Amask(i,j))*                        &
     &               Amask(i,j-1)*Amask(i,j+1)*Amask(i,j)
#else
          Awrk1(i,j)=0.25_r8*                                           &
     &               (A(i,j-1)+A(i,j+1)-2.0_r8*A(i,j))
#endif
        END DO
      END DO
!
!  Add the changes to the field.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
          Awrk2(i,j)=A(i,j)+Awrk1(i,j)
        END DO
      END DO
!
!  Initialize filter in the X-direction.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
#ifdef MASKING
          Awrk1(i,j)=0.25_r8*                                           &
     &               (Awrk2(i-1,j)*Amask(i-1,j)+                        &
     &                Awrk2(i+1,j)*Amask(i+1,j)-                        &
     &                2.0_r8*Awrk2(i,j)*Amask(i,j))*                    &
     &               Amask(i-1,j)*Amask(i+1,j)*Amask(i,j)
#else
          Awrk1(i,j)=0.25_r8*                                           &
     &               (Awrk2(i-1,j)+Awrk2(i+1,j)-2.0_r8*Awrk2(i,j))
#endif
        END DO
      END DO
!
!  Add changes to field.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          A(i,j)=Awrk2(i,j)+Awrk1(i,j)
        END DO
      END DO

      RETURN
      END SUBROUTINE shapiro2d_tile

#ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE shapiro3d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           IminS, ImaxS, JminS, JmaxS,            &
# ifdef MASKING
     &                           Amask,                                 &
# endif
     &                           A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Shapiro filter requested 3D field.
!-----------------------------------------------------------------------
!
!  This subroutine will apply a Shapiro filter of order 2 (defined
!  as twice the order in Shapiro (1970), with N even) to an array, A.
!  The order of the filter is reduced at the boundaries and at the
!  mask edges, if any.
!
!  Initialize filter in the Y-direction.
!
      DO k=LBk,UBk
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
# ifdef MASKING
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (A(i,j-1,k)*Amask(i,j-1)+                        &
     &                  A(i,j+1,k)*Amask(i,j+1)-                        &
     &                  2.0_r8*A(i,j,k)*Amask(i,j))*                    &
     &                 Amask(i,j-1)*Amask(i,j+1)*Amask(i,j)
# else
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (A(i,j-1,k)+A(i,j+1,k)-2.0_r8*A(i,j,k))
# endif
          END DO
        END DO
!
!  Add the changes to the field.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            Awrk2(i,j)=A(i,j,k)+Awrk1(i,j)
          END DO
        END DO
!
!  Initialize filter in the X-direction.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
# ifdef MASKING
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (Awrk2(i-1,j)*Amask(i-1,j)+                      &
     &                  Awrk2(i+1,j)*Amask(i+1,j)-                      &
     &                  2.0_r8*Awrk2(i,j)*Amask(i,j))*                  &
     &                 Amask(i-1,j)*Amask(i+1,j)*Amask(i,j)
# else
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (Awrk2(i-1,j)+Awrk2(i+1,j)-2.0_r8*Awrk2(i,j))
# endif
          END DO
        END DO
!
!  Add changes to field.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            A(i,j,k)=Awrk2(i,j)+Awrk1(i,j)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE shapiro3d_tile
#endif
      END MODULE shapiro_mod
