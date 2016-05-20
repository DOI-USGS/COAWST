      SUBROUTINE ana_mask (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets analytical Land/Sea masking.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam                    
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_mask_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &                      
     &                    GRID(ng) % pmask,                             &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(15)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_mask
!
!***********************************************************************
      SUBROUTINE ana_mask_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &    
     &                          pmask, rmask, umask, vmask)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE                       
      real(r8), intent(out) :: pmask(LBi:,LBj:)
      real(r8), intent(out) :: rmask(LBi:,LBj:)
      real(r8), intent(out) :: umask(LBi:,LBj:)
      real(r8), intent(out) :: vmask(LBi:,LBj:)
#else                      
      real(r8), intent(out) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      real(r8) :: mask(IminS:ImaxS,JminS:JmaxS)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set Land/Sea mask of RHO-points: Land=0, Sea=1.
!-----------------------------------------------------------------------
!
!  Notice that private scratch array "mask" is used to allow
!  computation within a parallel loop.
!
#ifdef DOUBLE_GYRE
      Imin=-2+(Lm(ng)+1)/2
      Imax=Imin+2
      Jmin=-2+(Mm(ng)+1)/2
      Jmax=Jmin+2
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          mask(i,j)=1.0_r8
          IF (((Imin.le.i).and.(i.le.Imax)).and.                        &
     &        ((Jmin.le.j).and.(j.le.Jmax))) THEN
            mask(i,j)=0.0_r8
          END IF
        END DO
      END DO
#elif defined FLT_TEST
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          mask(i,j)=1.0_r8
          IF (j.eq.1 ) mask(i,j)=0.0_r8
          IF (j.eq.Mm(ng)) mask(i,j)=0.0_r8
          IF ((i.ge.((Lm(ng)+1)/2)).and.                                &
     &        (i.le.((Lm(ng)+1)/2+1)).and.                              &
     &        (j.ge.((Mm(ng)+1)/2)).and.                                &
     &        (j.le.((Mm(ng)+1)/2+1))) mask(i,j)=0.0_r8
        END DO
      END DO
#elif defined LAKE_SIGNELL
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          mask(i,j)=1.0_r8
        END DO
      END DO
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstrm1,Jendp1
          mask(Istr-1,j)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstrm1,Jendp1
          mask(Iend+1,j)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istrm1,Iendp1
          mask(i,Jstr-1)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istrm1,Iendp1
          mask(i,Jend+1)=0.0_r8
        END DO
      END IF
#elif defined RIVERPLUME1
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          mask(i,j)=1.0_r8
        END DO
      END DO
      DO i=Istrm2,MIN(5,Iendp2)
        DO j=Jstrm2,MIN(Mm(ng)-18,Jendp2)
          mask(i,j)=0.0_r8
        END DO
        DO j=MAX(Jstrm2,Mm(ng)-16),Jendp2
          mask(i,j)=0.0_r8
        END DO
      END DO
#elif defined RIVERPLUME2
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          mask(i,j)=1.0_r8
        END DO
      END DO
      DO i=Istrm2,MIN(5,Iendp2)
        DO j=Jstrm2,MIN(Mm(ng)-11,Jendp2)
          mask(i,j)=0.0_r8
        END DO
        DO j=MAX(Jstrm2,Mm(ng)-9),Jendp2
          mask(i,j)=0.0_r8
        END DO
      END DO
#elif defined SHOREFACE
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          mask(i,j)=1.0_r8
        END DO
      END DO
#else
      ana_mask.h: no values provided for mask.
#endif
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rmask(i,j)=mask(i,j)
        END DO
      END DO                    
!
!-----------------------------------------------------------------------
!  Compute Land/Sea mask of U- and V-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          umask(i,j)=mask(i-1,j)*mask(i,j)
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vmask(i,j)=mask(i,j-1)*mask(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Land/Sea mask of PSI-points.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendT
        DO i=IstrP,IendT
          pmask(i,j)=mask(i-1,j-1)*mask(i,j-1)*                         &
     &               mask(i-1,j  )*mask(i,j  )
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask, pmask, umask, vmask)
#endif

      RETURN
      END SUBROUTINE ana_mask_tile
