!!
      SUBROUTINE ana_mask (ng, tile, model)
!
!! git $Id$
!! svn $Id: ana_mask.h 1151 2023-02-09 03:08:53Z arango $
!!======================================================================
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
!
! Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
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
        ANANAME(15)=MyFile
      END IF
!
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
      USE mod_parallel
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE stats_mod, ONLY : stats_2dfld
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
      logical, save :: first = .TRUE.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
!
      real(r8) :: mask(IminS:ImaxS,JminS:JmaxS)
!
      TYPE (T_STATS), save :: Stats(4)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize field statictics structure.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
        DO i=1,SIZE(Stats,1)
          Stats(i) % count=0.0_r8
          Stats(i) % min=Large
          Stats(i) % max=-Large
          Stats(i) % avg=0.0_r8
          Stats(i) % rms=0.0_r8
        END DO
      END IF
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
!  Report statitics.
!-----------------------------------------------------------------------
!
      CALL stats_2dfld (ng, tile, iNLM, p2dvar, Stats(1),               &
     &                  LBi, UBi, LBj, UBj, pmask)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'mask on PSI-points: mask_psi',               &
     &                    ng, Stats(1)%min, Stats(1)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(2),               &
     &                  LBi, UBi, LBj, UBj, rmask)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'mask on RHO-points: mask_rho',               &
     &                    ng, Stats(2)%min, Stats(2)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, u2dvar, Stats(3),               &
     &                  LBi, UBi, LBj, UBj, umask)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'mask on U-points: mask_u',                   &
     &                    ng, Stats(3)%min, Stats(3)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, v2dvar, Stats(4),               &
     &                  LBi, UBi, LBj, UBj, vmask)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'mask on V-points: mask_v',                   &
     &                    ng, Stats(4)%min, Stats(4)%max
      END IF
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
!
  10  FORMAT (3x,' ANA_MASK    - ',a,/,19x,                             &
     &        '(Grid = ',i2.2,', Min = ',1p,e15.8,0p,                   &
     &                         ' Max = ',1p,e15.8,0p,')')
!
      RETURN
      END SUBROUTINE ana_mask_tile
