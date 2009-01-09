      SUBROUTINE ana_tclima (ng, tile, model)
!
!! svn $Id: ana_tclima.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets analytical tracer climatology fields.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tclima_tile (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      CLIMA(ng) % tclm)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(33),'(a,a)') TRIM(Adir), '/ana_tclima.h'
      END IF

      RETURN
      END SUBROUTINE ana_tclima
!
!***********************************************************************
      SUBROUTINE ana_tclima_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tclm)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tclm(LBi:,LBj:,:,:)
#else
      real(r8), intent(out) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tclm(i,j,k,itemp)=???
            tclm(i,j,k,isalt)=???
          END DO
        END DO
      END DO
#else
      ana_tclima.h: No values provided for tclm.
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      DO itrc=1,NAT
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tclm(:,:,:,itrc))
      END DO
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NAT,         &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tclm)
#endif

      RETURN
      END SUBROUTINE ana_tclima_tile
