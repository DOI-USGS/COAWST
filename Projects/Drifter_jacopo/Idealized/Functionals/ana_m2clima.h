      SUBROUTINE ana_m2clima (ng, tile, model)
!
!! svn $Id: ana_m2clima.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets analytical 2D momentum climatology fields.        !
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
      CALL ana_m2clima_tile (ng, model, Istr, Iend, Jstr, Jend,         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng) % ubarclm,                       &
     &                       CLIMA(ng) % vbarclm)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(11),'(a,a)') TRIM(Adir), '/ana_m2clima.h'
      END IF

      RETURN
      END SUBROUTINE ana_m2clima
!
!***********************************************************************
      SUBROUTINE ana_m2clima_tile (ng, model, Istr, Iend, Jstr, Jend,   &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ubarclm, vbarclm)
!***********************************************************************
!
      USE mod_param
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ubarclm(LBi:,LBj:)
      real(r8), intent(out) :: vbarclm(LBi:,LBj:)
#else
      real(r8), intent(out) :: ubarclm(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vbarclm(LBi:UBi,LBj:UBj)
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
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubarclm(i,j)=???
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbarclm(i,j)=???
        END DO
      END DO
#else
      ana_m2clima.h: No values provided for ubarclm and vbarclm.
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubarclm)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbarclm)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ubarclm, vbarclm)
#endif

      RETURN
      END SUBROUTINE ana_m2clima_tile
