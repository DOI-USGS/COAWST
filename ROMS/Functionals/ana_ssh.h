      SUBROUTINE ana_ssh (ng, tile, model)
!
!! svn $Id: ana_ssh.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical sea surface height climatology.        !
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
      CALL ana_ssh_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   CLIMA(ng) % ssh)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(28)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_ssh
!
!***********************************************************************
      SUBROUTINE ana_ssh_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         ssh)
!***********************************************************************
!
      USE mod_param
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ssh(LBi:,LBj:)
#else
      real(r8), intent(out) :: ssh(LBi:UBi,LBj:UBj)
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
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set sea surface height (meters).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ssh(i,j)=???
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ssh)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ssh)
#endif
      RETURN
      END SUBROUTINE ana_ssh_tile
