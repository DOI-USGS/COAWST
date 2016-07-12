      SUBROUTINE ana_dqdsst (ng, tile, model)
!
!! svn $Id: ana_dqdsst.h 795 2016-05-11 01:42:43Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets surface net heat flux sensitivity dQdSTT to    !
!  SST using analytical expressions.  The forcing dQdSTT is usually    !
!  computed in units of (Watts/m2/degC).  It needs to be scaled to     !
!  (m/s/degC) by dividing by rho0*Cp. These forcing fields are used    !
!  when surface heat flux correction is activated:                     !
!                                                                      !
!       Q_model ~ Q + dQdSST * (T_model - SST)                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_dqdsst_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % Hz,                              &
     &                      FORCES(ng) % dqdt)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(38)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_dqdsst
!
!***********************************************************************
      SUBROUTINE ana_dqdsst_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Hz, dqdt)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
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
      real(r8), intent(in)  :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: dqdt(LBi:,LBj:)
#else
      real(r8), intent(in)  :: Hz(LBi:UBi,LBj:UBj,1:N(ng))
      real(r8), intent(out) :: dqdt(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set surface heat flux sensitivity to SST (m/s/degC).
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          dqdt(i,j)=???
        END DO
      END DO
#else
      ana_dqdsst.h: no values provided for dqdt.
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dqdt)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dqdt)
#endif

      RETURN
      END SUBROUTINE ana_dqdsst_tile
