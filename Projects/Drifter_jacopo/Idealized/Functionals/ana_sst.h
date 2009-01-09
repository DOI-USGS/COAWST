      SUBROUTINE ana_sst (ng, tile, model)
!
!! svn $Id: ana_sst.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This subroutine sets sea surface temperature SST  (Celsius)  and    !
!  surface net heat flux sensitivity dQdSTT to SST using analytical    !
!  expressions.  The forcing dQdSTT is usually computed in units of    !
!  (Watts/m2/degC).  It needs to be scaled to (m/s) by dividing by     !
!  rho0*Cp.  These forcing fields are used  when flux correction is    !
!  activated:                                                          !
!                                                                      !
!       Q_model ~ Q + dQdSST * (T_model - SST)                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_sst_tile (ng, model, Istr, Iend, Jstr, Jend,             &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   FORCES(ng) % sst,                              &
     &                   FORCES(ng) % dqdt)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(30),'(a,a)') TRIM(Adir), '/ana_sst.h'
      END IF

      RETURN
      END SUBROUTINE ana_sst
!
!***********************************************************************
      SUBROUTINE ana_sst_tile (ng, model, Istr, Iend, Jstr, Jend,       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         sst, dqdt)
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
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: sst(LBi:,LBj:)
      real(r8), intent(out) :: dqdt(LBi:,LBj:)
#else
      real(r8), intent(out) :: sst(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: dqdt(LBi:UBi,LBj:UBj)
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
!  Set sea surface temperature (Celsius) and heat flux sensitivity to
!  SST (Watts/m2).
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          sst(i,j)=???
          dqdt(i,j)=???
        END DO
      END DO
#else
      ana_sst.h: no values provided for sst and dqdt.
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sst)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        dqdt)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sst, dqdt)
#endif

      RETURN
      END SUBROUTINE ana_sst_tile
