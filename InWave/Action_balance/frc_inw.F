#include "cppdefs.h"
      MODULE frc_inw_mod
#if defined INWAVE_MODEL
# if defined WEC_MELLOR || defined WEC_VF
!
!svn $Id: celer_inw.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This routine computes the energy dissipation                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: frc_inw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE frc_inw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_ocean
      USE mod_coupling
      USE mod_stepping
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
!#  ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!#  endif

      CALL frc_inw_tile(ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 WAVEP(ng)% AC, WAVEP(ng)% Tr,                    &
     &                 WAVEP(ng)% Ta, WAVEP(ng)% kwc,                   &
     &                 WAVEG(ng) % WD,                                  &
     &                 WAVEP(ng) % h_tot,                               &
     &                 FORCES(ng) % Hwave,                              &
     &                 FORCES(ng) % Dwave,                              &
     &                 FORCES(ng) % Lwave                               &
#  if defined WAVES_UB
     &                 ,FORCES(ng) % Uwave_rms                          &
#  endif
#  if defined WAVES_TOP_PERIOD
     &                 ,FORCES(ng) % Pwave_top                          &
#  endif
#  if defined WAVES_BOT_PERIOD
     &                 ,FORCES(ng) % Pwave_bot                          &
#  endif
     &                 )
!#  ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!#  endif
      RETURN
      END SUBROUTINE frc_inw
!
!***********************************************************************
      SUBROUTINE frc_inw_tile(ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       AC, Tr, Ta, kwc, WD,                       &
     &                       h_tot,                                     &
     &                       Hwave, Dwave, Lwave                        &
#  if defined WAVES_UB
     &                       ,Uwave_rms                                 &
#  endif
#  if defined WAVES_TOP_PERIOD
     &                       ,Pwave_top                                 &
#  endif
#  if defined WAVES_BOT_PERIOD
     &                       ,Pwave_bot                                 &
#  endif
     &                       )
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      USE exchange_2d_mod

#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
      USE bc_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(in) :: Ta(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Hwave(LBi:,LBj:)
      real(r8), intent(inout) :: Dwave(LBi:,LBj:)
      real(r8), intent(inout) :: Lwave(LBi:,LBj:)
      real(r8), intent(inout) :: h_tot(LBi:,LBj:)
#  if defined WAVES_UB
      real(r8), intent(inout) :: Uwave_rms(LBi:,LBj:)
#  endif
#  if defined WAVES_TOP_PERIOD
      real(r8), intent(inout) :: Pwave_top(LBi:,LBj:)
#  endif
#  if defined WAVES_BOT_PERIOD
      real(r8), intent(inout) :: Pwave_bot(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: WD(:)
# else
      real(r8), intent(in) :: Tr(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: Ta(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
      real(r8), intent(inout) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Lwave(LBi:UBi,LBj:UBj)

#  if defined WAVES_UB
      real(r8), intent(inout) :: Uwave_rms(LBi:UBi,LBj:UBj)
#  endif
#  if defined WAVES_TOP_PERIOD
      real(r8), intent(inout) :: Pwave_top(LBi:UBi,LBj:UBj)
#  endif
#  if defined WAVES_BOT_PERIOD
      real(r8), intent(inout) :: Pwave_bot(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: WD(ND)
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d
      real(r8) :: EW, KM, TRM, H, DIRM, twopi
      real(r8) :: EN
      real(r8), parameter :: Lmin = 10.0_r8
      real(r8), parameter :: Pmin = 2.0_r8
!

# include "set_bounds.h"

!!!mai: BE CAREFUL WITH THE DIRECTIONS TALK TO JOHN!!!

      twopi=2.0_r8*pi

      DO j=Jstr,Jend
        DO i=Istr,Iend
          EW=0.000001_r8
          KM=0.0_r8
          TRM=0.0_r8
          DIRM=0.0_r8
          IF(h_tot(i,j).ge.Dcrit(ng))then
            DO d=1,ND
!=======================================================================
!  Compute the energy from the action balance and the wave height
!=======================================================================
              EN=AC(i,j,d,3)*twopi/Tr(i,j,d)
              IF(EN.lt.0.0_r8)EN=0.0_r8
!=======================================================================
!  Compute the total energy
!=======================================================================
              EW=EW+EN
!=======================================================================
!  Compute the mean wave number and intrinsic periods
!  What we do is give more importance to those wave 
!  numbers with more energy
!=======================================================================
             KM=KM+kwc(i,j,d)*EN
             TRM=TRM+Tr(i,j,d)*EN
             DIRM=DIRM+WD(d)*EN
            ENDDO
            IF(EW.ge.0.001_r8)then
              KM=KM/EW
              Lwave(i,j)=twopi/KM
              Dwave(i,j)=DIRM/EW
!=======================================================================
!  Compute the wave height
!=======================================================================
!             Hwave(i,j)=(8.0_r8*EW/(g*rho0))**0.5_r8
              Hwave(i,j)=(16.0_r8*EW/(g*rho0))**0.5_r8
#  if defined WAVES_UB
              Uwave_rms(i,j)=0.0_r8   !!!!We need to define this
#  endif
#  if defined WAVES_TOP_PERIOD
              Pwave_top(i,j)=TRM/EW
#  endif
#  if defined WAVES_BOT_PERIOD
              Pwave_bot(i,j)=TRM/EW
#  endif
            ELSE
              Hwave(i,j)=0.0_r8
              Dwave(i,j)=0.0_r8
              Lwave(i,j)=Lmin
#  if defined WAVES_UB
              Uwave_rms(i,j)=0.0_r8
#  endif
#  if defined WAVES_TOP_PERIOD
              Pwave_top(i,j)=Pmin
#  endif
#  if defined WAVES_BOT_PERIOD
              Pwave_bot(i,j)=Pmin
#  endif
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      
# if defined WDISS_INWAVE
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Hwave)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dwave)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Lwave)
#  if defined WAVES_UB
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Uwave_rms)
#  endif
#  if defined WAVES_TOP_PERIOD
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Pwave_top)
#  endif
#  if defined WAVES_BOT_PERIOD
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Pwave_bot)
#  endif
# endif

#  ifdef DISTRIBUTE
!
!   Exchange boundary data.
!
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Hwave, Dwave, Lwave)
#  if defined WAVES_UB
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Uwave_rms)
#  endif
#  if defined WAVES_TOP_PERIOD
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Pwave_top)
#  endif
#  if defined WAVES_BOT_PERIOD
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Pwave_bot)
#  endif
!
#  endif

      RETURN
      END SUBROUTINE frc_inw_tile
# endif
#endif
      END MODULE frc_inw_mod


