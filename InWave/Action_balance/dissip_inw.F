#include "cppdefs.h"
      MODULE dissip_inw_mod

#if defined INWAVE_MODEL
# if defined ENERGY_DISSIPATION
!
!=======================================================================
!                                                                      !
!  This routine computes the energy dissipation                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: dissip_inw_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE dissip_inw (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_ocean
      USE mod_coupling
      USE mod_stepping
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
!
!  Local variable declarations.
!
# include "tile.h"
!
!#  ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!#  endif

      CALL dissip_inw_tile(ng, tile,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 WAVEP(ng) % h_tot,                               &
     &                 WAVEP(ng) % AC, WAVEP(ng) % Tr, WAVEP(ng) % kwc  &
# ifdef MASKING
     &                 ,GRID(ng) % rmask                                &
# endif
#  if defined WDISS_INWAVE
     &                 ,FORCES(ng) % Dissip_break                       &
     &                 ,FORCES(ng) % Dissip_wcap                        &
#  endif
     &                 )
!#  ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!#  endif
      RETURN
      END SUBROUTINE dissip_inw
!
!***********************************************************************
      SUBROUTINE dissip_inw_tile(ng, tile,                              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       h_tot,                                     &
     &                       AC, Tr, kwc                                &
# ifdef MASKING
     &                       ,rmask                                     &
# endif
#  if defined WDISS_INWAVE
     &                       ,Dissip_break                              &
     &                       ,Dissip_wcap                               &
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
      integer, intent(in) :: nstp, nout
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   if defined WDISS_INWAVE
      real(r8), intent(inout) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcap(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   if defined WDISS_INWAVE
      real(r8), intent(inout) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcap(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Tr(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, d
      real(r8) :: EW, KM, TRM, H, Qb, Hmax_r, diff, Emax_r,twopi
#   if !defined WDISS_INWAVE
      real(r8) :: Db
#   endif
      real(r8), parameter :: alfa =1.0_r8
      real(r8), parameter :: break =0.45_r8
      real(r8), parameter :: n_r =15.0_r8
      real(r8), parameter :: EWlim = 0.00001_r8
      real(r8):: EN(ND)

# include "set_bounds.h"

      twopi=2.0_r8*pi
      
      DO j=Jstr,Jend
        DO i=Istr,Iend
          EW=0.0_r8
          KM=0.0_r8
          TRM=0.0_r8
          DO d=1,ND
!=======================================================================
!  Compute the energy from action balance and wave heigth
!=======================================================================
            EN(d)=AC(i,j,d,nout)*twopi/(Tr(i,j,d))
!=======================================================================
!  Compute the total energy
!=======================================================================
            EW=EW+EN(d)
!=======================================================================
!  Compute the mean wave number and intrinsic periods
!  What we do is give more importance to those wave 
!  numbers with more energy
!=======================================================================
            KM=KM+kwc(i,j,d)*EN(d)
            TRM=TRM+Tr(i,j,d)*EN(d)
          ENDDO
          KM=KM/(ND*max(EW,EWlim))
          TRM=TRM/(ND*max(EW,EWlim))
          IF (EW.lt.0.0_r8)EW=0.0_r8
!=======================================================================
!  Compute the wave height
!=======================================================================
!         H=(8.0_r8*EW/(g*rho0))**0.5_r8
          H=(16.0_r8*EW/(g*rho0))**0.5_r8
!=======================================================================
!  Compute the energy dissipation
!=======================================================================
#  ifdef ROELVINK
        if (h_tot(i,j).ge.Dcrit(ng))then
          Hmax_r=break*(max(h_tot(i,j),0.0_r8))
          Qb=min(1.0_r8,1.0_r8-exp(-(H/Hmax_r)**n_r))
#   if defined WDISS_INWAVE
          if(TRM.gt.0.0_r8)then
            Dissip_break(i,j)=2.0_r8*alfa*1./TRM*EW*Qb*dt(ng)
          endif
#    ifdef MASKING
          Dissip_break(i,j)=Dissip_break(i,j)*rmask(i,j)
#    endif
#   else
          if(TRM.gt.0.0_r8)then
            Db=2.0_r8*alfa*1./TRM*EW*Qb*dt(ng)
          endif
#   endif
         else
#   if defined WDISS_INWAVE
          Dissip_break(i,j)=0.0_r8
#   else
          Db=0.0_r8
#   endif
         endif
!=======================================================================
!  Distribute dissipation over directions and recompute Ac
!=======================================================================
          DO d=1,ND
#   if defined WDISS_INWAVE
            IF(Dissip_break(i,j).gt.0.0_r8.and.EW.gt.0.0_r8)THEN
#   else
            IF(Db.gt.0.0_r8.and.EW.gt.0.0_r8)THEN
#   endif
#   if defined WDISS_INWAVE
              EN(d)=EN(d)-Dissip_break(i,j)*EN(d)/EW
#   else
              EN(d)=EN(d)-Db*EN(d)/EW
#   endif
              AC(i,j,d,nout)=EN(d)*Tr(i,j,d)/(twopi)
            ENDIF
          ENDDO
#  else
          Hmax_r=break*(max(h_tot(i,j),0.0_r8))
!         Emax_r=0.125_r8*g*rho0*Hmax_r**2.0_r8
          Emax_r=0.0625_r8*g*rho0*Hmax_r**2.0_r8
          diff=EW-Emax_r
#   if defined WDISS_INWAVE
          Dissip_break(i,j)=max (0.0_r8,diff)
#    ifdef MASKING
          Dissip_break(i,j)=Dissip_break(i,j)*rmask(i,j)
#    endif
          Dissip_wcap(i,j)=0.0_r8
#   else
          Db=max (0.0_r8,diff)
#   endif
!=======================================================================
!  Distribute dissipation over directions and recompute Ac
!=======================================================================
          DO d=1,ND
            IF(H.gt.0.0_r8.and.EW.gt.0.0_r8)THEN
#   if defined WDISS_INWAVE
              EN(d)=EN(d)-Dissip_break(i,j)*EN(d)/EW
#   else
              EN(d)=EN(d)-Db*EN(d)/EW
#   endif
              AC(i,j,d,nout)=EN(d)*Tr(i,j,d)/(twopi)
            ENDIF
          ENDDO
#  endif

#   if defined WDISS_INWAVE
          Dissip_break(i,j)=Dissip_break(i,j)/(dt(ng)*rho0)
#   endif

        ENDDO
      ENDDO

# if defined WDISS_INWAVE
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_break)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Dissip_wcap)

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dissip_break, Dissip_wcap)
#  endif

# endif

      RETURN
      END SUBROUTINE dissip_inw_tile
# endif
#endif
      END MODULE dissip_inw_mod


