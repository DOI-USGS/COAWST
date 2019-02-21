#include "cppdefs.h"
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!
!BOP
!
! !MODULE: ice_coupling - message passing to and from the coupler
!
! !DESCRIPTION:
!
! Message passing to and from the coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!         Tony Craig, NCAR, Dec-30-2002, modified for cpl6
!
! !INTERFACE:
!
      module ice_fakecpl
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_broadcast
      use ice_constants
!      USE ice_constants, only: field_type_scalar, field_loc_center, &
!               field_loc_NEcorner
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_init
      use ice_history
      use ice_restart
      use ice_timers
      use ice_fileunits
      use ice_domain_size, only: max_blocks
      use ice_restart_shared, only: restart
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: nblocks, distrb_info

      implicit none

      contains

      !------------------------------------

      SUBROUTINE cice_fakecpl

      USE mod_param
      USE mod_parallel
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_stepping
!
      USE distribute_mod
      USE mod_grid
      USE CICE_RunMod
      USE mod_forces
      USE mod_ocean
      USE mod_coupling
      USE mod_ice
      USE mod_ncparam
      use ice_gather_scatter
      use ice_calendar, only: istep, istep1, stop_now, calendar
      use ice_calendar, only: ice_time=>time, ice_dt=>dt
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ng

!  KLUDGE - REVISIT THIS
      integer :: ii, jj, k
      integer ::                                                        &
     &   iblk           ,&! block indices
     &   ilo,ihi,jlo,jhi  ! beginning and end of physical domain
      type (block) ::                                                   &
     &   this_block       ! block information for current block

      real(r8) :: wrk_big(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: wrk_small(Lm(1),Mm(1))
      real(r8) :: urho(Lm(1),Mm(1))
      real(r8) :: vrho(Lm(1),Mm(1))
      real(r8) :: ua(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: va(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: zetaa(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: pma(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: pna(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: tltx(Lm(1),Mm(1))
      real(r8) :: tlty(Lm(1),Mm(1))
      real(r8) :: aicea(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: vicea(0:Lm(1)+1,0:Mm(1)+1)
      real(r8) :: fresh_small(Lm(1),Mm(1))
      real(r8) :: sss_small(Lm(1),Mm(1))
      integer :: LBi, UBi, LBj, UBj, tile
      integer :: gfactor, gtype, status
      real(r8) :: Hscale

!
!-----------------------------------------------------------------------
! Gather fields for the ice model
!-----------------------------------------------------------------------
!
! KLUDGE here! hard-coding ng=1 (and above)
        ng = 1
        Hscale = 1.0_r8/(rho0*Cp)
        DO tile=first_tile(ng),last_tile(ng),+1
          LBi=BOUNDS(ng)%LBi(tile)
          UBi=BOUNDS(ng)%UBi(tile)
          LBj=BOUNDS(ng)%LBj(tile)
          UBj=BOUNDS(ng)%UBj(tile)
          CALL ice_timer_start(timer_step)
! Ocean velocity
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, u2dvar,             &
     &                OCEAN(ng)%u(LBi,LBj,N(ng),nrhs(ng)),ua)
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, v2dvar,             &
     &                OCEAN(ng)%v(LBi,LBj,N(ng),nrhs(ng)),va)
          CALL roms_u2rho(ua, urho)
          CALL roms_v2rho(va, vrho)
          call scatter_global(uocn, urho, master_task, distrb_info,     &
     &             field_loc_center, field_type_scalar)
          call scatter_global(vocn, vrho, master_task, distrb_info,     &
     &             field_loc_center, field_type_scalar)

! SST
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                OCEAN(ng)%t(LBi,LBj,N(ng),nrhs(ng),itemp),wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(sst, wrk_small, master_task, distrb_info, &
     &             field_loc_center, field_type_scalar)

! SSS
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                OCEAN(ng)%t(LBi,LBj,N(ng),nrhs(ng),isalt),wrk_big)
          sss_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(sss, sss_small, master_task, distrb_info, &
     &             field_loc_center, field_type_scalar)

! Surface tilt
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                COUPLING(ng)%Zt_avg1(LBi,LBj),zetaa)
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                GRID(ng)%pm(LBi,LBj),pma)
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                GRID(ng)%pn(LBi,LBj),pna)
          CALL roms_zeta2tlt(zetaa, pma, pna, tltx, tlty)
          call scatter_global(ss_tltx, tltx, master_task, distrb_info,  &
     &             field_loc_NEcorner, field_type_scalar)
          call scatter_global(ss_tlty, tlty, master_task, distrb_info,  &
     &             field_loc_NEcorner, field_type_scalar)

! Air temperature
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%PotT, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
! CHECK THIS
          call scatter_global(PotT, wrk_small, master_task, distrb_info,&
     &             field_loc_center, field_type_scalar)
          PotT = PotT + 273.16
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%Tair, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(Tair, wrk_small, master_task, distrb_info,&
     &             field_loc_center, field_type_scalar)
          Tair = Tair + 273.16

! Humidity
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%Hair, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(Qa, wrk_small, master_task, distrb_info,  &
     &             field_loc_center, field_type_scalar)

! Winds
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%Uwind, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(uatm, wrk_small, master_task, distrb_info,&
     &             field_loc_center, field_type_scalar)
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%Vwind, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(vatm, wrk_small, master_task, distrb_info,&
     &             field_loc_center, field_type_scalar)
          wind(:,:,:) = sqrt(uatm(:,:,:)**2                             &
     &                           + vatm(:,:,:)**2)  ! wind speed, (m/s)

! Rain
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%rain, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(frain, wrk_small, master_task,            &
     &             distrb_info, field_loc_center, field_type_scalar)
#ifdef SNOWFALL
! Snow
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%snow, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(fsnow, wrk_small, master_task,            &
     &             distrb_info, field_loc_center, field_type_scalar)
#endif

! Ice growth potential
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                ICE(ng)%wfr, wrk_big)
! units already W/m2
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(frzmlt, wrk_small, master_task,           &
     &             distrb_info, field_loc_center, field_type_scalar)

! Shortwave
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%SW_down, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(swvdr, wrk_small, master_task,            &
     &             distrb_info, field_loc_center, field_type_scalar)
! HACK has to be this order:
          fsw = swvdr
          swvdf = 0.24*swvdr
          swidr = 0.31*swvdr
          swidf = 0.16*swvdr
          swvdr = 0.29*swvdr

! Longwave
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                FORCES(ng)%LW_down, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(flw, wrk_small, master_task, distrb_info, &
     &             field_loc_center, field_type_scalar)

! call CICE_MODEL
          zlvl = 2.0_r8   ! MERRA
          IF (iic(ng).eq.ntstart(ng) .and. .not. restart) THEN
#ifdef INI_GLORYS_ICE
! ai
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                ICE(ng)%ai, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(aice_init2, wrk_small, master_task,       &
     &             distrb_info, field_loc_center, field_type_scalar)
! hi
          call gather(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,             &
     &                ICE(ng)%hi, wrk_big)
          wrk_small = wrk_big(1:Lm(ng),1:Mm(Ng))
          call scatter_global(hice_init, wrk_small, master_task,        &
     &             distrb_info, field_loc_center, field_type_scalar)
#endif
            CALL init_state
          END IF
          CALL ice_step

          CALL ice_timer_stop(timer_step)

! Fields to go back to ROMS, starting with ice concentration
          CALL gather_global(wrk_small, aice, master_task, distrb_info)
          aicea(1:Lm(ng),1:Mm(Ng)) = wrk_small
          aicea(0,1:Mm(Ng)) = aicea(1,1:Mm(Ng))
          aicea(Lm(ng)+1,1:Mm(Ng)) = aicea(Lm(ng),1:Mm(Ng))
          aicea(:,0) = aicea(:,1)
          aicea(:,Mm(Ng)+1) = aicea(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,            &
     &             aicea, ICE(ng)%ai, ICE(ng)%ai, .true.)

! ice volume
          CALL gather_global(wrk_small, vice, master_task, distrb_info)
          vicea(1:Lm(ng),1:Mm(Ng)) = wrk_small
          vicea(0,1:Mm(Ng)) = vicea(1,1:Mm(Ng))
          vicea(Lm(ng)+1,1:Mm(Ng)) = vicea(Lm(ng),1:Mm(Ng))
          vicea(:,0) = vicea(:,1)
          vicea(:,Mm(Ng)+1) = vicea(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar,            &
     &             vicea, ICE(ng)%hi, ICE(ng)%ai, .true.)

! Surface stress
          CALL gather_global(wrk_small, strocnxT, master_task,          &
     &                     distrb_info)
! Ocean stress from ice opposite to ice stress from ocean, divide by
! rho0.
          wrk_small = -wrk_small/rho0
          wrk_big(1:Lm(ng),1:Mm(Ng)) = wrk_small
          wrk_big(0,1:Mm(Ng)) = wrk_big(1,1:Mm(Ng))
          wrk_big(Lm(ng)+1,1:Mm(Ng)) = wrk_big(Lm(ng),1:Mm(Ng))
          wrk_big(:,0) = wrk_big(:,1)
          wrk_big(:,Mm(Ng)+1) = wrk_big(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, u2dvar, wrk_big,   &
     &             FORCES(ng)%sustr, ICE(ng)%ai, .false.)
          CALL gather_global(wrk_small, strocnyT, master_task,          &
     &                     distrb_info)
          wrk_small = -wrk_small/rho0
          wrk_big(1:Lm(ng),1:Mm(Ng)) = wrk_small
          wrk_big(0,1:Mm(Ng)) = wrk_big(1,1:Mm(Ng))
          wrk_big(Lm(ng)+1,1:Mm(Ng)) = wrk_big(Lm(ng),1:Mm(Ng))
          wrk_big(:,0) = wrk_big(:,1)
          wrk_big(:,Mm(Ng)+1) = wrk_big(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, v2dvar, wrk_big,   &
     &             FORCES(ng)%svstr, ICE(ng)%ai, .false.)

! Shortwave
          CALL gather_global(wrk_small, fswthru, master_task,           &
     &                       distrb_info)
          wrk_small = Hscale*wrk_small
          wrk_big(1:Lm(ng),1:Mm(Ng)) = wrk_small
          wrk_big(0,1:Mm(Ng)) = wrk_big(1,1:Mm(Ng))
          wrk_big(Lm(ng)+1,1:Mm(Ng)) = wrk_big(Lm(ng),1:Mm(Ng))
          wrk_big(:,0) = wrk_big(:,1)
          wrk_big(:,Mm(Ng)+1) = wrk_big(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar, wrk_big,   &
     &             FORCES(ng)%srflx, ICE(ng)%ai, .false.)

! Surface tracer fluxes
          CALL gather_global(wrk_small, fhocn, master_task,             &
     &                     distrb_info)
! scale from W/m^2 to ROMS units
          wrk_small = Hscale*wrk_small
          wrk_big(1:Lm(ng),1:Mm(Ng)) = wrk_small
          wrk_big(0,1:Mm(Ng)) = wrk_big(1,1:Mm(Ng))
          wrk_big(Lm(ng)+1,1:Mm(Ng)) = wrk_big(Lm(ng),1:Mm(Ng))
          wrk_big(:,0) = wrk_big(:,1)
          wrk_big(:,Mm(Ng)+1) = wrk_big(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar, wrk_big,   &
     &             FORCES(ng)%stflx(LBi,LBj,1), ICE(ng)%ai, .false.)
          CALL gather_global(wrk_small, fsalt, master_task,             &
     &                     distrb_info)
          CALL gather_global(fresh_small, fresh, master_task,           &
     &                     distrb_info)
! add fresh water flux to salt flux
          wrk_small = wrk_small - sss_small*1.e-3_r8*fresh_small
          wrk_big(1:Lm(ng),1:Mm(Ng)) = wrk_small
          wrk_big(0,1:Mm(Ng)) = wrk_big(1,1:Mm(Ng))
          wrk_big(Lm(ng)+1,1:Mm(Ng)) = wrk_big(Lm(ng),1:Mm(Ng))
          wrk_big(:,0) = wrk_big(:,1)
          wrk_big(:,Mm(Ng)+1) = wrk_big(:,Mm(Ng))
          call scatter(ng, iNLM, LBi, UBi, LBj, UBj, r2dvar, wrk_big,   &
     &             FORCES(ng)%stflx(LBi,LBj,2), ICE(ng)%ai, .false.)

! Copying more from CICE_Run

          istep  = istep  + 1    ! update time step counters
          istep1 = istep1 + 1
          ice_time = ice_time + ice_dt       ! determine the time and date
          call calendar(ice_time)    ! at the end of the timestep
!          if (stop_now >= 1) exit timeLoop

          CALL init_flux_atm   ! reset atmosphere fluxes to zero for next step
          CALL init_flux_ocn   ! reset tracer fluxes to zero for next step

        END DO

      RETURN
      END SUBROUTINE cice_fakecpl
      !------------------------------------

      subroutine roms_u2rho(uvar,rhovar)

      use mod_param

      real (kind=dbl_kind), intent(in) :: uvar(0:Lm(1)+1,0:Mm(1)+1)
      real (kind=dbl_kind), intent(out) :: rhovar(Lm(1),Mm(1))

      integer(kind=int_kind) :: i,j   !loop index!

      !kca - changes made to uv2rho and rho2uv subroutines

      do j=1,Mm(1)
        do i=1,Lm(1)
!          if ((uvar(i+1,j).eq.c0) .or. (uvar(i,j).eq.c0)) then
!            rhovar(i,j) = uvar(i+1,j) + uvar(i,j)
!          else
            rhovar(i,j) = 0.5_r8 * (uvar(i+1,j) + uvar(i,j))
!          endif
!          rhovar(i,j) = 0.0_r8
        enddo
      enddo

      end subroutine roms_u2rho

      !-------------------------------------

      subroutine roms_v2rho(vvar,rhovar)

      use mod_param

      real (kind=dbl_kind), intent(in) :: vvar(0:Lm(1)+1,0:Mm(1)+1)
      real (kind=dbl_kind), intent(out) :: rhovar(Lm(1),Mm(1))

      integer(kind=int_kind) :: i,j   !loop index

      do j=1,Mm(1)
        do i=1,Lm(1)
!          if ((vvar(i,j+1).eq.c0) .or. (vvar(i,j).eq.c0)) then
!            rhovar(i,j) = vvar(i,j+1) + vvar(i,j)
!          else
            rhovar(i,j) = 0.5_r8 * (vvar(i,j+1) + vvar(i,j))
!          endif
!          rhovar(i,j) = 0.0_r8
        enddo
      enddo

      end subroutine roms_v2rho

      !--------------------------------------

      subroutine roms_rho2u(rhovar,uvar)

      use mod_param

      real (kind=dbl_kind), intent(in) :: rhovar(Lm(1)+2,Mm(1)+2)
      real (kind=dbl_kind), intent(out) :: uvar(Lm(1)+2,Mm(1)+2)

      integer(kind=int_kind) :: i,j   !loop index

      do j=1,Mm(1)+2
        do i=2,Lm(1)+2
          if ((rhovar(i,j).eq.c0) .or. (rhovar(i-1,j).eq.c0)) then
            uvar(i,j) = rhovar(i,j) + rhovar(i-1,j)
          else
            uvar(i,j) = 0.5_r8 * ( rhovar(i,j) + rhovar(i-1,j) )
          endif
        enddo
      enddo

      do j=1,Mm(1)+2
        !uvar(1,j) = rhovar(1,j)
         uvar(1,j) = c0
      enddo

      end subroutine roms_rho2u

      !--------------------------------------

      subroutine roms_rho2v(rhovar,vvar)

      use mod_param

      real (kind=dbl_kind), intent(in) :: rhovar(Lm(1)+2,Mm(1)+2)
      real (kind=dbl_kind), intent(out) :: vvar(Lm(1)+2,Mm(1)+2)

      integer(kind=int_kind) :: i,j   !loop index

      do i=1,Lm(1)+2
        do j=2,Mm(1)+2
          if ((rhovar(i,j).eq.c0) .or. (rhovar(i,j-1).eq.c0)) then
            vvar(i,j) = rhovar(i,j) + rhovar(i,j-1)
          else
            vvar(i,j) = 0.5_r8 * (rhovar(i,j) + rhovar(i,j-1))
          endif
        enddo
      enddo

      do i=1,Lm(1)+2
        !vvar(i,1) = rhovar(i,1)
         vvar(i,1) = c0
      enddo

      end subroutine roms_rho2v

      !--------------------------------------

      subroutine roms_zeta2tlt(zeta, pm, pn, tltx, tlty)

      use mod_param

      real (kind=dbl_kind), intent(in) :: zeta(0:Lm(1)+1,0:Mm(1)+1)
      real (kind=dbl_kind), intent(in) :: pm(0:Lm(1)+1,0:Mm(1)+1)
      real (kind=dbl_kind), intent(in) :: pn(0:Lm(1)+1,0:Mm(1)+1)
      real (kind=dbl_kind), intent(out) :: tltx(Lm(1),Mm(1))
      real (kind=dbl_kind), intent(out) :: tlty(Lm(1),Mm(1))

      integer(kind=int_kind) :: i,j   !loop index!

      real (kind=dbl_kind) :: tmp(0:Lm(1)+1,0:Mm(1)+1)

! at ROMS u-points
      do j=0,Mm(1)+1
        do i=1,Lm(1)+1
          if ((zeta(i-1,j).eq.c0) .or. (zeta(i,j).eq.c0)) then
            tmp(i,j) = c0
          else
            tmp(i,j) = 0.5*(pm(i-1,j)+pm(i,j))*(zeta(i,j)-zeta(i-1,j))
          endif
        enddo
      enddo

! at ROMS psi-points (NE_corner convention)
      do j=1,Mm(1)
        do i=1,Lm(1)
          tltx(i,j) = 0.5*(tmp(i+1,j)+tmp(i+1,j+1))
        enddo
      enddo

! at ROMS v-points
      do j=1,Mm(1)+1
        do i=0,Lm(1)+1
          if ((zeta(i,j-1).eq.c0) .or. (zeta(i,j).eq.c0)) then
            tmp(i,j) = c0
          else
            tmp(i,j) = 0.5*(pn(i,j-1)+pn(i,j))*(zeta(i,j)-zeta(i,j-1))
          endif
        enddo
      enddo

! at ROMS psi-points (NE_corner convention)
      do j=1,Mm(1)
        do i=1,Lm(1)
          tlty(i,j) = 0.5*(tmp(i,j+1)+tmp(i+1,j+1))
        enddo
      enddo

      end subroutine roms_zeta2tlt

      end module ice_fakecpl
