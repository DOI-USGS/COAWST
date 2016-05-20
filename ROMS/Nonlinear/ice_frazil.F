#include "cppdefs.h"
      MODULE ice_frazil_mod
#if (defined ICE_MODEL && defined ICE_THERMO && !defined TS_FIXED) \
      || defined CICE_MODEL
!
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the frazil ice growth in the water when the
!  water temperature gets below freezing. It adjusts the water
!  temperature and salinity accordingly.
!
!  Reference: Steele et al. (1989). JPO, 19, 139-147.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_frazil

      CONTAINS

      SUBROUTINE ice_frazil (ng, tile)

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_stepping
# ifdef CICE_MODEL
      USE mod_forces
# endif

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 51)
# endif
!
      CALL ice_frazil_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nnew(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
# ifdef CICE_MODEL
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % stflx,                         &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      GRID(ng) % z_w,                             &
# endif
     &                      ICE(ng) % wfr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 51)
# endif
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      subroutine ice_frazil_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nnew,                                 &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            Hz, z_r, rho, t,                      &
# ifdef CICE_MODEL
     &                            sustr, svstr, stflx, ai, hi, z_w,     &
# endif
     &                            wfr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef CICE_MODEL
      USE shr_const_mod
# endif
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod
      USE distribute_mod, ONLY : mp_reduce
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  ifdef CICE_MODEL
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:)
      real(r8), intent(in) :: hi(LBi:,LBj:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(out) :: wfr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef CICE_MODEL
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(out) :: wfr(LBi:UBi,LBj:UBj)
# endif
!
! Local variable definitions
!
      integer :: i, j, k, itrc

      real(r8), parameter :: Lhat = 79.2_r8
      real(r8), parameter :: r = 0.5_r8

      real(r8) :: t_freeze
      real(r8) :: s1
      real(r8) :: z1
      real(r8) :: gamma_k
      real(r8) :: t_fr
      real(r8) :: ice_dens
      real(r8) :: delta_wfr
# ifdef CICE_MODEL
      real(r8), parameter :: depressT = -0.054_r8
# else
      real(r8), parameter :: depressT = -0.0543_r8
# endif

# ifdef CICE_MODEL
#  ifdef ICE_LOG_LAYER
! Compute heat flux from ocean a la Mellor and Kantha, 1989
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: t0mk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: s0mk
      real(r8), parameter :: prt = 13._r8
      real(r8), parameter :: prs = 2432._r8
      real(r8), parameter :: tpr = 0.85_r8
      real(r8), parameter :: nu = 1.8E-6_r8
      real(r8), parameter :: z0ii = 0.02_r8
      real(r8), parameter :: kappa = 0.4_r8
      real(r8), parameter :: hfus = 3.347E+5_r8         ! [J kg-1]
      real(r8), parameter :: ykf = 3.14
      real(r8), parameter :: eps = 1.0e-20_r8
      real(r8) :: utau, dztop, rno, termt, terms, cht, chs, z0, zdz0
#  endif
      real(r8) :: dz, heat, tfr
      real(r8), parameter :: bl_thick = 5.0_r8
      real(r8) :: ice_thick, thickness, fract
      real(r8), parameter :: latent_heat_fusion = 3.337e5_r8   ! J/kg
      real(r8), parameter :: cp_sw              = 3996.0_r8    ! J/kg/K
      real(r8), parameter :: rho_sw             = 1026.0_r8    ! kg/m^3
      real(r8), parameter :: rho_fw             = 1000.0_r8    ! kg/m^3
      real(r8), parameter :: cp_over_lhfusion   =                       &
     &                 rho_sw*cp_sw/(latent_heat_fusion*rho_fw)
      real(r8), parameter :: rhocpr = 4093740._r8       ! [W s m-3 K-1]
      real(r8), parameter :: salice = 4.0_r8
      real(r8), parameter :: salref = 34.7_r8
      real(r8) :: potice
# endif

# ifdef DISTRIBUTE
      real(r8), allocatable           :: buffer(:)
      character (len=3), allocatable  :: op_handle(:)
# endif

!  Inline functions
!  Freezing temperature (Gill, 1982)
!     t_freeze(s1,z1) = -0.0575*s1 + 1.710523d-3*sqrt(s1)**3
!    &       - 2.154996d-4*s1*s1 + 0.000753*z1
!  Freezing temperature (Steele et al. 1989)
!      t_freeze(s1,z1) = -0.0543*s1 + 0.000759*z1
# ifdef CICE_MODEL
! Need to match CICE
      t_freeze(s1,z1) = -1.8_r8
!      t_freeze(s1,z1) = depressT*s1
# else
! Temperature or potential temperature??
      t_freeze(s1,z1) = depressT*s1
# endif

# include "set_bounds.h"

# ifdef DISTRIBUTE
      IF (.not. allocated(buffer))    allocate (buffer(1))
      IF (.not. allocated(op_handle)) allocate (op_handle(1))
# endif

# ifdef CICE_MODEL
      ice_dens = SHR_CONST_RHOICE
# else
      ice_dens = rhoice(ng)
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wfr(i,j) = 0.0_r8
        ENDDO
      ENDDO
# ifdef CICE_MODEL
      DO j=Jstr,Jend
        DO i=Istr,Iend
#  ifdef WET_DRY
          IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
#  elif defined MASKING
          IF (rmask(i,j) .ne. 0.0_r8) THEN
#  endif
            thickness = 0.0_r8
            DO k=N(ng),1,-1
              dz = Hz(i,j,k)
              t_fr = t_freeze(t(i,j,k,nnew,isalt),z_r(i,j,k))
              IF ((thickness + dz) > bl_thick) THEN
                fract = (bl_thick - thickness)/dz
                dz = fract*dz
              END IF
              thickness = thickness + dz
              potice = (t_fr - t(i,j,k,nnew,itemp)) * dz
!
!*** if potice < 0, compute again later
!*** if potice > 0, keep on freezing (wfr > 0)
!
              IF ((wfr(i,j) + potice) > 0) THEN
                wfr(i,j) = wfr(i,j) + potice
              ELSE
                potice = -wfr(i,j)
                wfr(i,j) = 0.0_r8
              ENDIF
              t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) +             &
     &                potice/Hz(i,j,k)
              t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt)               &
     &              + (salref-salice)*potice*cp_over_lhfusion/Hz(i,j,k)
!*** jump out if done with surface layer
              IF (z_w(i,j,N(ng))-z_w(i,j,k) > bl_thick) EXIT
            END DO
#  if defined WET_DRY || defined MASKING
          END IF
#  endif
        END DO
      END DO
# else
! Old ice model
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
#  ifdef WET_DRY
            IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
#  elif defined MASKING
            IF (rmask(i,j) .ne. 0.0_r8) THEN
#  endif
              t_fr = t_freeze(t(i,j,k,nnew,isalt),z_r(i,j,k))
              IF (t(i,j,k,nnew,itemp) .lt. t_fr) THEN
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                     (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r)     &
     &                         - depressT * t(i,j,k,nnew,isalt))
                IF (gamma_k .lt. 0.0_r8) THEN
                  print *, 'trouble in ice_frazil', i, j, k,            &
     &             t(i,j,k,nnew,itemp), t(i,j,k,nnew,isalt),            &
     &             t_fr, wfr(i,j), gamma_k, Hz(i,j,k)
                  exit_flag = 9
                END IF
                wfr(i,j) = wfr(i,j) + gamma_k * Hz(i,j,k) *             &
     &                    (rho0 + rho(i,j,k) ) / ice_dens
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) + gamma_k *   &
     &                 (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &                                  (1.0_r8 + gamma_k)
              ELSE IF (wfr(i,j) > 0 .and.                               &
     &                 t(i,j,k,nnew,itemp) .gt. t_fr) THEN
! Use heat at this level to melt some ice from below.
! gamma_k becomes negative here.
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                     (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r)     &
     &                         - depressT * t(i,j,k,nnew,isalt))
                delta_wfr = gamma_k * Hz(i,j,k) *                       &
     &                    (rho0 + rho(i,j,k) ) / rhoice(ng)
                IF ((wfr(i,j) + delta_wfr) > 0) THEN
                  wfr(i,j) = wfr(i,j) + delta_wfr
                ELSE
                  gamma_k = -wfr(i,j) * rhoice(ng) /                    &
     &                     (Hz(i,j,k)*(rho0+rho(i,j,k)))
                  wfr(i,j) = 0.0_r8
                ENDIF
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) + gamma_k *   &
     &                 (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &                                  (1.0_r8 + gamma_k)
              END IF
#  if defined WET_DRY || defined MASKING
            END IF
#  endif
          END DO
          wfr(i,j) = wfr(i,j)/dt(ng)
          IF (wfr(i,j) .lt. 0.0_r8) THEN
            print *, 'trouble in ice_frazil', i, j,                     &
     &         t(i,j,N(ng),nnew,itemp), t(i,j,N(ng),nnew,isalt),        &
     &         wfr(i,j), gamma_k, Hz(i,j,N(ng))
            exit_flag = 9
          END IF
        END DO
      END DO
# endif

# ifdef CICE_MODEL
! Now compute heat flux from ocean and subtract from wfr.
! Using old salt flux, solve for surface salinity s0mk.
#  ifdef ICE_LOG_LAYER
      DO j=Jstr,Jend
        DO i=Istr,Iend
          dztop=z_w(i,j,N(ng))-z_r(i,j,N(ng))
          ice_thick = 0.05_r8+hi(i,j)/MAX(ai(i,j),eps)
          utau = sqrt(sqrt(                                             &
     &             (0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2                &
     &           + (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2                &
     &                   )    )
          utau = max(utau,1.E-4_r8)
!  Need some roughness estimate here
          z0 = max(z0ii*ice_thick,0.01_r8)
          z0 = min(z0,0.1_r8)
!
!     *** Yaglom and Kader formulation for z0t and z0s
!
          zdz0 = dztop/z0   !WPB
          zdz0 = MAX(zdz0,3._r8)

!  Was this:
          rno = utau*0.09_r8/nu
! Mellor
!          rno = utau*z0/nu
          termt = ykf*sqrt(rno)*prt**0.666667_r8
          terms = ykf*sqrt(rno)*prs**0.666667_r8
          cht = utau/(tpr*log(zdz0)/kappa+termt)
          chs = utau/(tpr*log(zdz0)/kappa+terms)
          IF (ai(i,j) .le. min_a(ng)) THEN
            s0mk(i,j) = t(i,j,N(ng),nnew,isalt)
            t0mk(i,j) = t(i,j,N(ng),nnew,itemp)
          ELSE
            s0mk(i,j) = t(i,j,N(ng),nnew,isalt) -                       &
     &                         stflx(i,j,isalt)/chs
            s0mk(i,j) = max(s0mk(i,j),0._r8)
            s0mk(i,j) = min(s0mk(i,j),40._r8)
            t0mk(i,j) = t_freeze(s0mk(i,j),0.0_r8)
          END IF
! convert units to W/m^2
          wfr(i,j) = (wfr(i,j)                                          &
     &        -cht*(t(i,j,N(ng),nnew,itemp)-t0mk(i,j)))*rhocpr
        END DO
      END DO
#  else
      DO j=Jstr,Jend
        DO i=Istr,Iend
          heat = 0.0_r8
          thickness = 0.0_r8
          DO k=N(ng),1,-1
            dz = Hz(i,j,k)
            tfr = t_freeze(t(i,j,k,nnew,isalt),0.0_r8)
            IF ((thickness + dz) < bl_thick) THEN
              thickness = thickness + dz
              heat = heat + (t(i,j,k,nnew,itemp)-tfr)*dz
            ELSE
              fract = (bl_thick - thickness)/dz
!              thickness = thickness + dz*fract
              heat = heat + (t(i,j,k,nnew,itemp)-tfr)*dz*fract
            END IF
            IF (z_w(i,j,N(ng))-z_w(i,j,k) > bl_thick) EXIT
          END DO
! convert units to W/m^2
          wfr(i,j) = (wfr(i,j) - heat)*rhocpr/dt(ng)
        END DO
      END DO
#  endif
# endif

      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wfr)
# ifdef DISTRIBUTE
      buffer(1) = exit_flag
      op_handle(1) = 'MAX'
      CALL mp_reduce (ng, iNLM, 1, buffer, op_handle)
      exit_flag = int(buffer(1))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    wfr)
# endif
      RETURN
      END SUBROUTINE ice_frazil_tile

#endif
      END MODULE ice_frazil_mod

