#include "cppdefs.h"
      MODULE ice_frazil_mod
#if defined ICE_MODEL && defined ICE_THERMO
!
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
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
     &                            Hz, z_r, rho, t, wfr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
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
      real(r8) :: delta_wfr

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
! Temperature or potential temperature??
      t_freeze(s1,z1) = -0.0543*s1

# include "set_bounds.h"

# ifdef DISTRIBUTE
      IF (.not. allocated(buffer))    allocate (buffer(1))
      IF (.not. allocated(op_handle)) allocate (op_handle(1))
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wfr(i,j) = 0.0_r8
        ENDDO
      ENDDO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
# ifdef WET_DRY
            IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
# elif defined MASKING
            IF (rmask(i,j) .ne. 0.0_r8) THEN
# endif
              t_fr = t_freeze(t(i,j,k,nnew,isalt),z_r(i,j,k))
              IF (t(i,j,k,nnew,itemp) .lt. t_fr) THEN
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                     (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r)     &
     &                         + 0.0543_r8 * t(i,j,k,nnew,isalt))
                IF (gamma_k .lt. 0.0_r8) THEN
                  print *, 'trouble in ice_frazil', i, j, k,            &
     &             t(i,j,k,nnew,itemp), t(i,j,k,nnew,isalt),            &
     &             t_fr, wfr(i,j), gamma_k, Hz(i,j,k)
                  exit_flag = 9
                END IF
                wfr(i,j) = wfr(i,j) + gamma_k * Hz(i,j,k) *             &
     &                    (1000.0_r8 + rho(i,j,k) ) / rhoice(ng)
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
     &                         + 0.0543_r8 * t(i,j,k,nnew,isalt))
                delta_wfr = gamma_k * Hz(i,j,k) *                       &
     &                    (rho0 + rho(i,j,k) ) / rhoice(ng)
!                print *, 'ICE_FRAZIL before', t_fr, t(i,j,k,nnew,itemp),&
!     &               t(i,j,k,nnew,isalt)
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
!                print *, 'ICE_FRAZIL after', gamma_k,                   &
!     &               t(i,j,k,nnew,itemp), t(i,j,k,nnew,isalt)
              END IF
# if defined WET_DRY || defined MASKING
            END IF
# endif
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
# ifdef CICE
! Now compute heat flux from ocean and subtract from wfr
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

