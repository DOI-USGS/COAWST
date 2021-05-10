#include "cppdefs.h"

      MODULE marsh_sed_erosion_mod
!
#if defined SEDIMENT  && defined MARSH_DYNAMICS \
        && defined MARSH_SED_EROSION  
!
!svn $Id: marsh_sed_erosion.F 429 2019-04-20 17:30:26Z arango $
!======================================================================!
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==============================================Tarandeep S. Kalra======!
!==============================================Julia M. Moriarty=======!
!================================================Neil K. Ganju  =======!
!==============================================John C. Warner==========!
!                                                                      ! 
!  This routine calculates the export of sediment through bedload      !
!  exchange from the marsh cells when lateral wave thrust acts upon    !
!  them. For high resolution (of the order of 1m or less, one can      !
!  compute the lateral retreat of marsh cells and convert them to open !
!  water cells.                                                        !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC  :: marsh_sed_erosion

      CONTAINS
!
!***********************************************************************
      SUBROUTINE marsh_sed_erosion (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_vegetation 
      USE mod_vegarr
      USE mod_sedbed
      USE mod_stepping
      USE marsh_wave_thrust_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL marsh_sed_erosion_tile (ng, tile,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % h,                              &
     &                       nstp(ng), nnew(ng),                        &
     &                       VEG(ng) % umask_marsh,                     &
     &                       VEG(ng) % vmask_marsh,                     &
     &                       VEG(ng) % Thrust_xi,                       &
     &                       VEG(ng) % Thrust_eta,                      &
     &                       VEG(ng) % Thrust_total,                    &
     &                       VEG(ng) % marsh_mask,                      &
     &                       VEG(ng) % marsh_flux_out,                  &
# ifdef MARSH_RETREAT
     &                       VEG(ng) % marsh_retreat,                   &
# endif
# ifdef MARSH_STOCH
     &                       VEG(ng) % marsh_stoch,                     &
# endif
     &                    SEDBED(ng) % bed_frac,                        &
     &                    SEDBED(ng) % bed,                             &
     &                    SEDBED(ng) % bed_mass)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE marsh_sed_erosion
!
!***********************************************************************
      SUBROUTINE marsh_sed_erosion_tile (ng, tile,                      &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             om_r, on_r,                          &
     &                             om_u, om_v,                          &
     &                             on_u, on_v,                          &
     &                             pm, pn,                              &
     &                             h,                                   &
     &                             nstp, nnew,                          &
     &                             umask_marsh, vmask_marsh,            &
     &                             Thrust_xi, Thrust_eta,               &
     &                             Thrust_total,                        &
     &                             marsh_mask, marsh_flux_out,          &
# ifdef MARSH_RETREAT
     &                             marsh_retreat,                       &
# endif
# ifdef MARSH_STOCH
     &                             marsh_stoch,                         &
# endif
     &                             bed_frac, bed, bed_mass)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_vegetation 
      USE mod_vegarr
      USE mod_sediment
      USE marsh_wave_thrust_mod
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE bc_3d_mod, ONLY : bc_r3d_tile
!      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d,         &
     &                            mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
!
      real(r8), intent(in)    :: umask_marsh(LBi:,LBj:)
      real(r8), intent(in)    :: vmask_marsh(LBi:,LBj:)
      real(r8), intent(in)    :: Thrust_xi(LBi:,LBj:)
      real(r8), intent(in)    :: Thrust_eta(LBi:,LBj:)
      real(r8), intent(in)    :: Thrust_total(LBi:,LBj:)
      real(r8), intent(inout) :: marsh_mask(LBi:,LBj:)
      real(r8), intent(inout) :: marsh_flux_out(LBi:,LBj:,:)
# ifdef MARSH_RETREAT
      real(r8), intent(inout) :: marsh_retreat(LBi:,LBj:)
# endif
# ifdef MARSH_STOCH
      real(r8), intent(in) :: marsh_stoch(LBi:,LBj:)
# endif
!
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
!
# else
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
!
      real(r8), intent(in)    :: umask_marsh(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: vmask_marsh(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: Thrust_xi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: Thrust_eta(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: Thrust_total(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: marsh_mask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: marsh_flux_out(LBi:UBi,LBj:UBj,NST)
# ifdef MARSH_RETREAT
      real(r8), intent(inout) :: marsh_retreat(LBi:UBi,LBj:UBj)
# endif
# ifdef MARSH_STOCH
      real(r8), intent(in) :: marsh_stoch(LBi:UBi,LBj:UBj)
# endif
!
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
!
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ised
      real(r8), parameter :: half=0.5_r8, one=1.0_r8 
      real(r8), parameter :: eps = 1.0E-14_r8 !
      real(r8), parameter :: kN_N=1000.0_r8 ! convert Thrust from kN to N 
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff_retreat 
      real(r8) :: bed_thick_old, bathy
# if defined MARSH_RETREAT 
      real(r8) :: cff3_ceil
# endif
!
# include "set_bounds.h"
!
# if defined MARSH_RETREAT 
      cff3_ceil=0.0_r8
# endif
!
      cff1=kN_N*kfac_marsh(ng)*dt(ng)
!
      J_LOOP: DO j=Jstr,Jend
        DO i=Istr,Iend
          SED_LOOP: DO ised=1,NST
!
! Convert bed mass from kg/sq.m to kg in each cell
!
            cff=bed_mass(i,j,1,nnew,ised)*om_r(i,j)*on_r(i,j)
!
!# ifdef MARSH_STOCH 
!            cff1=kN_N*kfac_marsh(ng)*marsh_stoch(i,j)*dt(ng)
!# else
!# endif 
!
!---------------------------------------------------------------------
! Add bed mass to sea cell providing horizontal face thrust 
! you cannot add more than the available amount from adjacent cell 
! how to take into account the fraction of bed mass coming into the cell
! Rules-cff4 is mass and can never be negative 
!---------------------------------------------------------------------
!
            cff2=half*(one-SIGN(one,Thrust_xi(i,j  )))
            cff3=half*(one+SIGN(one,Thrust_xi(i+1,j)))
            cff4=cff1*                                                  &
     &                ( Thrust_xi(i+1,j  )*bed_frac(i+1,j,1,ised)*      &
     &                                         on_u(i+1,j  )*cff3       &
     &                 -Thrust_xi(i,j)*bed_frac(i-1,j,1,ised)*          &
     &                                         on_u(i,j)*cff2     )     
            cff=cff+cff4
! 
!---------------------------------------------------------------------
! Subtract bed mass from marsh cell if thrust acts on horizontal face 
! Marsh cell
!---------------------------------------------------------------------
!
            cff2=half*(one+SIGN(one,Thrust_xi(i,j  )))
            cff3=half*(one-SIGN(one,Thrust_xi(i+1,j)))
            cff4=cff1*bed_frac(i,j,1,ised)*                             &
     &                ( Thrust_xi(i,j  )*on_u(i,j  )*cff2               &
     &                 -Thrust_xi(i+1,j)*on_u(i+1,j)*cff3         )     
            cff=cff-cff4
!
!---------------------------------------------------------------------
! Ensure that cff is not negative 
!---------------------------------------------------------------------
!
            cff=MAX(cff,0.0_r8)
!
!---------------------------------------------------------------------
! Sediment flux for each sediment class out of marsh cells
!---------------------------------------------------------------------
!
            marsh_flux_out(i,j,ised)=cff4+marsh_flux_out(i,j,ised)
!
!---------------------------------------------------------------------
! Add bed mass to sea cell if it provides vertical face thrust 
!---------------------------------------------------------------------
!
            cff2=half*(one-SIGN(one,Thrust_eta(i,j  )))
            cff3=half*(one+SIGN(one,Thrust_eta(i,j+1)))
            cff4=cff1*                                                  &
     &                ( Thrust_eta(i,j+1)*bed_frac(i,j+1,1,ised)*       &
     &                                        om_v(i,j+1)*cff3          &
     &                   -Thrust_eta(i,j  )*bed_frac(i,j-1,1,ised)*     &
     &                                        om_v(i,j  )*cff2    )
            cff=cff+cff4
!
!---------------------------------------------------------------------
! Subtract bed mass from marsh cell if thrust acts on vertical face 
! Only can subtract from a cell if the sediment class exists  
!---------------------------------------------------------------------
!
            cff2=half*(one+SIGN(one,Thrust_eta(i,j  )))
            cff3=half*(one-SIGN(one,Thrust_eta(i,j+1)))
            cff4=cff1*bed_frac(i,j,1,ised)*                             &
     &                ( Thrust_eta(i,j  )*om_v(i,j  )*cff2              &
     &                 -Thrust_eta(i,j+1)*om_v(i,j+1)*cff3        )
            cff=cff-cff4
!
!---------------------------------------------------------------------
! Ensure that cff is not negative 
!---------------------------------------------------------------------
!
            cff=MAX(cff,0.0_r8)
!
!---------------------------------------------------------------------
! Sediment flux for each sediment class out of marsh cells
!---------------------------------------------------------------------
!
            marsh_flux_out(i,j,ised)=cff4+marsh_flux_out(i,j,ised)
!
!---------------------------------------------------------------------
! Update bed mass 
!---------------------------------------------------------------------
!
            bed_mass(i,j,1,nnew,ised)=cff*pm(i,j)*pn(i,j)
!        
          END DO SED_LOOP
!
# if defined MARSH_RETREAT
!
!---------------------------------------------------------------------
! recalculate thickness and fractions for all layers.
!---------------------------------------------------------------------
!
          cff3_ceil=MAX(CEILING(Thrust_total(i,j)),0)
!
          bed_thick_old=bed(i,j,1,ithck)
          cff3=0.0_r8
          DO ised=1,NST
            cff3=cff3+bed_mass(i,j,1,nnew,ised)
          END DO
          IF (cff3.eq.0.0_r8) THEN
            cff3=eps
          END IF
          bed(i,j,1,ithck)=0.0_r8
          DO ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
            bed(i,j,1,ithck)=MAX(bed(i,j,1,ithck)+                      &
     &                       bed_mass(i,j,1,nnew,ised)/                 &
     &                       (Srho(ised,ng)*                            &
     &                       (1.0_r8-bed(i,j,1,iporo))),0.0_r8)
          END DO
! 
! Determine the cell that is under marsh wave attack 
!
          cff1=(ABS(bed(i,j,1,ithck)-bed_thick_old))*cff3_ceil
          marsh_retreat(i,j)=marsh_retreat(i,j)+cff1
          IF(marsh_retreat(i,j).gt.SCARP_HGHT(ng))THEN
            marsh_mask(i,j)=0.0_r8
          END IF
# endif
!
        END DO
      END DO J_LOOP
!
!---------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!---------------------------------------------------------------------
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  marsh_mask)
      DO ised=1,NST
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    marsh_flux_out(:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
      END DO
      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        marsh_mask)
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, NST,                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    marsh_flux_out) 
!       
      CALL mp_exchange4d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:))
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed) 
# endif
!
# if defined MARSH_RETREAT
!
!---------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!---------------------------------------------------------------------
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  marsh_retreat)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        marsh_retreat)
#  endif
# endif 
!
      END SUBROUTINE marsh_sed_erosion_tile
#endif
      END MODULE marsh_sed_erosion_mod

