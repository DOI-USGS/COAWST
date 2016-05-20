#include "cppdefs.h"
      MODULE set_tides_mod
#if defined NONLINEAR && \
  ( defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES )
!
!svn $Id: set_tides.F 1382 2011-02-19 01:44:37Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        Robert Hetland   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine adds tidal elevation (m) and tidal currents (m/s) to   !
!  sea surface height and 2D momentum climatologies, respectively.     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_tides
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_tides (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_tides
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 11)
# endif
      CALL set_tides_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     NTC(ng),                                     &
     &                     GRID(ng) % angler,                           &
# ifdef POT_TIDES
     &                     TIDES(ng) % POT_Tamp,                        &
     &                     TIDES(ng) % POT_Tphase,                      &
     &                     TIDES(ng) % Ptide,                           &
# endif
# ifdef TIDES_ASTRO
     &                     GRID(ng) % latr,                             &
     &                     TIDES(ng) % Vu_sat,                          &
     &                     TIDES(ng) % f_sat,                           &
# endif
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef SSH_TIDES
     &                     TIDES(ng) % SSH_Tamp,                        &
     &                     TIDES(ng) % SSH_Tphase,                      &
# endif
# ifdef UV_TIDES
     &                     TIDES(ng) % UV_Tangle,                       &
     &                     TIDES(ng) % UV_Tphase,                       &
     &                     TIDES(ng) % UV_Tmajor,                       &
     &                     TIDES(ng) % UV_Tminor,                       &
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                     TIDES(ng) % SinOmega,                        &
     &                     TIDES(ng) % CosOmega,                        &
# endif
# ifdef UV_WAVEDRAG
     &                     TIDES(ng) % sum_ubot,                        &
     &                     TIDES(ng) % sum_vbot,                        &
     &                     TIDES(ng) % filt_ubot,                       &
     &                     TIDES(ng) % filt_vbot,                       &
# endif
     &                     TIDES(ng) % Tperiod)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 11)
# endif

      RETURN
      END SUBROUTINE set_tides
!
!***********************************************************************
      SUBROUTINE set_tides_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NTC,                                   &
     &                           angler,                                &
# ifdef POT_TIDES
     &                           POT_Tamp, POT_Tphase, Ptide,           &
# endif
# ifdef TIDES_ASTRO
     &                           latr, Vu_sat, f_sat,                   &
# endif
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef SSH_TIDES
     &                           SSH_Tamp, SSH_Tphase,                  &
# endif
# ifdef UV_TIDES
     &                           UV_Tangle, UV_Tphase,                  &
     &                           UV_Tmajor, UV_Tminor,                  &
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                           SinOmega, CosOmega,                    &
# endif
# ifdef UV_WAVEDRAG
     &                           sum_ubot, sum_vbot,                    &
     &                           filt_ubot, filt_vbot,                  &
# endif
     &                           Tperiod)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_ncparam
      USE mod_scalars
# ifdef TIDES_ASTRO
      USE mod_tides, ONLY : tide_astro
# endif
# ifdef UV_WAVEDRAG
      USE mod_tides, ONLY : nstep_hour, hour_steps
      USE mod_ocean
      USE mod_stepping, ONLY : nrhs
# endif
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_boundary
# endif
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variables declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: NTC
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  ifdef TIDES_ASTRO
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(inout) :: Vu_sat(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_sat(LBi:,LBj:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Tperiod(MTC)
#  ifdef SSH_TIDES
      real(r8), intent(inout) :: SSH_Tamp(LBi:,LBj:,:)
      real(r8), intent(inout) :: SSH_Tphase(LBi:,LBj:,:)
#  endif
#  ifdef UV_TIDES
      real(r8), intent(in) :: UV_Tangle(LBi:,LBj:,:)
      real(r8), intent(inout) :: UV_Tmajor(LBi:,LBj:,:)
      real(r8), intent(inout) :: UV_Tminor(LBi:,LBj:,:)
      real(r8), intent(inout) :: UV_Tphase(LBi:,LBj:,:)
#  endif
#  ifdef POT_TIDES
      real(r8), intent(inout) :: POT_Tamp(LBi:,LBj:,:)
      real(r8), intent(inout) :: POT_Tphase(LBi:,LBj:,:)
      real(r8), intent(out) :: Ptide(LBi:,LBj:)
#  endif
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(inout) :: SinOmega(:)
      real(r8), intent(inout) :: CosOmega(:)
#  endif
#  ifdef UV_WAVEDRAG
      real(r8), intent(inout) :: sum_ubot(LBi:,LBj:)
      real(r8), intent(inout) :: sum_vbot(LBi:,LBj:)
      real(r8), intent(inout) :: filt_ubot(LBi:,LBj:)
      real(r8), intent(inout) :: filt_vbot(LBi:,LBj:)
#  endif
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#  ifdef TIDES_ASTRO
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Vu_sat(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(inout) :: f_sat(LBi:UBi,LBj:UBj,MTC)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Tperiod(MTC)
#  ifdef SSH_TIDES
      real(r8), intent(in) :: SSH_Tamp(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(in) :: SSH_Tphase(LBi:UBi,LBj:UBj,MTC)
#  endif
#  ifdef UV_TIDES
      real(r8), intent(in) :: UV_Tangle(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(in) :: UV_Tmajor(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(in) :: UV_Tminor(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(in) :: UV_Tphase(LBi:UBi,LBj:UBj,MTC)
#  endif
#  ifdef POT_TIDES
      real(r8), intent(inout) :: POT_Tamp(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(inout) :: POT_Tphase(LBi:UBi,LBj:UBj,MTC)
      real(r8), intent(out) :: Ptide(LBi:UBi,LBj:UBj)
#  endif
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(inout) :: SinOmega(MTC)
      real(r8), intent(inout) :: CosOmega(MTC)
#  endif
#  ifdef UV_WAVEDRAG
      real(r8), intent(inout) :: sum_ubot(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: sum_vbot(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: filt_ubot(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: filt_vbot(LBi:UBi,LBj:UBj)
#  endif
# endif
!
!  Local variables declarations.
!
      logical :: update

# ifdef DISTRIBUTE
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, itide, j

      real(r8) :: Cangle, Cphase, Sangle, Sphase
      real(r8) :: angle, cff, phase, omega, ramp, tide_time
      real(r8) :: bry_cor, bry_pgr, bry_str, bry_val

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Etide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Utide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vtide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk

# include "set_bounds.h"

# ifdef DISTRIBUTE
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      ILB=BOUNDS(ng)%LBi(-1)
      IUB=BOUNDS(ng)%UBi(-1)
      JLB=BOUNDS(ng)%LBj(-1)
      JUB=BOUNDS(ng)%UBj(-1)
# endif
!
!=======================================================================
!  Process tidal forcing used at the lateral boundaries.
!=======================================================================
!
!  In refinement nesting applications, tidal forcing is only necessary
!  in the coarser large scale grid.
!
      NEEDED : IF (.not.(RefinedGrid(ng).and.                           &
     &                    RefineScale(ng).gt.0)) THEN
!
!  Set time-ramping parameter.
!
# ifdef RAMP_TIDES
        ramp=TANH((tdays(ng)-dstart)/1.0_r8)
# else
        ramp=1.0_r8
# endif
# ifdef TIDES_ASTRO
        tide_time = time(ng)-tide_start*day2sec
        call tide_astro(tide_time, Vu_sat, f_sat, latr,                 &
     &                  ng, tile, LBi, UBi, LBj, UBj)
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!-----------------------------------------------------------------------
!  Compute harmonic used to detide output fields.
!-----------------------------------------------------------------------
!
        cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
        DO itide=1,NTC
          IF (Tperiod(itide).gt.0.0_r8) THEN
            omega=cff/Tperiod(itide)
            SinOmega(itide)=SIN(omega)
            CosOmega(itide)=COS(omega)
          ELSE
            SinOmega(itide)=0.0_r8
            CosOmega(itide)=0.0_r8
          END IF
        END DO
# endif
# ifdef SSH_TIDES
!
!-----------------------------------------------------------------------
!  Add tidal elevation (m) to sea surface height climatology.
!-----------------------------------------------------------------------
!
        IF (LsshCLM(ng)) THEN
          Etide(:,:)=0.0_r8
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          DO itide=1,NTC
            IF (Tperiod(itide).gt.0.0_r8) THEN
              omega=cff/Tperiod(itide)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
#  ifdef TIDES_ASTRO
!  Convert Vu_sat from cycles to radians
                  Etide(i,j)=Etide(i,j)+                                &
     &                       ramp*f_sat(i,j,itide)*SSH_Tamp(i,j,itide)* &
     &                       COS(-SSH_Tphase(i,j,itide)+                &
     &                       2.0_r8*pi*Vu_sat(i,j,itide))
#  else
                    Etide(i,j)=Etide(i,j)+                              &
     &                         ramp*SSH_Tamp(i,j,itide)*                &
     &                         COS(omega-SSH_Tphase(i,j,itide))
#  endif
#  ifdef MASKING
                    Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
                END DO
              END DO
            END IF
          END DO
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              CLIMA(ng)%ssh(i,j)=CLIMA(ng)%ssh(i,j)+Etide(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              CLIMA(ng)%ssh)
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        CLIMA(ng)%ssh)
#  endif
        END IF
!
!  If appropriate, load tidal forcing into boundary arrays.  The "zeta"
!  boundary arrays are important for the Flather or reduced physics
!  boundary conditions for 2D momentum. To avoid having two boundary
!  points for these arrays, the values of "zeta_west" and "zeta_east"
!  are averaged at u-points.  Similarly, the values of "zeta_south"
!  and "zeta_north" is averaged at v-points. Noticed that these
!  arrays are also used for the clamped conditions for free-surface.
!  This averaging is less important for that type ob boundary
!  conditions.
!
        IF (LBC(iwest,isFsur,ng)%acquire.or.                            &
     &      LBC(iwest,isUbar,ng)%acquire.or.                            &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            Etide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
                omega=cff/Tperiod(itide)
                DO j=JstrR,JendR
                  DO i=Istr-1,Istr
#  ifdef TIDES_ASTRO
!  Convert Vu_sat from cycles to radians
                    Etide(i,j)=Etide(i,j)+                              &
     &                       ramp*f_sat(i,j,itide)*SSH_Tamp(i,j,itide)* &
     &                       COS(-SSH_Tphase(i,j,itide)+                &
     &                       2.0_r8*pi*Vu_sat(i,j,itide))
#  else
                    Etide(i,j)=Etide(i,j)+                              &
     &                         ramp*SSH_Tamp(i,j,itide)*                &
     &                         COS(omega-SSH_Tphase(i,j,itide))
#  endif
#  ifdef MASKING
                    Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
                  END DO
                END DO
              END IF
            END DO
            DO j=JstrR,JendR
#  ifdef ADD_FSOBC
              BOUNDARY(ng)%zeta_west(j)=BOUNDARY(ng)%zeta_west(j)+      &
     &                                  0.5_r8*(Etide(Istr-1,j)+        &
     &                                          Etide(Istr  ,j))
#  else
              BOUNDARY(ng)%zeta_west(j)=0.5_r8*(Etide(Istr-1,j)+        &
     &                                          Etide(Istr  ,j))
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_west)
#  endif
        END IF
!
        IF (LBC(ieast,isFsur,ng)%acquire.or.                            &
     &      LBC(ieast,isUbar,ng)%acquire.or.                            &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            Etide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
                omega=cff/Tperiod(itide)
                DO j=JstrR,JendR
                  DO i=Iend,Iend+1
#  ifdef TIDES_ASTRO
!  Convert Vu_sat from cycles to radians
                    Etide(i,j)=Etide(i,j)+                              &
     &                       ramp*f_sat(i,j,itide)*SSH_Tamp(i,j,itide)* &
     &                       COS(-SSH_Tphase(i,j,itide)+                &
     &                       2.0_r8*pi*Vu_sat(i,j,itide))
#  else
                    Etide(i,j)=Etide(i,j)+                              &
     &                         ramp*SSH_Tamp(i,j,itide)*                &
     &                         COS(omega-SSH_Tphase(i,j,itide))
#  endif
#  ifdef MASKING
                    Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
                  END DO
                END DO
              END IF
            END DO
            DO j=JstrR,JendR
#  ifdef ADD_FSOBC
              BOUNDARY(ng)%zeta_east(j)=BOUNDARY(ng)%zeta_east(j)+      &
     &                                  0.5_r8*(Etide(Iend  ,j)+        &
     &                                          Etide(Iend+1,j))
#  else
              BOUNDARY(ng)%zeta_east(j)=0.5_r8*(Etide(Iend  ,j)+        &
     &                                          Etide(Iend+1,j))
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_east)
#  endif
        END IF
!
        IF (LBC(isouth,isFsur,ng)%acquire.or.                           &
     &      LBC(isouth,isUbar,ng)%acquire.or.                           &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            Etide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
                omega=cff/Tperiod(itide)
                DO j=Jstr-1,Jstr
                  DO i=IstrR,IendR
#  ifdef TIDES_ASTRO
!  Convert Vu_sat from cycles to radians
                    Etide(i,j)=Etide(i,j)+                              &
     &                       ramp*f_sat(i,j,itide)*SSH_Tamp(i,j,itide)* &
     &                       COS(-SSH_Tphase(i,j,itide)+                &
     &                       2.0_r8*pi*Vu_sat(i,j,itide))
#  else
                    Etide(i,j)=Etide(i,j)+                              &
     &                         ramp*SSH_Tamp(i,j,itide)*                &
     &                         COS(omega-SSH_Tphase(i,j,itide))
#  endif
#  ifdef MASKING
                    Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
                  END DO
                END DO
              END IF
            END DO
            DO i=IstrR,IendR
#  ifdef ADD_FSOBC
              BOUNDARY(ng)%zeta_south(i)=BOUNDARY(ng)%zeta_south(i)+    &
     &                                   0.5_r8*(Etide(i,Jstr-1)+       &
     &                                           Etide(i,Jstr  ))
#  else
              BOUNDARY(ng)%zeta_south(i)=0.5_r8*(Etide(i,Jstr-1)+       &
     &                                           Etide(i,Jstr  ))
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_south)
#  endif
        END IF
!
        IF (LBC(inorth,isFsur,ng)%acquire.or.                           &
     &      LBC(inorth,isUbar,ng)%acquire.or.                           &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            Etide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
                omega=cff/Tperiod(itide)
                DO j=Jend,Jend+1
                  DO i=IstrR,IendR
#  ifdef TIDES_ASTRO
!  Convert Vu_sat from cycles to radians
                    Etide(i,j)=Etide(i,j)+                              &
     &                       ramp*f_sat(i,j,itide)*SSH_Tamp(i,j,itide)* &
     &                       COS(-SSH_Tphase(i,j,itide)+                &
     &                       2.0_r8*pi*Vu_sat(i,j,itide))
#  else
                    Etide(i,j)=Etide(i,j)+                              &
     &                         ramp*SSH_Tamp(i,j,itide)*                &
     &                         COS(omega-SSH_Tphase(i,j,itide))
#  endif
#  ifdef MASKING
                    Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
                  END DO
                END DO
              END IF
            END DO
            DO i=IstrR,IendR
#  ifdef ADD_FSOBC
              BOUNDARY(ng)%zeta_north(i)=BOUNDARY(ng)%zeta_north(i)+    &
     &                                   0.5_r8*(Etide(i,Jend  )+       &
     &                                           Etide(i,Jend+1))
#  else
              BOUNDARY(ng)%zeta_north(i)=0.5_r8*(Etide(i,Jend  )+       &
     &                                           Etide(i,Jend+1))
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%zeta_north)
#  endif
        END IF
# endif

# ifdef POT_TIDES
!
!-----------------------------------------------------------------------
!  Compute tidal potential (m) for forcing
!-----------------------------------------------------------------------
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Ptide(i,j)=0.0_r8
          END DO
        END DO
        cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
        DO itide=1,NTC
          IF (Tperiod(itide).gt.0.0_r8) THEN
#  ifdef TIDES_ASTRO
            tide_time = time(ng)-tide_start*day2sec
#  endif
            omega=cff/Tperiod(itide)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
#  ifdef TIDES_ASTRO
!  Convert Vu_sat from cycles to radians
                Ptide(i,j)=Ptide(i,j)+                                  &
     &                     ramp*POT_Tamp(i,j,itide)*                    &
     &                     COS(-POT_Tphase(i,j,itide)+                  &
     &                     2.0_r8*pi*Vu_sat(i,j,itide))
#  else
                Ptide(i,j)=Ptide(i,j)+                                  &
     &                     ramp*POT_Tamp(i,j,itide)*                    &
     &                     COS(omega-POT_Tphase(i,j,itide))
#  endif
#  ifdef MASKING
                Ptide(i,j)=Ptide(i,j)*rmask(i,j)
#  endif
              END DO
            END DO
          END IF
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Ptide)
        END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic(ng),               &
     &                      NSperiodic(ng),                             &
     &                      Ptide)
#  endif
# endif

# if defined UV_TIDES
#  ifdef ADD_M2OBC
!
!-----------------------------------------------------------------------
!  Add tidal currents (m/s) to 2D momentum climatologies.
!-----------------------------------------------------------------------
!
        IF (Lm2CLM(ng)) THEN
          Utide(:,:)=0.0_r8
          Vtide(:,:)=0.0_r8
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          DO itide=1,NTC
            IF (Tperiod(itide).gt.0.0_r8) THEN
#   ifdef TIDES_ASTRO
              tide_time = time(ng)-tide_start*day2sec
#   endif
              omega=cff/Tperiod(itide)
              DO j=MIN(JstrR,Jstr-1),JendR
                DO i=MIN(IstrR,Istr-1),IendR
                  angle=UV_Tangle(i,j,itide)-angler(i,j)
                  Cangle=COS(angle)
                  Sangle=SIN(angle)
#   ifdef TIDES_ASTRO
                  phase=-UV_Tphase(i,j,itide)
!  Convert Vu_sat from cycles to radians
                  phase = phase + 2.0_r8*pi*Vu_sat(i,j,itide)
#   else
                  phase=omega-UV_Tphase(i,j,itide)
#   endif
                  Cphase=COS(phase)
                  Sphase=SIN(phase)
                  Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-         &
     &                      UV_Tminor(i,j,itide)*Sangle*Sphase
                  Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+         &
     &                      UV_Tminor(i,j,itide)*Cangle*Sphase
#   ifdef TIDES_ASTRO
                  Uwrk(i,j) = Uwrk(i,j)*f_sat(i,j,itide)
                  Vwrk(i,j) = Vwrk(i,j)*f_sat(i,j,itide)
#   endif
                END DO
              END DO
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  Utide(i,j)=Utide(i,j)+                                &
     &                   ramp*0.5_r8*(Uwrk(i-1,j)+Uwrk(i,j))
#   ifdef MASKING
                  Utide(i,j)=Utide(i,j)*umask(i,j)
#   endif
                END DO
              END DO
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  Vtide(i,j)=(Vtide(i,j)+                               &
     &                    ramp*0.5_r8*(Vwrk(i,j-1)+Vwrk(i,j)))
#   ifdef MASKING
                  Vtide(i,j)=Vtide(i,j)*vmask(i,j)
#   endif
                END DO
              END DO
            END IF
          END DO
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
          DO j=JstrR,JendR
            DO i=Istr,IendR
              CLIMA(ng)%ubarclm(i,j)=CLIMA(ng)%ubarclm(i,j)+Utide(i,j)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              CLIMA(ng)%vbarclm(i,j)=CLIMA(ng)%vbarclm(i,j)+Vtide(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              CLIMA(ng)%ubarclm)
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              CLIMA(ng)%vbarclm)
          END IF
#   ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iNLM, 2,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        CLIMA(ng)%ubarclm,                        &
     &                        CLIMA(ng)%vbarclm)
#   endif
       END IF
#  endif
!
!  If appropriate, load tidal forcing into boundary arrays.
!
        IF (LBC(iwest,isUbar,ng)%acquire.and.                           &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            Utide(:,:)=0.0_r8
            Vtide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
#  ifdef TIDES_ASTRO
                tide_time = time(ng)-tide_start*day2sec
#  endif
                omega=cff/Tperiod(itide)
                DO j=MIN(JstrR,Jstr-1),JendR
                  DO i=Istr-1,Istr
                    angle=UV_Tangle(i,j,itide)-angler(i,j)
                    Cangle=COS(angle)
                    Sangle=SIN(angle)
#  ifdef TIDES_ASTRO
                    phase=-UV_Tphase(i,j,itide)
!  Convert Vu_sat from cycles to radians
                    phase = phase + 2.0_r8*pi*Vu_sat(i,j,itide)
#  else
                    phase=omega-UV_Tphase(i,j,itide)
#  endif
                    Cphase=COS(phase)
                    Sphase=SIN(phase)
                    Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-       &
     &                        UV_Tminor(i,j,itide)*Sangle*Sphase
                    Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+       &
     &                        UV_Tminor(i,j,itide)*Cangle*Sphase
#  ifdef TIDES_ASTRO
                    Uwrk(i,j) = Uwrk(i,j)*f_sat(i,j,itide)
                    Vwrk(i,j) = Vwrk(i,j)*f_sat(i,j,itide)
#  endif
                  END DO
                END DO
                DO j=JstrR,JendR
                  Utide(Istr,j)=Utide(Istr,j)+                          &
     &                       ramp*0.5_r8*(Uwrk(Istr-1,j)+Uwrk(Istr,j))
#  ifdef MASKING
                  Utide(Istr,j)=Utide(Istr,j)*umask(Istr,j)
#  endif
                END DO
                DO j=Jstr,JendR
                  Vtide(Istr-1,j)=Vtide(Istr-1,j)+                      &
     &                   ramp*0.5_r8*(Vwrk(Istr-1,j-1)+Vwrk(Istr-1,j))
#  ifdef MASKING
                  Vtide(Istr-1,j)=Vtide(Istr-1,j)*vmask(Istr-1,j)
#  endif
                END DO
              END IF
            END DO
            DO j=Jstr,JendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%ubar_west(j)=BOUNDARY(ng)%ubar_west(j)+      &
     &                                  Utide(Istr,j)
#  else
              BOUNDARY(ng)%ubar_west(j)=Utide(Istr,j)
#  endif
            END DO
            DO j=Jstr,JendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%vbar_west(j)=BOUNDARY(ng)%vbar_west(j)+      &
     &                                  Vtide(Istr-1,j)
#  else
              BOUNDARY(ng)%vbar_west(j)=Vtide(Istr-1,j)
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_west)
          CALL mp_boundary (ng, iNLM, Jstr,  JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_west)
#  endif
        END IF
        IF (LBC(ieast,isUbar,ng)%acquire.and.                           &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            Utide(:,:)=0.0_r8
            Vtide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
#  ifdef TIDES_ASTRO
                tide_time = time(ng)-tide_start*day2sec
#  endif
                omega=cff/Tperiod(itide)
                DO j=MIN(JstrR,Jstr-1),JendR
                  DO i=Iend,Iend+1
                    angle=UV_Tangle(i,j,itide)-angler(i,j)
                    Cangle=COS(angle)
                    Sangle=SIN(angle)
#  ifdef TIDES_ASTRO
                    phase=-UV_Tphase(i,j,itide)
!  Convert Vu_sat from cycles to radians
                    phase = phase + 2.0_r8*pi*Vu_sat(i,j,itide)
#  else
                    phase=omega-UV_Tphase(i,j,itide)
#  endif
                    Cphase=COS(phase)
                    Sphase=SIN(phase)
                    Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-       &
     &                        UV_Tminor(i,j,itide)*Sangle*Sphase
                    Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+       &
     &                        UV_Tminor(i,j,itide)*Cangle*Sphase
#  ifdef TIDES_ASTRO
                    Uwrk(i,j) = Uwrk(i,j)*f_sat(i,j,itide)
                    Vwrk(i,j) = Vwrk(i,j)*f_sat(i,j,itide)
#  endif
                  END DO
                END DO
                DO j=JstrR,JendR
                  Utide(Iend+1,j)=Utide(Iend+1,j)+                      &
     &                     ramp*0.5_r8*(Uwrk(Iend,j)+Uwrk(Iend+1,j))
#  ifdef MASKING
                  Utide(Iend+1,j)=Utide(Iend+1,j)*umask(Iend+1,j)
#  endif
                END DO
                DO j=Jstr,JendR
                  Vtide(Iend+1,j)=Vtide(Iend+1,j)+                      &
     &                   ramp*0.5_r8*(Vwrk(Iend+1,j-1)+Vwrk(Iend+1,j))
#  ifdef MASKING
                  Vtide(Iend+1,j)=Vtide(Iend+1,j)*vmask(Iend+1,j)
#  endif
                END DO
              END IF
            END DO
            DO j=JstrR,JendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%ubar_east(j)=BOUNDARY(ng)%ubar_east(j)+      &
     &                                  Utide(Iend+1,j)
#  else
              BOUNDARY(ng)%ubar_east(j)=Utide(Iend+1,j)
#  endif
            END DO
            DO j=Jstr,JendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%vbar_east(j)=BOUNDARY(ng)%vbar_east(j)+      &
     &                                  Vtide(Iend+1,j)
#  else
              BOUNDARY(ng)%vbar_east(j)=Vtide(Iend+1,j)
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, JstrR, JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_east)
          CALL mp_boundary (ng, iNLM, Jstr,  JendR,                     &
     &                      JLB, JUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_east)
#  endif
        END IF

        IF (LBC(isouth,isUbar,ng)%acquire.and.                          &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            Utide(:,:)=0.0_r8
            Vtide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
# ifdef TIDES_ASTRO
                tide_time = time(ng)-tide_start*day2sec
# endif
                omega=cff/Tperiod(itide)
                DO j=Jstr-1,Jstr
                  DO i=MIN(IstrR,Istr-1),IendR
                    angle=UV_Tangle(i,j,itide)-angler(i,j)
                    Cangle=COS(angle)
                    Sangle=SIN(angle)
# ifdef TIDES_ASTRO
                    phase=-UV_Tphase(i,j,itide)
!  Convert Vu_sat from cycles to radians
                    phase = phase + 2.0_r8*pi*Vu_sat(i,j,itide)
# else
                    phase=omega-UV_Tphase(i,j,itide)
# endif
                    Cphase=COS(phase)
                    Sphase=SIN(phase)
                    Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-       &
     &                        UV_Tminor(i,j,itide)*Sangle*Sphase
                    Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+       &
     &                        UV_Tminor(i,j,itide)*Cangle*Sphase
# ifdef TIDES_ASTRO
                    Uwrk(i,j) = Uwrk(i,j)*f_sat(i,j,itide)
                    Vwrk(i,j) = Vwrk(i,j)*f_sat(i,j,itide)
# endif
                  END DO
                END DO
                DO i=Istr,IendR
                  Utide(i,Jstr-1)=Utide(i,Jstr-1)+                      &
     &                 ramp*0.5_r8*(Uwrk(i-1,Jstr-1)+Uwrk(i,Jstr-1))
#  ifdef MASKING
                  Utide(i,Jstr-1)=Utide(i,Jstr-1)*umask(i,Jstr-1)
#  endif
                END DO
                DO i=IstrR,IendR
                  Vtide(i,Jstr)=Vtide(i,Jstr)+                          &
     &                      ramp*0.5_r8*(Vwrk(i,Jstr-1)+Vwrk(i,Jstr))
#  ifdef MASKING
                  Vtide(i,Jstr)=Vtide(i,Jstr)*vmask(i,Jstr)
#  endif
                END DO
              END IF
            END DO
            DO i=IstrR,IendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%ubar_south(i)=BOUNDARY(ng)%ubar_south(i)+    &
     &                                   Utide(i,Jstr-1)
#  else
              BOUNDARY(ng)%ubar_south(i)=Utide(i,Jstr-1)
#  endif
            END DO
            DO i=IstrR,IendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%vbar_south(i)=BOUNDARY(ng)%vbar_south(i)+    &
     &                                   Vtide(i,Jstr)
#  else
              BOUNDARY(ng)%vbar_south(i)=Vtide(i,Jstr)
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, Istr,  IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_south)
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_south)
#  endif
        END IF
        IF (LBC(inorth,isUbar,ng)%acquire.and.                          &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          update=.FALSE.
          cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            Utide(:,:)=0.0_r8
            Vtide(:,:)=0.0_r8
            DO itide=1,NTC
              IF (Tperiod(itide).gt.0.0_r8) THEN
# ifdef TIDES_ASTRO
                tide_time = time(ng)-tide_start*day2sec
# endif
                omega=cff/Tperiod(itide)
                DO j=Jend,Jend+1
                  DO i=MIN(IstrR,Istr-1),IendR
                    angle=UV_Tangle(i,j,itide)-angler(i,j)
                    Cangle=COS(angle)
                    Sangle=SIN(angle)
# ifdef TIDES_ASTRO
                    phase=-UV_Tphase(i,j,itide)
!  Convert Vu_sat from cycles to radians
                    phase = phase + 2.0_r8*pi*Vu_sat(i,j,itide)
# else
                    phase=omega-UV_Tphase(i,j,itide)
# endif
                    Cphase=COS(phase)
                    Sphase=SIN(phase)
                    Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-           &
     &                        UV_Tminor(i,j,itide)*Sangle*Sphase
                    Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+           &
     &                        UV_Tminor(i,j,itide)*Cangle*Sphase
# ifdef TIDES_ASTRO
                    Uwrk(i,j) = Uwrk(i,j)*f_sat(i,j,itide)
                    Vwrk(i,j) = Vwrk(i,j)*f_sat(i,j,itide)
# endif
                  END DO
                END DO
                DO i=Istr,IendR
                  Utide(i,Jend+1)=Utide(i,Jend+1)+                          &
     &                    ramp*0.5_r8*(Uwrk(i-1,Jend+1)+Uwrk(i,Jend+1))
#  ifdef MASKING
                  Utide(i,Jend+1)=Utide(i,Jend+1)*umask(i,Jend+1)
#  endif
                END DO
                DO i=IstrR,IendR
                  Vtide(i,Jend+1)=Vtide(i,Jend+1)+                          &
     &                      ramp*0.5_r8*(Vwrk(i,Jend)+Vwrk(i,Jend+1))
#  ifdef MASKING
                  Vtide(i,Jend+1)=Vtide(i,Jend+1)*vmask(i,Jend+1)
#  endif
                END DO
              END IF
            END DO

            DO i=Istr,IendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%ubar_north(i)=BOUNDARY(ng)%ubar_north(i)+    &
     &                                   Utide(i,Jend+1)
#  else
              BOUNDARY(ng)%ubar_north(i)=Utide(i,Jend+1)
#  endif
            END DO
            DO i=IstrR,IendR
#  ifdef ADD_M2OBC
              BOUNDARY(ng)%vbar_north(i)=BOUNDARY(ng)%vbar_north(i)+    &
     &                                   Vtide(i,Jend+1)
#  else
              BOUNDARY(ng)%vbar_north(i)=Vtide(i,Jend+1)
#  endif
            END DO
            update=.TRUE.
          END IF
#  ifdef DISTRIBUTE
          CALL mp_boundary (ng, iNLM, Istr,  IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%ubar_north)
          CALL mp_boundary (ng, iNLM, IstrR, IendR,                     &
     &                      ILB, IUB, 1, 1, update,                     &
     &                      BOUNDARY(ng)%vbar_north)
#  endif
        END IF
# endif
# ifdef UV_WAVEDRAG
        IF ((MOD(iic(ng)-ntstart(ng),nstep_hour).eq.0).and.             &
     &      (iic(ng).ne.ntstart(ng))) THEN
!          IF (MOD(iic(ng)-ntstart(ng),nstep_hour*24).eq.0) THEN
          IF (hour_steps == 24 .or. hour_steps == 0) THEN
! At end of day, add in the 25th value, then compute new filt_ubot value,
! and finally add in the first value for the next day.
            hour_steps = 1
            DO j=JstrR,JendR
              DO i=Istr,Iend+1
                sum_ubot(i,j)=sum_ubot(i,j)+OCEAN(ng)%u(i,j,1,nrhs(ng))
                filt_ubot(i,j)=0.04_r8*sum_ubot(i,j)
                sum_ubot(i,j)=OCEAN(ng)%u(i,j,1,nrhs(ng))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=IstrR,IendR
                sum_vbot(i,j)=sum_vbot(i,j)+OCEAN(ng)%v(i,j,1,nrhs(ng))
                filt_vbot(i,j)=0.04_r8*sum_vbot(i,j)
                sum_vbot(i,j)=OCEAN(ng)%v(i,j,1,nrhs(ng))
              END DO
            END DO
          ELSE
            hour_steps = hour_steps + 1
            DO j=JstrR,JendR
              DO i=Istr,Iend+1
                sum_ubot(i,j)=sum_ubot(i,j)+OCEAN(ng)%u(i,j,1,nrhs(ng))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=IstrR,IendR
                sum_vbot(i,j)=sum_vbot(i,j)+OCEAN(ng)%v(i,j,1,nrhs(ng))
              END DO
            END DO
          END IF
        END IF
# endif
      END IF NEEDED

      RETURN
      END SUBROUTINE set_tides_tile
#endif
      END MODULE set_tides_mod
