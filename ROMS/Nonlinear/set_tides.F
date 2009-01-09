#include "cppdefs.h"
      MODULE set_tides_mod
#if defined NONLINEAR && (defined SSH_TIDES || defined UV_TIDES)
!
!svn $Id: set_tides.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group        Robert Hetland   !
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

      PRIVATE
      PUBLIC  :: set_tides

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_tides (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
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
     &                     NTC(ng),                                     &
     &                     GRID(ng) % angler,                           &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
     &                     TIDES(ng) % Tperiod,                         &
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
# ifdef ZCLIMATOLOGY
     &                     CLIMA(ng) % ssh,                             &
# endif
# ifdef M2CLIMATOLOGY
     &                     CLIMA(ng) % ubarclm,                         &
     &                     CLIMA(ng) % vbarclm,                         &
# endif
     &                     BOUNDARY(ng) % zeta_east,                    &
     &                     BOUNDARY(ng) % zeta_west,                    &
     &                     BOUNDARY(ng) % zeta_south,                   &
     &                     BOUNDARY(ng) % zeta_north,                   &
     &                     BOUNDARY(ng) % ubar_east,                    &
     &                     BOUNDARY(ng) % vbar_east,                    &
     &                     BOUNDARY(ng) % ubar_west,                    &
     &                     BOUNDARY(ng) % vbar_west,                    &
     &                     BOUNDARY(ng) % ubar_south,                   &
     &                     BOUNDARY(ng) % vbar_south,                   &
     &                     BOUNDARY(ng) % ubar_north,                   &
     &                     BOUNDARY(ng) % vbar_north)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 11)
# endif
      RETURN
      END SUBROUTINE set_tides
!
!***********************************************************************
      SUBROUTINE set_tides_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NTC,                                   &
     &                           angler,                                &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
     &                           Tperiod,                               &
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
# ifdef ZCLIMATOLOGY
     &                           ssh,                                   &
# endif
# ifdef M2CLIMATOLOGY
     &                           ubarclm, vbarclm,                      &
# endif
     &                           zeta_east, zeta_west,                  &
     &                           zeta_south, zeta_north,                &
     &                           ubar_east, vbar_east,                  &
     &                           ubar_west, vbar_west,                  &
     &                           ubar_south, vbar_south,                &
     &                           ubar_north, vbar_north)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_boundary
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variables declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: NTC
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Tperiod(MTC)
#  ifdef SSH_TIDES
      real(r8), intent(in) :: SSH_Tamp(LBi:,LBj:,:)
      real(r8), intent(in) :: SSH_Tphase(LBi:,LBj:,:)
#  endif
#  ifdef UV_TIDES
      real(r8), intent(in) :: UV_Tangle(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tmajor(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tminor(LBi:,LBj:,:)
      real(r8), intent(in) :: UV_Tphase(LBi:,LBj:,:)
#  endif
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(inout) :: SinOmega(:)
      real(r8), intent(inout) :: CosOmega(:)
#  endif
#  ifdef ZCLIMATOLOGY
      real(r8), intent(inout) :: ssh(LBi:,LBj:)
#  endif
#  ifdef M2CLIMATOLOGY
      real(r8), intent(inout) :: ubarclm(LBi:,LBj:)
      real(r8), intent(inout) :: vbarclm(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: zeta_east(LOWER_BOUND_J:)
      real(r8), intent(inout) :: zeta_west(LOWER_BOUND_J:)
      real(r8), intent(inout) :: zeta_south(LOWER_BOUND_I:)
      real(r8), intent(inout) :: zeta_north(LOWER_BOUND_I:)
      real(r8), intent(inout) :: ubar_east(LOWER_BOUND_J:)
      real(r8), intent(inout) :: vbar_east(LOWER_BOUND_J:)
      real(r8), intent(inout) :: ubar_west(LOWER_BOUND_J:)
      real(r8), intent(inout) :: vbar_west(LOWER_BOUND_J:)
      real(r8), intent(inout) :: ubar_south(LOWER_BOUND_I:)
      real(r8), intent(inout) :: vbar_south(LOWER_BOUND_I:)
      real(r8), intent(inout) :: ubar_north(LOWER_BOUND_I:)
      real(r8), intent(inout) :: vbar_north(LOWER_BOUND_I:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
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
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(inout) :: SinOmega(MTC)
      real(r8), intent(inout) :: CosOmega(MTC)
#  endif
#  ifdef ZCLIMATOLOGY
      real(r8), intent(inout) :: ssh(LBi:UBi,LBj:UBj)
#  endif
#  ifdef M2CLIMATOLOGY
      real(r8), intent(inout) :: ubarclm(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vbarclm(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: zeta_east(ETA_DIM)
      real(r8), intent(inout) :: zeta_west(ETA_DIM)
      real(r8), intent(inout) :: zeta_south(XI_DIM)
      real(r8), intent(inout) :: zeta_north(XI_DIM)
      real(r8), intent(inout) :: ubar_east(ETA_DIM)
      real(r8), intent(inout) :: vbar_east(ETA_DIM)
      real(r8), intent(inout) :: ubar_west(ETA_DIM)
      real(r8), intent(inout) :: vbar_west(ETA_DIM)
      real(r8), intent(inout) :: ubar_south(XI_DIM)
      real(r8), intent(inout) :: vbar_south(XI_DIM)
      real(r8), intent(inout) :: ubar_north(XI_DIM)
      real(r8), intent(inout) :: vbar_north(XI_DIM)
# endif
!
!  Local variables declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      logical :: update

# ifdef DISTRIBUTE
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: i, itide, j

      real(r8) :: Cangle, Cphase, Sangle, Sphase
      real(r8) :: angle, cff, phase, omega, ramp
      real(r8) :: bry_cor, bry_pgr, bry_str, bry_val

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Etide
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Utide
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Uwrk
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Vtide
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Vwrk

# include "set_bounds.h"

# ifdef DISTRIBUTE
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif
!
!  Set time-ramping parameter.
!
# ifdef RAMP_TIDES
      ramp=TANH((tdays(ng)-dstart)/1.0_r8)
# else
      ramp=1.0_r8
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
      Etide(:,:)=0.0_r8
      cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
      DO itide=1,NTC
        IF (Tperiod(itide).gt.0.0_r8) THEN
          omega=cff/Tperiod(itide)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Etide(i,j)=Etide(i,j)+                                    &
     &                   ramp*SSH_Tamp(i,j,itide)*                      &
     &                   COS(omega-SSH_Tphase(i,j,itide))
#  ifdef MASKING
              Etide(i,j)=Etide(i,j)*rmask(i,j)
#  endif
            END DO
          END DO
        END IF
      END DO
#  if defined ZCLIMATOLOGY && defined ADD_FSOBC
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ssh(i,j)=ssh(i,j)+Etide(i,j)
        END DO
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ssh)
#   endif
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ssh)
#   endif
#  endif
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
#  if defined WEST_FSOBC || defined WEST_M2OBC
      update=.FALSE.
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
#   ifdef ADD_FSOBC
          zeta_west(j)=zeta_west(j)+                                    &
     &                 0.5_r8*(Etide(Istr-1,j)+Etide(Istr,j))
#   else
          zeta_west(j)=0.5_r8*(Etide(Istr-1,j)+Etide(Istr,j))
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, JstrR, JendR, JLB, JUB, 1, 1, update, &
     &                  zeta_west)
#   endif
#  endif
#  if defined EAST_FSOBC || defined EAST_M2OBC
      update=.FALSE.
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
#   ifdef ADD_FSOBC
          zeta_east(j)=zeta_east(j)+                                    &
     &                 0.5_r8*(Etide(Iend,j)+Etide(Iend+1,j))
#   else
          zeta_east(j)=0.5_r8*(Etide(Iend,j)+Etide(Iend+1,j))
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, JstrR, JendR, JLB, JUB, 1, 1, update, &
     &                  zeta_east)
#   endif
#  endif
#  if defined SOUTH_FSOBC || defined SOUTH_M2OBC
      update=.FALSE.
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
#   ifdef ADD_FSOBC
          zeta_south(i)=zeta_south(i)+                                  &
     &                  0.5_r8*(Etide(i,Jstr-1)+Etide(i,Jstr))
#   else
          zeta_south(i)=0.5_r8*(Etide(i,Jstr-1)+Etide(i,Jstr))
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, IstrR, IendR, ILB, IUB, 1, 1, update, &
     &                  zeta_south)
#   endif
#  endif
#  if defined NORTH_FSOBC || defined NORTH_M2OBC
      update=.FALSE.
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
#   ifdef ADD_FSOBC
          zeta_north(i)=zeta_north(i)+                                  &
     &                  0.5_r8*(Etide(i,Jend)+Etide(i,Jend+1))
#   else
          zeta_north(i)=0.5_r8*(Etide(i,Jend)+Etide(i,Jend+1))
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, IstrR, IendR, ILB, IUB, 1, 1, update, &
     &                  zeta_north)
#   endif
#  endif
# endif
# if defined UV_TIDES
!
!-----------------------------------------------------------------------
!  Add tidal currents (m/s) to 2D momentum climatologies.
!-----------------------------------------------------------------------
!
      Utide(:,:)=0.0_r8
      Vtide(:,:)=0.0_r8
      cff=2.0_r8*pi*(time(ng)-tide_start*day2sec)
      DO itide=1,NTC
        IF (Tperiod(itide).gt.0.0_r8) THEN
          omega=cff/Tperiod(itide)
          DO j=MIN(JstrR,Jstr-1),JendR
            DO i=MIN(IstrR,Istr-1),IendR
              angle=UV_Tangle(i,j,itide)-angler(i,j)
              Cangle=COS(angle)
              Sangle=SIN(angle)
              phase=omega-UV_Tphase(i,j,itide)
              Cphase=COS(phase)
              Sphase=SIN(phase)
              Uwrk(i,j)=UV_Tmajor(i,j,itide)*Cangle*Cphase-             &
     &                  UV_Tminor(i,j,itide)*Sangle*Sphase
              Vwrk(i,j)=UV_Tmajor(i,j,itide)*Sangle*Cphase+             &
     &                  UV_Tminor(i,j,itide)*Cangle*Sphase
            END DO
          END DO
          DO j=JstrR,JendR
            DO i=Istr,IendR
              Utide(i,j)=Utide(i,j)+                                    &
     &                   ramp*0.5_r8*(Uwrk(i-1,j)+Uwrk(i,j))
#  ifdef MASKING
              Utide(i,j)=Utide(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              Vtide(i,j)=(Vtide(i,j)+                                   &
     &                    ramp*0.5_r8*(Vwrk(i,j-1)+Vwrk(i,j)))
#  ifdef MASKING
              Vtide(i,j)=Vtide(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
        END IF
      END DO
#  if defined M2CLIMATOLOGY && defined ADD_M2OBC
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubarclm(i,j)=ubarclm(i,j)+Utide(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbarclm(i,j)=vbarclm(i,j)+Vtide(i,j)
        END DO
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubarclm)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbarclm)
#   endif
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ubarclm, vbarclm)
#   endif
#  endif
!
!  If appropriate, load tidal forcing into boundary arrays.
!
#  ifdef WEST_M2OBC
      update=.FALSE.
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
#   ifdef ADD_M2OBC
          ubar_west(j)=ubar_west(j)+Utide(Istr,j)
#   else
          ubar_west(j)=Utide(Istr,j)
#   endif
        END DO
        DO j=Jstr,JendR
#   ifdef ADD_M2OBC
          vbar_west(j)=vbar_west(j)+Vtide(Istr-1,j)
#   else
          vbar_west(j)=Vtide(Istr-1,j)
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, JstrR, JendR, JLB, JUB, 1, 1, update, &
     &                  ubar_west)
      CALL mp_boundary (ng, iNLM, Jstr,  JendR, JLB, JUB, 1, 1, update, &
     &                  vbar_west)
#   endif
#  endif
#  ifdef EAST_M2OBC
      update=.FALSE.
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
#   ifdef ADD_M2OBC
          ubar_east(j)=ubar_east(j)+Utide(Iend+1,j)
#   else
          ubar_east(j)=Utide(Iend+1,j)
#   endif
        END DO
        DO j=Jstr,JendR
#   ifdef ADD_M2OBC
          vbar_east(j)=vbar_east(j)+Vtide(Iend+1,j)
#   else
          vbar_east(j)=Vtide(Iend+1,j)
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, JstrR, JendR, JLB, JUB, 1, 1, update, &
     &                  ubar_east)
      CALL mp_boundary (ng, iNLM, Jstr,  JendR, JLB, JUB, 1, 1, update, &
     &                  vbar_east)
#   endif
#  endif
#  ifdef SOUTH_M2OBC
      update=.FALSE.
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,IendR
#   ifdef ADD_M2OBC
          ubar_south(i)=ubar_south(i)+Utide(i,Jstr-1)
#   else
          ubar_south(i)=Utide(i,Jstr-1)
#   endif
        END DO
        DO i=IstrR,IendR
#   ifdef ADD_M2OBC
          vbar_south(i)=vbar_south(i)+Vtide(i,Jstr)
#   else
          vbar_south(i)=Vtide(i,Jstr)
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, Istr,  IendR, ILB, IUB, 1, 1, update, &
     &                  ubar_south)
      CALL mp_boundary (ng, iNLM, IstrR, IendR, ILB, IUB, 1, 1, update, &
     &                  vbar_south)
#   endif
#  endif
#  ifdef NORTH_M2OBC
      update=.FALSE.
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,IendR
#   ifdef ADD_M2OBC
          ubar_north(i)=ubar_north(i)+Utide(i,Jend+1)
#   else
          ubar_north(i)=Utide(i,Jend+1)
#   endif
        END DO
        DO i=IstrR,IendR
#   ifdef ADD_M2OBC
          vbar_north(i)=vbar_north(i)+Vtide(i,Jend+1)
#   else
          vbar_north(i)=Vtide(i,Jend+1)
#   endif
        END DO
        update=.TRUE.
      END IF
#   ifdef DISTRIBUTE
      CALL mp_boundary (ng, iNLM, Istr,  IendR, ILB, IUB, 1, 1, update, &
     &                  ubar_north)
      CALL mp_boundary (ng, iNLM, IstrR, IendR, ILB, IUB, 1, 1, update, &
     &                  vbar_north)
#   endif
#  endif
# endif
      RETURN
      END SUBROUTINE set_tides_tile
#endif
      END MODULE set_tides_mod
