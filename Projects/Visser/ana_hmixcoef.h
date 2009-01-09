      SUBROUTINE ana_hmixcoef (ng, tile, model)
!
!! svn $Id: ana_hmixcoef.h 999 2007-08-13 14:17:47Z jcwarner $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!                                                                      !
!  This routine rescales horizontal mixing coefficients according      !
!  to the grid size.  Also,  if applicable,  increases horizontal      !
!  in sponge areas.                                                    !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_scalars
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"

      CALL ana_hmixcoef_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                        MIXING(ng) % diff2,                       &
# endif
# ifdef TS_DIF4
     &                        MIXING(ng) % diff4,                       &
# endif
#endif
#ifdef UV_VIS2
     &                        MIXING(ng) % visc2_p,                     &
     &                        MIXING(ng) % visc2_r,                     &
#endif
#ifdef UV_VIS4
     &                        MIXING(ng) % visc4_p,                     &
     &                        MIXING(ng) % visc4_r,                     &
#endif
# ifdef TKE_WAVEDISS
     &                        FORCES(ng) % Wave_dissip,                 &
# endif
     &                        GRID(ng) % grdscl,                        &
     &                        GRID(ng) % xr,                            &
     &                        GRID(ng) % yr)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME( 8),'(a,a)') TRIM(Adir), '/ana_hmixcoef.h'
      END IF

      RETURN
      END SUBROUTINE ana_hmixcoef
!
!***********************************************************************
      SUBROUTINE ana_hmixcoef_tile (ng, model, Istr, Iend, Jstr, Jend,  &
     &                              LBi, UBi, LBj, UBj,                 &
#ifdef SOLVE3D
# ifdef TS_DIF2
     &                              diff2,                              &
# endif
# ifdef TS_DIF4
     &                              diff4,                              &
# endif
#endif
#ifdef UV_VIS2
     &                              visc2_p,                            &
     &                              visc2_r,                            &
#endif
#ifdef UV_VIS4
     &                              visc4_p,                            &
     &                              visc4_r,                            &
#endif
# ifdef TKE_WAVEDISS
     &                              Wave_dissip,                        &
# endif
     &                              grdscl, xr, yr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:,LBj:,:)
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:,LBj:,:)
#  endif
# endif
# ifdef UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc2_r(LBi:,LBj:)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc4_r(LBi:,LBj:)
# endif
# ifdef TKE_WAVEDISS
      real(r8), intent(in) :: Wave_dissip(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
#  ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
#  endif
# endif
# ifdef UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc2_r(LBi:UBi,LBj:UBj)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc4_r(LBi:UBi,LBj:UBj)
# endif
# ifdef TKE_WAVEDISS
      real(r8), intent(in) :: Wave_dissip(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
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
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Iwrk, i, j, itrc
      real(r8) :: cff, cff1, cff2, fac

#include "set_bounds.h"

#ifdef VISC_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal viscosity according to the grid size.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!            
# ifdef UV_VIS2
      cff=visc2(ng)/grdmax(ng)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc2_r(i,j)=cff*grdscl(i,j)
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc2_p(i,j)=cff*(grdscl(i,j  )+grdscl(i-1,j  )+              &
     &                      grdscl(i,j-1)+grdscl(i-1,j-1))
        END DO
      END DO
# endif
# ifdef UV_VIS4
      cff=visc4(ng)/(grdmax(ng)**3)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc4_r(i,j)=cff*grdscl(i,j)**3
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc4_p(i,j)=cff*(grdscl(i,j  )**3+grdscl(i-1,j  )**3+        &
     &                      grdscl(i,j-1)**3+grdscl(i-1,j-1)**3)
        END DO
      END DO
# endif
#endif
#ifdef DIFF_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal diffusion according to the grid size.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!            
# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        cff=tnu2(itrc,ng)/grdmax(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff2(i,j,itrc)=cff*grdscl(i,j)
          END DO
        END DO
      END DO
# endif
# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        cff=tnu4(itrc,ng)/(grdmax(ng)**3)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff4(i,j,itrc)=cff*grdscl(i,j)**3
          END DO
        END DO
      END DO
# endif
#endif
#ifdef SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------
!
!! User modifiable section.  Please specify the appropiate sponge area
!! by increasing its horizontal mixing coefficients.
!!
# if defined VISSER
      DO i=IstrR,IendR
        DO j=JstrR,JendR
          cff=visc2(ng)+Wave_dissip(i,j)*100000.0_r8
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
#  if defined no_use_this
      DO j=JstrR,JendR

        cff=visc2(ng)*1.0_r8
        visc2_r(54,j)=cff
        visc2_p(54,j)=cff

        cff=visc2(ng)*6.0_r8
        visc2_r(55,j)=cff
        visc2_p(55,j)=cff

        cff=visc2(ng)*6.0_r8
        visc2_r(56,j)=cff
        visc2_p(56,j)=cff

        cff=visc2(ng)*6.0_r8
        visc2_r(57,j)=cff
        visc2_p(57,j)=cff

        cff=visc2(ng)*8.0_r8
        visc2_r(58,j)=cff
        visc2_p(58,j)=cff

        cff=visc2(ng)*8.0_r8
        visc2_r(59,j)=cff
        visc2_p(59,j)=cff

        cff=visc2(ng)*8.0_r8
        visc2_r(60,j)=cff
        visc2_p(60,j)=cff

        cff=visc2(ng)*6._r8
        visc2_r(61,j)=cff
        visc2_p(61,j)=cff

        cff=visc2(ng)*4._r8
        visc2_r(62,j)=cff
        visc2_p(62,j)=cff

        cff=visc2(ng)*4._r8
        visc2_r(63,j)=cff
        visc2_p(63,j)=cff

        cff=visc2(ng)*2._r8
        visc2_r(64,j)=cff
        visc2_p(64,j)=cff

        cff=visc2(ng)*2._r8
        visc2_r(65,j)=cff
        visc2_p(65,j)=cff
      END DO
#  endif
# elif defined VISSER_notsure
      DO i=IstrR,IendR
        DO j=JstrR,JendR
          cff=visc2(ng)
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
      DO j=JstrR,JendR

        cff=visc2(ng)*1.0_r8
        visc2_r(54,j)=cff
        visc2_p(54,j)=cff

        cff=visc2(ng)*1.0_r8
        visc2_r(55,j)=cff
        visc2_p(55,j)=cff

        cff=visc2(ng)*1.0_r8
        visc2_r(56,j)=cff
        visc2_p(56,j)=cff

        cff=visc2(ng)*1.0_r8
        visc2_r(57,j)=cff
        visc2_p(57,j)=cff

        cff=visc2(ng)*1.0_r8
        visc2_r(58,j)=cff
        visc2_p(58,j)=cff

        cff=visc2(ng)*1.0_r8
        visc2_r(59,j)=cff
        visc2_p(59,j)=cff

        cff=visc2(ng)*1.0_r8
        visc2_r(60,j)=cff
        visc2_p(60,j)=cff

        cff=visc2(ng)*1._r8
        visc2_r(61,j)=cff
        visc2_p(61,j)=cff

        cff=visc2(ng)*1._r8
        visc2_r(62,j)=cff
        visc2_p(62,j)=cff

        cff=visc2(ng)*1._r8
        visc2_r(63,j)=cff
        visc2_p(63,j)=cff

        cff=visc2(ng)*1._r8
        visc2_r(64,j)=cff
        visc2_p(64,j)=cff

        cff=visc2(ng)*1._r8
        visc2_r(65,j)=cff
        visc2_p(65,j)=cff

      END DO
# endif
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!            
# if defined EW_PERIODIC || defined NS_PERIODIC
#  ifdef UV_VIS2
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        visc2_r)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        visc2_p)
#  endif
#  ifdef UV_VIS4
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        visc4_r)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        visc4_p)
#  endif
#  ifdef SOLVE3D
#   ifdef TS_DIF2
      DO itrc=1,NT(ng)
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          diff2(:,:,itrc))
      END DO
#   endif
#   ifdef TS_DIF4
      DO itrc=1,NT(ng)
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          diff4(:,:,itrc))
      END DO
#   endif
#  endif
# endif
# ifdef DISTRIBUTE
#  ifdef UV_VIS2
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    visc2_r, visc2_p)
#  endif
#  ifdef UV_VIS4
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    visc4_r, visc4_p)
#  endif
#  ifdef SOLVE3D
#   ifdef TS_DIF2
      CALL mp_exchange3d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    diff2)
#   endif
#   ifdef TS_DIF4
      CALL mp_exchange3d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    diff4)
#   endif
#  endif
# endif

#endif
      RETURN
      END SUBROUTINE ana_hmixcoef_tile
