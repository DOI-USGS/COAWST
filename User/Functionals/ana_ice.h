      SUBROUTINE ana_ice (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for ice fields                 !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ice
      USE mod_ocean
      USE mod_ncparam
!
      implicit none


      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_ice_tile (ng, tile, model,                               &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ICE(ng) % ui,                              &
     &                       ICE(ng) % vi,                              &
     &                       ICE(ng) % uie,                             &
     &                       ICE(ng) % vie,                             &
     &                       ICE(ng) % ai,                              &
     &                       ICE(ng) % hi,                              &
     &                       ICE(ng) % hsn,                             &
     &                       ICE(ng) % ti,                              &
     &                       ICE(ng) % ageice,                          &
#ifdef MELT_PONDS
     &                       ICE(ng) % apond,                           &
     &                       ICE(ng) % hpond,                           &
#endif
     &                       ICE(ng) % sig11,                           &
     &                       ICE(ng) % sig22,                           &
     &                       ICE(ng) % sig12,                           &
#ifdef NCEP_FLUXES
     &                       FORCES(ng) % wg2_d,                        &
     &                       FORCES(ng) % cd_d,                         &
     &                       FORCES(ng) % ch_d,                         &
     &                       FORCES(ng) % ce_d,                         &
     &                       FORCES(ng) % wg2_m,                        &
     &                       FORCES(ng) % cd_m,                         &
     &                       FORCES(ng) % ch_m,                         &
     &                       FORCES(ng) % ce_m,                         &
     &                       FORCES(ng) % rhoa_n,                       &
#endif
     &                       ICE(ng) % tis,                             &
     &                       ICE(ng) % s0mk,                            &
     &                       ICE(ng) % t0mk,                            &
     &                       ICE(ng) % utau_iw,                         &
     &                       ICE(ng) % chu_iw,                          &
#ifdef ICE_BIO
     &                       ICE(ng) % IcePhL,                          &
     &                       ICE(ng) % IceNO3,                          &
     &                       ICE(ng) % IceNH4,                          &
     &                       ICE(ng) % IceLog,                          &
#endif
     &                       OCEAN(ng) % t )
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(46)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_ice
!
!***********************************************************************
      SUBROUTINE ana_ice_tile (ng, tile, model,                         &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ui, vi, uie, vie, ai, hi, hsn,       &
     &                             ti, ageice,                          &
#ifdef MELT_PONDS
     &                             apond, hpond,                        &
#endif
     &                             sig11, sig22, sig12,                 &
#ifdef NCEP_FLUXES
     &                             wg2_d, cd_d, ch_d, ce_d,             &
     &                             wg2_m, cd_m, ch_m, ce_m,             &
     &                             rhoa_n,                              &
#endif
     &                             tis, s0mk, t0mk, utau_iw, chu_iw,    &
#ifdef ICE_BIO
     &                             IcePhL, IceNO3,                      &
     &                             IceNH4, IceLog,                      &
#endif
     &                             t )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: uie(LBi:,LBj:,:)
      real(r8), intent(inout) :: vie(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
# ifdef MELT_PONDS
      real(r8), intent(inout) :: apond(LBi:,LBj:,:)
      real(r8), intent(inout) :: hpond(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: sig11(LBi:,LBj:,:)
      real(r8), intent(inout) :: sig22(LBi:,LBj:,:)
      real(r8), intent(inout) :: sig12(LBi:,LBj:,:)
# ifdef NCEP_FLUXES
      real(r8), intent(inout) :: wg2_d(LBi:,LBj:)
      real(r8), intent(inout) :: cd_d(LBi:,LBj:)
      real(r8), intent(inout) :: ch_d(LBi:,LBj:)
      real(r8), intent(inout) :: ce_d(LBi:,LBj:)
      real(r8), intent(inout) :: wg2_m(LBi:,LBj:)
      real(r8), intent(inout) :: cd_m(LBi:,LBj:)
      real(r8), intent(inout) :: ch_m(LBi:,LBj:)
      real(r8), intent(inout) :: ce_m(LBi:,LBj:)
      real(r8), intent(inout) :: rhoa_n(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: tis(LBi:,LBj:)
      real(r8), intent(inout) :: s0mk(LBi:,LBj:)
      real(r8), intent(inout) :: t0mk(LBi:,LBj:)
      real(r8), intent(inout) :: utau_iw(LBi:,LBj:)
      real(r8), intent(inout) :: chu_iw(LBi:,LBj:)
# ifdef ICE_BIO
      real(r8), intent(inout) :: IcePhL(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNO3(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNH4(LBi:,LBj:,:)
      integer, intent(inout) :: IceLog(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: uie(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: vie(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
# ifdef MELT_PONDS
      real(r8), intent(inout) :: apond(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hpond(LBi:UBi,LBj:UBj,2)
# endif
      real(r8), intent(inout) :: sig11(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sig22(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sig12(LBi:UBi,LBj:UBj,2)
# ifdef NCEP_FLUXES
      real(r8), intent(inout) :: wg2_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cd_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ch_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ce_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wg2_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cd_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ch_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ce_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rhoa_n(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: tis(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: s0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: utau_iw(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: chu_iw(LBi:UBi,LBj:UBj)
# ifdef ICE_BIO
      real(r8), intent(inout) :: IcePhL(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNO3(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNH4(LBi:UBi,LBj:UBj,2)
      integer, intent(inout) :: IceLog(LBi:UBi,LBj:UBj,2)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: r2

#include "set_bounds.h"

#ifdef ICE_BASIN
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ui(i,j,1) = 0._r8
          uie(i,j,1) = 0._r8
          ui(i,j,2) = ui(i,j,1)
          uie(i,j,2) = uie(i,j,1)
        ENDDO
      ENDDO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vi(i,j,1) = 0._r8
          vie(i,j,1) = 0._r8
          vi(i,j,2) = vi(i,j,1)
          vie(i,j,2) = vie(i,j,1)
        ENDDO
      ENDDO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ai(i,j,1) = 1._r8
          hi(i,j,1) = 2._r8
#ifdef NO_SNOW
          hsn(i,j,1) = 0.0_r8
#else
          hsn(i,j,1) = 0.2_r8
#endif
          ti(i,j,1) = -5._r8
# ifdef MELT_PONDS
          apond(i,j,1) = 0._r8
          hpond(i,j,1) = 0._r8
# endif
          ageice(i,j,1) = 0._r8
          sig11(i,j,1) = 0._r8
          sig22(i,j,1) = 0._r8
          sig12(i,j,1) = 0._r8
          ai(i,j,2) = ai(i,j,1)
          hi(i,j,2) = hi(i,j,1)
          hsn(i,j,2) = hsn(i,j,1)
          ti(i,j,2) = ti(i,j,1)
# ifdef MELT_PONDS
          apond(i,j,2) = apond(i,j,1)
          hpond(i,j,2) = hpond(i,j,1)
# endif
          ageice(i,j,2) = ageice(i,j,1)
          sig11(i,j,2) = sig11(i,j,1)
          sig22(i,j,2) = sig22(i,j,1)
          sig12(i,j,2) = sig12(i,j,1)
# ifdef NCEP_FLUXES
          wg2_d(i,j) = 1._r8
          cd_d(i,j) = 0.00319_r8
          ch_d(i,j) = 1.0E-4_r8
          ce_d(i,j) = 1.0E-4_r8
          wg2_m(i,j) = 1._r8
          cd_m(i,j) = 0.00319_r8
          ch_m(i,j) = 1.0E-4_r8
          ce_m(i,j) = 1.0E-4_r8
          rhoa_n(i,j) = 1.4_r8
# endif
          tis(i,j) = -10._r8
          s0mk(i,j) = t(i,j,N(ng),1,isalt)
          t0mk(i,j) = t(i,j,N(ng),1,itemp)
          utau_iw(i,j) = 0.001_r8
          chu_iw(i,j) = 0.001125_r8
#elif defined ICE_OCEAN_1D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ui(i,j,1) = 0.0_r8
          uie(i,j,1) = 0.0_r8
          ui(i,j,2) = ui(i,j,1)
          uie(i,j,2) = uie(i,j,1)
        ENDDO
      ENDDO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vi(i,j,1) = 0.0_r8
          vie(i,j,1) = 0.0_r8
          vi(i,j,2) = vi(i,j,1)
          vie(i,j,2) = vie(i,j,1)
        ENDDO
      ENDDO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ai(i,j,1) = 0._r8
          hi(i,j,1) = 0._r8
          hsn(i,j,1) = 0.2_r8
          ti(i,j,1) = -5._r8
# ifdef MELT_PONDS
          apond(i,j,1) = 0._r8
          hpond(i,j,1) = 0._r8
# endif
          ageice(i,j,1) = 0._r8
          sig11(i,j,1) = 0._r8
          sig22(i,j,1) = 0._r8
          sig12(i,j,1) = 0._r8
          ai(i,j,2) = ai(i,j,1)
          hi(i,j,2) = hi(i,j,1)
          hsn(i,j,2) = hsn(i,j,1)
          ti(i,j,2) = ti(i,j,1)
# ifdef MELT_PONDS
          apond(i,j,2) = apond(i,j,1)
          hpond(i,j,2) = hpond(i,j,1)
# endif
          ageice(i,j,2) = ageice(i,j,1)
          sig11(i,j,2) = sig11(i,j,1)
          sig22(i,j,2) = sig22(i,j,1)
          sig12(i,j,2) = sig12(i,j,1)
# ifdef NCEP_FLUXES
          wg2_d(i,j) = 1._r8
          cd_d(i,j) = 0.00319_r8
          ch_d(i,j) = 1.0E-4_r8
          ce_d(i,j) = 1.0E-4_r8
          wg2_m(i,j) = 1._r8
          cd_m(i,j) = 0.00319_r8
          ch_m(i,j) = 1.0E-4_r8
          ce_m(i,j) = 1.0E-4_r8
          rhoa_n(i,j) = 1.4_r8
# endif
          tis(i,j) = -10._r8
          s0mk(i,j) = t(i,j,N(ng),1,isalt)
          t0mk(i,j) = t(i,j,N(ng),1,itemp)
          utau_iw(i,j) = 0.001_r8
          chu_iw(i,j) = 0.001125_r8
# ifdef ICE_BIO
          IcePhL(i,j,1) = 0._r8
          IceNO3(i,j,1) = 0._r8
          IceNH4(i,j,1) = 0._r8
          IceLog(i,j,1) = -1
          IcePhL(i,j,2) = IcePhL(i,j,1)
          IceNO3(i,j,2) = IceNO3(i,j,1)
          IceNH4(i,j,2) = IceNH4(i,j,1)
          IceLog(i,j,2) = IceLog(i,j,1)
# endif
#else
        Must define a case for ice initialization.
#endif
        ENDDO
      ENDDO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO i=1,2
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ui(:,:,i))
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            uie(:,:,i))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vi(:,:,i))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vie(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ai(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hi(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hsn(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ti(:,:,i))
#ifdef MELT_PONDS
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            apond(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hpond(:,:,i))
#endif
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ageice(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            sig11(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            sig22(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            sig12(:,:,i))
#ifdef ICE_BIO
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IcePhL(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IceNO3(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IceNH4(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IceLog(:,:,i))
#endif
        END DO
#ifdef NCEP_FLUXES
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wg2_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          cd_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ch_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ce_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wg2_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          cd_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ch_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ce_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoa_n)
#endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tis)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          s0mk)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          t0mk)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          utau_iw)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          chu_iw)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ui, uie, vi, vie)
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ai, hi, hsn, ti)
      CALL mp_exchange3d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sig11, sig12, sig22)
# ifdef MELT_PONDS
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    apond, hpond)
# endif
# ifdef ICE_BIO
      CALL mp_exchange3d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    IcePhL, IceNO3, IceNH4)
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    IceLog)
# endif
# ifdef NCEP_FLUXES
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    wg2_d, cd_d, ch_d, ce_d)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    wg2_m, cd_m, ch_m, ce_m)
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoa_n)
# endif
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tis, s0mk, t0mk, utau_iw)
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    chu_iw)
#endif

      RETURN
      END SUBROUTINE ana_ice_tile
