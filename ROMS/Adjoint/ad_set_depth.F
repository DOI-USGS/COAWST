#include "cppdefs.h"
      MODULE ad_set_depth_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_set_depth.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!  BASIC STATE variables needed: NONE                                  !
!  Independent Variables: ad_Hz, ad_z_r, ad_z_w                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_set_depth, ad_set_depth_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
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
      CALL ad_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % ad_h,                          &
# ifdef ICESHELF
     &                        GRID(ng) % zice,                          &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                        OCEAN(ng) % ad_bed,                       &
     &                        GRID(ng) % ad_bed_thick0,                 &
# endif
     &                        COUPLING(ng) % Zt_avg1,                   &
     &                        COUPLING(ng) % ad_Zt_avg1,                &
     &                        GRID(ng) % ad_Hz,                         &
     &                        GRID(ng) % ad_z_r,                        &
     &                        GRID(ng) % ad_z_w)
      RETURN
      END SUBROUTINE ad_set_depth

!
!***********************************************************************
      SUBROUTINE ad_set_depth_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              h, ad_h,                            &
# ifdef ICESHELF
     &                              zice,                               &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                              ad_bed, bed_thick0,                 &
# endif
     &                              Zt_avg1, ad_Zt_avg1,                &
     &                              ad_Hz, ad_z_r, ad_z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
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
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout):: ad_bed_thick0(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBi)
#  endif
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
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
      integer :: i, j, k, kbed

      real(r8) :: cff, cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: adfac, ad_hinv, ad_hwater, ad_z_r0, ad_z_w0

# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_BedThick
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_z_r0=0.0_r8
      ad_z_w0=0.0_r8
      ad_hinv=0.0_r8
      ad_hwater=0.0_r8
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_BedThick(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute time evolving depths and vertical thicknesses.
!-----------------------------------------------------------------------
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_z_r, tl_Hz)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_z_r, ad_Hz)
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_z_w)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_z_w)
!>    CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_h)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_h)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC

!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_Hz)
!>
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_Hz)
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_z_r)
!>
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_z_r)
!>    CALL exchange_w3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
!>   &                        tl_z_w)
!>
      CALL ad_exchange_w3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 0, N(ng),          &
     &                           ad_z_w)
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_h)
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_h)
#  endif
# endif
!
!  Compute vertical depths (meters, negative) at RHO- and W-points,
!  and vertical grid thicknesses (Hz=d(z_w)/d(s), meters).
!
      DO j=JstrR,JendR
        DO k=N(ng),1,-1
          cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
          cff1_w=SCALARS(ng)%Cs_w(k)
          cff2_w=SCALARS(ng)%sc_w(k)+1.0_r8
          cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
          cff1_r=SCALARS(ng)%Cs_r(k)
          cff2_r=SCALARS(ng)%sc_r(k)+1.0_r8
          DO i=IstrR,IendR
            hwater=h(i,j)
# ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            hwater=hwater-BedThick(i,j)
# endif
            hinv=1.0_r8/hwater
            z_w0=cff_w+cff1_w*hwater
            z_r0=cff_r+cff1_r*hwater

!>          tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
!>
            ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+ad_Hz(i,j,k)
            ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-ad_Hz(i,j,k)
            ad_Hz(i,j,k)=0.0_r8
!>          tl_z_r(i,j,k)=tl_z_r0+                                      &
!>   &                    tl_Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)+           &
!>   &                    Zt_avg1(i,j)*(tl_z_r0*hinv+z_r0*tl_hinv)
!>
            adfac=Zt_avg1(i,j)*ad_z_r(i,j,k)
            ad_z_r0=ad_z_r0+hinv*adfac+ad_z_r(i,j,k)
            ad_hinv=ad_hinv+z_r0*adfac
            ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                            &
     &                      (1.0_r8+z_r0*hinv)*ad_z_r(i,j,k)
            ad_z_r(i,j,k)=0.0_r8
!>          tl_z_r0=cff1_r*tl_hwater
!>
            ad_hwater=ad_hwater+cff1_r*ad_z_r0
            ad_z_r0=0.0_r8
!>          tl_z_w(i,j,k)=tl_z_w0+                                      &
!>   &                    tl_Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)+           &
!>   &                    Zt_avg1(i,j)*(tl_z_w0*hinv+z_w0*tl_hinv)
!>
            adfac=Zt_avg1(i,j)*ad_z_w(i,j,k)
            ad_z_w0=ad_z_w0+hinv*adfac+ad_z_w(i,j,k)
            ad_hinv=ad_hinv+z_w0*adfac
            ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                            &
     &                      (1.0_r8+z_w0*hinv)*ad_z_w(i,j,k)
            ad_z_w(i,j,k)=0.0_r8
!>          tl_z_w0=cff1_w*tl_hwater
!>
            ad_hwater=ad_hwater+cff1_w*ad_z_w0
            ad_z_w0=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>          tl_hwater=tl_hwater-tl_BedThick(i,j)
!>
            ad_BedThick(i,j)=ad_BedThick(i,j)-ad_hwater
# endif
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END DO

        DO i=IstrR,IendR
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>        tl_z_w(i,j,0)=tl_z_w(i,j,0)+tl_BedThick(i,j)
!>
          ad_BedThick(i,j)=ad_BedThick(i,j)+ad_z_w(i,j,0)
!>        tl_BedThick(i,j)=tl_BedThick(i,j)-tl_bed_thick0(i,j)
!>
          ad_bed_thick0(i,j)=ad_bed_thick0(i,j)-ad_BedThick(i,j)
          DO kbed=1,Nbed
!>          tl_BedThick(i,j)=tl_BedThick(i,j)+tl_bed(i,j,kbed,ithck)
!>
            ad_bed(i,j,kbed,ithck)=ad_bed(i,j,kbed,ithck)+
     &                             ad_BedThick(i,j)
          END DO
          ad_BedThick(i,j)=0.0_r8
# endif
!>        tl_z_w(i,j,0)=-tl_h(i,j)
!>
          ad_h(i,j)=ad_h(i,j)-ad_z_w(i,j,0)
          ad_z_w(i,j,0)=0.0_r8
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_set_depth_tile
#endif
      END MODULE ad_set_depth_mod
