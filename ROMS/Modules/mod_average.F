#include "cppdefs.h"
      MODULE mod_average
#ifdef AVERAGES
!
!svn $Id: mod_average.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  2D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgu2d     2D velocity component (m/s) in the XI-direction.         !
!  avgv2d     2D velocity component (m/s) in the ETA-direction.        !
!  avgzeta    Free surface (m).                                        !
!                                                                      !
!  3D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgAKs     Vertical diffusion of Salinity (m2/s).                   !
!  avgAKt     Vertical diffusion of temperature (m2/s).                !
!  avgAKv     Vertical viscosity (m2/s).                               !
!  avglhf     Latent heat flux (W/m2).                                 !
!  avglrf     Longwave radiation flux (W/m2).                          !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                         !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).             !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                         !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).             !
!  avgbus     Bottom u-momentum stress (N/m2).                         !
!  avgbvs     Bottom v-momentum stress (N/m2).                         !
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgshf     Sensible heat flux (W/m2).                               !
!  avgsrf     Shortwave radiation flux (W/m2).                         !
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net salt flux (kg/m2/s).                         !
!  avgevap    Surface net evaporation (kg/m2/s).                       !
!  avgrain    Surface net rain fall (kg/m2/s).                         !
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  avgUT      Quadratic term <u*t> for potential temperature and       !
!               salinity at U-points.                                  !
!  avgVT      Quadratic term <v*t> for potential temperature and       !
!               salinity at V-points.                                  !
!  avgTT      Quadratic term <t*t> for tracers.                        !
!  avgUU      Quadratic term <u*u> for 3D momentum at U-points.        !
!  avgUV      Quadratic term <u*v> for 3D momentum at RHO-points.      !
!  avgVV      Quadratic term <v*v> for 3D momentum at V-points.        !
!  avgU2      Quadratic term <ubar*ubar> for 2D momentum at U-points.  !
!  avgV2      Quadratic term <vbar*vbar> for 2D momentum at V-points.  !
!  avgZZ      Quadratic term <zeta*zeta> for free-surface.             !
# if defined AVERAGES_UV_KIRBY
!  avguWave   Kirby and Chen velocity in the XI-direction.             !
!  avgvWave   Kirby and Chen velocity in the ETA- direction.           !
# endif
!  avgu3d     3D velocity component (m/s) in the XI-direction.         !
!  avgv3d     3D velocity component (m/s) in the ETA-direction.        !
!  avgw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!  avgwvel    3D "true" vertical velocity (m/s).                       !
!                                                                      !
# if defined FORWARD_WRITE && defined SOLVE3D
!  avgDU_avg1 time-averaged u-flux for 3D momentum coupling.           !
!  avgDU_avg2 time-averaged u-flux for 3D momentum coupling.           !
!  avgDV_avg1 time-averaged v-flux for 3D momentum coupling.           !
!  avgDV_avg2 time-averaged v-flux for 3D momentum coupling.           !
!                                                                      !
# endif
# ifdef AVERAGES_NEARSHORE
!  Time-averaged radiation stresses.                                   !
!                                                                      !
!  avgu2Sd    2D stokes velocity component (m/s) in the XI-direction.  !
!  avgv2Sd    2D stokes velocity component (m/s) in the ETA-direction. !
!  avgu2RS    2D radiation stress tensor in the XI-direction.          !
!  avgv2RS    2D radiation stress tensor in the ETA-direction.         !
!  avgSxx2d   2D radiation stress, xx-component.                       !
!  avgSxy2d   2D radiation stress, xy-component.                       !
!  avgSyy2d   2D radiation stress, yy-component.                       !
!                                                                      !
!  avgu3Sd    3D stokes velocity component (m/s) in the XI-direction.  !
!  avgv3Sd    3D stokes velocity component (m/s) in the ETA-direction. !
!  avgu3RS    3D radiation stress tensor in the XI-direction.          !
!  avgv3RS    3D radiation stress tensor in the ETA-direction.         !
!  avgSxx3d   3D radiation stress, xx-component.                       !
!  avgSxy3d   3D radiation stress, xy-component.                       !
!  avgSyy3d   3D radiation stress, yy-component.                       !
!  avgSzx3d   3D radiation stress, zx-component.                       !
!  avgSzy3d   3D radiation stress, zy-component.                       !
!                                                                      !
# endif
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_AVERAGE

          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)
# ifdef AVERAGES_UV_KIRBY
          real(r8), pointer :: avguWave(:,:)
          real(r8), pointer :: avgvWave(:,:)
# endif
# ifdef AVERAGES_NEARSHORE
          real(r8), pointer :: avgu2Sd(:,:)
          real(r8), pointer :: avgv2Sd(:,:)
          real(r8), pointer :: avgu2RS(:,:)
          real(r8), pointer :: avgv2RS(:,:)
          real(r8), pointer :: avgSxx2d(:,:)
          real(r8), pointer :: avgSxy2d(:,:)
          real(r8), pointer :: avgSyy2d(:,:)
# endif
          real(r8), pointer :: avgzeta(:,:)
# ifdef AVERAGES_QUADRATIC
          real(r8), pointer :: avgU2(:,:)
          real(r8), pointer :: avgV2(:,:)
          real(r8), pointer :: avgZZ(:,:)
# endif
# ifdef SOLVE3D
          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)
          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)
          real(r8), pointer :: avgw3d(:,:,:)
          real(r8), pointer :: avgwvel(:,:,:)
#  ifdef FORWARD_WRITE
          real(r8), pointer :: avgDU_avg1(:,:)
          real(r8), pointer :: avgDU_avg2(:,:)
          real(r8), pointer :: avgDV_avg1(:,:)
          real(r8), pointer :: avgDV_avg2(:,:)
#  endif
#  ifdef AVERAGES_NEARSHORE
          real(r8), pointer :: avgu3Sd(:,:,:)
          real(r8), pointer :: avgv3Sd(:,:,:)
          real(r8), pointer :: avgu3RS(:,:,:)
          real(r8), pointer :: avgv3RS(:,:,:)
          real(r8), pointer :: avgSxx3d(:,:,:)
          real(r8), pointer :: avgSxy3d(:,:,:)
          real(r8), pointer :: avgSyy3d(:,:,:)
          real(r8), pointer :: avgSzx3d(:,:,:)
          real(r8), pointer :: avgSzy3d(:,:,:)
#  endif
#  ifdef AVERAGES_QUADRATIC
#   ifdef SALINITY
          real(r8), pointer :: avgSS(:,:,:)
#   endif
          real(r8), pointer :: avgTT(:,:,:,:)
          real(r8), pointer :: avgHuonT(:,:,:,:)
          real(r8), pointer :: avgHvomT(:,:,:,:)
          real(r8), pointer :: avgUT(:,:,:,:)
          real(r8), pointer :: avgVT(:,:,:,:)
          real(r8), pointer :: avgHuon(:,:,:)
          real(r8), pointer :: avgHvom(:,:,:)
          real(r8), pointer :: avgUU(:,:,:)
          real(r8), pointer :: avgUV(:,:,:)
          real(r8), pointer :: avgVV(:,:,:)
#  endif
#  ifdef AVERAGES_AKS
          real(r8), pointer :: avgAKs(:,:,:)
#  endif
#  ifdef AVERAGES_AKT
          real(r8), pointer :: avgAKt(:,:,:)
#  endif
#  ifdef AVERAGES_AKV
          real(r8), pointer :: avgAKv(:,:,:)
#  endif
#  ifdef AVERAGES_FLUXES
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
#   ifdef BULK_FLUXES
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
#    ifdef EMINUSP
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
#    endif
#   endif
#   ifdef SHORTWAVE
          real(r8), pointer :: avgsrf(:,:)
#   endif
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
# endif
# ifdef AVERAGES_FLUXES
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)
          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)
# endif
        END TYPE T_AVERAGE

        TYPE (T_AVERAGE), allocatable :: AVERAGE(:)

      CONTAINS

      SUBROUTINE allocate_average (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1 ) allocate ( AVERAGE(Ngrids) )
!
      allocate ( AVERAGE(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgzeta(LBi:UBi,LBj:UBj) )
# ifdef AVERAGES_UV_KIRBY
      allocate ( AVERAGE(ng) % avguWave(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgvWave(LBi:UBi,LBj:UBj) )
# endif
# ifdef AVERAGES_NEARSHORE
      allocate ( AVERAGE(ng) % avgu2Sd(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgv2Sd(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgu2RS(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgv2RS(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgSxx2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgSxy2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgSyy2d(LBi:UBi,LBj:UBj) )
# endif
# ifdef AVERAGES_QUADRATIC
      allocate ( AVERAGE(ng) % avgU2(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgV2(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgZZ(LBi:UBi,LBj:UBj) )
# endif
# ifdef SOLVE3D
      allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
#  ifdef FORWARD_WRITE
      allocate ( AVERAGE(ng) % avgDU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgDU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgDV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgDV_avg2(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef AVERAGES_NEARSHORE
      allocate ( AVERAGE(ng) % avgu3Sd(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgv3Sd(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgu3RS(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgv3RS(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgSxx3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgSxy3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgSyy3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgSzx3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgSzy3d(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
      allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( AVERAGE(ng) % avgwvel(LBi:UBi,LBj:UBj,0:N(ng)) )
#  ifdef AVERAGES_QUADRATIC
      allocate ( AVERAGE(ng) % avgTT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgHuonT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgHvomT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgUT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgVT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgHuon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgHvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgUU(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgUV(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgVV(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
#  ifdef AVERAGES_AKS
      allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_AKT
      allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_AKV
      allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_FLUXES
      allocate ( AVERAGE(ng) % avgstf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgswf(LBi:UBi,LBj:UBj) )
#   ifdef BULK_FLUXES
      allocate ( AVERAGE(ng) % avglhf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avglrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgshf(LBi:UBi,LBj:UBj) )
#    ifdef EMINUSP
      allocate ( AVERAGE(ng) % avgevap(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgrain(LBi:UBi,LBj:UBj) )
#    endif
#   endif
#   ifdef SHORTWAVE
      allocate ( AVERAGE(ng) % avgsrf(LBi:UBi,LBj:UBj) )
#   endif
#  endif

#  ifdef LMD_BKPP
      allocate ( AVERAGE(ng) % avghbbl(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef LMD_SKPP
      allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
#  endif
# endif
# ifdef AVERAGES_FLUXES
      allocate ( AVERAGE(ng) % avgsus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgbus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgbvs(LBi:UBi,LBj:UBj) )
# endif

      RETURN
      END SUBROUTINE allocate_average

      SUBROUTINE initialize_average (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
# else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          AVERAGE(ng) % avgu2d(i,j) = IniVal
          AVERAGE(ng) % avgv2d(i,j) = IniVal
          AVERAGE(ng) % avgzeta(i,j) = IniVal
# ifdef AVERAGES_UV_KIRBY
          AVERAGE(ng) % avguWave(i,j) = IniVal
          AVERAGE(ng) % avgvWave(i,j) = IniVal
# endif
# ifdef AVERAGES_NEARSHORE
          AVERAGE(ng) % avgu2Sd(i,j) = IniVal
          AVERAGE(ng) % avgv2Sd(i,j) = IniVal
          AVERAGE(ng) % avgu2RS(i,j) = IniVal
          AVERAGE(ng) % avgv2RS(i,j) = IniVal
          AVERAGE(ng) % avgSxx2d(i,j) = IniVal
          AVERAGE(ng) % avgSxy2d(i,j) = IniVal
          AVERAGE(ng) % avgSyy2d(i,j) = IniVal
# endif
# ifdef AVERAGES_QUADRATIC
          AVERAGE(ng) % avgU2(i,j) = IniVal
          AVERAGE(ng) % avgV2(i,j) = IniVal
          AVERAGE(ng) % avgZZ(i,j) = IniVal
# endif
# ifdef AVERAGES_FLUXES
          AVERAGE(ng) % avgsus(i,j) = IniVal
          AVERAGE(ng) % avgsvs(i,j) = IniVal
          AVERAGE(ng) % avgbus(i,j) = IniVal
          AVERAGE(ng) % avgbvs(i,j) = IniVal
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
          AVERAGE(ng) % avgDU_avg1(i,j) = IniVal
          AVERAGE(ng) % avgDU_avg2(i,j) = IniVal
          AVERAGE(ng) % avgDV_avg1(i,j) = IniVal
          AVERAGE(ng) % avgDV_avg2(i,j) = IniVal
#  endif
#  ifdef AVERAGES_FLUXES
          AVERAGE(ng) % avgstf(i,j) = IniVal
          AVERAGE(ng) % avgswf(i,j) = IniVal
#   ifdef BULK_FLUXES
          AVERAGE(ng) % avglhf(i,j) = IniVal
          AVERAGE(ng) % avglrf(i,j) = IniVal
          AVERAGE(ng) % avgshf(i,j) = IniVal
#    ifdef EMINUSP
          AVERAGE(ng) % avgevap(i,j) = IniVal
          AVERAGE(ng) % avgrain(i,j) = IniVal
#    endif
#   endif
#   ifdef SHORTWAVE
          AVERAGE(ng) % avgsrf(i,j) = IniVal
#   endif
#  endif
#  ifdef LMD_BKPP
          AVERAGE(ng) % avghbbl(i,j) = IniVal
#  endif
#  ifdef LMD_SKPP
          AVERAGE(ng) % avghsbl(i,j) = IniVal
#  endif
# endif
        END DO
# ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgrho(i,j,k) = IniVal
            AVERAGE(ng) % avgu3d(i,j,k) = IniVal
            AVERAGE(ng) % avgv3d(i,j,k) = IniVal
#  ifdef AVERAGES_NEARSHORE
            AVERAGE(ng) % avgu3Sd(i,j,k) = IniVal
            AVERAGE(ng) % avgv3Sd(i,j,k) = IniVal
            AVERAGE(ng) % avgu3RS(i,j,k) = IniVal
            AVERAGE(ng) % avgv3RS(i,j,k) = IniVal
            AVERAGE(ng) % avgSxx3d(i,j,k) = IniVal
            AVERAGE(ng) % avgSxy3d(i,j,k) = IniVal
            AVERAGE(ng) % avgSyy3d(i,j,k) = IniVal
            AVERAGE(ng) % avgSzx3d(i,j,k) = IniVal
            AVERAGE(ng) % avgSzy3d(i,j,k) = IniVal
#  endif
#  ifdef AVERAGES_QUADRATIC
            AVERAGE(ng) % avgHuon(i,j,k) = IniVal
            AVERAGE(ng) % avgHvom(i,j,k) = IniVal
            AVERAGE(ng) % avgUU(i,j,k) = IniVal
            AVERAGE(ng) % avgUV(i,j,k) = IniVal
            AVERAGE(ng) % avgVV(i,j,k) = IniVal
#  endif
          END DO
        END DO
        DO k=0,N(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgw3d(i,j,k) = IniVal
            AVERAGE(ng) % avgwvel(i,j,k) = IniVal
#  ifdef AVERAGES_AKS
            AVERAGE(ng) % avgAKs(i,j,k) = IniVal
#  endif
#  ifdef AVERAGES_AKT
            AVERAGE(ng) % avgAKt(i,j,k) = IniVal
#  endif
#  ifdef AVERAGES_AKV
            AVERAGE(ng) % avgAKv(i,j,k) = IniVal
#  endif
          END DO
        END DO
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgt(i,j,k,itrc) = IniVal
#  ifdef AVERAGES_QUADRATIC
              AVERAGE(ng) % avgTT(i,j,k,itrc) = IniVal
              AVERAGE(ng) % avgHuonT(i,j,k,itrc) = IniVal
              AVERAGE(ng) % avgHvomT(i,j,k,itrc) = IniVal
              AVERAGE(ng) % avgUT(i,j,k,itrc) = IniVal
              AVERAGE(ng) % avgVT(i,j,k,itrc) = IniVal
#  endif
            END DO
          END DO
        END DO
# endif
      END DO

      RETURN
      END SUBROUTINE initialize_average
#endif
      END MODULE mod_average
