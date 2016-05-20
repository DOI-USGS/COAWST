#include "cppdefs.h"
      MODULE mod_average
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
!
!svn $Id: mod_average.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  The strategy here is  to define all possible pointers in the        !
!  time-averaged structure and allocate only those requested by        !
!  the user. This will facilitate a better management of memory.       !
!                                                                      !
!  Time-averaged state variables for output purposes.                  !
!                                                                      !
!  avgu2d     2D velocity  component (m/s) in the XI-direction.        !
!  avgv2d     2D velocity  component (m/s) in the ETA-direction.       !
!  avgu2dE    2D Eastward  component (m/s) at RHO-points.              !
!  avgv2dN    2D Northward component (m/s) at RHO-points.              !
!  avgzeta    Free surface (m).                                        !

# ifdef SOLVE3D
!  avgUwind   2D wind velocity component (m/s) in the XI-direction.    !
!  avgVwind   2D wind velocity component (m/s) in the ETA-direction.   !
!  avgUwindE  2D wind velocity component (m/s) to the east.            !
!  avgVwindN  2D wind velocity component (m/s) to the north.           !
!  avgu3d     3D velocity  component (m/s) in the XI-direction.        !
!  avgv3d     3D velocity  component (m/s) in the ETA-direction.       !
!  avgu3dE    3D Eastward  component (m/s) at RHO-points.              !
!  avgv3dN    3D Northward component (m/s) at RHO-points.              !
!  avgw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!  avgwvel    3D "true" vertical velocity (m/s).                       !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!  avgAKt     Vertical diffusion of temperature (m2/s).                !
!  avgAKv     Vertical viscosity (m2/s).                               !
!  avgAKs     Vertical diffusion of Salinity (m2/s).                   !
#  endif
#  ifdef LMD_SKPP
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
#  endif
#  ifdef LMD_BKPP
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
#  endif
#  if defined UV_KIRBY
!  avguWave   Kirby and Chen velocity in the XI-direction.             !
!  avgvWave   Kirby and Chen velocity in the ETA- direction.           !
#  endif
# endif
# if defined FORWARD_WRITE && defined SOLVE3D
!                                                                      !
!  Time-averaged 2D/3D coupling terms.                                 !
!                                                                      !
!  avgDU_avg1 time-averaged u-flux for 3D momentum coupling.           !
!  avgDU_avg2 time-averaged u-flux for 3D momentum coupling.           !
!  avgDV_avg1 time-averaged v-flux for 3D momentum coupling.           !
!  avgDV_avg2 time-averaged v-flux for 3D momentum coupling.           !
# endif
!                                                                      !
!  Time-averaged surface and bottom fluxes.                            !
!                                                                      !
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!  avgbus     Bottom u-momentum stress (N/m2).                         !
!  avgbvs     Bottom v-momentum stress (N/m2).                         !

# ifdef SOLVE3D
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net freshwater flux (kg/m2/s).                   !
#  ifdef SHORTWAVE
!  avgsrf     Shortwave radiation flux (W/m2).                         !
#  endif
#  ifdef BULK_FLUXES
!  avglhf     Latent heat flux (W/m2).                                 !
!  avglrf     Longwave radiation flux (W/m2).                          !
!  avgshf     Sensible heat flux (W/m2).                               !
#   ifdef EMINUSP
!  avgevap    Surface net evaporation (kg/m2/s).                       !
!  avgrain    Surface net rain fall (kg/m2/s).                         !
#   endif
#  endif
# endif
# ifdef WEC
!  Time-averaged wave efect on currents.                               !
!                                                                      !
!  avgu2Sd    2D stokes velocity component (m/s) in the XI-direction.  !
!  avgv2Sd    2D stokes velocity component (m/s) in the ETA-direction. !
!  avgu2rs    2D radiation stress tensor in the XI-direction.          !
!  avgv2rs    2D radiation stress tensor in the ETA-direction.         !
#  ifdef WEC_MELLOR
!  avgSxx2d   2D radiation stress, xx-component.                       !
!  avgSxy2d   2D radiation stress, xy-component.                       !
!  avgSyy2d   2D radiation stress, yy-component.                       !
#  endif
#  ifdef SOLVE3D
!  avgu3Sd    3D stokes velocity component (m/s) in the XI-direction.  !
!  avgv3Sd    3D stokes velocity component (m/s) in the ETA-direction. !
!  avgu3rs    3D radiation stress tensor in the XI-direction.          !
!  avgv3rs    3D radiation stress tensor in the ETA-direction.         !
#   ifdef WEC_MELLOR
!  avgSxx3d   3D radiation stress, xx-component.                       !
!  avgSxy3d   3D radiation stress, xy-component.                       !
!  avgSyy3d   3D radiation stress, yy-component.                       !
!  avgSzx3d   3D radiation stress, zx-component.                       !
!  avgSzy3d   3D radiation stress, zy-component.                       !
#   endif
#  endif
# endif
!                                                                      !
!  Time-averaged quadratic fields.                                     !
!                                                                      !
!  avgZZ      Quadratic term <zeta*zeta> for free-surface.             !
!  avgU2      Quadratic term <ubar*ubar> for 2D momentum at U-points.  !
!  avgV2      Quadratic term <vbar*vbar> for 2D momentum at V-points.  !

# ifdef SOLVE3D
!  avgUU      Quadratic term <u*u> for 3D momentum at U-points.        !
!  avgVV      Quadratic term <v*v> for 3D momentum at V-points.        !
!  avgUV      Quadratic term <u*v> for 3D momentum at RHO-points.      !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                         !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                         !
!  avgTT      Quadratic term <t*t> for tracers.                        !
!  avgUT      Quadratic term <u*t> for potential temperature and       !
!               salinity at U-points.                                  !
!  avgVT      Quadratic term <v*t> for potential temperature and       !
!               salinity at V-points.                                  !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).             !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).             !
# endif
!                                                                      !
!  Time-averages vorticity fields.                                     !
!                                                                      !
!  avgpvor2d  2D, vertically integrated, potential vorticity.          !
!  avgrvor2d  2D, vertically integrated, relative vorticity.           !
!  avgpvor3d  3D potential vorticity.                                  !
!  rvorvor2d  3D relative vorticity.                                   !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_AVERAGE
!
!  Time-averaged state variables.
!
          real(r8), pointer :: avgzeta(:,:)

          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)

          real(r8), pointer :: avgu2dE(:,:)
          real(r8), pointer :: avgv2dN(:,:)

# ifdef SOLVE3D

#  ifdef BEST_NPZ
#   ifdef STATIONARY
          real(r8), pointer :: avgst(:,:,:,:)
#   endif
#   ifdef STATIONARY2
          real(r8), pointer :: avgst2(:,:,:)
#   endif
#   ifdef PROD3
          real(r8), pointer :: avgpt3(:,:,:,:)
#   endif
#   ifdef PROD2
          real(r8), pointer :: avgpt2(:,:,:)
#   endif
#   ifdef BENTHIC
          real(r8), pointer :: avgbt(:,:,:,:)
#   endif
#   ifdef ICE_BIO
#     ifdef CLIM_ICE_1D   
          real(r8), pointer :: avgit(:,:,:)
#     elif defined BERING_10K
          real(r8), pointer :: avgIPhL(:,:)
          real(r8),pointer:: avgINO3(:,:)
          real(r8), pointer :: avgINH4(:,:)
          real(r8), pointer :: avgILog(:,:)
#     endif
#   endif
#   if defined BIOFLUX
          real(r8), pointer :: avgbflx(:,:)
#   endif
#  endif

          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)

          real(r8), pointer :: avgu3dE(:,:,:)
          real(r8), pointer :: avgv3dN(:,:,:)

          real(r8), pointer :: avgw3d(:,:,:)
          real(r8), pointer :: avgwvel(:,:,:)

          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)

#  ifdef NEMURO_SED1
          real(r8), pointer :: avgPONsed(:,:)
          real(r8), pointer :: avgOPALsed(:,:)
          real(r8), pointer :: avgDENITsed(:,:)
          real(r8), pointer :: avgPONbur(:,:)
          real(r8), pointer :: avgOPALbur(:,:)
#  endif
#  ifdef PRIMARY_PROD
          real(r8), pointer :: avgNPP(:,:)
#  endif
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          real(r8), pointer :: avgAKv(:,:,:)
          real(r8), pointer :: avgAKt(:,:,:)
          real(r8), pointer :: avgAKs(:,:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
          real(r8), pointer :: avgDU_avg1(:,:)
          real(r8), pointer :: avgDU_avg2(:,:)

          real(r8), pointer :: avgDV_avg1(:,:)
          real(r8), pointer :: avgDV_avg2(:,:)
# endif
!
!  Time-averaged surface and bottom fluxes.
!
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)

          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)

# ifdef BBL
          real(r8), pointer :: avgUbrs(:,:)
          real(r8), pointer :: avgVbrs(:,:)
          real(r8), pointer :: avgUbws(:,:)
          real(r8), pointer :: avgVbws(:,:)
          real(r8), pointer :: avgUbcs(:,:)
          real(r8), pointer :: avgVbcs(:,:)
          real(r8), pointer :: avgUVwc(:,:)

          real(r8), pointer :: avgUbot(:,:)
          real(r8), pointer :: avgVbot(:,:)
          real(r8), pointer :: avgUbur(:,:)
          real(r8), pointer :: avgVbvr(:,:)
# endif

# ifdef SOLVE3D
#  if defined RUNOFF
          real(r8), pointer :: avgRunoff(:,:)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          real(r8), pointer :: avgPair(:,:)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          real(r8), pointer :: avgUwind(:,:)
          real(r8), pointer :: avgVwind(:,:)
          real(r8), pointer :: avgUwindE(:,:)
          real(r8), pointer :: avgVwindN(:,:)
#  endif
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
#  ifdef SHORTWAVE
          real(r8), pointer :: avgsrf(:,:)
#  endif
#  ifdef BULK_FLUXES
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
          real(r8), pointer :: avgsssflx(:,:)
#   ifdef EMINUSP
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
#   endif
#  endif
# endif
# ifdef WEC
          real(r8), pointer :: avgu2Sd(:,:)
          real(r8), pointer :: avgv2Sd(:,:)
          real(r8), pointer :: avgu2rs(:,:)
          real(r8), pointer :: avgv2rs(:,:)

          real(r8), pointer :: avgWztw(:,:)
          real(r8), pointer :: avgWqsp(:,:)
          real(r8), pointer :: avgWbeh(:,:)
#  ifdef WEC_MELLOR
          real(r8), pointer :: avgSxx2d(:,:)
          real(r8), pointer :: avgSxy2d(:,:)
          real(r8), pointer :: avgSyy2d(:,:)
#  endif
#  ifdef SOLVE3D
          real(r8), pointer :: avgu3Sd(:,:,:)
          real(r8), pointer :: avgv3Sd(:,:,:)
          real(r8), pointer :: avgw3Sd(:,:,:)
          real(r8), pointer :: avgw3St(:,:,:)
          real(r8), pointer :: avgu3rs(:,:,:)
          real(r8), pointer :: avgv3rs(:,:,:)
#   ifdef WEC_MELLOR
          real(r8), pointer :: avgSxx3d(:,:,:)
          real(r8), pointer :: avgSxy3d(:,:,:)
          real(r8), pointer :: avgSyy3d(:,:,:)
          real(r8), pointer :: avgSzx3d(:,:,:)
          real(r8), pointer :: avgSzy3d(:,:,:)
#   endif
#  endif
# endif
# ifdef WAVES_HEIGHT
          real(r8), pointer :: avgWamp(:,:)
# endif
# ifdef WAVES_LENGTH
          real(r8), pointer :: avgWlen(:,:)
# endif
# ifdef WAVES_LENGTHP
          real(r8), pointer :: avgWlep(:,:)
# endif
# ifdef WAVES_DIR
          real(r8), pointer :: avgWdir(:,:)
# endif
# ifdef WAVES_TOP_PERIOD
          real(r8), pointer :: avgWptp(:,:)
# endif
# ifdef WAVES_BOT_PERIOD
          real(r8), pointer :: avgWpbt(:,:)
# endif
# ifdef BBL_MODEL
          real(r8), pointer :: avgWorb(:,:)
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          real(r8), pointer :: avgWdif(:,:)
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          real(r8), pointer :: avgWdib(:,:)
          real(r8), pointer :: avgWdiw(:,:)
# endif
# ifdef ROLLER_SVENDSEN
          real(r8), pointer :: avgWbrk(:,:)
# endif
# ifdef WEC_ROLLER
          real(r8), pointer :: avgWdis(:,:)
          real(r8), pointer :: avgWrol(:,:)
# endif
# ifdef UV_KIRBY
          real(r8), pointer :: avguWav(:,:)
          real(r8), pointer :: avgvWav(:,:)
# endif
!
!  Time-averaged quadratic fields.
!
          real(r8), pointer :: avgZZ(:,:)
          real(r8), pointer :: avgU2(:,:)
          real(r8), pointer :: avgV2(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: avgUU(:,:,:)
          real(r8), pointer :: avgUV(:,:,:)
          real(r8), pointer :: avgVV(:,:,:)

          real(r8), pointer :: avgHuon(:,:,:)
          real(r8), pointer :: avgHvom(:,:,:)

          real(r8), pointer :: avgTT(:,:,:,:)
          real(r8), pointer :: avgUT(:,:,:,:)
          real(r8), pointer :: avgVT(:,:,:,:)

          real(r8), pointer :: avgHuonT(:,:,:,:)
          real(r8), pointer :: avgHvomT(:,:,:,:)
#  ifdef ICE_MODEL
          real(r8), pointer :: avguice(:,:)
          real(r8), pointer :: avgvice(:,:)
          real(r8), pointer :: avguiceE(:,:)
          real(r8), pointer :: avgviceN(:,:)
          real(r8), pointer :: avgaice(:,:)
          real(r8), pointer :: avghice(:,:)
          real(r8), pointer :: avgtice(:,:)
          real(r8), pointer :: avgtimid(:,:)
          real(r8), pointer :: avghsno(:,:)
#   ifdef MELT_PONDS
          real(r8), pointer :: avgapond(:,:)
          real(r8), pointer :: avghpond(:,:)
#   endif
          real(r8), pointer :: avgiomflx(:,:)
          real(r8), pointer :: avgageice(:,:)
          real(r8), pointer :: avgsig11(:,:)
          real(r8), pointer :: avgsig12(:,:)
          real(r8), pointer :: avgsig22(:,:)
          real(r8), pointer :: avgT0mk(:,:)
          real(r8), pointer :: avgS0mk(:,:)
          real(r8), pointer :: avgWfr(:,:)
          real(r8), pointer :: avgWai(:,:)
          real(r8), pointer :: avgWao(:,:)
          real(r8), pointer :: avgWio(:,:)
          real(r8), pointer :: avgWro(:,:)
          real(r8), pointer :: avgWdiv(:,:)
          real(r8), pointer :: avgchu_iw(:,:)
          real(r8), pointer :: avgutau_iw(:,:)
#  endif
# endif
!
!  Time-averaged vorticity fields.
!
          real(r8), pointer :: avgpvor2d(:,:)
          real(r8), pointer :: avgrvor2d(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: avgpvor3d(:,:,:)
          real(r8), pointer :: avgrvor3d(:,:,:)
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
      USE mod_ncparam
      USE mod_scalars
# if defined BEST_NPZ
      USE mod_biology
# endif
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
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        allocate ( AVERAGE(ng) % avgzeta(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idUbar,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbar,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idu2dE,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2dE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idv2dN,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2dN(LBi:UBi,LBj:UBj) )
      END IF

# ifdef SOLVE3D

#  ifdef BEST_NPZ
#   ifdef STATIONARY
      allocate ( AVERAGE(ng) % avgst(LBi:UBi,LBj:UBj,N(ng),NTS(ng)) )
#   endif
#   ifdef STATIONARY2
      allocate ( AVERAGE(ng) % avgst2(LBi:UBi,LBj:UBj,NTS2(ng)) )
#   endif
#   ifdef PROD3
      allocate ( AVERAGE(ng) % avgpt3(LBi:UBi,LBj:UBj,N(ng),NPT3(ng)) )
#   endif
#   ifdef PROD2
      allocate ( AVERAGE(ng) % avgpt2(LBi:UBi,LBj:UBj,NPT2(ng)) )
#   endif
#   ifdef BENTHIC
       allocate ( AVERAGE(ng) % avgbt(LBi:UBi,LBj:UBj,NBL(ng),NBeT(ng)) )
#   endif
#   ifdef ICE_BIO
#    ifdef CLIM_ICE_1D
        allocate ( AVERAGE(ng) % avgit(LBi:UBi,LBj:UBj,NIceT(ng)) )
#    elif defined BERING_10K
      IF (Aout(idIceNO3,ng)) THEN
        allocate ( AVERAGE(ng) % avgINO3(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idIcePhL,ng)) THEN
        allocate ( AVERAGE(ng) % avgIPhL(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idIceNH4,ng)) THEN
        allocate ( AVERAGE(ng) % avgINH4(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idIceLog,ng)) THEN
        allocate ( AVERAGE(ng) % avgILog(LBi:UBi,LBj:UBj) )
      END IF
#    endif
#   endif

#   ifdef BIOFLUX
        allocate ( AVERAGE(ng) % avgbflx(NAT+NBT,NAT+NBT) )
#   endif
#  endif

#  ifdef BIO_GOANPZ
        allocate ( AVERAGE(ng) % avgst(LBi:UBi,LBj:UBj,N(ng),NTS(ng)) )
#  endif
      IF (Aout(idUvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idVvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF

      IF (Aout(idu3dE,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3dE(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idv3dN,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3dN(LBi:UBi,LBj:UBj,N(ng)) )
      END IF

      IF (Aout(idOvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idWvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgwvel(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF

      IF (Aout(idDano,ng)) THEN
        allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF

#  ifdef NEMURO_SED1
      IF (Aout(idPONsed,ng)) THEN
        allocate ( AVERAGE(ng) % avgPONsed(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idOPALsed,ng)) THEN
        allocate ( AVERAGE(ng) % avgOPALsed(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idDENITsed,ng)) THEN
        allocate ( AVERAGE(ng) % avgDENITsed(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idPONbur,ng)) THEN
        allocate ( AVERAGE(ng) % avgPONbur(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idOPALbur,ng)) THEN
        allocate ( AVERAGE(ng) % avgOPALbur(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  ifdef PRIMARY_PROD
      IF (Aout(idNPP,ng)) THEN
        allocate ( AVERAGE(ng) % avgNPP(LBi:UBi,LBj:UBj) )
      END IF
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idTdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idSdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
#  endif
#  ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        allocate ( AVERAGE(ng) % avghbbl(LBi:UBi,LBj:UBj) )
      END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
      IF (Aout(idUfx1,ng)) THEN
        allocate ( AVERAGE(ng) % avgDU_avg1(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUfx2,ng)) THEN
        allocate ( AVERAGE(ng) % avgDU_avg2(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idVfx1,ng)) THEN
        allocate ( AVERAGE(ng) % avgDV_avg1(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVfx2,ng)) THEN
        allocate ( AVERAGE(ng) % avgDV_avg2(LBi:UBi,LBj:UBj) )
      END IF
# endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        allocate ( AVERAGE(ng) % avgsus(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVsms,ng)) THEN
        allocate ( AVERAGE(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idUbms,ng)) THEN
        allocate ( AVERAGE(ng) % avgbus(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbms,ng)) THEN
        allocate ( AVERAGE(ng) % avgbvs(LBi:UBi,LBj:UBj) )
      END IF

# ifdef BBL
      IF (Aout(idUbrs,ng)) THEN
        allocate ( AVERAGE(ng) % avgUbrs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbrs,ng)) THEN
        allocate ( AVERAGE(ng) % avgVbrs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUbws,ng)) THEN
        allocate ( AVERAGE(ng) % avgUbws(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbws,ng)) THEN
        allocate ( AVERAGE(ng) % avgVbws(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUbcs,ng)) THEN
        allocate ( AVERAGE(ng) % avgUbcs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbcs,ng)) THEN
        allocate ( AVERAGE(ng) % avgVbcs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUVwc,ng)) THEN
        allocate ( AVERAGE(ng) % avgUVwc(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idUbot,ng)) THEN
        allocate ( AVERAGE(ng) % avgUbot(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbot,ng)) THEN
        allocate ( AVERAGE(ng) % avgVbot(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idUbur,ng)) THEN
        allocate ( AVERAGE(ng) % avgUbur(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVbvr,ng)) THEN
        allocate ( AVERAGE(ng) % avgVbvr(LBi:UBi,LBj:UBj) )
      END IF
# endif

# ifdef SOLVE3D
#  if defined RUNOFF
      IF (Aout(idRunoff,ng)) THEN
        allocate ( AVERAGE(ng) % avgRunoff(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      IF (Aout(idPair,ng)) THEN
        allocate ( AVERAGE(ng) % avgPair(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
      IF (Aout(idUair,ng)) THEN
        allocate ( AVERAGE(ng) % avgUwind(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVair,ng)) THEN
        allocate ( AVERAGE(ng) % avgVwind(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idUairE,ng)) THEN
        allocate ( AVERAGE(ng) % avgUwindE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVairN,ng)) THEN
        allocate ( AVERAGE(ng) % avgVwindN(LBi:UBi,LBj:UBj) )
      END IF
#  endif
      IF (Aout(idTsur(itemp),ng)) THEN
        allocate ( AVERAGE(ng) % avgstf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTsur(isalt),ng)) THEN
        allocate ( AVERAGE(ng) % avgswf(LBi:UBi,LBj:UBj) )
      END IF
#  ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        allocate ( AVERAGE(ng) % avgsrf(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  ifdef BULK_FLUXES
      IF (Aout(idLhea,ng)) THEN
        allocate ( AVERAGE(ng) % avglhf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idLrad,ng)) THEN
        allocate ( AVERAGE(ng) % avglrf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idShea,ng)) THEN
        allocate ( AVERAGE(ng) % avgshf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idSSSf,ng)) THEN
        allocate ( AVERAGE(ng) % avgsssflx(LBi:UBi,LBj:UBj) )
      END IF
#   ifdef EMINUSP
      IF (Aout(idevap,ng)) THEN
        allocate ( AVERAGE(ng) % avgevap(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idrain,ng)) THEN
        allocate ( AVERAGE(ng) % avgrain(LBi:UBi,LBj:UBj) )
      END IF
#   endif
#  endif
# endif

# ifdef WEC
!
!  Time-averaged wave effects on currents.
!
      IF (Aout(idU2Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2Sd(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idV2Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2Sd(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idU2rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2rs(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idV2rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2rs(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idWztw,ng)) THEN
        allocate ( AVERAGE(ng) % avgWztw(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWqsp,ng)) THEN
        allocate ( AVERAGE(ng) % avgWqsp(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWbeh,ng)) THEN
        allocate ( AVERAGE(ng) % avgWbeh(LBi:UBi,LBj:UBj) )
      END IF

#  ifdef WEC_MELLOR
      IF (Aout(idW2xx,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxx2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idW2xy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxy2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idW2yy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSyy2d(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  ifdef SOLVE3D
      IF (Aout(idU3Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3Sd(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idV3Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3Sd(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idW3Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgW3Sd(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idW3St,ng)) THEN
        allocate ( AVERAGE(ng) % avgW3St(LBi:UBi,LBj:UBj,0:N(ng)) )
      END IF
      IF (Aout(idU3rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3rs(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idV3rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3rs(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
#   ifdef WEC_MELLOR
      IF (Aout(idW3xx,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxx3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idW3xy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxy3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idW3yy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSyy3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idW3zx,ng)) THEN
        allocate ( AVERAGE(ng) % avgSzx3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idW3zy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSzy3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
#   endif
#  endif
# endif
# ifdef WAVES_HEIGHT
      IF (Aout(idWamp,ng)) THEN
        allocate ( AVERAGE(ng) % avgWamp(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef WAVES_LENGTH
      IF (Aout(idWlen,ng)) THEN
        allocate ( AVERAGE(ng) % avgWlen(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef WAVES_LENGTHP
      IF (Aout(idWlep,ng)) THEN
        allocate ( AVERAGE(ng) % avgWlep(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef WAVES_DIR
      IF (Aout(idWdir,ng)) THEN
        allocate ( AVERAGE(ng) % avgWdir(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef WAVES_TOP_PERIOD
      IF (Aout(idWptp,ng)) THEN
        allocate ( AVERAGE(ng) % avgWptp(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef WAVES_BOT_PERIOD
      IF (Aout(idWpbt,ng)) THEN
        allocate ( AVERAGE(ng) % avgWpbt(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef BBL_MODEL
      IF (Aout(idWorb,ng)) THEN
        allocate ( AVERAGE(ng) % avgWorb(LBi:UBi,LBj:UBj) )
      END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      IF (Aout(idWdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgWdif(LBi:UBi,LBj:UBj) )
      END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
      IF (Aout(idWdib,ng)) THEN
        allocate ( AVERAGE(ng) % avgWdib(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWdiw,ng)) THEN
        allocate ( AVERAGE(ng) % avgWdiw(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef ROLLER_SVENDSEN
      IF (Aout(idWbrk,ng)) THEN
        allocate ( AVERAGE(ng) % avgWbrk(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef WEC_ROLLER
      IF (Aout(idWdis,ng)) THEN
        allocate ( AVERAGE(ng) % avgWdis(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWrol,ng)) THEN
        allocate ( AVERAGE(ng) % avgWrol(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef UV_KIRBY
      IF (Aout(idUwav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUwav(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVwav,ng)) THEN
        allocate ( AVERAGE(ng) % avgVwav(LBi:UBi,LBj:UBj) )
      END IF
# endif
# ifdef ICE_MODEL
      IF (Aout(idUice,ng)) THEN
        allocate ( AVERAGE(ng) % avguice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idVice,ng)) THEN
        allocate ( AVERAGE(ng) % avgvice(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idUiceE,ng)) THEN
        allocate ( AVERAGE(ng) % avguiceE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idViceN,ng)) THEN
        allocate ( AVERAGE(ng) % avgviceN(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout(idAice,ng)) THEN
        allocate ( AVERAGE(ng) % avgaice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idHice,ng)) THEN
        allocate ( AVERAGE(ng) % avghice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTice,ng)) THEN
        allocate ( AVERAGE(ng) % avgtice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTimid,ng)) THEN
        allocate ( AVERAGE(ng) % avgtimid(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idHsno,ng)) THEN
        allocate ( AVERAGE(ng) % avghsno(LBi:UBi,LBj:UBj) )
      END IF
#  ifdef MELT_PONDS
      IF (Aout(idApond,ng)) THEN
        allocate ( AVERAGE(ng) % avgapond(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idHpond,ng)) THEN
        allocate ( AVERAGE(ng) % avghpond(LBi:UBi,LBj:UBj) )
      END IF
#  endif
      IF (Aout(idIomflx,ng)) THEN
        allocate ( AVERAGE(ng) % avgiomflx(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idAgeice,ng)) THEN
        allocate ( AVERAGE(ng) % avgageice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idSig11,ng)) THEN
        allocate ( AVERAGE(ng) % avgsig11(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idSig12,ng)) THEN
        allocate ( AVERAGE(ng) % avgsig12(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idSig22,ng)) THEN
        allocate ( AVERAGE(ng) % avgsig22(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idT0mk,ng)) THEN
        allocate ( AVERAGE(ng) % avgT0mk(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idS0mk,ng)) THEN
        allocate ( AVERAGE(ng) % avgS0mk(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWfr,ng)) THEN
        allocate ( AVERAGE(ng) % avgWfr(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWai,ng)) THEN
        allocate ( AVERAGE(ng) % avgWai(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWao,ng)) THEN
        allocate ( AVERAGE(ng) % avgWao(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWio,ng)) THEN
        allocate ( AVERAGE(ng) % avgWio(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWro,ng)) THEN
        allocate ( AVERAGE(ng) % avgWro(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idWdiv,ng)) THEN
        allocate ( AVERAGE(ng) % avgWdiv(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idTauiw,ng)) THEN
        allocate ( AVERAGE(ng) % avgutau_iw(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idChuiw,ng)) THEN
        allocate ( AVERAGE(ng) % avgchu_iw(LBi:UBi,LBj:UBj) )
      END IF
# endif
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        allocate ( AVERAGE(ng) % avgU2(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idU2av,ng)) THEN
        allocate ( AVERAGE(ng) % avgV2(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(idV2av,ng)) THEN
        allocate ( AVERAGE(ng) % avgZZ(LBi:UBi,LBj:UBj) )
      END IF

# ifdef SOLVE3D
      IF (Aout(idUUav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUU(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idVVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgVV(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idUVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUV(LBi:UBi,LBj:UBj,N(ng)) )
      END IF

      IF (Aout(idHUav,ng)) THEN
        allocate ( AVERAGE(ng) % avgHuon(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(idHVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgHvom(LBi:UBi,LBj:UBj,N(ng)) )
      END IF

      IF (ANY(Aout(idTTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgTT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(idUTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgUT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(idVTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgVT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF

      IF (ANY(Aout(iHUTav(:),ng))) THEN
        allocate ( AVERAGE(ng)% avgHuonT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
      IF (ANY(Aout(iHVTav(:),ng))) THEN
        allocate ( AVERAGE(ng)% avgHvomT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      END IF
# endif
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        allocate ( AVERAGE(ng) % avgpvor2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout(id2dRV,ng)) THEN
        allocate ( AVERAGE(ng) % avgrvor2d(LBi:UBi,LBj:UBj) )
      END IF
# ifdef SOLVE3D
      IF (Aout(id3dPV,ng)) THEN
        allocate ( AVERAGE(ng) % avgpvor3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
      IF (Aout(id3dRV,ng)) THEN
        allocate ( AVERAGE(ng) % avgrvor3d(LBi:UBi,LBj:UBj,N(ng)) )
      END IF
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
      USE mod_ncparam
      USE mod_scalars
# if defined BEST_NPZ
      USE mod_biology
# endif
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
      integer :: itrc, itrc2, k
# endif

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgzeta(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2d(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idu2dE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2dE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idv2dN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2dN(i,j) = IniVal
          END DO
        END DO
      END IF

# ifdef SOLVE3D
      IF (Aout(idUvel,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVvel,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idu3dE,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3dE(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idv3dN,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3dN(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idOvel,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgw3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWvel,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgwvel(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idDano,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgrho(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgt(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

#  ifdef NEMURO_SED1
      IF (Aout(idPONsed,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgPONsed(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idOPALsed,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgOPALsed(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idDENITsed,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDENITsed(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idPONbur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgPONbur(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idOPALbur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgOPALbur(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  ifdef PRIMARY_PROD
      IF (Aout(idNPP,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgNPP(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKv(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTdif,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKt(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSdif,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKs(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif

#  ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghsbl(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghbbl(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
      IF (Aout(idUfx1,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDU_avg1(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUfx2,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDU_avg2(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idVfx1,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDV_avg1(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVfx2,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDV_avg2(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsvs(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbvs(i,j) = IniVal
          END DO
        END DO
      END IF

# ifdef BBL
      IF (Aout(idUbrs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbrs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbrs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbrs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUbws,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbws(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbws,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbws(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUbcs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbcs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbcs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbcs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUVwc,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUVwc(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUbot,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUbot(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbot,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgVbot(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUbur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUbur(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbvr,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgVbvr(i,j) = IniVal
          END DO
        END DO
      END IF
# endif

# ifdef SOLVE3D
#  if defined RUNOFF
      IF (Aout(idRunoff,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgRunoff(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      IF (Aout(idPair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgPair(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
      IF (Aout(idUair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUwind(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgVwind(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUairE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUwindE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVairN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgVwindN(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
      IF (Aout(idTsur(itemp),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgstf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTsur(isalt),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgswf(i,j) = IniVal
          END DO
        END DO
      END IF

#  ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsrf(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  ifdef BULK_FLUXES
      IF (Aout(idLhea,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avglhf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idLrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avglrf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idShea,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgshf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idSSSf,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsssflx(i,j) = IniVal
          END DO
        END DO
      END IF

#   ifdef EMINUSP
      IF (Aout(idevap,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgevap(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idrain,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgrain(i,j) = IniVal
          END DO
        END DO
      END IF
#   endif
#  endif
# endif

# ifdef WEC
!
!  Time-averaged WEC.
!
      IF (Aout(idU2Sd,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2Sd(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2Sd,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2Sd(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idU2rs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2rs(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2rs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2rs(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idWztw,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWztw(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWqsp,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWqsp(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWbeh,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWbeh(i,j) = IniVal
          END DO
        END DO
      END IF

#  ifdef WEC_MELLOR
      IF (Aout(idW2xx,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgSxx2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idW2xy,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgSxy2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idW2yy,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgSyy2d(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  ifdef SOLVE3D
      IF (Aout(idU3Sd,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3Sd(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idV3Sd,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3Sd(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3Sd,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgW3Sd(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3St,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgW3St(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idU3rs,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3rs(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idV3rs,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3rs(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
#   ifdef WEC_MELLOR
      IF (Aout(idW3xx,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSxx3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3xy,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSxy3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3yy,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSyy3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3zx,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSzx3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3zy,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSzy3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
#   endif
#  endif
# endif
# ifdef WAVES_HEIGHT
      IF (Aout(idWamp,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWamp(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef WAVES_LENGTH
      IF (Aout(idWlen,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWlen(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef WAVES_LENGTHP
      IF (Aout(idWlep,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWlep(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef WAVES_DIR
      IF (Aout(idWdir,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWdir(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef WAVES_TOP_PERIOD
      IF (Aout(idWptp,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWptp(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef WAVES_BOT_PERIOD
      IF (Aout(idWpbt,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWpbt(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef BBL_MODEL
      IF (Aout(idWorb,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWorb(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      IF (Aout(idWdif,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWdif(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
      IF (Aout(idWdib,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWdib(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWdiw,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWdiw(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef ROLLER_SVENDSEN
      IF (Aout(idWbrk,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWbrk(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef WEC_ROLLER
      IF (Aout(idWdis,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWdis(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWrol,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWrol(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
# ifdef UV_KIRBY
      IF (Aout(idUwav,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUwav(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVwav,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgVwav(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgU2(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idU2av,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgV2(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2av,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgZZ(i,j) = IniVal
          END DO
        END DO
      END IF

# ifdef SOLVE3D
      IF (Aout(idUUav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgUU(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgVV(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idUVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgUV(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idHUav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgHuon(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgHvom(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (ANY(Aout(idTTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgTT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idUTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgUT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idVTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgVT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (ANY(Aout(iHUTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgHuonT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(iHVTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgHvomT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
# endif
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgpvor2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(id2dRV,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgrvor2d(i,j) = IniVal
          END DO
        END DO
      END IF
# ifdef SOLVE3D

#ifdef BEST_NPZ
# ifdef STATIONARY
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgst(i,j,k,itrc) = IniVal
            END DO
          END DO
        END DO
# endif
# ifdef STATIONARY2
        DO itrc=1,NTS2(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgst2(i,j,itrc) = IniVal
          END DO
        END DO
# endif
# ifdef PROD3
        DO itrc=1,NPT3(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgpt3(i,j,k,itrc) = IniVal
            END DO
          END DO
        END DO
# endif
# ifdef PROD2
        DO itrc=1,NPT2(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgpt2(i,j,itrc) = IniVal
          END DO
        END DO
# endif
# ifdef BENTHIC
        DO itrc=1,NBeT(ng)
          DO k=1,NBL(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgbt(i,j,k,itrc) = IniVal
              END DO
          END DO
        END DO
# endif
# ifdef ICE_BIO
#  ifdef CLIM_ICE_1D
        DO itrc=1,NIceT(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgit(i,j,itrc) = IniVal
          END DO
        END DO
#  elif defined BERING_10K       
      IF (Aout(idIceNO3,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgINO3(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idIceNH4,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgINH4(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idIcePhL,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgIPhL(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idIceLog,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgILog(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
# endif
# ifdef BIOFLUX
        DO itrc=1,NT(ng)
          DO itrc2=1,NT(ng)
            AVERAGE(ng) % avgbflx(itrc,itrc2) = IniVal
          END DO
        END DO
# endif

#endif

# ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgst(i,j,k,itrc) = IniVal
            END DO
          END DO
        END DO
# endif
      IF (Aout(id3dPV,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgpvor3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(id3dRV,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgrvor3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
# endif
# ifdef ICE_MODEL
      IF (Aout(idUice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avguice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgvice(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUiceE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avguiceE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idViceN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgviceN(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idAice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgaice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idHice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgtice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTimid,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgtimid(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idHsno,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghsno(i,j) = IniVal
          END DO
        END DO
      END IF
#  ifdef MELT_PONDS
      IF (Aout(idApond,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgapond(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idHpond,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghpond(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
      IF (Aout(idIomflx,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgiomflx(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idAgeice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgageice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idSig11,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsig11(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idSig12,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsig12(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idSig22,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsig22(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idT0mk,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgT0mk(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idS0mk,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgS0mk(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWfr,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWfr(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWai,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWai(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWao,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWao(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWio,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWio(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWro,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWro(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idWdiv,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgWdiv(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idTauiw,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgutau_iw(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idChuiw,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgchu_iw(i,j) = IniVal
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE initialize_average
#endif
      END MODULE mod_average
