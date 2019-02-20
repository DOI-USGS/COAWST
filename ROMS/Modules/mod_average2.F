#include "cppdefs.h"
      MODULE mod_average2
#ifdef AVERAGES2
!
!svn $Id: mod_average.F 707 2008-08-19 17:58:55Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  FOR the AVERAGE2 option, these become surface only fields.          !
!  The strategy here is  to define all possible pointers in the        !
!  time-averaged structure and allocate only those requested by        !
!  the user. This will facilitate a better management of memory.       !
!                                                                      !
!  Time-averaged fields for output purposes.                           !
!                                                                      !
!  avgu2d     2D velocity component (m/s) in the XI-direction.         !
!  avgv2d     2D velocity component (m/s) in the ETA-direction.        !
!  avgu2dE    2D Eastward component (m/s) at RHO-points.               !
!  avgv2dN    2D Northward component (m/s) at RHO-points.              !
!  avgzeta    Free surface (m).                                        !

# ifdef SOLVE3D
!  avguwind   2D wind velocity component (m/s) in the XI-direction.    !
!  avgvwind   2D wind velocity component (m/s) in the ETA-direction.   !
!  avguwindE  2D wind velocity component (m/s) to the East.            !
!  avgvwindN  2D wind velocity component (m/s) to the North.           !
!  avgu3d     3D velocity component (m/s) in the XI-direction.         !
!  avgv3d     3D velocity component (m/s) in the ETA-direction.        !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgt       Tracer type variables (usually, potential temperature    !
!              and salinity).                                          !

#  if defined LMD_SKPP
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
#  endif
#  if defined LMD_BKPP
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
#  endif
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                         !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).             !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                         !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).             !
!  avgbus     Bottom u-momentum stress (N/m2).                         !
!  avgbvs     Bottom v-momentum stress (N/m2).                         !
!  avgsssflx  Sea surface salinity flux correction.                    !
!  avgshf     Sensible heat flux (W/m2).                               !
#  ifdef SHORTWAVE
!  avgsrf     Shortwave radiation flux (W/m2).                         !
#  endif
#  ifdef BULK_FLUXES
!  avglhf     Latent heat flux (W/m2).                                 !
!  avglrf     Longwave radiation flux (W/m2).                          !
#  endif
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net salt flux (kg/m2/s).                         !
!  avgevap    Surface net evaporation (kg/m2/s).                       !
!  avgrain    Surface net rain fall (kg/m2/s).                         !
# endif
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_AVERAGE2

          real(r8), pointer :: avgzeta(:,:)

          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)

          real(r8), pointer :: avgu2dE(:,:)
          real(r8), pointer :: avgv2dN(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: avgrho(:,:)
          real(r8), pointer :: avgt(:,:,:)
#  ifdef BIO_GOANPZ
          real(r8), pointer :: avgst(:,:,:)
#  endif
          real(r8), pointer :: avgu3d(:,:)
          real(r8), pointer :: avgv3d(:,:)
          real(r8), pointer :: avgu3dE(:,:)
          real(r8), pointer :: avgv3dN(:,:)

          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
#  ifdef BULK_FLUXES
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
          real(r8), pointer :: avguwind(:,:)
          real(r8), pointer :: avgvwind(:,:)
          real(r8), pointer :: avguwindE(:,:)
          real(r8), pointer :: avgvwindN(:,:)
          real(r8), pointer :: avgsssflx(:,:)
#   ifdef EMINUSP
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
#   endif
#  endif
#  ifdef SHORTWAVE
          real(r8), pointer :: avgsrf(:,:)
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
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
          real(r8), pointer :: avgchu_iw(:,:)
          real(r8), pointer :: avgutau_iw(:,:)
#  endif
# endif
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)
          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)

        END TYPE T_AVERAGE2

        TYPE (T_AVERAGE2), allocatable :: AVERAGE2(:)

      CONTAINS

      SUBROUTINE allocate_average2 (ng, LBi, UBi, LBj, UBj)
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
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1 ) allocate ( AVERAGE2(Ngrids) )
!
!  Time-averaged state variables.
!
      IF (Aout2(idFsur,ng)) THEN
        allocate ( AVERAGE2(ng) % avgzeta(LBi:UBi,LBj:UBj) )
      ENDIF

      IF (Aout2(idUbar,ng)) THEN
        allocate ( AVERAGE2(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVbar,ng)) THEN
        allocate ( AVERAGE2(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout2(idu2dE,ng)) THEN
        allocate ( AVERAGE2(ng) % avgu2dE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idv2dN,ng)) THEN
        allocate ( AVERAGE2(ng) % avgv2dN(LBi:UBi,LBj:UBj) )
      END IF
# ifdef SOLVE3D
      IF (Aout2(idUvel,ng)) THEN
        allocate ( AVERAGE2(ng) % avgu3d(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVvel,ng)) THEN
        allocate ( AVERAGE2(ng) % avgv3d(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout2(idu3dE,ng)) THEN
        allocate ( AVERAGE2(ng) % avgu3dE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idv3dN,ng)) THEN
        allocate ( AVERAGE2(ng) % avgv3dN(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout2(idDano,ng)) THEN
        allocate ( AVERAGE2(ng) % avgrho(LBi:UBi,LBj:UBj) )
      END IF
      IF (ANY(Aout2(idTvar(:),ng))) THEN
        allocate ( AVERAGE2(ng) % avgt(LBi:UBi,LBj:UBj,NT(ng)) )
      END IF

#  ifdef LMD_BKPP
      IF (Aout2(idHbbl,ng)) THEN
        allocate ( AVERAGE2(ng) % avghbbl(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  ifdef LMD_SKPP
      IF (Aout2(idHsbl,ng)) THEN
        allocate ( AVERAGE2(ng) % avghsbl(LBi:UBi,LBj:UBj) )
      END IF
#  endif
# endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout2(idUsms,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsus(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVsms,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout2(idUbms,ng)) THEN
        allocate ( AVERAGE2(ng) % avgbus(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVbms,ng)) THEN
        allocate ( AVERAGE2(ng) % avgbvs(LBi:UBi,LBj:UBj) )
      END IF
# ifdef SOLVE3D
#  ifdef BIO_GOANPZ
        allocate ( AVERAGE2(ng) % avgst(LBi:UBi,LBj:UBj,NTS(ng)) )
#  endif
      IF (Aout2(idTsur(itemp),ng)) THEN
        allocate ( AVERAGE2(ng) % avgstf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idTsur(isalt),ng)) THEN
        allocate ( AVERAGE2(ng) % avgswf(LBi:UBi,LBj:UBj) )
      END IF
#  ifdef SHORTWAVE
      IF (Aout2(idSrad,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsrf(LBi:UBi,LBj:UBj) )
      END IF
#  endif
#  ifdef BULK_FLUXES
      IF (Aout2(idLhea,ng)) THEN
        allocate ( AVERAGE2(ng) % avglhf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idLrad,ng)) THEN
        allocate ( AVERAGE2(ng) % avglrf(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idShea,ng)) THEN
        allocate ( AVERAGE2(ng) % avgshf(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout2(idUair,ng)) THEN
        allocate ( AVERAGE2(ng) % avguwind(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVair,ng)) THEN
        allocate ( AVERAGE2(ng) % avgvwind(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idUairE,ng)) THEN
        allocate ( AVERAGE2(ng) % avguwindE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVairN,ng)) THEN
        allocate ( AVERAGE2(ng) % avgvwindN(LBi:UBi,LBj:UBj) )
      END IF

      IF (Aout2(idSSSf,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsssflx(LBi:UBi,LBj:UBj) )
      END IF
#   ifdef EMINUSP
      IF (Aout2(idevap,ng)) THEN
        allocate ( AVERAGE2(ng) % avgevap(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idrain,ng)) THEN
        allocate ( AVERAGE2(ng) % avgrain(LBi:UBi,LBj:UBj) )
      END IF
#   endif
#  endif

#  ifdef ICE_MODEL
      IF (Aout2(idUice,ng)) THEN
        allocate ( AVERAGE2(ng) % avguice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idVice,ng)) THEN
        allocate ( AVERAGE2(ng) % avgvice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idUiceE,ng)) THEN
        allocate ( AVERAGE2(ng) % avguiceE(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idViceN,ng)) THEN
        allocate ( AVERAGE2(ng) % avgviceN(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idAice,ng)) THEN
        allocate ( AVERAGE2(ng) % avgaice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idHice,ng)) THEN
        allocate ( AVERAGE2(ng) % avghice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idTice,ng)) THEN
        allocate ( AVERAGE2(ng) % avgtice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idTimid,ng)) THEN
        allocate ( AVERAGE2(ng) % avgtimid(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idHsno,ng)) THEN
        allocate ( AVERAGE2(ng) % avghsno(LBi:UBi,LBj:UBj) )
      END IF
#   ifdef MELT_PONDS
      IF (Aout2(idApond,ng)) THEN
        allocate ( AVERAGE2(ng) % avgapond(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idHpond,ng)) THEN
        allocate ( AVERAGE2(ng) % avghpond(LBi:UBi,LBj:UBj) )
      END IF
#   endif
      IF (Aout2(idIomflx,ng)) THEN
        allocate ( AVERAGE2(ng) % avgiomflx(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idAgeice,ng)) THEN
        allocate ( AVERAGE2(ng) % avgageice(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idSig11,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsig11(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idSig12,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsig12(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idSig22,ng)) THEN
        allocate ( AVERAGE2(ng) % avgsig22(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idT0mk,ng)) THEN
        allocate ( AVERAGE2(ng) % avgT0mk(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idS0mk,ng)) THEN
        allocate ( AVERAGE2(ng) % avgS0mk(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idWfr,ng)) THEN
        allocate ( AVERAGE2(ng) % avgWfr(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idWai,ng)) THEN
        allocate ( AVERAGE2(ng) % avgWai(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idWao,ng)) THEN
        allocate ( AVERAGE2(ng) % avgWao(LBi:UBi,LBj:UBj) ) 
      END IF
      IF (Aout2(idWio,ng)) THEN
        allocate ( AVERAGE2(ng) % avgWio(LBi:UBi,LBj:UBj) ) 
      END IF
      IF (Aout2(idWro,ng)) THEN
        allocate ( AVERAGE2(ng) % avgWro(LBi:UBi,LBj:UBj) ) 
      END IF
      IF (Aout2(idTauiw,ng)) THEN
        allocate ( AVERAGE2(ng) % avgutau_iw(LBi:UBi,LBj:UBj) )
      END IF
      IF (Aout2(idChuiw,ng)) THEN
        allocate ( AVERAGE2(ng) % avgchu_iw(LBi:UBi,LBj:UBj) )
      END IF
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_average2

      SUBROUTINE initialize_average2 (ng, tile)
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
      IF (Aout2(idFsur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgzeta(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idUbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgu2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgv2d(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idu2dE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgu2dE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idv2dN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgv2dN(i,j) = IniVal
          END DO
        END DO
      END IF
# ifdef SOLVE3D
      IF (Aout2(idUvel,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgu3d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVvel,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgv3d(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idu3dE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgu3dE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idV3dN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgV3dN(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idDano,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgrho(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (ANY(Aout2(idTvar(:),ng))) THEN
        DO itrc=1,NT(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE2(ng) % avgt(i,j,itrc) = IniVal
            END DO
          END DO
        END DO
      END IF
#  ifdef LMD_BKPP
      IF (Aout2(idHbbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avghbbl(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  ifdef LMD_SKPP
      IF (Aout2(idHsbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avghsbl(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
# endif
      IF (Aout2(idUsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsvs(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idUbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgbus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgbvs(i,j) = IniVal
          END DO
        END DO
      END IF
# ifdef SOLVE3D
      IF (Aout2(idTsur(itemp),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgstf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idTsur(isalt),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgswf(i,j) = IniVal
          END DO
        END DO
      END IF
#  ifdef SHORTWAVE
      IF (Aout2(idSrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsrf(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  ifdef BULK_FLUXES
      IF (Aout2(idLhea,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avglhf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idLrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avglrf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idShea,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgshf(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idUair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avguwind(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgvwind(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idUairE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avguwindE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVairN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgvwindN(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout2(idSSSf,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsssflx(i,j) = IniVal
          END DO
        END DO
      END IF
#   ifdef EMINUSP
      IF (Aout2(idevap,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgevap(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idrain,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgrain(i,j) = IniVal
          END DO
        END DO
      END IF
#   endif
#  endif
#  ifdef ICE_MODEL
      IF (Aout2(idUice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avguice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idVice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgvice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idAice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgaice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idHice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avghice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idTice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgtice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idTimid,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgtimid(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idHsno,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avghsno(i,j) = IniVal
          END DO
        END DO
      END IF
#   ifdef MELT_PONDS
      IF (Aout2(idApond,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgapond(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idHpond,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avghpond(i,j) = IniVal
          END DO
        END DO
      END IF
#   endif
      IF (Aout2(idIomflx,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgiomflx(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idAgeice,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgageice(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idSig11,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsig11(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idSig12,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsig12(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idSig22,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgsig22(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idT0mk,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgT0mk(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idS0mk,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgS0mk(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idWfr,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgWfr(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idWai,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgWai(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idWao,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgWao(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idWio,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgWio(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idWro,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgWro(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idTauiw,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgutau_iw(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout2(idChuiw,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgchu_iw(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  ifdef BIO_GOANPZ
      DO itrc=1,NTS(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE2(ng) % avgst(i,j,itrc) = IniVal
          END DO
        END DO
      END DO
#  endif
# endif


      RETURN
      END SUBROUTINE initialize_average2
#endif
      END MODULE mod_average2
