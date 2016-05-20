#include "cppdefs.h"
      MODULE mod_filter
#ifdef FILTERED
!
!=======================================================================
!  Copyright (c) 2002-2014 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  2D Time-filtered fields for output purposes.                        !
!                                                                      !
!  filu2d     2D velocity component (m/s) in the XI-direction.         !
!  filv2d     2D velocity component (m/s) in the ETA-direction.        !
!  filzeta    Free surface (m).                                        !
!                                                                      !
!  3D Time-filtered fields for output purposes.                        !
!                                                                      !
!  filrho     Density anomaly (kg/m3).                                 !
!  filt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  filu3d     3D velocity component (m/s) in the XI-direction.         !
!  filv3d     3D velocity component (m/s) in the ETA-direction.        !
!  filomega   S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!  filw3d     Vertical W-velocity (m/s).                               !
!                                                                      !
!  2D Time-filtered ice fields for output purposes.                    !
!                                                                      !
!  filrho     Density anomaly (kg/m3).                                 !
!                                                                      !
!  Filter characteristics                                              !
!                                                                      !
!  FILN       Number of concurrent filter sums.                        !
!  dtfil      Filter sampling interval in seconds.                     !
!  nw         Number of filter weights.                                !
!  nfil       Filter length in timesteps.                              !
!  fcount     Tracks how many filtering intervals have been summed.    !
!  fw         Filter weights.                                          !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_FILTER

          real(r8), pointer :: filzeta(:,:,:)
          real(r8), pointer :: filu2d(:,:,:)
          real(r8), pointer :: filv2d(:,:,:)
# ifdef SOLVE3D
          real(r8), pointer :: filu3d(:,:,:,:)
          real(r8), pointer :: filv3d(:,:,:,:)
          real(r8), pointer :: filrho(:,:,:,:)
          real(r8), pointer :: filt(:,:,:,:,:)
          real(r8), pointer :: filomega(:,:,:,:)
          real(r8), pointer :: filw3d(:,:,:,:)
#  ifdef LMD_SKPP
          real(r8), pointer :: filhsbl(:,:,:)
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: filhbbl(:,:,:)
#  endif
          real(r8), pointer :: filAKv(:,:,:,:)
          real(r8), pointer :: filAKt(:,:,:,:)
          real(r8), pointer :: filAKs(:,:,:,:)
          real(r8), pointer :: filstf(:,:,:)
          real(r8), pointer :: filswf(:,:,:)
#  ifdef BULK_FLUXES
          real(r8), pointer :: fillhf(:,:,:)
          real(r8), pointer :: filshf(:,:,:)
          real(r8), pointer :: fillrf(:,:,:)
#  endif
#  ifdef SHORTWAVE
          real(r8), pointer :: filsrf(:,:,:)
#  endif
#  ifdef ICE_MODEL
          real(r8), pointer :: filui(:,:,:)
          real(r8), pointer :: filvi(:,:,:)
          real(r8), pointer :: filai(:,:,:)
          real(r8), pointer :: filhi(:,:,:)
          real(r8), pointer :: filhsn(:,:,:)
#   ifdef MELT_PONDS
          real(r8), pointer :: filapond(:,:,:)
          real(r8), pointer :: filhpond(:,:,:)
#   endif
          real(r8), pointer :: filtis(:,:,:)
          real(r8), pointer :: filti(:,:,:)
          real(r8), pointer :: filAgeice(:,:,:)
          real(r8), pointer :: fils0mk(:,:,:)
          real(r8), pointer :: filt0mk(:,:,:)
          real(r8), pointer :: filWfr(:,:,:)
          real(r8), pointer :: filWai(:,:,:)
          real(r8), pointer :: filWao(:,:,:)
          real(r8), pointer :: filWio(:,:,:)
          real(r8), pointer :: filWro(:,:,:)
          real(r8), pointer :: filWdiv(:,:,:)
          real(r8), pointer :: filsig11(:,:,:)
          real(r8), pointer :: filsig22(:,:,:)
          real(r8), pointer :: filsig12(:,:,:)
          real(r8), pointer :: filio_mflux(:,:,:)
          real(r8), pointer :: filutau_iw(:,:,:)
          real(r8), pointer :: filchu_iw(:,:,:)
#  endif
# endif
          real(r8), pointer :: filsus(:,:,:)
          real(r8), pointer :: filsvs(:,:,:)
          real(r8), pointer :: filbus(:,:,:)
          real(r8), pointer :: filbvs(:,:,:)

        END TYPE T_FILTER

        TYPE (T_FILTER), allocatable :: FILTER(:)

!    Hourly filter weights for pl33 filter.
#undef PL33
!    Hourly filter weights for 5 day lolopass2 filter
#define LOLOPASS2
# ifdef PL33
        integer, parameter :: FILN = 3
        integer, parameter :: dtfil = 3600
        integer, parameter :: nw = 72
        integer :: nfil
        integer, dimension(FILN) :: fcount
        real(r8), dimension(nw) :: fw =                                 &
     &            (/                                                    &
     &             0.0,      0.0,      0.0,      0.0,      0.0,         &
     &            -0.00027, -0.00114, -0.00211, -0.00317, -0.00427,     &
     &            -0.00537, -0.00641, -0.00735, -0.00811, -0.00864,     &
     &            -0.00887, -0.00872, -0.00816, -0.00714, -0.00560,     &
     &            -0.00355, -0.00097,  0.00213,  0.00574,  0.00980,     &
     &             0.01425,  0.01902,  0.02400,  0.02911,  0.03423,     &
     &             0.03923,  0.04399,  0.04842,  0.05237,  0.05576,     &
     &             0.05850,  0.06051,  0.06174,  0.06215,               &
     &             0.06174,  0.06051,  0.05850,                         &
     &             0.05576,  0.05237,  0.04842,  0.04399,  0.03923,     &
     &             0.03423,  0.02911,  0.02400,  0.01902,  0.01425,     &
     &             0.00980,  0.00574,  0.00213, -0.00097, -0.00355,     &
     &            -0.00560, -0.00714, -0.00816, -0.00872, -0.00887,     &
     &            -0.00864, -0.00811, -0.00735, -0.00641, -0.00537,     &
     &            -0.00427, -0.00317, -0.00211, -0.00114, -0.00027      &
     &              /)
# elif defined LOLOPASS2
        integer, parameter :: FILN = 5
        integer, parameter :: dtfil = 3600
        integer, parameter :: nw = 120
        integer :: nfil
        integer, dimension(FILN) :: fcount
        real(r8), dimension(nw) :: fw =                                 &
     &   (/                                                             &
     &     0.00000e+0, 2.24482e-6, 1.01382e-5, 2.69017e-5, 5.52832e-5,  &
     &     9.72726e-5, 1.53843e-4, 2.24728e-4, 3.08265e-4, 4.01287e-4,  &
     &     4.99108e-4, 5.95586e-4, 6.83272e-4, 7.53646e-4, 7.97444e-4,  &
     &     8.05050e-4, 7.66951e-4, 6.74250e-4, 5.19196e-4, 2.95727e-4,  &
     &     0.00000e+0,-3.69122e-4,-8.09640e-4,-1.31607e-3,-1.87920e-3,  &
     &    -2.48599e-3,-3.11952e-3,-3.75924e-3,-4.38118e-3,-4.95844e-3,  &
     &    -5.46179e-3,-5.86035e-3,-6.12246e-3,-6.21655e-3,-6.11216e-3,  &
     &    -5.78090e-3,-5.19754e-3,-4.34092e-3,-3.19491e-3,-1.74920e-3,  &
     &     0.00000e+0, 2.04941e-3, 4.38837e-3, 6.99864e-3, 9.85452e-3,  &
     &     1.29230e-2, 1.61646e-2, 1.95334e-2, 2.29786e-2, 2.64454e-2,  &
     &     2.98760e-2, 3.32112e-2, 3.63914e-2, 3.93586e-2, 4.20574e-2,  &
     &     4.44366e-2, 4.64503e-2, 4.80592e-2, 4.92317e-2, 4.99447e-2,  &
     &     5.01840e-2,                                                  &
     &     4.99447e-2, 4.92317e-2, 4.80592e-2, 4.64503e-2, 4.44366e-2,  &
     &     4.20574e-2, 3.93586e-2, 3.63914e-2, 3.32112e-2, 2.98760e-2,  &
     &     2.64454e-2, 2.29786e-2, 1.95334e-2, 1.61646e-2, 1.29230e-2,  &
     &     9.85452e-3, 6.99864e-3, 4.38837e-3, 2.04941e-3, 0.00000e+0,  &
     &    -1.74920e-3,-3.19491e-3,-4.34092e-3,-5.19754e-3,-5.78090e-3,  &
     &    -6.11216e-3,-6.21655e-3,-6.12246e-3,-5.86035e-3,-5.46179e-3,  &
     &    -4.95844e-3,-4.38118e-3,-3.75924e-3,-3.11952e-3,-2.48599e-3,  &
     &    -1.87920e-3,-1.31607e-3,-8.09640e-4,-3.69122e-4, 0.00000e+0,  &
     &     2.95727e-4, 5.19196e-4, 6.74250e-4, 7.66951e-4, 8.05050e-4,  &
     &     7.97444e-4, 7.53646e-4, 6.83272e-4, 5.95586e-4, 4.99108e-4,  &
     &     4.01287e-4, 3.08265e-4, 2.24728e-4, 1.53843e-4, 9.72726e-5,  &
     &     5.52832e-5, 2.69017e-5, 1.01382e-5, 2.24482e-6               &
     &              /)
# else
        integer, parameter :: FILN = 2
        integer, parameter :: dtfil = 900
        integer, parameter :: nw = 4
        integer :: nfil
        integer, dimension(FILN) :: fcount
        real(r8), dimension(nw) :: fw =                                 &
     &            (/                                                    &
     &             0.25, 0.25, 0.25, 0.25                               &
     &              /)
# endif

      CONTAINS

      SUBROUTINE allocate_filter (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
      IF (ng.eq.1 ) allocate ( FILTER(Ngrids) )
!
      IF (Aout(idFsur,ng)) THEN
        allocate ( FILTER(ng) % filzeta(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
        allocate ( FILTER(ng) % filu2d(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
        allocate ( FILTER(ng) % filv2d(LBi:UBi,LBj:UBj,FILN) )
      END IF
# ifdef SOLVE3D
      IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
        allocate ( FILTER(ng) % filu3d(LBi:UBi,LBj:UBj,N(ng),FILN) )
      END IF
      IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
        allocate ( FILTER(ng) % filv3d(LBi:UBi,LBj:UBj,N(ng),FILN) )
      END IF
      IF (Aout(idDano,ng)) THEN
        allocate ( FILTER(ng) % filrho(LBi:UBi,LBj:UBj,N(ng),FILN) )
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        allocate ( FILTER(ng) % filt(LBi:UBi,LBj:UBj,N(ng),NAT,FILN) )
      END IF
      IF (Aout(idOvel,ng)) THEN
        allocate ( FILTER(ng) % filomega(LBi:UBi,LBj:UBj,0:N(ng),FILN) )
      END IF
      IF (Aout(idWvel,ng)) THEN
        allocate ( FILTER(ng) % filw3d(LBi:UBi,LBj:UBj,0:N(ng),FILN) )
      END IF
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        allocate ( FILTER(ng) % filAKv(LBi:UBi,LBj:UBj,0:N(ng),FILN) )
      END IF
      IF (Aout(idTdif,ng)) THEN
        allocate ( FILTER(ng) % filAKt(LBi:UBi,LBj:UBj,0:N(ng),FILN) )
      END IF
      IF (Aout(idSdif,ng)) THEN
        allocate ( FILTER(ng) % filAKs(LBi:UBi,LBj:UBj,0:N(ng),FILN) )
      END IF
#  endif
#  ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        allocate ( FILTER(ng) % filhsbl(LBi:UBi,LBj:UBj,FILN) )
      END IF
#  endif
#  ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        allocate ( FILTER(ng) % filhbbl(LBi:UBi,LBj:UBj,FILN) )
      END IF
#  endif
# endif
      IF (Aout(idUsms,ng)) THEN
        allocate ( FILTER(ng) % filsus(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idVsms,ng)) THEN
        allocate ( FILTER(ng) % filsvs(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idUbms,ng)) THEN
        allocate ( FILTER(ng) % filbus(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idVbms,ng)) THEN
        allocate ( FILTER(ng) % filbvs(LBi:UBi,LBj:UBj,FILN) )
      END IF
# ifdef SOLVE3D
      IF (Aout(idTsur(itemp),ng)) THEN
        allocate ( FILTER(ng) % filstf(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idTsur(isalt),ng)) THEN
        allocate ( FILTER(ng) % filswf(LBi:UBi,LBj:UBj,FILN) )
      END IF
#  ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        allocate ( FILTER(ng) % filsrf(LBi:UBi,LBj:UBj,FILN) )
      END IF
#  endif
#  ifdef BULK_FLUXES
      IF (Aout(idLhea,ng)) THEN
        allocate ( FILTER(ng) % fillhf(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idLrad,ng)) THEN
        allocate ( FILTER(ng) % filshf(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idShea,ng)) THEN
        allocate ( FILTER(ng) % fillrf(LBi:UBi,LBj:UBj,FILN) )
      END IF
#  endif
#  ifdef ICE_MODEL
      IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
        allocate ( FILTER(NG) % filui(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
        allocate ( FILTER(NG) % filvi(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idAice,ng)) THEN
        allocate ( FILTER(NG) % filai(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idHice,ng)) THEN
        allocate ( FILTER(NG) % filhi(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idTice,ng)) THEN
        allocate ( FILTER(NG) % filtis(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idTimid,ng)) THEN
        allocate ( FILTER(NG) % filti(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idHsno,ng)) THEN
        allocate ( FILTER(NG) % filhsn(LBi:UBi,LBj:UBj,FILN) )
      END IF
#   ifdef MELT_PONDS
      IF (Aout(idApond,ng)) THEN
        allocate ( FILTER(NG) % filapond(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idHpond,ng)) THEN
        allocate ( FILTER(NG) % filhpond(LBi:UBi,LBj:UBj,FILN) )
      END IF
#   endif
      IF (Aout(idIomflx,ng)) THEN
        allocate ( FILTER(NG) % filio_mflux(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idAgeice,ng)) THEN
        allocate ( FILTER(NG) % filAgeice(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idSig11,ng)) THEN
        allocate ( FILTER(NG) % filsig11(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idSig12,ng)) THEN
        allocate ( FILTER(NG) % filsig22(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idSig22,ng)) THEN
        allocate ( FILTER(NG) % filsig12(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idS0mk,ng)) THEN
        allocate ( FILTER(NG) % fils0mk(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idT0mk,ng)) THEN
        allocate ( FILTER(NG) % filt0mk(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idWfr,ng)) THEN
        allocate ( FILTER(NG) % filWfr(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idWai,ng)) THEN
        allocate ( FILTER(NG) % filWai(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idWao,ng)) THEN
        allocate ( FILTER(NG) % filWao(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idWio,ng)) THEN
        allocate ( FILTER(NG) % filWio(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idWro,ng)) THEN
        allocate ( FILTER(NG) % filWro(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idWdiv,ng)) THEN
        allocate ( FILTER(NG) % filWdiv(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idTauiw,ng)) THEN
        allocate ( FILTER(NG) % filutau_iw(LBi:UBi,LBj:UBj,FILN) )
      END IF
      IF (Aout(idChuiw,ng)) THEN
        allocate ( FILTER(NG) % filchu_iw(LBi:UBi,LBj:UBj,FILN) )
      END IF
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_filter

      SUBROUTINE initialize_filter (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
      integer :: i, j, if, k, itrc
      integer :: Imin, Imax, Jmin, Jmax
      real(r8), parameter :: IniVal = 0.0_r8
!
# include "set_bounds.h"

#ifdef _OPENMP 
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
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (Aout(idFsur,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filzeta(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filu2d(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filv2d(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
# ifdef SOLVE3D
      IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
        DO if=1,FILN
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filu3d(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
        DO if=1,FILN
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filv3d(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWvel,ng)) THEN
        DO if=1,FILN
          DO k=0,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filw3d(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idOvel,ng)) THEN
        DO if=1,FILN
          DO k=0,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filomega(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idDano,ng)) THEN
        DO if=1,FILN
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filrho(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        DO if=1,FILN
          DO itrc=1,NAT
            DO k=1,N(ng)
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  FILTER(ng) % filt(i,j,k,itrc,if)= IniVal
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        DO if=1,FILN
          DO k=0,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filAKv(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTdif,ng)) THEN
        DO if=1,FILN
          DO k=0,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filAKt(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSdif,ng)) THEN
        DO if=1,FILN
          DO k=0,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                FILTER(ng) % filAKs(i,j,k,if) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
#  endif
#  ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filhbbl(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif
#  ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filhsbl(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif
# endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filsus(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVsms,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filsvs(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idUbms,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filbus(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVbms,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filbvs(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF

# ifdef SOLVE3D
      IF (Aout(idTsur(itemp),ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filstf(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTsur(isalt),ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filswf(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF

#  ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filsrf(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif

#  ifdef BULK_FLUXES
      IF (Aout(idLhea,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % fillhf(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idLrad,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filshf(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF 
      IF (Aout(idShea,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % fillrf(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif
#  ifdef ICE_MODEL
      IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filui(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filvi(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idAice,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filai(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHice,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filhi(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHsno,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filhsn(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#   ifdef MELT_PONDS
      IF (Aout(idApond,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filapond(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHpond,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filhpond(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#   endif
      IF (Aout(idIomflx,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filio_mflux(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTice,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filtis(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTimid,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filti(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idAgeice,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filAgeice(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSig11,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filsig11(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSig22,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filsig22(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSig12,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filsig12(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idS0mk,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % fils0mk(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idT0mk,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filt0mk(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWfr,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filWfr(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWai,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filWai(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWao,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filWao(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWio,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filWio(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWro,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filWro(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWdiv,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filWdiv(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTauiw,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filutau_iw(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idChuiw,ng)) THEN
        DO if=1,FILN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filchu_iw(i,j,if) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif
# endif
      RETURN
      END SUBROUTINE initialize_filter
#endif
      END MODULE mod_filter

