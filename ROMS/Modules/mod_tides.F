#include "cppdefs.h"
      MODULE mod_tides
#if defined SSH_TIDES || defined UV_TIDES
!
!svn $Id: mod_tides.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Tidal Components:                                                   !
!                                                                      !
!  Each of the following arrays has a dimension in tidal components    !
!  classified by period:                                               !
!                                                                      !
!    semi-diurnal:  M2, S2, N2, K2  (12.42, 12.00, 12.66, 11.97h)      !
!         diurnal:  K1, O1, P1, Q1  (23.93, 25.82, 24.07, 26.87h)      !
!                                                                      !
!  and other longer periods. The order of these tidal components is    !
!  irrelevant here.  The number of components to USE is depends on     !
!  the regional application.                                           !
!                                                                      !
!  CosOmega     Cosine tidal harmonics for current omega(t).           !
!  SinOmega     Sine tidal harmonics for current omega(t).             !
!  SSH_Tamp     Tidal elevation amplitude (m) at RHO-points.           !
!  SSH_Tphase   Tidal elevation phase (degrees/360) at RHO-points.     !
!  Tperiod      Tidal period (s).                                      !
!  UV_Tangle    Tidal current angle (radians; counterclockwise         !
!                 from EAST and rotated to curvilinear grid) at        !
!                 RHO-points.                                          !
!  UV_Tmajor    Maximum tidal current: tidal ellipse major axis        !
!                 (m/s) at RHO-points.                                 !
!  UV_Tminor    Minimum tidal current: tidal ellipse minor axis        !
!                 (m/s) at RHO-points.                                 !
!  UV_Tphase    Tidal current phase (degrees/360) at RHO-points.       !
!                                                                      !
# if defined AVERAGES_DETIDE && defined AVERAGES
!                                                                      !
!  Detided time-averaged fields via least-squares fitting. Notice that !
!  the harmonics for the state variable have an extra dimension of     !
!  size (0:2*NTC) to store several terms:                              !
!                                                                      !
!               index 0               mean term (accumulated sum)      !
!                     1:NTC           accumulated sine terms           !
!                     NTC+1:2*NTC     accumulated cosine terms         !
!                                                                      !
!  CosW_avg     Current time-average window COS(omega(k)*t).           !
!  CosW_sum     Time-accumulated COS(omega(k)*t).                      !
!  SinW_avg     Current time-average window SIN(omega(k)*t).           !
!  SinW_sum     Time-accumulated SIN(omega(k)*t).                      !
!  CosWCosW     Time-accumulated COS(omega(k)*t)*COS(omega(l)*t).      !
!  SinWSinW     Time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t).      !
!  SinWCosW     Time-accumulated SIN(omega(k)*t)*COS(omega(l)*t).      !
!                                                                      !
!  ubar_detided Time-averaged and detided 2D u-momentum.               !
!  ubar_tide    Time-accumulated 2D u-momentum tide harmonics.         !
!  vbar_detided Time-averaged and detided 2D v-momentum.               !
!  vbar_tide    Time-accumulated 2D v-momentum tide harmonics.         !
!  zeta_detided Time-averaged and detided free-surface.                !
!  zeta_tide    Time-accumulated free-surface tide harmonics.          !
#  ifdef SOLVE3D
!  u_detided    Time-averaged and detided 3D u-momentum.               !
!  u_tide       Time-accumulated 3D u-momentum tide harmonics.         !
!  v_detided    Time-averaged and detided 3D v-momentum.               !
!  v_tide       Time-accumulated 3D v-momentum tide harmonics.         !
#  endif
!
# endif
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_TIDES

          real(r8), pointer :: Tperiod(:)
# if defined AVERAGES_DETIDE && defined AVERAGES
          real(r8), pointer :: CosOmega(:)
          real(r8), pointer :: SinOmega(:)
          real(r8), pointer :: CosW_avg(:)
          real(r8), pointer :: CosW_sum(:)
          real(r8), pointer :: SinW_avg(:)
          real(r8), pointer :: SinW_sum(:)
          real(r8), pointer :: CosWCosW(:,:)
          real(r8), pointer :: SinWSinW(:,:)
          real(r8), pointer :: SinWCosW(:,:)
# endif
# if defined SSH_TIDES
          real(r8), pointer :: SSH_Tamp(:,:,:)
          real(r8), pointer :: SSH_Tphase(:,:,:)
# endif
# if defined UV_TIDES
          real(r8), pointer :: UV_Tangle(:,:,:)
          real(r8), pointer :: UV_Tmajor(:,:,:)
          real(r8), pointer :: UV_Tminor(:,:,:)
          real(r8), pointer :: UV_Tphase(:,:,:)
# endif
# if defined AVERAGES_DETIDE && defined AVERAGES
          real(r8), pointer :: ubar_detided(:,:)
          real(r8), pointer :: ubar_tide(:,:,:)

          real(r8), pointer :: vbar_detided(:,:)
          real(r8), pointer :: vbar_tide(:,:,:)

          real(r8), pointer :: zeta_detided(:,:)
          real(r8), pointer :: zeta_tide(:,:,:)
#  ifdef SOLVE3D
          real(r8), pointer :: u_detided(:,:,:)
          real(r8), pointer :: u_tide(:,:,:,:)

          real(r8), pointer :: v_detided(:,:,:)
          real(r8), pointer :: v_tide(:,:,:,:)
#  endif
# endif

        END TYPE T_TIDES

        TYPE (T_TIDES), allocatable :: TIDES(:)

      CONTAINS

      SUBROUTINE allocate_tides (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#  if defined AVERAGES_DETIDE && defined AVERAGES
      USE distribute_mod, ONLY : mp_bcastl
#  endif
# endif
!
! Inported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, Nrec, nvd, Vsize(4)
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( TIDES(Ngrids) )
!
!  Inquire about the maximum number of tidal components.
!
      IF (ng.eq.1) THEN
        IF (InpThread) THEN
          CALL opencdf (ng, nFfiles(ng), FRCname(1,ng), ncfile, N(ng),  &
     &                  idTper, Nrec, nvd, Vsize)
          MTC=Vsize(1)
# if defined AVERAGES_DETIDE && defined AVERAGES
!
!  Check if detiding variables are available and set definition switch.
!
          DO i=1,nvars
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsuH))) THEN
              LdefTIDE(ng)=.FALSE.
            END IF
          END DO
# endif
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, MTC, 1)
#  if defined AVERAGES_DETIDE && defined AVERAGES
        CALL mp_bcastl (ng, iNLM, LdefTIDE, Ngrids)
#  endif
# endif
      END IF
!
!  Allocate tidal forcing variables.
!
      allocate ( TIDES(ng) % Tperiod(MTC)  )

# if defined AVERAGES_DETIDE && defined AVERAGES
      allocate ( TIDES(ng) % CosOmega(MTC) )
      allocate ( TIDES(ng) % SinOmega(MTC) )
      allocate ( TIDES(ng) % CosW_avg(MTC) )
      allocate ( TIDES(ng) % CosW_sum(MTC) )
      allocate ( TIDES(ng) % SinW_avg(MTC) )
      allocate ( TIDES(ng) % SinW_sum(MTC) )
      allocate ( TIDES(ng) % CosWCosW(MTC,MTC) )
      allocate ( TIDES(ng) % SinWSinW(MTC,MTC) )
      allocate ( TIDES(ng) % SinWCosW(MTC,MTC) )
# endif

# if defined SSH_TIDES
      allocate ( TIDES(ng) % SSH_Tamp(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % SSH_Tphase(LBi:UBi,LBj:UBj,MTC) )
# endif

# if defined UV_TIDES
      allocate ( TIDES(ng) % UV_Tangle(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % UV_Tmajor(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % UV_Tminor(LBi:UBi,LBj:UBj,MTC) )
      allocate ( TIDES(ng) % UV_Tphase(LBi:UBi,LBj:UBj,MTC) )
# endif

# if defined AVERAGES_DETIDE && defined AVERAGES
      allocate ( TIDES(ng) % ubar_detided(LBi:UBi,LBj:UBj) )
      allocate ( TIDES(ng) % ubar_tide(LBi:UBi,LBj:UBj,0:2*MTC) )

      allocate ( TIDES(ng) % vbar_detided(LBi:UBi,LBj:UBj) )
      allocate ( TIDES(ng) % vbar_tide(LBi:UBi,LBj:UBj,0:2*MTC) )

      allocate ( TIDES(ng) % zeta_detided(LBi:UBi,LBj:UBj) )
      allocate ( TIDES(ng) % zeta_tide(LBi:UBi,LBj:UBj,0:2*MTC) )
#  ifdef SOLVE3D
      allocate ( TIDES(ng) % u_detided(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( TIDES(ng) % u_tide(LBi:UBi,LBj:UBj,N(ng),0:2*MTC) )

      allocate ( TIDES(ng) % v_detided(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( TIDES(ng) % v_tide(LBi:UBi,LBj:UBj,N(ng),0:2*MTC) )
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_tides

      SUBROUTINE initialize_tides (ng, tile)
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
      integer :: i, itide, j, jtide, k

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_TEST) THEN
        DO itide=1,MTC
          TIDES(ng) % Tperiod(itide) = IniVal
        END DO
# if defined AVERAGES_DETIDE && defined AVERAGES
        DO jtide=1,MTC
          TIDES(ng) % CosOmega(jtide) = IniVal
          TIDES(ng) % SinOmega(jtide) = IniVal
          TIDES(ng) % CosW_avg(jtide) = IniVal
          TIDES(ng) % CosW_sum(jtide) = IniVal
          TIDES(ng) % SinW_avg(jtide) = IniVal
          TIDES(ng) % SinW_sum(jtide) = IniVal
          DO itide=1,MTC
            TIDES(ng) % CosWCosW(itide,jtide) = IniVal
            TIDES(ng) % SinWSinW(itide,jtide) = IniVal
            TIDES(ng) % SinWCosW(itide,jtide) = IniVal
          END DO
        END DO
# endif
      END IF

      DO itide=1,MTC
# if defined SSH_TIDES
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % SSH_Tamp(i,j,itide) = IniVal
            TIDES(ng) % SSH_Tphase(i,j,itide) = IniVal
          END DO
        END DO
# endif
# if defined UV_TIDES
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            TIDES(ng) % UV_Tangle(i,j,itide) = IniVal
            TIDES(ng) % UV_Tmajor(i,j,itide) = IniVal
            TIDES(ng) % UV_Tminor(i,j,itide) = IniVal
            TIDES(ng) % UV_Tphase(i,j,itide) = IniVal
          END DO
        END DO
# endif
      END DO

# if defined AVERAGES_DETIDE && defined AVERAGES
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          TIDES(ng) % ubar_detided(i,j) = IniVal
          TIDES(ng) % vbar_detided(i,j) = IniVal
          TIDES(ng) % zeta_detided(i,j) = IniVal
          DO itide=0,2*MTC
            TIDES(ng) % ubar_tide(i,j,itide) = IniVal
            TIDES(ng) % vbar_tide(i,j,itide) = IniVal
            TIDES(ng) % zeta_tide(i,j,itide) = IniVal
          END DO
        END DO
#  ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=Imin,Imax
            TIDES(ng) % u_detided(i,j,k) = IniVal
            TIDES(ng) % v_detided(i,j,k) = IniVal
            DO itide=0,2*MTC
              TIDES(ng) % u_tide(i,j,k,itide) = IniVal
              TIDES(ng) % v_tide(i,j,k,itide) = IniVal
            END DO
          END DO
        END DO
#  endif
      END DO
# endif

      RETURN
      END SUBROUTINE initialize_tides
#endif
      END MODULE mod_tides
