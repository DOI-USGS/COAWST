#include "cppdefs.h"
      MODULE mod_clima
#if defined AD_SENSITIVITY  || defined CLIMATOLOGY      || \
    defined OBS_SENSITIVITY || defined OPT_OBSERVATIONS || \
    defined SO_SEMI
!
!svn $Id: mod_clima.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sea surface height fields.                                          !
!                                                                      !
!   ssh         Climatology for sea surface height (m).                !
!   sshG        Latest two-time snapshots of input "ssh" grided        !
!                 data used for interpolation.                         !
!   zeta_ads    Sensitivity functional for sea surface height.         !
!   zeta_adsF   Latest two-time snapshots of input "zeta_ads" grided   !
!                 data used fot interpolation.                         !
!                                                                      !
!  2D momentum fields.                                                 !
!                                                                      !
!   ubarclm     Vertically integrated U-momentum climatology (m/s).    !
!   ubarclmG    Latest two-time snapshots of input "ubarclm" grided    !
!                 data used for interpolation.                         !
!   ubar_ads    Sensitivity functional for vertically integrated       !
!                 U-momentum.                                          !
!   ubar_adsG   Latest two-time snapshots of input "ubar_ads" grided   !
!                 data used for interpolation.                         !
!   vbarclm     Vertically integrated V-momentum climatology (m/s).    !
!   vbarclmG    Latest two-time snapshots of input "vbarclm" grided    !
!                 data used for interpolation.                         !
!   vbar_ads    Sensitivity functional for vertically integrated       !
!                 V-momentum.                                          !
!   vbar_adsG   Latest two-time snapshots of input "vbar_ads" grided   !
!                 data used for interpolation.                         !
!                                                                      !
!  Tracer fields.                                                      !
!                                                                      !
!   tclm        Climatology for tracer type variables (usually,        !
!                 temperature: degC; salinity: PSU).                   !
!   tclmG       Latest two-time snapshots of input "tclm" grided       !
!                 data used for interpolation.                         !
!   t_ads       Sensitivity functional for tracer type variables.      !
!   t_adsG      Latest two-time snapshots of input "t_ads" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  3D momentum climatology.                                            !
!                                                                      !
!   uclm        3D U-momentum climatology (m/s).                       !
!   uclmG       Latest two-time snapshots of input "uclm" grided       !
!                 data used for interpolation.                         !
!   u_ads       Sensitivity functional for 3D U-momentum.              !
!   u_adsG      Latest two-time snapshots of input "u_ads" grided      !
!                 data used for interpolation.                         !
!   vclm        3D V-momentum climatology (m/s).                       !
!   vclmG       Latest two-time snapshots of input "vclm" grided       !
!                 data used for interpolation.                         !
!   v_ads       Sensitivity functional for 3D V-momentum.              !
!   v_adsG      Latest two-time snapshots of input "v_ads" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Nudging variables.                                                  !
!                                                                      !
!   M2nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 2D momentum data.                                    !
!   M3nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 3D momentum data.                                    !
!   Tnudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 tracer data.                                         !
!   Znudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 sea surface height data.                             !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_CLIMA

# ifdef ZCLIMATOLOGY
          real(r8), pointer :: ssh(:,:)
#  ifndef ANA_SSH
          real(r8), pointer :: sshG(:,:,:)
#  endif
# endif
# ifdef ZCLM_NUDGING
          real(r8), pointer :: Znudgcof(:,:)
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
          real(r8), pointer :: zeta_ads(:,:)
          real(r8), pointer :: zeta_adsG(:,:,:)
# endif
# ifdef M2CLIMATOLOGY
          real(r8), pointer :: ubarclm(:,:)
          real(r8), pointer :: vbarclm(:,:)
#  ifndef ANA_M2CLIMA
          real(r8), pointer :: ubarclmG(:,:,:)
          real(r8), pointer :: vbarclmG(:,:,:)
#  endif
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
          real(r8), pointer :: ubar_ads(:,:)
          real(r8), pointer :: vbar_ads(:,:)
          real(r8), pointer :: ubar_adsG(:,:,:)
          real(r8), pointer :: vbar_adsG(:,:,:)
# endif
# ifdef M2CLM_NUDGING
          real(r8), pointer :: M2nudgcof(:,:)
# endif
# ifdef SOLVE3D
#  ifdef TCLIMATOLOGY
          real(r8), pointer :: tclm(:,:,:,:)
#   ifndef ANA_TCLIMA
          real(r8), pointer :: tclmG(:,:,:,:,:)
#   endif
#  endif
#  ifdef TCLM_NUDGING
          real(r8), pointer :: Tnudgcof(:,:,:)
#  endif
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
          real(r8), pointer :: t_ads(:,:,:,:)
          real(r8), pointer :: t_adsG(:,:,:,:,:)
#  endif
#  ifdef M3CLIMATOLOGY
          real(r8), pointer :: uclm(:,:,:)
          real(r8), pointer :: vclm(:,:,:)
#   ifndef ANA_M3CLIMA
          real(r8), pointer :: uclmG(:,:,:,:)
          real(r8), pointer :: vclmG(:,:,:,:)
#   endif
#  endif
#  ifdef M3CLM_NUDGING
          real(r8), pointer :: M3nudgcof(:,:)
#  endif
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
          real(r8), pointer :: u_ads(:,:,:)
          real(r8), pointer :: v_ads(:,:,:)
          real(r8), pointer :: u_adsG(:,:,:,:)
          real(r8), pointer :: v_adsG(:,:,:,:)
#  endif

# endif

        END TYPE T_CLIMA

        TYPE (T_CLIMA), allocatable :: CLIMA(:)

      CONTAINS

      SUBROUTINE allocate_clima (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1) allocate ( CLIMA(Ngrids) )
!
# ifdef ZCLIMATOLOGY
      allocate ( CLIMA(ng) % ssh(LBi:UBi,LBj:UBj) )
#  ifndef ANA_SSH
      allocate ( CLIMA(ng) % sshG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef ZCLM_NUDGING
      allocate ( CLIMA(ng) % Znudgcof(LBi:UBi,LBj:UBj) )
# endif

# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
      allocate ( CLIMA(ng) % zeta_ads(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % zeta_adsG(LBi:UBi,LBj:UBj,2) )
# endif

# ifdef M2CLIMATOLOGY
      allocate ( CLIMA(ng) % ubarclm(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % vbarclm(LBi:UBi,LBj:UBj) )
#  ifndef ANA_M2CLIMA
      allocate ( CLIMA(ng) % ubarclmG(LBi:UBi,LBj:UBj,2) )
      allocate ( CLIMA(ng) % vbarclmG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef M2CLM_NUDGING
      allocate ( CLIMA(ng) % M2nudgcof(LBi:UBi,LBj:UBj) )
# endif

# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
      allocate ( CLIMA(ng) % ubar_ads(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % vbar_ads(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % ubar_adsG(LBi:UBi,LBj:UBj,2) )
      allocate ( CLIMA(ng) % vbar_adsG(LBi:UBi,LBj:UBj,2) )
# endif

# ifdef SOLVE3D
#  ifdef TCLIMATOLOGY
      allocate ( CLIMA(ng) % tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
#   ifndef ANA_TCLIMA
      allocate ( CLIMA(ng) % tclmG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
#   endif
#  endif

#  ifdef TCLM_NUDGING
      allocate ( CLIMA(ng) % Tnudgcof(LBi:UBi,LBj:UBj,NT(ng)) )
#  endif

#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
      allocate ( CLIMA(ng) % t_ads(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( CLIMA(ng) % t_adsG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
#  endif

#  ifdef M3CLIMATOLOGY
      allocate ( CLIMA(ng) % uclm(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( CLIMA(ng) % vclm(LBi:UBi,LBj:UBj,N(ng)) )
#   ifndef ANA_M3CLIMA
      allocate ( CLIMA(ng) % uclmG(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( CLIMA(ng) % vclmG(LBi:UBi,LBj:UBj,N(ng),2) )
#   endif
#  endif

#  ifdef M3CLM_NUDGING
      allocate ( CLIMA(ng) % M3nudgcof(LBi:UBi,LBj:UBj) )
#  endif

#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
      allocate ( CLIMA(ng) % u_ads(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( CLIMA(ng) % v_ads(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( CLIMA(ng) % u_adsG(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( CLIMA(ng) % v_adsG(LBi:UBi,LBj:UBj,N(ng),2) )
#  endif

# endif

      RETURN
      END SUBROUTINE allocate_clima

      SUBROUTINE initialize_clima (ng, tile)
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
      DO j=Jmin,Jmax
        DO i=Imin,Imax
# ifdef ZCLIMATOLOGY
          CLIMA(ng) % ssh(i,j) = IniVal
#  ifndef ANA_SSH
          CLIMA(ng) % sshG(i,j,1) = IniVal
          CLIMA(ng) % sshG(i,j,2) = IniVal
#  endif
# endif
# ifdef ZCLM_NUDGING
          CLIMA(ng) % Znudgcof(i,j) = IniVal
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
          CLIMA(ng) % zeta_ads(i,j) = IniVal
          CLIMA(ng) % zeta_adsG(i,j,1) = IniVal
          CLIMA(ng) % zeta_adsG(i,j,2) = IniVal
# endif
# ifdef M2CLIMATOLOGY
          CLIMA(ng) % ubarclm(i,j) = IniVal
          CLIMA(ng) % vbarclm(i,j) = IniVal
#  ifndef ANA_M2CLIMA
          CLIMA(ng) % ubarclmG(i,j,1) = IniVal
          CLIMA(ng) % ubarclmG(i,j,2) = IniVal
          CLIMA(ng) % vbarclmG(i,j,1) = IniVal
          CLIMA(ng) % vbarclmG(i,j,2) = IniVal
#  endif
# endif
# ifdef M2CLM_NUDGING
          CLIMA(ng) % M2nudgcof(i,j) = IniVal
# endif
# if defined M3CLM_NUDGING && defined SOLVE3D
          CLIMA(ng) % M3nudgcof(i,j) = IniVal
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
          CLIMA(ng) % ubar_ads(i,j) = IniVal
          CLIMA(ng) % vbar_ads(i,j) = IniVal
          CLIMA(ng) % ubar_adsG(i,j,1) = IniVal
          CLIMA(ng) % ubar_adsG(i,j,2) = IniVal
          CLIMA(ng) % vbar_adsG(i,j,1) = IniVal
          CLIMA(ng) % vbar_adsG(i,j,2) = IniVal
# endif
        END DO
# ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=Imin,Imax
#  ifdef M3CLIMATOLOGY
            CLIMA(ng) % uclm(i,j,k) = IniVal
            CLIMA(ng) % vclm(i,j,k) = IniVal
#   ifndef ANA_M3CLIMA
            CLIMA(ng) % uclmG(i,j,k,1) = IniVal
            CLIMA(ng) % uclmG(i,j,k,2) = IniVal
            CLIMA(ng) % vclmG(i,j,k,1) = IniVal
            CLIMA(ng) % vclmG(i,j,k,2) = IniVal
#   endif
#  endif
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
            CLIMA(ng) % u_ads(i,j,k) = IniVal
            CLIMA(ng) % v_ads(i,j,k) = IniVal
            CLIMA(ng) % u_adsG(i,j,k,1) = IniVal
            CLIMA(ng) % u_adsG(i,j,k,2) = IniVal
            CLIMA(ng) % v_adsG(i,j,k,1) = IniVal
            CLIMA(ng) % v_adsG(i,j,k,2) = IniVal
#  endif
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
#  ifdef TCLIMATOLOGY
              CLIMA(ng) % tclm(i,j,k,itrc) = IniVal
#   ifndef ANA_TCLIMA
              CLIMA(ng) % tclmG(i,j,k,1,itrc) = IniVal
              CLIMA(ng) % tclmG(i,j,k,2,itrc) = IniVal
#   endif
#  endif
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI
              CLIMA(ng) % t_ads(i,j,k,itrc) = IniVal
              CLIMA(ng) % t_adsG(i,j,k,1,itrc) = IniVal
              CLIMA(ng) % t_adsG(i,j,k,2,itrc) = IniVal
#  endif
            END DO
          END DO
#  ifdef TCLM_NUDGING
          DO i=Imin,Imax
            CLIMA(ng) % Tnudgcof(i,j,itrc) = IniVal
          END DO
#  endif
        END DO
# endif
      END DO

      RETURN
      END SUBROUTINE initialize_clima
#endif
      END MODULE mod_clima
