#include "cppdefs.h"
      MODULE mod_clima
!
!svn $Id: mod_clima.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_CLIMA
!
!  Climatology/Nudging arrays.
!
          real(r8), pointer :: ssh(:,:)
#ifndef ANA_SSH
          real(r8), pointer :: sshG(:,:,:)
#endif

          real(r8), pointer :: ubarclm(:,:)
          real(r8), pointer :: vbarclm(:,:)
#ifndef ANA_M2CLIMA
          real(r8), pointer :: ubarclmG(:,:,:)
          real(r8), pointer :: vbarclmG(:,:,:)
#endif

#ifdef SOLVE3D
          real(r8), pointer :: uclm(:,:,:)
          real(r8), pointer :: vclm(:,:,:)
# ifndef ANA_M3CLIMA
          real(r8), pointer :: uclmG(:,:,:,:)
          real(r8), pointer :: vclmG(:,:,:,:)
# endif

          real(r8), pointer :: tclm(:,:,:,:)
# ifndef ANA_TCLIMA
          real(r8), pointer :: tclmG(:,:,:,:,:)
# endif
#endif
!
!  Nudging coefficient arrays.
!
          real(r8), pointer :: M2nudgcof(:,:)
#ifdef SOLVE3D
          real(r8), pointer :: M3nudgcof(:,:,:)

          real(r8), pointer :: Tnudgcof(:,:,:,:)
#endif

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint-based algorithms arrays.
!
          real(r8), pointer :: zeta_ads(:,:)
          real(r8), pointer :: zeta_adsG(:,:,:)

          real(r8), pointer :: ubar_ads(:,:)
          real(r8), pointer :: vbar_ads(:,:)
          real(r8), pointer :: ubar_adsG(:,:,:)
          real(r8), pointer :: vbar_adsG(:,:,:)

# ifdef SOLVE3D
          real(r8), pointer :: u_ads(:,:,:)
          real(r8), pointer :: v_ads(:,:,:)
          real(r8), pointer :: u_adsG(:,:,:,:)
          real(r8), pointer :: v_adsG(:,:,:,:)
          real(r8), pointer :: wvel_ads(:,:,:)
          real(r8), pointer :: wvel_adsG(:,:,:,:)

          real(r8), pointer :: t_ads(:,:,:,:)
          real(r8), pointer :: t_adsG(:,:,:,:,:)
# endif
#endif

#ifdef AKTCLIMATOLOGY
          real(r8), pointer :: Aclm(:,:,:)
          real(r8), pointer :: AclmG(:,:,:,:)
#endif
#ifdef OCLIMATOLOGY
          real(r8), pointer :: oclm(:,:,:)
# ifndef ANA_OCLIMA
          real(r8), pointer :: oclmG(:,:,:,:)
# endif
#endif
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
      USE mod_scalars
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( CLIMA(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Climatology/Nudging arrays.
!
      IF (LsshCLM(ng)) THEN
        allocate ( CLIMA(ng) % ssh(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d

#ifndef ANA_SSH
        allocate ( CLIMA(ng) % sshG(LBi:UBi,LBj:UBj,2) )
        Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#endif
      END IF
!
      IF (Lm2CLM(ng)) THEN
        allocate ( CLIMA(ng) % ubarclm(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d

        allocate ( CLIMA(ng) % vbarclm(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d

#ifndef ANA_M2CLIMA
        allocate ( CLIMA(ng) % ubarclmG(LBi:UBi,LBj:UBj,2) )
        Dmem(ng)=Dmem(ng)+2.0_r8*size2d

        allocate ( CLIMA(ng) % vbarclmG(LBi:UBi,LBj:UBj,2) )
        Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#endif
      END IF

#ifdef SOLVE3D
!
      IF (Lm3CLM(ng)) THEN
        allocate ( CLIMA(ng) % uclm(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

        allocate ( CLIMA(ng) % vclm(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

# ifndef ANA_M3CLIMA
        allocate ( CLIMA(ng) % uclmG(LBi:UBi,LBj:UBj,N(ng),2) )
        Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng),r8)*size2d

        allocate ( CLIMA(ng) % vclmG(LBi:UBi,LBj:UBj,N(ng),2) )
        Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng),r8)*size2d
# endif
      END IF
!
      IF (ANY(LtracerCLM(:,ng)).or.ANY(LnudgeTCLM(:,ng))) THEN
        allocate ( CLIMA(ng) % tclm(LBi:UBi,LBj:UBj,N(ng),NTCLM(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NTCLM(ng),r8)*size2d

# ifndef ANA_TCLIMA
        allocate ( CLIMA(ng) % tclmG(LBi:UBi,LBj:UBj,N(ng),2,           &
     &                               NTCLM(ng)) )
        Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng)*NTCLM(ng),r8)*size2d
# endif
      END IF
#endif
!
!  Nudging coefficient arrays.
!
      IF (LnudgeM2CLM(ng)) THEN
        allocate ( CLIMA(ng) % M2nudgcof(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

#ifdef SOLVE3D
      IF (LnudgeM3CLM(ng)) THEN
        allocate ( CLIMA(ng) % M3nudgcof(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (ANY(LnudgeTCLM(:,ng))) THEN
        allocate ( CLIMA(ng) % Tnudgcof(LBi:UBi,LBj:UBj,N(ng),          &
     &                                  NTCLM(ng)) )
        Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng)*NTCLM(ng),r8)*size2d
      END IF
#endif

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint-based algorithms arrays.
!
      allocate ( CLIMA(ng) % zeta_ads(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( CLIMA(ng) % zeta_adsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( CLIMA(ng) % ubar_ads(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( CLIMA(ng) % vbar_ads(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( CLIMA(ng) % ubar_adsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( CLIMA(ng) % vbar_adsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

# ifdef SOLVE3D
      allocate ( CLIMA(ng) % u_ads(LBi:UBi,LBj:UBj,N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

      allocate ( CLIMA(ng) % v_ads(LBi:UBi,LBj:UBj,N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d

      allocate ( CLIMA(ng) % u_adsG(LBi:UBi,LBj:UBj,N(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng),r8)*size2d

      allocate ( CLIMA(ng) % v_adsG(LBi:UBi,LBj:UBj,N(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng),r8)*size2d

      allocate ( CLIMA(ng) % wvel_ads(LBi:UBi,LBj:UBj,0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d

      allocate ( CLIMA(ng) % wvel_adsG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng)+1,r8)*size2d

      allocate ( CLIMA(ng) % t_ads(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d

      allocate ( CLIMA(ng) % t_adsG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng)*NT(ng),r8)*size2d
# endif
#endif

#if defined AKTCLIMATOLOGY && defined OFFLINE_FLOATS
      allocate ( CLIMA(ng) % Aclm(LBi:UBi,LBj:UBj,0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
# ifndef ANA_AKTCLIMA
      allocate ( CLIMA(ng) % AclmG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng)+1,r8)*size2d
# endif
#endif
#ifdef OCLIMATOLOGY
      allocate ( CLIMA(ng) % oclm(LBi:UBi,LBj:UBj,0:N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
# ifndef ANA_OCLIMA
      allocate ( CLIMA(ng) % oclmG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(N(ng)+1,r8)*size2d
# endif
#endif
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
#ifdef SOLVE3D
      integer :: itrc, k
#endif

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#else
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
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Climatology/Nudging arrays.
!
      IF (LsshCLM(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            CLIMA(ng) % ssh(i,j) = IniVal
#ifndef ANA_SSH
            CLIMA(ng) % sshG(i,j,1) = IniVal
            CLIMA(ng) % sshG(i,j,2) = IniVal
#endif
          END DO
        END DO
      END IF
!
      IF (Lm2CLM(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            CLIMA(ng) % ubarclm(i,j) = IniVal
            CLIMA(ng) % vbarclm(i,j) = IniVal
#ifndef ANA_M2CLIMA
            CLIMA(ng) % ubarclmG(i,j,1) = IniVal
            CLIMA(ng) % ubarclmG(i,j,2) = IniVal
            CLIMA(ng) % vbarclmG(i,j,1) = IniVal
            CLIMA(ng) % vbarclmG(i,j,2) = IniVal
#endif
          END DO
        END DO
      END IF

#ifdef SOLVE3D
!
      IF (Lm3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              CLIMA(ng) % uclm(i,j,k) = IniVal
              CLIMA(ng) % vclm(i,j,k) = IniVal
# ifndef ANA_M3CLIMA
              CLIMA(ng) % uclmG(i,j,k,1) = IniVal
              CLIMA(ng) % uclmG(i,j,k,2) = IniVal
              CLIMA(ng) % vclmG(i,j,k,1) = IniVal
              CLIMA(ng) % vclmG(i,j,k,2) = IniVal
# endif
            END DO
          END DO
        END DO
      END IF
!
      IF (ANY(LtracerCLM(:,ng)).or.ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                CLIMA(ng) % tclm(i,j,k,itrc) = IniVal
# ifndef ANA_TCLIMA
                CLIMA(ng) % tclmG(i,j,k,1,itrc) = IniVal
                CLIMA(ng) % tclmG(i,j,k,2,itrc) = IniVal
# endif
              END DO
            END DO
          END DO
        END DO
      END IF
#endif
!
!  Nudging coefficient arrays.
!
      IF (LnudgeM2CLM(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            CLIMA(ng) % M2nudgcof(i,j) = IniVal
          END DO
        END DO
      END IF

#ifdef SOLVE3D
!
      IF (LnudgeM3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              CLIMA(ng) % M3nudgcof(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                CLIMA(ng) % Tnudgcof(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
#endif

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint-based algorithms arrays.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          CLIMA(ng) % zeta_ads(i,j) = IniVal
          CLIMA(ng) % zeta_adsG(i,j,1) = IniVal
          CLIMA(ng) % zeta_adsG(i,j,2) = IniVal
!
          CLIMA(ng) % ubar_ads(i,j) = IniVal
          CLIMA(ng) % vbar_ads(i,j) = IniVal
          CLIMA(ng) % ubar_adsG(i,j,1) = IniVal
          CLIMA(ng) % ubar_adsG(i,j,2) = IniVal
          CLIMA(ng) % vbar_adsG(i,j,1) = IniVal
          CLIMA(ng) % vbar_adsG(i,j,2) = IniVal
        END DO

# ifdef SOLVE3D
!
        DO k=1,N(ng)
          DO i=Imin,Imax
            CLIMA(ng) % u_ads(i,j,k) = IniVal
            CLIMA(ng) % v_ads(i,j,k) = IniVal
            CLIMA(ng) % u_adsG(i,j,k,1) = IniVal
            CLIMA(ng) % u_adsG(i,j,k,2) = IniVal
            CLIMA(ng) % v_adsG(i,j,k,1) = IniVal
            CLIMA(ng) % v_adsG(i,j,k,2) = IniVal
          END DO
        END DO
!
        DO k=0,N(ng)
          DO i=Imin,Imax
            CLIMA(ng) % wvel_ads(i,j,k) = IniVal
            CLIMA(ng) % wvel_adsG(i,j,k,1) = IniVal
            CLIMA(ng) % wvel_adsG(i,j,k,2) = IniVal
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              CLIMA(ng) % t_ads(i,j,k,itrc) = IniVal
              CLIMA(ng) % t_adsG(i,j,k,1,itrc) = IniVal
              CLIMA(ng) % t_adsG(i,j,k,2,itrc) = IniVal
            END DO
          END DO
        END DO
# endif
      END DO
#endif

#ifdef AKTCLIMATOLOGY
      DO j=Jmin,Jmax
        DO k=0,N(ng)
          DO i=Imin,Imax
            CLIMA(ng) % Aclm(i,j,k) = IniVal
# ifndef ANA_AKTCLIMA
            CLIMA(ng) % AclmG(i,j,k,1) = IniVal
            CLIMA(ng) % AclmG(i,j,k,2) = IniVal
# endif
          END DO
        END DO
      END DO
#endif
#ifdef OCLIMATOLOGY
      DO j=Jmin,Jmax
        DO k=0,N(ng)
          DO i=Imin,Imax
            CLIMA(ng) % oclm(i,j,k) = IniVal
# ifndef ANA_OCLIMA
            CLIMA(ng) % oclmG(i,j,k,1) = IniVal
            CLIMA(ng) % oclmG(i,j,k,2) = IniVal
# endif
          END DO
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE initialize_clima
      END MODULE mod_clima
