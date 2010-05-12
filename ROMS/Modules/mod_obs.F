#include "cppdefs.h"
      MODULE mod_obs
#if defined ASSIMILATION || defined NUDGING
!
!svn $Id: mod_obs.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Assimilation of sea surface height data, SSH.                       !
!                                                                      !
!  EdatSSH         Latest two-time snapshots of SSH observation        !
!                    error variance data used for interpolation.       !
!  EmodSSH         SSH model error variance.                           !
!  EobsSSH         SSH observation error variance (non-dimensional).   !
!  SSHdat          Latest two-time snapshots of SSH observation        !
!                    data used for interpolation.                      !
!  SSHobs          Current SSH observations (m).                       !
!  tSSHobs         Time (days) of current SSH data snapshots.          !
!                                                                      !
!  Assimilation of sea surface temperature, SST.                       !
!                                                                      !
!  EdatSST         Latest two-time snapshots of SST observation        !
!                    error variance data used for interpolation.       !
!  EmodSST         SST model error variance.                           !
!  EobsSST         SST observation error variance (non-dimensional).   !
!  SSTdat          Latest two-time snapshots of SST observation        !
!                    data used for interpolation.                      !
!  SSTobs          Current SST observations (Celsius).                 !
!                                                                      !
!  Assimilation of tracers data.                                       !
!                                                                      !
!  EdatT           Latest two-time snapshots of tracers observation    !
!                    error variance data used for interpolation.       !
!  EmodT           Tracers model error variance.                       !
!  EobsT           Tracers observation error variance                  !
!                    (non-dimensional).                                !
!  Tdat            Latest two-time snapshots of tracers observation    !
!                    data used for interpolation.                      !
!  Tobs            Current Tracers observations (m).                   !
!                                                                      !
!  Horizontal currents observations for assimilation.                  !
!                                                                      !
!  EdatUV          Latest two-time snapshots of horizontal currents    !
!                    observations error variance data.                 !
!  EdatVsur        Latest two-time snapshots of surface currents       !
!                    observations error variance data.                 !
!  EmodU           U-velocity model error variance.                    !
!  EmodV           V-velocity model error variance.                    !
!  EobsUV          Horizontal currents observations error variance     !
!                    at RHO-points (non-dimensional).                  !
!  EobsVsur        Surface currents observations error variance        !
!                    error variance.                                   !
!  Udat            Latest two-time snapshots of U-velocity data.       !
!  Uobs            Current U-velocity observations.                    !
!  Usur            Current surface U-velocity observations.            !
!  Usurdat         Latest two-time snapshots of surface U-velocity     !
!                    data.                                             !
!  Vdat            Latest two-time snapshots of V-velocity data.       !
!  Vobs            Current V-velocity observations.                    !
!  Vsur            Current surface V-velocity observations.            !
!  Vsurdat         Latest two-time snapshots of surface U-velocity     !
!                    data.                                             !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_OBS

# ifdef ASSIMILATION_SSH
          real(r8), pointer :: EmodSSH(:,:)
          real(r8), pointer :: EobsSSH(:,:)
          real(r8), pointer :: SSHobs (:,:)
# endif

# if defined ASSIMILATION_SST || defined NUDGING_SST
          real(r8), pointer :: EobsSST(:,:)
          real(r8), pointer :: SSTobs (:,:)
#  ifdef ASSIMILATION_SST
          real(r8), pointer :: EmodSST(:,:)
#  endif
#  ifdef NUDGING_SST
          real(r8), pointer :: EdatSST(:,:,:)
          real(r8), pointer :: SSTdat (:,:,:)
#  endif
# endif

# if defined ASSIMILATION_T   || defined NUDGING_T   || \
     defined ASSIMILATION_SST || defined NUDGING_SST
          real(r8), pointer :: EobsT(:,:,:,:)
          real(r8), pointer :: Tobs (:,:,:,:)
#  if defined ASSIMILATION_T || defined ASSIMILATION_SST
          real(r8), pointer :: EmodT(:,:,:,:)
#  endif
#  ifdef NUDGING_T
          real(r8), pointer :: EdatT(:,:,:,:,:)
          real(r8), pointer :: Tdat (:,:,:,:,:)
#  endif
# endif

# if defined ASSIMILATION_UVsur || defined NUDGING_UVsur  || \
     defined ASSIMILATION_UV    || defined NUDGING_UV
          real(r8), pointer :: EobsUV(:,:,:)
          real(r8), pointer :: Uobs  (:,:,:)
          real(r8), pointer :: Vobs  (:,:,:)
#  if defined ASSIMILATION_UV || defined ASSIMILATION_UVsur
          real(r8), pointer :: EmodU (:,:,:)
          real(r8), pointer :: EmodV (:,:,:)
#  endif
#  ifdef NUDGING_UV
          real(r8), pointer :: Udat  (:,:,:,:)
          real(r8), pointer :: Vdat  (:,:,:,:)
          real(r8), pointer :: EdatUV(:,:,:,:)
#  endif
#  if defined ASSIMILATION_UVsur || defined NUDGING_UVsur
          real(r8), pointer :: Usur  (:,:)
          real(r8), pointer :: Vsur  (:,:)
          real(r8), pointer :: EobsVsur(:,:)
#  endif
#  ifdef NUDGING_UVsur
          real(r8), pointer :: Usurdat (:,:,:)
          real(r8), pointer :: Vsurdat (:,:,:)
          real(r8), pointer :: EdatVsur(:,:,:)
#  endif
# endif

        END TYPE T_OBS

        TYPE (T_OBS), allocatable :: OBS(:)

      CONTAINS

      SUBROUTINE allocate_obs (ng, LBi, UBi, LBj, UBj)
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
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OBS(Ngrids) )
!
# ifdef ASSIMILATION_SSH
      allocate ( OBS(ng) % EmodSSH(LBi:UBi,LBj:UBj) )
      allocate ( OBS(ng) % EobsSSH(LBi:UBi,LBj:UBj) )
      allocate ( OBS(ng) % SSHobs(LBi:UBi,LBj:UBj) )
# endif

# if defined ASSIMILATION_SST || defined NUDGING_SST
      allocate ( OBS(ng) % EobsSST(LBi:UBi,LBj:UBj) )
      allocate ( OBS(ng) % SSTobs(LBi:UBi,LBj:UBj) )
#  ifdef ASSIMILATION_SST
      allocate ( OBS(ng) % EmodSST(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef NUDGING_SST
      allocate ( OBS(ng) % EdatSST(LBi:UBi,LBj:UBj,2) )
      allocate ( OBS(ng) % SSTdat(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# if defined ASSIMILATION_T   || defined NUDGING_T   || \
     defined ASSIMILATION_SST || defined NUDGING_SST
      allocate ( OBS(ng) % EobsT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( OBS(ng) % Tobs(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
#  if defined ASSIMILATION_T || defined ASSIMILATION_SST
      allocate ( OBS(ng) % EmodT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
#  endif
#  ifdef NUDGING_T
      allocate ( OBS(ng) % EdatT(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
      allocate ( OBS(ng) % Tdat(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
#  endif
# endif

# if defined ASSIMILATION_UVsur || defined NUDGING_UVsur  || \
     defined ASSIMILATION_UV    || defined NUDGING_UV
      allocate ( OBS(ng) % EobsUV(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OBS(ng) % Uobs(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OBS(ng) % Vobs(LBi:UBi,LBj:UBj,N(ng)) )
#  if defined ASSIMILATION_UV || defined ASSIMILATION_UVsur
      allocate ( OBS(ng) % EmodU(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OBS(ng) % EmodV(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
#  ifdef NUDGING_UV
      allocate ( OBS(ng) % Udat(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OBS(ng) % Vdat(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OBS(ng) % EdatUV(LBi:UBi,LBj:UBj,N(ng),2) )
#  endif
#  if defined ASSIMILATION_UVsur || defined NUDGING_UVsur
      allocate ( OBS(ng) % Usur(LBi:UBi,LBj:UBj) )
      allocate ( OBS(ng) % Vsur(LBi:UBi,LBj:UBj) )
      allocate ( OBS(ng) % EobsVsur(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef NUDGING_UVsur
      allocate ( OBS(ng) % Usurdat(LBi:UBi,LBj:UBj,2) )
      allocate ( OBS(ng) % Vsurdat(LBi:UBi,LBj:UBj,2) )
      allocate ( OBS(ng) % EdatVsur(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_obs

      SUBROUTINE initialize_obs (ng, tile)
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
# ifdef ASSIMILATION_SSH
        DO i=Imin,Imax
          OBS(ng) % EmodSSH(i,j) = IniVal
          OBS(ng) % EobsSSH(i,j) = IniVal
          OBS(ng) % SSHobs(i,j) = IniVal
        END DO
# endif
# if defined ASSIMILATION_SST || defined NUDGING_SST
        DO i=Imin,Imax
          OBS(ng) % EobsSST(i,j) = IniVal
          OBS(ng) % SSTobs(i,j) = IniVal
#  ifdef ASSIMILATION_SST
          OBS(ng) % EmodSST(i,j) = IniVal
#  endif
#  ifdef NUDGING_SST
          OBS(ng) % EdatSST(i,j,1) = IniVal
          OBS(ng) % EdatSST(i,j,2) = IniVal
          OBS(ng) % SSTdat(i,j,1) = IniVal
          OBS(ng) % SSTdat(i,j,2) = IniVal
#  endif
        END DO
# endif
# if defined ASSIMILATION_T   || defined NUDGING_T   || \
     defined ASSIMILATION_SST || defined NUDGING_SST
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              OBS(ng) % EobsT(i,j,k,itrc) = IniVal
              OBS(ng) % Tobs(i,j,k,itrc) = IniVal
#  if defined ASSIMILATION_T || defined ASSIMILATION_SST
              OBS(ng) % EmodT(i,j,k,itrc) = IniVal
#  endif
#  ifdef NUDGING_T
              OBS(ng) % EdatT(i,j,k,1,itrc) = IniVal
              OBS(ng) % EdatT(i,j,k,2,itrc) = IniVal
              OBS(ng) % Tdat(i,j,k,1,itrc) = IniVal
              OBS(ng) % Tdat(i,j,k,2,itrc) = IniVal
#  endif
            END DO
          END DO
        END DO
# endif
# if defined ASSIMILATION_UVsur || defined NUDGING_UVsur  || \
     defined ASSIMILATION_UV    || defined NUDGING_UV
        DO k=1,N(ng)
          DO i=Imin,Imax
            OBS(ng) % EobsUV(i,j,k) = IniVal
            OBS(ng) % Uobs(i,j,k) = IniVal
            OBS(ng) % Vobs(i,j,k) = IniVal
#  if defined ASSIMILATION_UV || defined ASSIMILATION_UVsur
            OBS(ng) % EmodU(i,j,k) = IniVal
            OBS(ng) % EmodV(i,j,k) = IniVal
#  endif
#  ifdef NUDGING_UV
            OBS(ng) % Udat(i,j,k,1) = IniVal
            OBS(ng) % Udat(i,j,k,2) = IniVal
            OBS(ng) % Vdat(i,j,k,1) = IniVal
            OBS(ng) % Vdat(i,j,k,2) = IniVal
            OBS(ng) % EdatUV(i,j,k,1) = IniVal
            OBS(ng) % EdatUV(i,j,k,2) = IniVal
#  endif
          END DO
        END DO
#  if defined ASSIMILATION_UVsur || defined NUDGING_UVsur
        DO i=Imin,Imax
          OBS(ng) % Usur(i,j) = IniVal
          OBS(ng) % Vsur(i,j) = IniVal
          OBS(ng) % EobsVsur(i,j) = IniVal
        END DO
#  endif
#  ifdef NUDGING_UVsur
        DO i=Imin,Imax
          OBS(ng) % Usurdat(i,j,1) = IniVal
          OBS(ng) % Usurdat(i,j,2) = IniVal
          OBS(ng) % Vsurdat(i,j,1) = IniVal
          OBS(ng) % Vsurdat(i,j,2) = IniVal
          OBS(ng) % EdatVsur(i,j,1) = IniVal
          OBS(ng) % EdatVsur(i,j,2) = IniVal
        END DO
#  endif
# endif
      END DO

      RETURN
      END SUBROUTINE initialize_obs
#endif
      END MODULE mod_obs
