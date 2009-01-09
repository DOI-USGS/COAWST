#include "cppdefs.h"
      MODULE mod_sources
#if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
!
!svn $Id: mod_sources.F 1338 2008-01-24 02:47:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
# ifdef ANA_PSOURCE
!  Msrc       Maximum number of analytical point Sources/Sinks.        !
# endif
!  Nsrc       Number of point Sources/Sinks.                           !
!  Dsrc       Direction of point Sources/Sinks:                        !
!               Dsrc(:) = 0,  Along XI-direction.                      !
!               Dsrc(:) > 0,  Along ETA-direction.                     !
!  Fsrc       Point Source/Sinks identification flag:                  !
!               Fsrc(:) = 0,  All Tracer source/sink are off.          !
!               Fsrc(:) = 1,  Only temperature is on.                  !
!               Fsrc(:) = 2,  Only salinity is on.                     !
!               Fsrc(:) = 3,  Both temperature and salinity are on.    !
!               Fsrc(:) = 4,  Both nitrate and salinity are on.        !
!               Fsrc(:) = ... And other combinations.                  !
!                             (We need a more robust logic here)       !
!  Isrc       I-grid location of point Sources/Sinks,                  !
!               1 =< Isrc =< Lm(ng).                                   !
!  Jsrc       J-grid location of point Sources/Sinks,                  !
!               1 =< Jsrc =< Mm(ng).                                   !
!  Lsrc       Logical switch for each source point source/sink datum   !
!               indicating which tracer(s) is (are) activated          !
!               according to the "Fsrc" indentification flag.          !
!  Ltracer    Logical switch indicating which tracer field need to be  !
!               processed for Sources/Sinks terms.                     !
!  Qbar       Vertically integrated mass transport (m3/s) of point     !
!               Sources/Sinks at U- or V-points:                       !
!               Qbar -> positive, if the mass transport is in the      !
!                       positive U- or V-direction.                    !
!               Qbar -> negative, if the mass transport is in the      !
!                       negative U- or V-direction.                    !
!  QbarG      Latest two-time snapshots of vertically integrated       !
!               mass transport (m3/s) of point Sources/Sinks.          !
!  Qshape     Nondimensional shape function to distribute mass         !
!               mass point Sources/Sinks vertically.                   !
!  Qsrc       Mass transport profile (m3/s) of point Sources/Sinks.    !
!  Tsrc       Tracer (tracer units) point Sources/Sinks.               !
!  TsrcG      Latest two-time snapshots of tracer (tracer units)       !
!               point Sources/Sinks.                                   !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        USE mod_param

        implicit none

        integer, dimension(Ngrids) :: Nsrc
# ifdef ANA_PSOURCE
        integer, parameter :: Msrc = 200
# endif

        TYPE T_SOURCES

          logical, pointer :: Lsrc(:,:)
          logical, pointer :: Ltracer(:)

          integer, pointer :: Isrc(:)
          integer, pointer :: Jsrc(:)

          real(r8), pointer :: Dsrc(:)
          real(r8), pointer :: Fsrc(:)
          real(r8), pointer :: Qbar(:)
          real(r8), pointer :: Qshape(:,:)
          real(r8), pointer :: Qsrc(:,:)
          real(r8), pointer :: Tsrc(:,:,:)
          real(r8), pointer :: Xsrc(:)
          real(r8), pointer :: Ysrc(:)

# ifndef ANA_PSOURCE
          real(r8), pointer :: QbarG(:,:)
          real(r8), pointer :: TsrcG(:,:,:,:)
# endif

        END TYPE T_SOURCES

        TYPE (T_SOURCES), allocatable :: SOURCES(:)

      CONTAINS

      SUBROUTINE allocate_sources (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in the module  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifndef ANA_PSOURCE
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#  endif
# endif
!
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
! 
# ifndef ANA_PSOURCE
      integer :: Nrec, nvd, Vsize(4)
# endif
      integer :: is, itrc, k

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( SOURCES(Ngrids) )

# if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                              defined TS_PSOURCE || defined Q_PSOURCE)
!
!  Inquire about the number of point sources.
!
      IF (InpThread) THEN
        CALL opencdf (ng, nFfiles(ng), FRCname(1,ng), ncfile, N(ng),    &
     &                idRxpo, Nrec, nvd, Vsize)
        Nsrc(ng)=Vsize(1)
      END IF
#  ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iNLM, Nsrc(ng), 1)
#  endif
# else
!
!  Set number of point sources to maximum number of analytical sources.
!
      Nsrc(ng)=Msrc
# endif
!
!  Allocate point Sources/Sinks variables.
!
      allocate ( SOURCES(ng) % Isrc(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Jsrc(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Lsrc(Nsrc(ng),NT(ng)) )

      allocate ( SOURCES(ng) % Ltracer(NT(ng)) )

      allocate ( SOURCES(ng) % Dsrc(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Fsrc(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Qbar(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Qshape(Nsrc(ng),N(ng)) )

      allocate ( SOURCES(ng) % Qsrc(Nsrc(ng),N(ng)) )

      allocate ( SOURCES(ng) % Tsrc(Nsrc(ng),N(ng),NT(ng)) )

      allocate ( SOURCES(ng) % Xsrc(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Ysrc(Nsrc(ng)) )

# ifndef ANA_PSOURCE
      allocate ( SOURCES(ng) % QbarG(Nsrc(ng),2) )

      allocate ( SOURCES(ng) % TsrcG(Nsrc(ng),N(ng),2,NT(ng)) )
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO is=1,Nsrc(ng)
        SOURCES(ng) % Isrc(is) = 0
        SOURCES(ng) % Jsrc(is) = 0
        SOURCES(ng) % Dsrc(is) = IniVal
        SOURCES(ng) % Fsrc(is) = IniVal
        SOURCES(ng) % Xsrc(is) = IniVal
        SOURCES(ng) % Ysrc(is) = IniVal
        SOURCES(ng) % Qbar(is) = IniVal
# ifndef ANA_PSOURCE
        SOURCES(ng) % QbarG(is,1) = IniVal
        SOURCES(ng) % QbarG(is,2) = IniVal
# endif
      END DO
      DO k=1,N(ng)
        DO is=1,Nsrc(ng)
          SOURCES(ng) % Qshape(is,k) = IniVal
          SOURCES(ng) % Qsrc(is,k) = IniVal
        END DO
      END DO 
      DO itrc=1,NT(ng)
        SOURCES(ng) % Ltracer(itrc) = .FALSE.
        DO is=1,Nsrc(ng)
          SOURCES(ng) % Lsrc(is,itrc) = .FALSE.
        END DO
        DO k=1,N(ng)
          DO is=1,Nsrc(ng)
            SOURCES(ng) % Tsrc(is,k,itrc) = IniVal
# ifndef ANA_PSOURCE
            SOURCES(ng) % TsrcG(is,k,1,itrc) = IniVal
            SOURCES(ng) % TsrcG(is,k,2,itrc) = IniVal
# endif
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_sources
#endif
      END MODULE mod_sources
