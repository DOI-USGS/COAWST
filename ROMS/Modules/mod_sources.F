#include "cppdefs.h"
      MODULE mod_sources
#if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
!
!svn $Id: mod_sources.F 1338 2008-01-24 02:47:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Msrc       Maximum number of analytical point Sources/Sinks.        !
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

        integer, dimension(Ngrids) :: Msrc
        integer, dimension(Ngrids) :: Nsrc

        TYPE T_SOURCES

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
      USE mod_netcdf
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
!
# ifndef ANA_PSOURCE
      logical :: foundit

      integer :: Vid, ifile, nvatt, nvdim
# endif
      integer :: is, itrc, k, mg

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------

# if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                              defined TS_PSOURCE || defined Q_PSOURCE)
!
!  Inquire about the number of point sources.
!
      IF (ng.eq.1) THEN
        DO mg=1,Ngrids
          foundit=.FALSE.

          QUERY : DO ifile=1,nFfiles(mg)
            CALL netcdf_inq_var (ng, iNLM, FRCname(ifile,mg),           &
     &                           MyVarName = TRIM(Vname(1,idRxpo)),     &
     &                           SearchVar = foundit,                   &
     &                           VarID = Vid,                           &
     &                           nVardim = nvdim,                       &
     &                           nVarAtt = nvatt)
            IF (exit_flag.ne.NoError) RETURN

            IF (foundit) THEN
              Nsrc(mg)=var_Dsize(1)         ! first dimension
              Msrc(mg)=Nsrc(mg)
              EXIT QUERY
            END IF
          END DO QUERY
        END DO
      END IF
# else
!
!  Set number of point sources to maximum number of analytical sources.
!  Notice that a maximum of 200 analytical sources are set-up here.
!
      Msrc(ng)=200
      Nsrc(ng)=Msrc(ng)
# endif
!
!  Allocate structure.
!
      IF (ng.eq.1) allocate ( SOURCES(Ngrids) )
!
!  Allocate point Sources/Sinks variables.
!
      allocate ( SOURCES(ng) % Isrc(Nsrc(ng)) )

      allocate ( SOURCES(ng) % Jsrc(Nsrc(ng)) )

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
