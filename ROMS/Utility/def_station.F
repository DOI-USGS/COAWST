#include "cppdefs.h"
      MODULE def_station_mod
#ifdef STATIONS
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module creates output STATION data using either the standard   !
!  NetCDF library or the Parallel-IO (PIO) library.  It defines its    !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BIOLOGY
      USE mod_biology
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
!
# if defined BBL_MODEL || defined WAVES_OUTPUT
      USE bbl_output_mod,      ONLY : bbl_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,      ONLY : bbl_def_station_pio
#  endif
# endif
      USE def_dim_mod,         ONLY : def_dim
      USE def_info_mod,        ONLY : def_info
      USE def_var_mod,         ONLY : def_var
# ifdef ICE_MODEL
      USE ice_output_mod,      ONLY : ice_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,      ONLY : ice_def_station_pio
#  endif
# endif
# ifdef SEDIMENT
      USE sediment_output_mod, ONLY : sediment_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod, ONLY : sediment_def_station_pio
#  endif
# endif
      USE strings_mod,         ONLY : FoundError
# if defined WEC || defined WEC_VF
      USE wec_output_mod,      ONLY : wec_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,      ONLY : wec_def_station_pio
#  endif
# endif
      USE wrt_info_mod,        ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_station
      PRIVATE :: def_station_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_station_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_station (ng, ldef)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Create a new history file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (STA(ng)%IOtype)
        CASE (io_nf90)
          CALL def_station_nf90 (ng, iNLM, ldef)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_station_pio (ng, iNLM, ldef)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) STA(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_STATION - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_station
!
!***********************************************************************
      SUBROUTINE def_station_nf90 (ng, model, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, recdim, stadim
      integer :: status

      integer :: DimIDs(nDimID), pgrd(2)

# ifdef SOLVE3D
      integer :: itrc

      integer :: bgrd(3), rgrd(3), wgrd(3)
# endif
# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
      integer :: igrd(3), cgrd(3)
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=STA(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new station data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), STA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
# ifdef SOLVE3D
        status=def_dim(ng, model, STA(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined INWAVE_MODEL
        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'ND',           &
                       ND, DimIDs(37))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'NDp1',         &
                       ND+1, DimIDs(38))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SEDIMENT
        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'NST',          &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef ECOSIM
        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nbands',       &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'station' ,     &
     &                 Nstation(ng), DimIDs(13))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef MULTI_SCALE_B
        status=def_dim(ng, model, STA(ng)%ncid, ncname, 'Nscale',       &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, model, STA(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
        stadim=DimIDs(13)
!
!  Define dimension vector for point variables.
!
        pgrd(1)=DimIDs(13)
        pgrd(2)=DimIDs(12)
# ifdef SOLVE3D
!
!  Define dimension vector for cast variables at vertical RHO-points.
!
        rgrd(1)=DimIDs( 9)
        rgrd(2)=DimIDs(13)
        rgrd(3)=DimIDs(12)
!
!  Define dimension vector for cast variables at vertical W-points.
!
        wgrd(1)=DimIDs(10)
        wgrd(2)=DimIDs(13)
        wgrd(3)=DimIDs(12)
!
!  Define dimension vector for sediment bed layer type variables.
!
        bgrd(1)=DimIDs(16)
        bgrd(2)=DimIDs(13)
        bgrd(3)=DimIDs(12)
# endif
#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Define dimension vector for energy wave layer type variables.
!
        igrd(1)=DimIDs(36)
        igrd(2)=DimIDs(13)
        igrd(3)=DimIDs(12)

        cgrd(1)=DimIDs(37)
        cgrd(2)=DimIDs(13)
        cgrd(3)=DimIDs(12)
#endif
!
!  Initialize unlimited time record dimension.
!
        STA(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, STA(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idtime),    &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface.
!
        IF (Sout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idFsur),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
# ifdef WET_DRY
     &                   SetFillVal = .FALSE.,                          &
# else
     &                   SetFillVal = .TRUE.,                           &
# endif
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Sout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUbar),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Sout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVbar),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward/Northward momentum component at RHO-points.
!
        IF (Sout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idu2dE),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Sout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idv2dN),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Sout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUvel),  &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Sout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVvel),  &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward/Northward momentum component at RHO-points.
!
        IF (Sout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idu3dE),  &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Sout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idv3dN),  &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define "true" vertical velocity.
!
        IF (Sout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idWvel),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" velocity (m3/s).
!
        IF (Sout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter3 second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idOvel),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Sout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
            status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Tid(itrc),  &
     &                     NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,       &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Sout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idDano),  &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef NEMURO_SED1
!
!  Define PON in sediment.
!
        IF (Sout(idPONsed,ng)) THEN
          Vinfo( 1)=Vname(1,idPONsed)
          Vinfo( 2)=Vname(2,idPONsed)
          Vinfo( 3)=Vname(3,idPONsed)
          Vinfo(14)=Vname(4,idPONsed)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idPONsed), &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define OPAL in sediment.
!
        IF (Sout(idOPALsed,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALsed)
          Vinfo( 2)=Vname(2,idOPALsed)
          Vinfo( 3)=Vname(3,idOPALsed)
          Vinfo(14)=Vname(4,idOPALsed)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idOPALsed),&
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define DENIT in sediment.
!
        IF (Sout(idDENITsed,ng)) THEN
          Vinfo( 1)=Vname(1,idDENITsed)
          Vinfo( 2)=Vname(2,idDENITsed)
          Vinfo( 3)=Vname(3,idDENITsed)
          Vinfo(14)=Vname(4,idDENITsed)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid,                        &
     &                   STA(ng)%Vid(idDENITsed),                       &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define buried PON in sediment.
!
        IF (Sout(idPONbur,ng)) THEN
          Vinfo( 1)=Vname(1,idPONbur)
          Vinfo( 2)=Vname(2,idPONbur)
          Vinfo( 3)=Vname(3,idPONbur)
          Vinfo(14)=Vname(4,idPONbur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idPONbur), &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define buried OPAL in sediment.
!
        IF (Sout(idOPALbur,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALbur)
          Vinfo( 2)=Vname(2,idOPALbur)
          Vinfo( 3)=Vname(3,idOPALbur)
          Vinfo(14)=Vname(4,idOPALbur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idOPALbur),&
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef PRIMARY_PROD
!
!  Define Net primary productivity.
!
        IF (Sout(idNPP,ng)) THEN
          Vinfo( 1)=Vname(1,idNPP)
          Vinfo( 2)=Vname(2,idNPP)
          Vinfo( 3)=Vname(3,idNPP)
          Vinfo(14)=Vname(4,idNPP)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idNPP),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Sout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idHsbl),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Sout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idHbbl),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define vertical viscosity coefficient.
!
        IF (Sout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVvis),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Sout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idTdif),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Sout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idSdif),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Sout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idMtke),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Sout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idMtls),  &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Sout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idPair),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Sout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUair),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Sout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVair),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface Eastward/Nortward winds.
!
        IF (Sout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          Vinfo( 2)=Vname(2,idUaiE)
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUaiE),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Sout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          Vinfo( 2)=Vname(2,idVaiN)
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVaiN),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define surface net heat flux.
!
        IF (Sout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          Vinfo( 2)=Vname(2,idTsur(itemp))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid,                       &
     &                   STA(ng)%Vid(idTsur(itemp)), NF_FOUT,           &
     &                   2, pgrd, Aval, Vinfo, ncname,                  &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SALINITY
!
!  Define surface net salt flux.
!
        IF (Sout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          Vinfo( 2)=Vname(2,idTsur(isalt))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid,                       &
     &                   STA(ng)%Vid(idTsur(isalt)), NF_FOUT,           &
     &                   2, pgrd, Aval, Vinfo, ncname,                  &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Sout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idLhea),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Sout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idShea),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Sout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idLrad),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Sout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idSrad),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined EMINUSP && defined BULK_FLUXES
!
!  Define E-P flux (as computed by bulk_flux.F).
!
        IF (Sout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idEmPf),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Sout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idevap),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Sout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idrain),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
!
!  Define surface U-momentum stress.
!
        IF (Sout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUsms),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Sout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVsms),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Sout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUbms),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Sout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVbms),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
        IF (Sout(idRwet,ng)) THEN
          Vinfo( 1)=Vname(1,idRwet)
          Vinfo( 2)=Vname(2,idRwet)
          Vinfo( 3)=Vname(3,idRwet)
          Vinfo(14)=Vname(4,idRwet)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idRwet),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out wet/dry mask at U-points.
!
        IF (Sout(idUwet,ng)) THEN
          Vinfo( 1)=Vname(1,idUwet)
          Vinfo( 2)=Vname(2,idUwet)
          Vinfo( 3)=Vname(3,idUwet)
          Vinfo(14)=Vname(4,idUwet)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idUwet),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out wet/dry mask at V-points.
!
        IF (Sout(idUwet,ng)) THEN
          Vinfo( 1)=Vname(1,idVwet)
          Vinfo( 2)=Vname(2,idVwet)
          Vinfo( 3)=Vname(3,idVwet)
          Vinfo(14)=Vname(4,idVwet)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, STA(ng)%ncid, STA(ng)%Vid(idVwet),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# if defined BBL_MODEL || defined WAVES_OUTPUT
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_station_nf90 (ng, model, ldef, Sout, STA,          &
     &                             pgrd, rgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Define the ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_station_nf90 (ng, model, ldef, Sout, STA,          &
     &                             pgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_station_nf90 (ng, model, ldef, Sout, STA,     &
     &                                  bgrd, pgrd, rgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Define the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
        CALL wec_def_station_nf90 (ng, model, ldef, Sout, STA,          &
     &                             pgrd, rgrd, igrd, cgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, STA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, STA(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=STA(ng)%name
!
!  Open stations file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, STA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = STA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = STA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  stations variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            STA(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            STA(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            STA(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            STA(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            STA(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            STA(ng)%Vid(idv2dN)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            STA(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            STA(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            STA(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            STA(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            STA(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            STA(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            STA(ng)%Vid(idDano)=var_id(i)
#  ifdef NEMURO_SED1
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONsed))) THEN
            got_var(idPONsed)=.TRUE.
            STA(ng)%Vid(idPONsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALsed))) THEN
            got_var(idOPALsed)=.TRUE.
            STA(ng)%Vid(idOPALsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDENITsed))) THEN
            got_var(idDENITsed)=.TRUE.
            STA(ng)%Vid(idDENITsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONbur))) THEN
            got_var(idPONbur)=.TRUE.
            STA(ng)%Vid(idPONbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALbur))) THEN
            got_var(idOPALbur)=.TRUE.
            STA(ng)%Vid(idOPALbur)=var_id(i)
#  endif
#  ifdef PRIMARY_PROD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNPP))) THEN
            got_var(idNPP)=.TRUE.
            STA(ng)%Vid(idNPP)=var_id(i)
#  endif
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            STA(ng)%Vid(idHsbl)=var_id(i)
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            STA(ng)%Vid(idHbbl)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            STA(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            STA(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            STA(ng)%Vid(idSdif)=var_id(i)
#  if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            STA(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            STA(ng)%Vid(idMtls)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            STA(ng)%Vid(idPair)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            STA(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            STA(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            STA(ng)%Vid(idUaiE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            STA(ng)%Vid(idVaiN)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            STA(ng)%Vid(idTsur(itemp))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            STA(ng)%Vid(idTsur(isalt))=var_id(i)
#  ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            STA(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            STA(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            STA(ng)%Vid(idLrad)=var_id(i)
#  endif
#  ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            STA(ng)%Vid(idSrad)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            STA(ng)%Vid(idEmPf)=var_id(i)
#  if defined EMINUSP && defined BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            STA(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            STA(ng)%Vid(idrain)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            STA(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            STA(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            STA(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            STA(ng)%Vid(idVbms)=var_id(i)
# ifdef WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            STA(ng)%Vid(idRwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            STA(ng)%Vid(idUwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            STA(ng)%Vid(idVwet)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              STA(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Sout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Sout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Sout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Sout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Sout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Sout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Sout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Sout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Sout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Sout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Sout(idOvel,ng)) THEN
          IF (Master) WRITE(stdout,60) TRIM(Vname(1,idOvel)),           &
     &                                 TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Sout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef NEMURO_SED1
        IF (.not.got_var(idPONsed).and.Sout(idPONsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPONsed)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALsed).and.Sout(idOPALsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOPALsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDENITsed).and.Sout(idDENITsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDENITsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPONbur).and.Sout(idPONbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPONbur)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALbur).and.Sout(idOPALbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOPALbur)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (.not.got_var(idNPP).and.Sout(idNPP,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idNPP)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Sout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Sout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idVvis).and.Sout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Sout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SALINITY
        IF (.not.got_var(idSdif).and.Sout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Sout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Sout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Sout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Sout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Sout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Sout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUaiE)) ,         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Sout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idTsur(itemp)).and.Sout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Sout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Sout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Sout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Sout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Sout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idEmPf).and.Sout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined EMINUSP && defined BULK_FLUXES
        IF (.not.got_var(idevap).and.Sout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Sout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Sout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Sout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Sout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Sout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef WET_DRY
        IF (.not.got_var(idRwet).and.Sout(idRwet,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet).and.Sout(idUwet,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet).and.Sout(idVwet,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Sout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif

# if defined BBL_MODEL || defined WAVES_OUTPUT
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for station variables. Get variable IDs.
!
        CALL bbl_def_station_nf90 (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# endif

# ifdef ICE_MODEL
!
!  Scan ice model variables from input NetCDF and activate switches for
!  station variables. Get variable IDs.
!
        CALL ice_def_station_nf90 (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SEDIMENT
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for station variables. Get variable IDs.
!
        CALL sediment_def_station_nf90 (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC || defined WEC_VF
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for station variables. Get variable IDs.
!
        CALL wec_def_station_nf90 (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      STA(ng)%Rindex=rec_size
        STA(ng)%Rindex=(ntstart(ng)-1)/nSTA(ng)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_STATION_NF90 - creating stations file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_STATION_NF90 - inquiring stations file,',t56,     &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_STATION_NF90 - unable to create stations NetCDF', &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_STATION_NF90 - unable to open stations NetCDF',   &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_STATION_NF90 - unable to find variable: ',a,2x,   &
     &        ' in stations NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_station_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_station_pio (ng, model, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, recdim, stadim
      integer :: status

      integer :: DimIDs(nDimID), pgrd(2)

#  ifdef SOLVE3D
      integer :: itrc

      integer :: bgrd(3), rgrd(3), wgrd(3)
#  endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=STA(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new station data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname), STA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
#  ifdef SOLVE3D
        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 's_rho',     &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 's_w',       &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'tracer',    &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SEDIMENT
        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'NST',       &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nbed',      &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef ECOSIM
        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nbands',    &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nphy',      &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nbac',      &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Ndom',      &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nfec',      &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'station' ,  &
     &                 Nstation(ng), DimIDs(13))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'boundary',  &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FOUR_DVAR
        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nstate',    &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef MULTI_SCALE_B
        status=def_dim(ng, model, STA(ng)%pioFile, ncname, 'Nscale',    &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, model, STA(ng)%pioFile, ncname,              &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
        stadim=DimIDs(13)
!
!  Define dimension vector for point variables.
!
        pgrd(1)=DimIDs(13)
        pgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
!
!  Define dimension vector for cast variables at vertical RHO-points.
!
        rgrd(1)=DimIDs( 9)
        rgrd(2)=DimIDs(13)
        rgrd(3)=DimIDs(12)
!
!  Define dimension vector for cast variables at vertical W-points.
!
        wgrd(1)=DimIDs(10)
        wgrd(2)=DimIDs(13)
        wgrd(3)=DimIDs(12)
!
!  Define dimension vector for sediment bed layer type variables.
!
        bgrd(1)=DimIDs(16)
        bgrd(2)=DimIDs(13)
        bgrd(3)=DimIDs(12)
#  endif
!
!  Initialize unlimited time record dimension.
!
        STA(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, STA(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        STA(ng)%pioVar(idtime)%dkind=PIO_TOUT
        STA(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, model, STA(ng)%pioFile,                      &
     &                 STA(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface.
!
        IF (Sout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idFsur)%dkind=PIO_FOUT
          STA(ng)%pioVar(idFsur)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idFsur)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
#  ifdef WET_DRY
     &                   SetFillVal = .FALSE.,                          &
#  else
     &                   SetFillVal = .TRUE.,                           &
#  endif
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Sout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idUbar)%dkind=PIO_FOUT
          STA(ng)%pioVar(idUbar)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idUbar)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Sout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVbar)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVbar)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVbar)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Sout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
          STA(ng)%pioVar(idu2dE)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idu2dE)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Sout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
          STA(ng)%pioVar(idv2dN)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idv2dN)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Sout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idUvel)%dkind=PIO_FOUT
          STA(ng)%pioVar(idUvel)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idUvel)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Sout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVvel)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVvel)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVvel)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Sout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
          STA(ng)%pioVar(idu3dE)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idu3dE)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Sout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
          STA(ng)%pioVar(idv3dN)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idv3dN)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Sout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWvel)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWvel)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idWvel)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component (m3/s).
!
        IF (Sout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter3 second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idOvel)%dkind=PIO_FOUT
          STA(ng)%pioVar(idOvel)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idOvel)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Sout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            STA(ng)%pioTrc(itrc)%dkind=PIO_FOUT
            STA(ng)%pioTrc(itrc)%gtype=0

#   ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#   endif
!
            status=def_var(ng, model, STA(ng)%pioFile,                  &
     &                     STA(ng)%pioTrc(itrc)%vd,                     &
     &                     PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,      &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Sout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idDano)%dkind=PIO_FOUT
          STA(ng)%pioVar(idDano)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idDano)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Sout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
          STA(ng)%pioVar(idHsbl)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idHsbl)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Sout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
          STA(ng)%pioVar(idHbbl)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idHbbl)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define vertical viscosity coefficient.
!
        IF (Sout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVvis)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVvis)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVvis)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Sout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idTdif)%dkind=PIO_FOUT
          STA(ng)%pioVar(idTdif)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idTdif)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Sout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idSdif)%dkind=PIO_FOUT
          STA(ng)%pioVar(idSdif)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idSdif)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Sout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idMtke)%dkind=PIO_FOUT
          STA(ng)%pioVar(idMtke)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idMtke)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Sout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idMtls)%dkind=PIO_FOUT
          STA(ng)%pioVar(idMtls)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idMtls)%vd,                     &
     &                   PIO_FOUT, 3, wgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Sout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idPair)%dkind=PIO_FOUT
          STA(ng)%pioVar(idPair)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idPair)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Sout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idUair)%dkind=PIO_FOUT
          STA(ng)%pioVar(idUair)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idUair)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Sout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVair)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVair)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVair)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define surface net heat flux.
!
        IF (Sout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          Vinfo( 2)=Vname(2,idTsur(itemp))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idTsur(itemp))%dkind=PIO_FOUT
          STA(ng)%pioVar(idTsur(itemp))%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idTsur(itemp))%vd,              &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SALINITY
!
!  Define surface net salt flux.
!
        IF (Sout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          Vinfo( 2)=Vname(2,idTsur(isalt))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idTsur(isalt))%dkind=PIO_FOUT
          STA(ng)%pioVar(idTsur(isalt))%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idTsur(isalt))%vd,              &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Sout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idLhea)%dkind=PIO_FOUT
          STA(ng)%pioVar(idLhea)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idLhea)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Sout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idShea)%dkind=PIO_FOUT
          STA(ng)%pioVar(idShea)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idShea)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Sout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idLrad)%dkind=PIO_FOUT
          STA(ng)%pioVar(idLrad)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idLrad)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Sout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idSrad)%dkind=PIO_FOUT
          STA(ng)%pioVar(idSrad)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idSrad)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define E-P flux (as computed by bulk_flux.F).
!
        IF (Sout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
          STA(ng)%pioVar(idEmPf)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idEmPf)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   if defined EMINUSP && defined BULK_FLUXES
!
!  Define evaporation rate.
!
        IF (Sout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idevap)%dkind=PIO_FOUT
          STA(ng)%pioVar(idevap)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idevap)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Sout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idrain)%dkind=PIO_FOUT
          STA(ng)%pioVar(idrain)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idrain)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
!
!  Define surface U-momentum stress.
!
        IF (Sout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          STA(ng)%pioVar(idUsms)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idUsms)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Sout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVsms)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVsms)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Sout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idUbms)%dkind=PIO_FOUT
          STA(ng)%pioVar(idUbms)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idUbms)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Sout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVbms)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVbms)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVbms)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
        IF (Sout(idRwet,ng)) THEN
          Vinfo( 1)=Vname(1,idRwet)
          Vinfo( 2)=Vname(2,idRwet)
          Vinfo( 3)=Vname(3,idRwet)
          Vinfo(14)=Vname(4,idRwet)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idRwet)%dkind=PIO_FOUT
          STA(ng)%pioVar(idRwet)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idRwet)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out wet/dry mask at U-points.
!
        IF (Sout(idUwet,ng)) THEN
          Vinfo( 1)=Vname(1,idUwet)
          Vinfo( 2)=Vname(2,idUwet)
          Vinfo( 3)=Vname(3,idUwet)
          Vinfo(14)=Vname(4,idUwet)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idUwet)%dkind=PIO_FOUT
          STA(ng)%pioVar(idUwet)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idUwet)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out wet/dry mask at V-points.
!
        IF (Sout(idUwet,ng)) THEN
          Vinfo( 1)=Vname(1,idVwet)
          Vinfo( 2)=Vname(2,idVwet)
          Vinfo( 3)=Vname(3,idVwet)
          Vinfo(14)=Vname(4,idVwet)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idVwet)%dkind=PIO_FOUT
          STA(ng)%pioVar(idVwet)%gtype=0
!
          status=def_var(ng, model, STA(ng)%pioFile,                    &
     &                   STA(ng)%pioVar(idVwet)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined BBL_MODEL || defined WAVES_OUTPUT
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_station_pio (ng, model, ldef, Sout, STA,           &
     &                            pgrd, rgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Define the ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_station_pio (ng, model, ldef, Sout, STA,           &
     &                            pgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_station_pio (ng, model, ldef, Sout, STA,      &
     &                                 bgrd, pgrd, rgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Define the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
        CALL wec_def_station_pio (ng, model, ldef, Sout, STA,           &
     &                            pgrd, rgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, STA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, STA(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=STA(ng)%name
!
!  Open stations file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1, STA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = STA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = STA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  stations variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            STA(ng)%pioVar(idtime)%vd=var_desc(i)
            STA(ng)%pioVar(idtime)%dkind=PIO_TOUT
            STA(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            STA(ng)%pioVar(idFsur)%vd=var_desc(i)
            STA(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            STA(ng)%pioVar(idFsur)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            STA(ng)%pioVar(idUbar)%vd=var_desc(i)
            STA(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            STA(ng)%pioVar(idUbar)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            STA(ng)%pioVar(idVbar)%vd=var_desc(i)
            STA(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVbar)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            STA(ng)%pioVar(idu2dE)%vd=var_desc(i)
            STA(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
            STA(ng)%pioVar(idu2dE)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            STA(ng)%pioVar(idv2dN)%vd=var_desc(i)
            STA(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
            STA(ng)%pioVar(idv2dN)%gtype=0
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            STA(ng)%pioVar(idUvel)%vd=var_desc(i)
            STA(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            STA(ng)%pioVar(idUvel)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            STA(ng)%pioVar(idVvel)%vd=var_desc(i)
            STA(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVvel)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            STA(ng)%pioVar(idu3dE)%vd=var_desc(i)
            STA(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
            STA(ng)%pioVar(idu3dE)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            STA(ng)%pioVar(idv3dN)%vd=var_desc(i)
            STA(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
            STA(ng)%pioVar(idv3dN)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            STA(ng)%pioVar(idWvel)%vd=var_desc(i)
            STA(ng)%pioVar(idWvel)%dkind=PIO_FOUT
            STA(ng)%pioVar(idWvel)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            STA(ng)%pioVar(idOvel)%vd=var_desc(i)
            STA(ng)%pioVar(idOvel)%dkind=PIO_FOUT
            STA(ng)%pioVar(idOvel)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            STA(ng)%pioVar(idDano)%vd=var_desc(i)
            STA(ng)%pioVar(idDano)%dkind=PIO_FOUT
            STA(ng)%pioVar(idDano)%gtype=0
#   ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            STA(ng)%pioVar(idHsbl)%vd=var_desc(i)
            STA(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
            STA(ng)%pioVar(idHsbl)%gtype=0
#   endif
#   ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            STA(ng)%pioVar(idHbbl)%vd=var_desc(i)
            STA(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
            STA(ng)%pioVar(idHbbl)%gtype=0
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            STA(ng)%pioVar(idVvis)%vd=var_desc(i)
            STA(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVvis)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            STA(ng)%pioVar(idTdif)%vd=var_desc(i)
            STA(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            STA(ng)%pioVar(idTdif)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            STA(ng)%pioVar(idSdif)%vd=var_desc(i)
            STA(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            STA(ng)%pioVar(idSdif)%gtype=0
#   if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            STA(ng)%pioVar(idMtke)%vd=var_desc(i)
            STA(ng)%pioVar(idMtke)%dkind=PIO_FOUT
            STA(ng)%pioVar(idMtke)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            STA(ng)%pioVar(idMtls)%vd=var_desc(i)
            STA(ng)%pioVar(idMtls)%dkind=PIO_FOUT
            STA(ng)%pioVar(idMtls)%gtype=0
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            STA(ng)%pioVar(idPair)%vd=var_desc(i)
            STA(ng)%pioVar(idPair)%dkind=PIO_FOUT
            STA(ng)%pioVar(idPair)%gtype=0
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            STA(ng)%pioVar(idUair)%vd=var_desc(i)
            STA(ng)%pioVar(idUair)%dkind=PIO_FOUT
            STA(ng)%pioVar(idUair)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            STA(ng)%pioVar(idVair)%vd=var_desc(i)
            STA(ng)%pioVar(idVair)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVair)%gtype=0
#   endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            STA(ng)%pioVar(idTsur(itemp))%vd=var_desc(i)
            STA(ng)%pioVar(idTsur(itemp))%dkind=PIO_FOUT
            STA(ng)%pioVar(idTsur(itemp))%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            STA(ng)%pioVar(idTsur(isalt))%vd=var_desc(i)
            STA(ng)%pioVar(idTsur(isalt))%dkind=PIO_FOUT
            STA(ng)%pioVar(idTsur(isalt))%gtype=0
#   ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            STA(ng)%pioVar(idLhea)%vd=var_desc(i)
            STA(ng)%pioVar(idLhea)%dkind=PIO_FOUT
            STA(ng)%pioVar(idLhea)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            STA(ng)%pioVar(idShea)%vd=var_desc(i)
            STA(ng)%pioVar(idShea)%dkind=PIO_FOUT
            STA(ng)%pioVar(idShea)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            STA(ng)%pioVar(idLrad)%vd=var_desc(i)
            STA(ng)%pioVar(idLrad)%dkind=PIO_FOUT
            STA(ng)%pioVar(idLrad)%gtype=0
#   endif
#   ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            STA(ng)%pioVar(idSrad)%vd=var_desc(i)
            STA(ng)%pioVar(idSrad)%dkind=PIO_FOUT
            STA(ng)%pioVar(idSrad)%gtype=0
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            STA(ng)%pioVar(idEmPf)%vd=var_desc(i)
            STA(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
            STA(ng)%pioVar(idEmPf)%gtype=0
#   if defined EMINUSP && defined BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            STA(ng)%pioVar(idevap)%vd=var_desc(i)
            STA(ng)%pioVar(idevap)%dkind=PIO_FOUT
            STA(ng)%pioVar(idevap)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            STA(ng)%pioVar(idrain)%vd=var_desc(i)
            STA(ng)%pioVar(idrain)%dkind=PIO_FOUT
            STA(ng)%pioVar(idrain)%gtype=0
#   endif
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            STA(ng)%pioVar(idUsms)%vd=var_desc(i)
            STA(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            STA(ng)%pioVar(idUsms)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            STA(ng)%pioVar(idVsms)%vd=var_desc(i)
            STA(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVsms)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            STA(ng)%pioVar(idUbms)%vd=var_desc(i)
            STA(ng)%pioVar(idUbms)%dkind=PIO_FOUT
            STA(ng)%pioVar(idUbms)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            STA(ng)%pioVar(idVbms)%vd=var_desc(i)
            STA(ng)%pioVar(idVbms)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVbms)%gtype=0
#  ifdef WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            STA(ng)%pioVar(idRwet)%vd=var_desc(i)
            STA(ng)%pioVar(idRwet)%dkind=PIO_FOUT
            STA(ng)%pioVar(idRwet)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            STA(ng)%pioVar(idUwet)%vd=var_desc(i)
            STA(ng)%pioVar(idUwet)%dkind=PIO_FOUT
            STA(ng)%pioVar(idUwet)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            STA(ng)%pioVar(idVwet)%vd=var_desc(i)
            STA(ng)%pioVar(idVwet)%dkind=PIO_FOUT
            STA(ng)%pioVar(idVwet)%gtype=0
#  endif
          END IF
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              STA(ng)%pioTrc(itrc)%vd=var_desc(i)
              STA(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              STA(ng)%pioTrc(itrc)%gtype=0
            END IF
          END DO
#  endif
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Sout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        IF (.not.got_var(idUbar).and.Sout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Sout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Sout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Sout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Sout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Sout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Sout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Sout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Sout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Sout(idOvel,ng)) THEN
          IF (Master) WRITE(stdout,60) TRIM(Vname(1,idOvel)),           &
     &                                 TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Sout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Sout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Sout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idVvis).and.Sout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Sout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SALINITY
        IF (.not.got_var(idSdif).and.Sout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Sout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Sout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Sout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Sout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Sout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idTsur(itemp)).and.Sout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Sout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Sout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Sout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Sout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Sout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idEmPf).and.Sout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   if defined EMINUSP && defined BULK_FLUXES
        IF (.not.got_var(idevap).and.Sout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Sout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
        IF (.not.got_var(idUsms).and.Sout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Sout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Sout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Sout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef WET_DRY
        IF (.not.got_var(idRwet).and.Sout(idRwet,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet).and.Sout(idUwet,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet).and.Sout(idVwet,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Sout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif

#  if defined BBL_MODEL || defined WAVES_OUTPUT
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for station variables. Get variable IDs.
!
        CALL bbl_def_station_pio (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ICE_MODEL
!
!  Scan ice model variables from input NetCDF and activate switches for
!  station variables. Get variable IDs.
!
        CALL ice_def_station_pio (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SEDIMENT
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for station variables. Get variable IDs.
!
        CALL sediment_def_station_pio (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WEC || defined WEC_VF
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for stations variables. Get variable IDs.
!
        CALL wec_def_station_pio (ng, model, ldef, Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      STA(ng)%Rindex=rec_size
        STA(ng)%Rindex=(ntstart(ng)-1)/nSTA(ng)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_STATION_PIO  - creating stations file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_STATION_PIO  - inquiring stations file',t56,      &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_STATION_PIO - unable to create stations NetCDF',  &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_STATION_PIO - unable to open stations NetCDF',    &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_STATION_PIO - unable to find variable: ',a,2x,    &
     &        ' in stations NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_station_pio
# endif
#endif
      END MODULE def_station_mod
