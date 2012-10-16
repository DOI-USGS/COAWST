#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE def_station (ng, ldef)
!
!svn $Id: def_station.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates station data NetCDF file, it defines its       !
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
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25

      logical :: got_var(NV)

      integer :: i, j, recdim, stadim
      integer :: status

      integer :: DimIDs(32), pgrd(2)
      integer :: Vsize(4)

      integer :: def_dim

# ifdef SOLVE3D
      integer :: itrc

      integer :: bgrd(3), rgrd(3), wgrd(3)
# endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
      SourceFile='def_station.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=STAname(ng)
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new station data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), ncSTAid(ng))
        IF (exit_flag.ne.NoError) THEN
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
        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 's_rho',          &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 's_w',            &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'tracer',         &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'NST',            &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'Nbed',           &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'Nphy',           &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'Nbac',           &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'Ndom',           &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'Nfec',           &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'station' ,       &
     &                 Nstation(ng), DimIDs(13))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'boundary',       &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, ncSTAid(ng), ncname, 'Nstate',         &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, ncSTAid(ng), ncname,                   &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

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
!
!  Initialize unlimited time record dimension.
!
        tSTAindx(ng)=0
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
        CALL def_info (ng, iNLM, ncSTAid(ng), ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, ncSTAid(ng), staVid(idtime,ng),        &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        IF (Sout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idFsur,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
        IF (Sout(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          Vinfo( 2)=Vname(2,idbath)
          Vinfo( 3)=Vname(3,idbath)
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idbath,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define 2D momentum in the XI-direction.
!
        IF (Sout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbar,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbar,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUvel,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVvel,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWvel,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idOvel,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
            status=def_var(ng, iNLM, ncSTAid(ng), staTid(itrc,ng),      &
     &                     NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,       &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idDano,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idHsbl,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idHbbl,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVvis,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idTdif,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSdif,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idMtke,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idMtls,ng),      &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng),                         &
     &                   staVid(idTsur(itemp),ng), NF_FOUT,             &
     &                   2, pgrd, Aval, Vinfo, ncname,                  &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
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
          status=def_var(ng, iNLM, ncSTAid(ng),                         &
     &                   staVid(idTsur(isalt),ng), NF_FOUT,             &
     &                   2, pgrd, Aval, Vinfo, ncname,                  &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idLhea,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idShea,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idLrad,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSrad,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  if defined EMINUSP && (defined BULK_FLUXES || defined CCSM_FLUX)
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idEmPf,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idevap,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idrain,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUsms,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVsms,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbms,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbms,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifdef SOLVE3D
#  ifdef BBL_MODEL
!
!  Define bottom U-current stress.
!
        IF (Sout(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          Vinfo( 2)=Vname(2,idUbrs)
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbrs,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-current stress.
!
        IF (Sout(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          Vinfo( 2)=Vname(2,idVbrs)
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbrs,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (Sout(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          Vinfo( 2)=Vname(2,idUbws)
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbws,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (Sout(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          Vinfo( 2)=Vname(2,idVbws)
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbws,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (Sout(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          Vinfo( 2)=Vname(2,idUbcs)
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbcs,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (Sout(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          Vinfo( 2)=Vname(2,idVbcs)
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbcs,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (Sout(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          Vinfo( 2)=Vname(2,idUbot)
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbot,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (Sout(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          Vinfo( 2)=Vname(2,idVbot)
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbot,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (Sout(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          Vinfo( 2)=Vname(2,idUbur)
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUbur,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (Sout(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          Vinfo( 2)=Vname(2,idVbvr)
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVbvr,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  ifdef SEDIMENT
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (Sout(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            Vinfo( 2)=Vname(2,idfrac(i))
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
            status=def_var(ng, iNLM, ncSTAid(ng), staVid(idfrac(i),ng), &
     &                     NF_FOUT, 3, bgrd, Aval, Vinfo, ncname,       &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Define sediment mass of each size class in each bed layer.
!
          IF (Sout(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            Vinfo( 2)=Vname(2,idBmas(i))
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
            status=def_var(ng, iNLM, ncSTAid(ng), staVid(idBmas(i),ng), &
     &                     NF_FOUT, 3, bgrd, Aval, Vinfo, ncname,       &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (Sout(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            Vinfo( 2)=Vname(2,idSbed(i))
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSbed(i),ng), &
     &                     NF_FOUT, 3, bgrd, Aval, Vinfo, ncname,       &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (Sout(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            Vinfo( 2)=Vname(2,idBott(i))
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, iNLM, ncSTAid(ng), staVid(idBott(i),ng), &
     &                     NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,       &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
# endif
# ifdef WEC_MELLOR
!
!  Define 2D radiation stress, Sxx-component.
!
        IF (Sout(idW2xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xx)
          Vinfo( 2)=Vname(2,idW2xx)
          Vinfo( 3)=Vname(3,idW2xx)
          Vinfo(14)=Vname(4,idW2xx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW2xx,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress, Sxy-component.
!
        IF (Sout(idW2xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xy)
          Vinfo( 2)=Vname(2,idW2xy)
          Vinfo( 3)=Vname(3,idW2xy)
          Vinfo(14)=Vname(4,idW2xy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW2xy,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress, Syy-component.
!
        IF (Sout(idW2yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2yy)
          Vinfo( 2)=Vname(2,idW2yy)
          Vinfo( 3)=Vname(3,idW2yy)
          Vinfo(14)=Vname(4,idW2yy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW2yy,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define 3D radiation stress, Sxx-component.
!
        IF (Sout(idW3xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xx)
          Vinfo( 2)=Vname(2,idW3xx)
          Vinfo( 3)=Vname(3,idW3xx)
          Vinfo(14)=Vname(4,idW3xx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3xx,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Sxy-component.
!
        IF (Sout(idW3xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xy)
          Vinfo( 2)=Vname(2,idW3xy)
          Vinfo( 3)=Vname(3,idW3xy)
          Vinfo(14)=Vname(4,idW3xy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3xy,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Syy-component.
!
        IF (Sout(idW3yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3yy)
          Vinfo( 2)=Vname(2,idW3yy)
          Vinfo( 3)=Vname(3,idW3yy)
          Vinfo(14)=Vname(4,idW3yy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3yy,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Spx-component.
!
        IF (Sout(idW3zx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zx)
          Vinfo( 2)=Vname(2,idW3zx)
          Vinfo( 3)=Vname(3,idW3zx)
          Vinfo(14)=Vname(4,idW3zx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3zx,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Spy-component.
!
        IF (Sout(idW3zy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zy)
          Vinfo( 2)=Vname(2,idW3zy)
          Vinfo( 3)=Vname(3,idW3zy)
          Vinfo(14)=Vname(4,idW3zy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3zy,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
# ifdef WEC
!
!  Define 2D u-Stokes drift velocity.
!
        IF (Sout(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          Vinfo( 2)=Vname(2,idU2Sd)
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idU2Sd,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D v-Stokes drift velocity.
!
        IF (Sout(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          Vinfo( 2)=Vname(2,idV2Sd)
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idV2Sd,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D total WEC u-stress.
!
        IF (Sout(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          Vinfo( 2)=Vname(2,idU2rs)
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idU2rs,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D total WEC v-stress.
!
        IF (Sout(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          Vinfo( 2)=Vname(2,idV2rs)
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idV2rs,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef SOLVE3D
!
!  Define 3D u-Stokes velocity.
!
        IF (Sout(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          Vinfo( 2)=Vname(2,idU3Sd)
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idU3Sd,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D v-Stokes velocity.
!
        IF (Sout(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          Vinfo( 2)=Vname(2,idV3Sd)
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idV3Sd,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D omega-Stokes velocity.
!
        IF (Sout(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          Vinfo( 2)=Vname(2,idW3Sd)
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3Sd,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!
!  Define 3D w-Stokes velocity test
!
        IF (Sout(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          Vinfo( 2)=Vname(2,idW3St)
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idW3St,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!
!  Define 3D total WEC u-stress.
!
        IF (Sout(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          Vinfo( 2)=Vname(2,idU3rs)
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idU3rs,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D total WEC v-stress.
!
        IF (Sout(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          Vinfo( 2)=Vname(2,idV3rs)
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idV3rs,ng),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
# ifdef WAVES_HEIGHT
!
!  Define wind-induced significant wave height.
!
        IF (Sout(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          Vinfo( 2)=Vname(2,idWamp)
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWamp,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_LENGTH
!
!  Define wind-induced mean wavelenght.
!
        IF (Sout(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          Vinfo( 2)=Vname(2,idWlen)
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWlen,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_LENGTHP
!
!  Define wind-induced mean wavelenght.
!
        IF (Sout(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          Vinfo( 2)=Vname(2,idWlep)
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWlep,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_DIR
!
!  Define wind-induced wave direction.
!
        IF (Sout(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          Vinfo( 2)=Vname(2,idWdir)
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWdir,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
!
!  Define wind-induced surface wave period.
!
        IF (Sout(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          Vinfo( 2)=Vname(2,idWptp)
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWptp,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# endif
# ifdef WAVES_BOT_PERIOD
!
!  Define wind-induced bottom wave period.
!
        IF (Sout(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          Vinfo( 2)=Vname(2,idWpbt)
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWpbt,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# endif
# ifdef WAVES_UB
!
!  Define wind-induced wave bottom orbital velocity.
!
        IF (Sout(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          Vinfo( 2)=Vname(2,idWorb)
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWorb,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
!
!  Define wave dissipation due to bottom friction.
!
        IF (Sout(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          Vinfo( 2)=Vname(2,idWdif)
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWdif,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR
!
!  Define wave dissipation due to breaking.
!
        IF (Sout(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdib)
          Vinfo( 2)=Vname(2,idWdib)
          Vinfo( 3)=Vname(3,idWdib)
          Vinfo(14)=Vname(4,idWdib)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWdib,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wave dissipation due to white capping.
!
        IF (Sout(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          Vinfo( 2)=Vname(2,idWdiw)
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWdiw,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined WEC_ROLLER
!
!  Define wave roller dissipation.
!
        IF (Sout(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          Vinfo( 2)=Vname(2,idWdis)
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWdis,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined ROLLER_RENIERS
!
!  Define wave roller action density.
!
        IF (Sout(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          Vinfo( 2)=Vname(2,idWrol)
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWrol,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined ROLLER_SVENDSEN
!
!  Define percent wave breaking.
!
        IF (Sout(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          Vinfo( 2)=Vname(2,idWbrk)
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWbrk,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined WEC_VF
!
!  Define WEC quasi-static sea level adjustment.
!
        IF (Sout(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          Vinfo( 2)=Vname(2,idWztw)
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWztw,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define WEC quasi-static pressure.
!
        IF (Sout(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          Vinfo( 2)=Vname(2,idWqsp)
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWqsp,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define WEC quasi-static sea level adjustment.
!
        IF (Sout(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          Vinfo( 2)=Vname(2,idWztw)
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWztw,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        IF (Sout(idUice,ng)) THEN
          Vinfo( 1)=Vname(1,idUice)
          Vinfo( 2)=Vname(2,idUice)
          Vinfo( 3)=Vname(3,idUice)
          Vinfo(14)=Vname(4,idUice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idUice,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Sout(idVice,ng)) THEN
          Vinfo( 1)=Vname(1,idVice)
          Vinfo( 2)=Vname(2,idVice)
          Vinfo( 3)=Vname(3,idVice)
          Vinfo(14)=Vname(4,idVice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idVice,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice concentration.
!
        IF (Sout(idAice,ng)) THEN
          Vinfo( 1)=Vname(1,idAice)
          Vinfo( 2)=Vname(2,idAice)
          Vinfo( 3)=Vname(3,idAice)
          Vinfo(14)=Vname(4,idAice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idAice,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice average thickness.
!
        IF (Sout(idHice,ng)) THEN
          Vinfo( 1)=Vname(1,idHice)
          Vinfo( 2)=Vname(2,idHice)
          Vinfo( 3)=Vname(3,idHice)
          Vinfo(14)=Vname(4,idHice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idHice,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice/snow surface temperature.
!
        IF (Sout(idTice,ng)) THEN
          Vinfo( 1)=Vname(1,idTice)
          Vinfo( 2)=Vname(2,idTice)
          Vinfo( 3)=Vname(3,idTice)
          Vinfo(14)=Vname(4,idTice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idTice,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define snow thickness.
!
        IF (Sout(idHsno,ng)) THEN
          Vinfo( 1)=Vname(1,idHsno)
          Vinfo( 2)=Vname(2,idHsno)
          Vinfo( 3)=Vname(3,idHsno)
          Vinfo(14)=Vname(4,idHsno)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idHsno,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface water thickness (on ice).
!
        IF (Sout(idSfwat,ng)) THEN
          Vinfo( 1)=Vname(1,idSfwat)
          Vinfo( 2)=Vname(2,idSfwat)
          Vinfo( 3)=Vname(3,idSfwat)
          Vinfo(14)=Vname(4,idSfwat)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSfwat,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice age
!
        IF (Sout(idAgeice,ng)) THEN
          Vinfo( 1)=Vname(1,idAgeice)
          Vinfo( 2)=Vname(2,idAgeice)
          Vinfo( 3)=Vname(3,idAgeice)
          Vinfo(14)=Vname(4,idAgeice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idAgeice,ng),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-ocean mass flux
!
        IF (Sout(idIomflx,ng)) THEN
          Vinfo( 1)=Vname(1,idIomflx)
          Vinfo( 2)=Vname(2,idIomflx)
          Vinfo( 3)=Vname(3,idIomflx)
          Vinfo(14)=Vname(4,idIomflx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idIomflx,ng),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice internal temperature.
!
        IF (Sout(idTimid,ng)) THEN
          Vinfo( 1)=Vname(1,idTimid)
          Vinfo( 2)=Vname(2,idTimid)
          Vinfo( 3)=Vname(3,idTimid)
          Vinfo(14)=Vname(4,idTimid)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idTimid,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water friction velocity.
!
        IF (Sout(idTauiw,ng)) THEN
          Vinfo( 1)=Vname(1,idTauiw)
          Vinfo( 2)=Vname(2,idTauiw)
          Vinfo( 3)=Vname(3,idTauiw)
          Vinfo(14)=Vname(4,idTauiw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idTauiw,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water momentum transfer coefficient.
!
        IF (Sout(idChuiw,ng)) THEN
          Vinfo( 1)=Vname(1,idChuiw)
          Vinfo( 2)=Vname(2,idChuiw)
          Vinfo( 3)=Vname(3,idChuiw)
          Vinfo(14)=Vname(4,idChuiw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idChuiw,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define under-ice temperature
!
        IF (Sout(idT0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idT0mk)
          Vinfo( 2)=Vname(2,idT0mk)
          Vinfo( 3)=Vname(3,idT0mk)
          Vinfo(14)=Vname(4,idT0mk)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idT0mk,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define under-ice salinity
!
        IF (Sout(idS0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idS0mk)
          Vinfo( 2)=Vname(2,idS0mk)
          Vinfo( 3)=Vname(3,idS0mk)
          Vinfo(14)=Vname(4,idS0mk)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idS0mk,ng),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define frazil ice growth rate
!
        IF (Sout(idWfr,ng)) THEN
          Vinfo( 1)=Vname(1,idWfr)
          Vinfo( 2)=Vname(2,idWfr)
          Vinfo( 3)=Vname(3,idWfr)
          Vinfo(14)=Vname(4,idWfr)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWfr,ng),       &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate
!
        IF (Sout(idWai,ng)) THEN
          Vinfo( 1)=Vname(1,idWai)
          Vinfo( 2)=Vname(2,idWai)
          Vinfo( 3)=Vname(3,idWai)
          Vinfo(14)=Vname(4,idWai)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWai,ng),       &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!
!  Define ice growth/melt rate
!
        IF (Sout(idWao,ng)) THEN
          Vinfo( 1)=Vname(1,idWao)
          Vinfo( 2)=Vname(2,idWao)
          Vinfo( 3)=Vname(3,idWao)
          Vinfo(14)=Vname(4,idWao)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWao,ng),       &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate
!
        IF (Sout(idWio,ng)) THEN
          Vinfo( 1)=Vname(1,idWio)
          Vinfo( 2)=Vname(2,idWio)
          Vinfo( 3)=Vname(3,idWio)
          Vinfo(14)=Vname(4,idWio)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWio,ng),       &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice melt runoff rate
!
        IF (Sout(idWro,ng)) THEN
          Vinfo( 1)=Vname(1,idWro)
          Vinfo( 2)=Vname(2,idWro)
          Vinfo( 3)=Vname(3,idWro)
          Vinfo(14)=Vname(4,idWro)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idWro,ng),       &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress sig11
!
        IF (Sout(idSig11,ng)) THEN
          Vinfo( 1)=Vname(1,idSig11)
          Vinfo( 2)=Vname(2,idSig11)
          Vinfo( 3)=Vname(3,idSig11)
          Vinfo(14)=Vname(4,idSig11)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSig11,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress sig12
!
        IF (Sout(idSig12,ng)) THEN
          Vinfo( 1)=Vname(1,idSig12)
          Vinfo( 2)=Vname(2,idSig12)
          Vinfo( 3)=Vname(3,idSig12)
          Vinfo(14)=Vname(4,idSig12)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSig12,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress sig22
!
        IF (Sout(idSig22,ng)) THEN
          Vinfo( 1)=Vname(1,idSig22)
          Vinfo( 2)=Vname(2,idSig22)
          Vinfo( 3)=Vname(3,idSig22)
          Vinfo(14)=Vname(4,idSig22)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, ncSTAid(ng), staVid(idSig22,ng),     &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, ncSTAid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, ncSTAid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=STAname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open stations file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, ncSTAid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
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
            staVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            staVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            staVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            staVid(idVbar,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            staVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            staVid(idVvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            staVid(idWvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            staVid(idOvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            staVid(idDano,ng)=var_id(i)
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            staVid(idHsbl,ng)=var_id(i)
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            staVid(idHbbl,ng)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            staVid(idVvis,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            staVid(idTdif,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            staVid(idSdif,ng)=var_id(i)
#  if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            staVid(idMtke,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            staVid(idMtls,ng)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            staVid(idTsur(itemp),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            staVid(idTsur(isalt),ng)=var_id(i)
#  ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            staVid(idLhea,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            staVid(idShea,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            staVid(idLrad,ng)=var_id(i)
#  endif
#  ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            staVid(idSrad,ng)=var_id(i)
#  endif
#  if defined EMINUSP && defined BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            staVid(idEmPf,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            staVid(idevap,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            staVid(idrain,ng)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            staVid(idUsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            staVid(idVsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            staVid(idUbms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            staVid(idVbms,ng)=var_id(i)
# ifdef SOLVE3D
#  ifdef BBL_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            staVid(idUbrs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            staVid(idVbrs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            staVid(idUbws,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            staVid(idVbws,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            staVid(idUbcs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            staVid(idVbcs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            staVid(idUbot,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            staVid(idVbot,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            staVid(idUbur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            staVid(idVbvr,ng)=var_id(i)
#  endif
#  ifdef ICE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.TRUE.
            staVid(idUice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.TRUE.
            staVid(idVice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.TRUE.
            staVid(idAice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.TRUE.
            staVid(idHice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.TRUE.
            staVid(idTice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.TRUE.
            staVid(idHsno,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSfwat))) THEN
            got_var(idSfwat)=.TRUE.
            staVid(idSfwat,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.TRUE.
            staVid(idAgeice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
            got_var(idIomflx)=.TRUE.
            staVid(idIomflx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.TRUE.
            staVid(idTimid,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.TRUE.
            staVid(idTauiw,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.TRUE.
            staVid(idChuiw,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.TRUE.
            staVid(idT0mk,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.TRUE.
            staVid(idS0mk,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
            got_var(idWfr)=.TRUE.
            staVid(idWfr,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
            got_var(idWai)=.TRUE.
            staVid(idWai,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
            got_var(idWao)=.TRUE.
            staVid(idWao,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
            got_var(idWio)=.TRUE.
            staVid(idWio,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
            got_var(idWro)=.TRUE.
            staVid(idWro,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.TRUE.
            staVid(idSig11,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.TRUE.
            staVid(idSig12,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.TRUE.
            staVid(idSig22,ng)=var_id(i)
# endif
# endif
# ifdef WEC_MELLOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xx))) THEN
            got_var(idW2xx)=.TRUE.
            staVid(idW2xx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xy))) THEN
            got_var(idW2xy)=.TRUE.
            staVid(idW2xy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2yy))) THEN
            got_var(idW2yy)=.TRUE.
            staVid(idW2yy,ng)=var_id(i)
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xx))) THEN
            got_var(idW3xx)=.TRUE.
            staVid(idW3xx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xy))) THEN
            got_var(idW3xy)=.TRUE.
            staVid(idW3xy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3yy))) THEN
            got_var(idW3yy)=.TRUE.
            staVid(idW3yy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zx))) THEN
            got_var(idW3zx)=.TRUE.
            staVid(idW3zx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zy))) THEN
            got_var(idW3zy)=.TRUE.
            staVid(idW3zy,ng)=var_id(i)
#  endif
# endif
# ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            staVid(idU2Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            staVid(idV2Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            staVid(idU2rs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            staVid(idV2rs,ng)=var_id(i)
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            staVid(idU3Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            staVid(idV3Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            staVid(idW3Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            staVid(idW3St,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            staVid(idU3rs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            staVid(idV3rs,ng)=var_id(i)
#  endif
# endif
# ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            staVid(idWamp,ng)=var_id(i)
# endif
# ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            staVid(idWlen,ng)=var_id(i)
# endif
# ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            staVid(idWlep,ng)=var_id(i)
# endif
# ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            staVid(idWdir,ng)=var_id(i)
# endif
# ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            staVid(idWptp,ng)=var_id(i)
# endif
# ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            staVid(idWpbt,ng)=var_id(i)
# endif
# ifdef WAVES_UB
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            staVid(idWorb,ng)=var_id(i)
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            staVid(idWdif,ng)=var_id(i)
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            staVid(idWdib,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            staVid(idWdiw,ng)=var_id(i)
# endif
# if defined WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            staVid(idWdis,ng)=var_id(i)
# endif
# if defined ROLLER_RENIERS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            staVid(idWrol,ng)=var_id(i)
# endif
# if defined ROLLER_SVENDSEN
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            staVid(idWbrk,ng)=var_id(i)
# endif
# if defined WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            staVid(idWztw,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            staVid(idWqsp,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            staVid(idWbeh,ng)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              staTid(itrc,ng)=var_id(i)
            END IF
          END DO
#  ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              staVid(idfrac(itrc),ng)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              staVid(idBmas(itrc),ng)=var_id(i)
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              staVid(idSbed(itrc),ng)=var_id(i)
            END IF
          END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              staVid(idBott(itrc),ng)=var_id(i)
            END IF
          END DO
#  endif
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
#  if defined EMINUSP && defined BULK_FLUXES
        IF (.not.got_var(idEmPf).and.Sout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
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
# ifdef SOLVE3D
#  ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.Sout(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.Sout(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.Sout(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.Sout(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.Sout(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.Sout(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.Sout(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.Sout(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.Sout(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.Sout(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef WEC_MELLOR
        IF (.not.got_var(idW2xx).and.Sout(idW2xx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2xy).and.Sout(idW2xy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2yy).and.Sout(idW2yy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idW3xx).and.Sout(idW3xx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3xy).and.Sout(idW3xy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3yy).and.Sout(idW3yy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zx).and.Sout(idW3zx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3zx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zy).and.Sout(idW3zy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3zy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef WEC
        IF (.not.got_var(idU2Sd).and.Sout(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Sout(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Sout(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Sout(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idU3Sd).and.Sout(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Sout(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.Sout(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.Sout(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Sout(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Sout(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.Sout(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.Sout(idWlen,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWlen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.Sout(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.Sout(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.Sout(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.Sout(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_UB
        IF (.not.got_var(idWorb).and.Sout(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (.not.got_var(idWdif).and.Sout(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (.not.got_var(idWdib).and.Sout(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.Sout(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined WEC_ROLLER
        IF (.not.got_var(idWdis).and.Sout(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined ROLLER_RENIERS
        IF (.not.got_var(idWrol).and.Sout(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined ROLLER_SVENDSEN
        IF (.not.got_var(idWbrk).and.Sout(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined WEC_VF
        IF (.not.got_var(idWztw).and.Sout(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.Sout(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.Sout(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWbeh)),          &
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
#  ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.Sout(idfrac(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idBmas(i)).and.Sout(idBmas(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.Sout(idSbed(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.Sout(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
# ifdef ICE_MODEL
        IF (.not.got_var(idUice).and.Sout(idUice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Sout(idVice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Sout(idAice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Sout(idHice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Sout(idTice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Sout(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSfwat).and.Sout(idSfwat,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSfwat)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAgeice).and.Sout(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Sout(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Sout(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Sout(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Sout(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Sout(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Sout(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Sout(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Sout(idWai,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Sout(idWao,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Sout(idWio,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Sout(idWro,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Sout(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Sout(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Sout(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      tSTAindx(ng)=rec_size
        tSTAindx(ng)=(ntstart(ng)-1)/nSTA(ng)
      END IF QUERY
!
  10  FORMAT (6x,'DEF_STATION - creating stations file: ',a)
  20  FORMAT (6x,'DEF_STATION - inquiring stations file: ',a)
  30  FORMAT (/,' DEF_STATION - unable to create stations NetCDF ',     &
     &        'file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_STATION - unable to open stations NetCDF file: ', &
     &        a)
  60  FORMAT (/,' DEF_STATION - unable to find variable: ',a,2x,        &
     &        ' in stations NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_station
#else
      SUBROUTINE def_station
      RETURN
      END SUBROUTINE def_station
#endif
