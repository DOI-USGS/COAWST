#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE def_station (ng,ldef)
!
!svn $Id: def_station.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
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
      integer, parameter :: Natt = 24

      logical :: got_var(NV)

      integer :: i, j, nrec, nvd, recdim, stadim
      integer :: status

      integer :: DimIDs(29), pgrd(2)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

# ifdef SOLVE3D
      integer :: itrc

      integer :: bgrd(3), rgrd(3), wgrd(3)
# endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
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
      IF (ldef.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncSTAid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iNLM, ncSTAid(ng), 1)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
# ifdef SOLVE3D
        status=nf90_def_dim(ncSTAid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncSTAid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncSTAid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
#  ifdef SEDIMENT
        status=nf90_def_dim(ncSTAid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#  endif
#  ifdef ECOSIM
        status=nf90_def_dim(ncSTAid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncSTAid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncSTAid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncSTAid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
#  endif
# endif
        status=nf90_def_dim(ncSTAid(ng),'station' ,Nstation(ng),        &
     &                      DimIDs(13))
        status=nf90_def_dim(ncSTAid(ng),'boundary',4,                   &
     &                      DimIDs(14))
# ifdef FOUR_DVAR
        status=nf90_def_dim(ncSTAid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
# endif
        status=nf90_def_dim(ncSTAid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
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
          DO j=1,80
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
        status=def_info(ng,ncSTAid(ng),ncname,DimIDs)
        if (exit_flag.ne.NoError) RETURN
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
        status=def_var(ncSTAid(ng),staVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        IF (Sout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ncSTAid(ng),staVid(idFsur,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbar,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbar,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUvel,ng),NF_FOUT,         &
     &                   3,rgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVvel,ng),NF_FOUT,         &
     &                   3,rgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idWvel,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idOvel,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncSTAid(ng),staTid(itrc,ng),NF_FOUT,         &
     &                     3,rgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idDano,ng),NF_FOUT,         &
     &                   3,rgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idHsbl,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idHbbl,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVvis,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idTdif,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idSdif,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idMtke,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idMtls,ng),NF_FOUT,         &
     &                   3,wgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idTsur(itemp),ng),NF_FOUT,  &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idTsur(isalt),ng),NF_FOUT,  &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idLhea,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idShea,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idLrad,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idSrad,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUsms,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVsms,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbms,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbms,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbrs,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbrs,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbws,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbws,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbcs,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbcs,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbot,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbot,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idUbur,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncSTAid(ng),staVid(idVbvr,ng),NF_FOUT,         &
     &                   2,pgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncSTAid(ng),staVid(idfrac(i),ng),NF_FOUT,    &
     &                     3,bgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncSTAid(ng),staVid(idBmas(i),ng),NF_FOUT,    &
     &                     3,bgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncSTAid(ng),staVid(idSbed(i),ng),NF_FOUT,    &
     &                     3,bgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncSTAid(ng),staVid(idBott(i),ng),NF_FOUT,    &
     &                     2,pgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncSTAid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iTLM, ncSTAid(ng), Master, ncname)
        if (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of stations NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        if (exit_flag.ne.NoError) RETURN
!
!  Open stations file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncSTAid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(ncname)
          exit_flag=3
          ioerror=status
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
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idtime)),    &
     &                            staVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idFsur)),    &
     &                            staVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbar)),    &
     &                            staVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbar)),    &
     &                            staVid(idVbar,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUvel)),    &
     &                            staVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVvel)),    &
     &                            staVid(idVvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idWvel)),    &
     &                            staVid(idWvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idOvel)),    &
     &                            staVid(idOvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              status=nf90_inq_varid(ncSTAid(ng),                        &
     &                              TRIM(Vname(1,idTvar(itrc))),        &
     &                              staTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idDano)),    &
     &                            staVid(idDano,ng))
          END IF
#  ifdef LMD_SKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idHsbl)),    &
     &                            staVid(idHsbl,ng))
          END IF
#  endif
#  ifdef LMD_BKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idHbbl)),    &
     &                            staVid(idHbbl,ng))
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVvis)),    &
     &                            staVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idTdif)),    &
     &                            staVid(idTdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idSdif)),    &
     &                            staVid(idSdif,ng))
          END IF
#  if defined GLS_MIXING || defined MY25_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idMtke)),    &
     &                            staVid(idMtke,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idMtls)),    &
     &                            staVid(idMtls,ng))
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),                          &
     &                            TRIM(Vname(1,idTsur(itemp))),         &
     &                            staVid(idTsur(itemp),ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),                          &
     &                            TRIM(Vname(1,idTsur(isalt))),         &
     &                            staVid(idTsur(isalt),ng))
          END IF
#  ifdef BULK_FLUXES
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idLhea)),    &
     &                            staVid(idLhea,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idShea)),    &
     &                            staVid(idShea,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idLrad)),    &
     &                            staVid(idLrad,ng))
          END IF
#  endif
#  ifdef SHORTWAVE
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idSrad)),    &
     &                            staVid(idSrad,ng))
          END IF
#  endif
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUsms)),    &
     &                            staVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVsms)),    &
     &                            staVid(idVsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbms)),    &
     &                            staVid(idUbms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbms)),    &
     &                            staVid(idVbms,ng))
          END IF
# ifdef SOLVE3D
#  ifdef BBL_MODEL
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbrs)),    &
     &                            staVid(idUbrs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbrs)),    &
     &                            staVid(idVbrs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbws)),    &
     &                            staVid(idUbws,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbws)),    &
     &                            staVid(idVbws,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbcs)),    &
     &                            staVid(idUbcs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbcs)),    &
     &                            staVid(idVbcs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbot)),    &
     &                            staVid(idUbot,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbot)),    &
     &                            staVid(idVbot,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idUbur)),    &
     &                            staVid(idUbur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            status=nf90_inq_varid(ncSTAid(ng),TRIM(Vname(1,idVbvr)),    &
     &                            staVid(idVbvr,ng))
          END IF
#  endif
#  ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              status=nf90_inq_varid(ncSTAid(ng),                        &
     &                              TRIM(Vname(1,idfrac(itrc))),        &
     &                              staVid(idfrac(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              status=nf90_inq_varid(ncSTAid(ng),                        &
     &                              TRIM(Vname(1,idBmas(itrc))),        &
     &                              staVid(idBmas(itrc),ng))
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              status=nf90_inq_varid(ncSTAid(ng),                        &
     &                              TRIM(Vname(1,idSbed(itrc))),        &
     &                              staVid(idSbed(itrc),ng))
            END IF
          END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              status=nf90_inq_varid(ncSTAid(ng),                        &
     &                              TRIM(Vname(1,idBott(itrc))),        &
     &                              staVid(idBott(itrc),ng))
            END IF
          END DO
#  endif
# endif
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Sout(idFsur,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Sout(idUbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Sout(idVbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Sout(idUvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Sout(idVvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Sout(idWvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Sout(idOvel,ng)) THEN
          WRITE(stdout,60) TRIM(Vname(1,idOvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Sout(idTvar(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano).and.Sout(idDano,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Sout(idHsbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHsbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Sout(idHbbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHbbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idVvis).and.Sout(idVvis,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvis)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Sout(idTdif,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SALINITY
        IF (.not.got_var(idSdif).and.Sout(idSdif,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Sout(idMtke,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idMtke)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Sout(idMtls,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idMtls)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idTsur(itemp)).and.Sout(idTsur(itemp),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Sout(idTsur(isalt),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Sout(idLhea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLhea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Sout(idShea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idShea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Sout(idLrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Sout(idSrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Sout(idUsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Sout(idVsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Sout(idUbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Sout(idVbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
#  ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.Sout(idUbrs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbrs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.Sout(idVbrs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbrs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.Sout(idUbws,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbws)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.Sout(idVbws,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbws)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.Sout(idUbcs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbcs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.Sout(idVbcs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbcs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.Sout(idUbot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.Sout(idVbot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.Sout(idUbur,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.Sout(idVbvr,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbvr)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.Sout(idfrac(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idfrac(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idBmas(i)).and.Sout(idBmas(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idBmas(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.Sout(idSbed(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idSbed(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.Sout(idBott(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idBott(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      tSTAindx(ng)=tsize
        tSTAindx(ng)=(ntstart(ng)-1)/nSTA(ng)
      END IF
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
