#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE def_filt (ng)
!
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates filter NetCDF file, it defines its             !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_filter
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25

      logical :: Ldefine, got_var(NV)

      integer :: i, j, ifield, nvd3, nvd4
      integer :: recdim, status, ifile, varid

      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)

# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=80 ) :: type
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
!=======================================================================
!  Create a new filter NetCDF file.
!=======================================================================
!
!  Activate creation of filter NetCDF file.  Create a new filter
!  file if during a restart run, the filter filename "FIL%name"
!  is different than the restart filename "RST%name".  I would prefer to
!  not create a new filter file during a restart run, but simply write to 
!  one assumed to have been created before, but this doesn't work yet.
!
      IF (exit_flag.ne.NoError) RETURN
      Ldefine=.FALSE.
      IF (((nrrec(ng).eq.0).and.(iic(ng).eq.ntstart(ng))).or.           &
     &      ((nrrec(ng).ne.0))) THEN
        Ldefine=.TRUE.
      END IF
!
      DO ifile = 1,nfile
        ncname=FIL(ifile,ng)%name
        IF (Ldefine) THEN
          CALL netcdf_create(ng, iNLM, TRIM(ncname), FIL(ifile,ng)%ncid)
          IF (exit_flag.ne.NoError) THEN
            IF (Master) WRITE (stdout,30) TRIM(ncname)
            RETURN
          END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xi_rho', IOBOUNDS(ng)%xi_rho, DimIDs( 1))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xi_u', IOBOUNDS(ng)%xi_u, DimIDs( 2))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xi_v', IOBOUNDS(ng)%xi_v, DimIDs( 3))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xi_psi', IOBOUNDS(ng)%xi_psi, DimIDs( 4))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'eta_rho', IOBOUNDS(ng)%eta_rho, DimIDs( 5))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'eta_u', IOBOUNDS(ng)%eta_u, DimIDs( 6))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'eta_v', IOBOUNDS(ng)%eta_v, DimIDs( 7))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'eta_psi', IOBOUNDS(ng)%eta_psi, DimIDs( 8))
          IF (exit_flag.ne.NoError) RETURN
# if defined WRITE_WATER && defined MASKING
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xy_rho', IOBOUNDS(ng)%xy_rho, DimIDs(17))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xy_u', IOBOUNDS(ng)%xy_u, DimIDs(18))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xy_v', IOBOUNDS(ng)%xy_v, DimIDs(19))
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xyz_rho', IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xyz_u', IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xyz_v', IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'xyz_w', IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
          IF (exit_flag.ne.NoError) RETURN
#  endif
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'N', N(ng), DimIDs( 9))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             's_rho', N(ng), DimIDs( 9))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             's_w', N(ng)+1, DimIDs(10))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'tracer', NT(ng), DimIDs(11))
          IF (exit_flag.ne.NoError) RETURN
# endif
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'boundary', 4, DimIDs(14))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
     &             'filterLevel', nf90_unlimited, DimIDs(32))
          IF (exit_flag.ne.NoError) RETURN
!         status = def_dim(ng, iNLM, FIL(ifile,ng)%ncid, ncname,        &
!     &            TRIM(ADJUSTL(Vname(5,idtime))),                      &
!     &            nf90_unlimited, DimIDs(12))
!         IF (exit_flag.ne.NoError) RETURN
          recdim=DimIDs(32)
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
          nvd3=2
          nvd4=2
# else
          nvd3=3
          nvd4=4
# endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
          t2dgrd(1)=DimIDs(17)
          t2dgrd(2)=DimIDs(32)
#  ifdef SOLVE3D
          t3dgrd(1)=DimIDs(20)
          t3dgrd(2)=DimIDs(32)
#  endif
# else
          t2dgrd(1)=DimIDs( 1)
          t2dgrd(2)=DimIDs( 5)
          t2dgrd(3)=DimIDs(32)
#  ifdef SOLVE3D
          t3dgrd(1)=DimIDs( 1)
          t3dgrd(2)=DimIDs( 5)
          t3dgrd(3)=DimIDs( 9)
          t3dgrd(4)=DimIDs(32)
#  endif
# endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
          u2dgrd(1)=DimIDs(18)
          u2dgrd(2)=DimIDs(32)
#  ifdef SOLVE3D
          u3dgrd(1)=DimIDs(21)
          u3dgrd(2)=DimIDs(32)
#  endif
# else
          u2dgrd(1)=DimIDs( 2)
          u2dgrd(2)=DimIDs( 6)
          u2dgrd(3)=DimIDs(32)
#  ifdef SOLVE3D
          u3dgrd(1)=DimIDs( 2)
          u3dgrd(2)=DimIDs( 6)
          u3dgrd(3)=DimIDs( 9)
          u3dgrd(4)=DimIDs(32)
#  endif
# endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
          v2dgrd(1)=DimIDs(19)
          v2dgrd(2)=DimIDs(32)
#  ifdef SOLVE3D
          v3dgrd(1)=DimIDs(22)
          v3dgrd(2)=DimIDs(32)
#  endif
# else
          v2dgrd(1)=DimIDs( 3)
          v2dgrd(2)=DimIDs( 7)
          v2dgrd(3)=DimIDs(32)
#  ifdef SOLVE3D
          v3dgrd(1)=DimIDs( 3)
          v3dgrd(2)=DimIDs( 7)
          v3dgrd(3)=DimIDs( 9)
          v3dgrd(4)=DimIDs(32)
#  endif
# endif
# ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
#  if defined WRITE_WATER && defined MASKING
          w3dgrd(1)=DimIDs(23)
          w3dgrd(2)=DimIDs(32)
#  else
          w3dgrd(1)=DimIDs( 1)
          w3dgrd(2)=DimIDs( 5)
          w3dgrd(3)=DimIDs(10)
          w3dgrd(4)=DimIDs(32)
#  endif
# endif
!
!  Initialize local information variable arrays.
!
          DO i=1,Natt
            DO j=1,LEN(Vinfo(1))
              Vinfo(i)=' '
            END DO
          END DO
          DO i=1,6
            Aval(i)=0.0_r8
          END DO
!
!-----------------------------------------------------------------------
!  Define common variables and their attributes.
!-----------------------------------------------------------------------
          IF (Master) WRITE (stdout,10) TRIM(ncname)
!
!  Define file type, in lieu of whole def_info.
!
          type='ROMS/TOMS filter status file for restarts'
          IF (exit_flag.eq.NoError) THEN
            status = nf90_put_att(FIL(ifile,ng)%ncid, nf90_global,      &
     &               'type', TRIM(type))
            IF (status.ne.nf90_noerr) THEN 
              IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  Define model time as filter levels.
!
          Vinfo( 1)='filterLevel'
          Vinfo( 2)='Number of days needed in the time filtering'
          status = def_var(ng, iNLM, FIL(ifile,ng)%ncid,                &
     &           varid, NF90_INT, 1, (/recdim/), Aval,                  &
     &           Vinfo, ncname, SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
!
!  Define ocean_time as scalar.
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
          status = def_var(ng, iNLM, FIL(ifile,ng)%ncid,                &
     &             FIL(ifile,ng)%Vid(idtime), NF_TYPE, 1, (/ recdim /), &
     &             Aval, Vinfo, ncname, SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
!
!  Define filter counts.
!
          Vinfo( 1)='fcount'
          Vinfo( 2)=Vname(2,idtime)
          status = def_var(ng, iNLM, FIL(ifile,ng)%ncid,                &
     &             varid, NF90_INT, 1, (/recdim/),                      &
     &             Aval, Vinfo, ncname, SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Define variables and their attributes to each filter file separately.
!-----------------------------------------------------------------------
      IF (Ldefine) THEN
!
!  Define free-surface.
!
        ncname=FIL(1,ng)%name
        IF (Aout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idFsur),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idUbar),                           &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idVbar),                           &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
          ncname=FIL(2,ng)%name
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, iNLM, FIL(2,ng)%ncid,                      &
     &                 FIL(2,ng)%Vid(idUvel),                           &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
          ncname=FIL(3,ng)%name
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, iNLM, FIL(3,ng)%ncid,                      &
     &                 FIL(3,ng)%Vid(idVvel),                           &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical (omega) velocity.
!
        IF (Aout(idOvel,ng)) THEN
          ncname=FIL(9,ng)%name
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, iNLM, FIL(9,ng)%ncid,                      &
     &                 FIL(9,ng)%Vid(idOvel),                           &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical (W) velocity.
!
        IF (Aout(idWvel,ng)) THEN
          ncname=FIL(9,ng)%name
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, iNLM, FIL(9,ng)%ncid,                      &
     &                 FIL(9,ng)%Vid(idWvel),                           &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define density anomaly.
!
        IF (Aout(idDano,ng)) THEN
          ncname=FIL(7,ng)%name
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, iNLM, FIL(7,ng)%ncid,                      &
     &                 FIL(7,ng)%Vid(idDano),                           &
     &                 NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NAT
          IF (Aout(idTvar(itrc),ng)) THEN
            ncname=FIL(8,ng)%name
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, FIL(8,ng)%ncid,                    &
     &                 FIL(8,ng)%Tid(itrc),                             &
     &                 NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!
!  Define vertical viscosity coefficient.
!
        IF (Aout(idVvis,ng)) THEN
          ncname=FIL(6,ng)%name
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, iNLM, FIL(6,ng)%ncid,                      &
     &                 FIL(6,ng)%Vid(idVvis),                           &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Aout(idTdif,ng)) THEN
          ncname=FIL(4,ng)%name
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, iNLM, FIL(4,ng)%ncid,                      &
     &                 FIL(4,ng)%Vid(idTdif),                           &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Aout(idSdif,ng)) THEN
          ncname=FIL(5,ng)%name
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, iNLM, FIL(5,ng)%ncid,                      &
     &                 FIL(5,ng)%Vid(idSdif),                           &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif

#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Aout(idHsbl,ng)) THEN
          ncname=FIL(1,ng)%name
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idHsbl),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Aout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                   FIL(1,ng)%Vid(idHbbl),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
!
!  Define surface net heat flux.
!
        ncname=FIL(1,ng)%name
        IF (Aout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          Vinfo( 2)=Vname(2,idTsur(itemp))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idTsur(itemp)),                    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Aout(idTsur(isalt),ng)) THEN
          ncname=FIL(1,ng)%name
          Vinfo( 1)=Vname(1,idTsur(isalt))
          Vinfo( 2)=Vname(2,idTsur(isalt))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, rain'
          Vinfo(12)='downward flux, evap'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idTsur(isalt)),                    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Aout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idLhea),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Aout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idShea),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Aout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idLrad),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Aout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idSrad),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
!
!  Define surface u-momentum stress.
!
        IF (Aout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idUsms),                           &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Aout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idVsms),                           &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Aout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idUbms),                           &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Aout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, iNLM, FIL(1,ng)%ncid,                      &
     &                 FIL(1,ng)%Vid(idVbms),                           &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        ncname=FIL(10,ng)%name
        IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
          Vinfo( 1)=Vname(1,idUice)
          Vinfo( 2)=Vname(2,idUice)
          Vinfo( 3)=Vname(3,idUice)
          Vinfo(14)=Vname(4,idUice)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUice,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idUice),                          &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
          Vinfo( 1)=Vname(1,idVice)
          Vinfo( 2)=Vname(2,idVice)
          Vinfo( 3)=Vname(3,idVice)
          Vinfo(14)=Vname(4,idVice)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVice,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idVice),                          &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice concentration
!
        IF (Aout(idAice,ng)) THEN
          Vinfo( 1)=Vname(1,idAice)
          Vinfo( 2)=Vname(2,idAice)
          Vinfo( 3)=Vname(3,idAice)
          Vinfo(14)=Vname(4,idAice)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAice,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idAice),                          &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice thickness
!
        IF (Aout(idHice,ng)) THEN
          Vinfo( 1)=Vname(1,idHice)
          Vinfo( 2)=Vname(2,idHice)
          Vinfo( 3)=Vname(3,idHice)
          Vinfo(14)=Vname(4,idHice)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHice,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idHice),                          &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define snow thickness
!
        IF (Aout(idHsno,ng)) THEN
          Vinfo( 1)=Vname(1,idHsno)
          Vinfo( 2)=Vname(2,idHsno)
          Vinfo( 3)=Vname(3,idHsno)
          Vinfo(14)=Vname(4,idHsno)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsno,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idHsno),                          &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#ifdef MELT_PONDS
!
!  Define surface melt ponds
!
        IF (Aout(idApond,ng)) THEN
          Vinfo( 1)=Vname(1,idApond)
          Vinfo( 2)=Vname(2,idApond)
          Vinfo( 3)=Vname(3,idApond)
          Vinfo(14)=Vname(4,idApond)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idApond,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idApond),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface melt ponds
!
        IF (Aout(idHpond,ng)) THEN
          Vinfo( 1)=Vname(1,idHpond)
          Vinfo( 2)=Vname(2,idHpond)
          Vinfo( 3)=Vname(3,idHpond)
          Vinfo(14)=Vname(4,idHpond)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHpond,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idHpond),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Define ice surface temperature
!
        IF (Aout(idTice,ng)) THEN
          Vinfo( 1)=Vname(1,idTice)
          Vinfo( 2)=Vname(2,idTice)
          Vinfo( 3)=Vname(3,idTice)
          Vinfo(14)=Vname(4,idTice)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTice,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idTice),                          &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice temperature
!
        IF (Aout(idTimid,ng)) THEN
          Vinfo( 1)=Vname(1,idTimid)
          Vinfo( 2)=Vname(2,idTimid)
          Vinfo( 3)=Vname(3,idTimid)
          Vinfo(14)=Vname(4,idTimid)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTimid,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idTimid),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define age of ice
!
        IF (Aout(idAgeice,ng)) THEN
          Vinfo( 1)=Vname(1,idAgeice)
          Vinfo( 2)=Vname(2,idAgeice)
          Vinfo( 3)=Vname(3,idAgeice)
          Vinfo(14)=Vname(4,idAgeice)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAgeice,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idAgeice),                        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define under-ice salinity
!
        IF (Aout(idS0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idS0mk)
          Vinfo( 2)=Vname(2,idS0mk)
          Vinfo( 3)=Vname(3,idS0mk)
          Vinfo(14)=Vname(4,idS0mk)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idS0mk,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idS0mk),                          &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define under-ice temperature
!
        IF (Aout(idT0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idT0mk)
          Vinfo( 2)=Vname(2,idT0mk)
          Vinfo( 3)=Vname(3,idT0mk)
          Vinfo(14)=Vname(4,idT0mk)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idT0mk,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idT0mk),                          &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sig11 component of internal ice stress
!
        IF (Aout(idSig11,ng)) THEN
          Vinfo( 1)=Vname(1,idSig11)
          Vinfo( 2)=Vname(2,idSig11)
          Vinfo( 3)=Vname(3,idSig11)
          Vinfo(14)=Vname(4,idSig11)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig11,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idSig11),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sig12 component of internal ice stress
!
        IF (Aout(idSig12,ng)) THEN
          Vinfo( 1)=Vname(1,idSig12)
          Vinfo( 2)=Vname(2,idSig12)
          Vinfo( 3)=Vname(3,idSig12)
          Vinfo(14)=Vname(4,idSig12)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig12,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idSig12),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sig22 component of internal ice stress
!
        IF (Aout(idSig22,ng)) THEN
          Vinfo( 1)=Vname(1,idSig22)
          Vinfo( 2)=Vname(2,idSig22)
          Vinfo( 3)=Vname(3,idSig22)
          Vinfo(14)=Vname(4,idSig22)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig22,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idSig22),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-ocean mass flux
!
        IF (Aout(idIomflx,ng)) THEN
          Vinfo( 1)=Vname(1,idIomflx)
          Vinfo( 2)=Vname(2,idIomflx)
          Vinfo( 3)=Vname(3,idIomflx)
          Vinfo(14)=Vname(4,idIomflx)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIomflx,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idIomflx),                        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water friction velocity
!
        IF (Aout(idTauiw,ng)) THEN
          Vinfo( 1)=Vname(1,idTauiw)
          Vinfo( 2)=Vname(2,idTauiw)
          Vinfo( 3)=Vname(3,idTauiw)
          Vinfo(14)=Vname(4,idTauiw)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTauiw,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idTauiw),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water momentum transfer coefficient
!
        IF (Aout(idChuiw,ng)) THEN
          Vinfo( 1)=Vname(1,idChuiw)
          Vinfo( 2)=Vname(2,idChuiw)
          Vinfo( 3)=Vname(3,idChuiw)
          Vinfo(14)=Vname(4,idChuiw)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idChuiw,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idChuiw),                         &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define frazil ice growth rate.
!
        IF (Aout(idWfr,ng)) THEN
          Vinfo( 1)=Vname(1,idWfr)
          Vinfo( 2)=Vname(2,idWfr)
          Vinfo( 3)=Vname(3,idWfr)
          Vinfo(14)=Vname(4,idWfr)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWfr,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idWfr),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-air melt rate.
!
        IF (Aout(idWai,ng)) THEN
          Vinfo( 1)=Vname(1,idWai)
          Vinfo( 2)=Vname(2,idWai)
          Vinfo( 3)=Vname(3,idWai)
          Vinfo(14)=Vname(4,idWai)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWai,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idWai),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ocean-air ice growth rate.
!
        IF (Aout(idWao,ng)) THEN
          Vinfo( 1)=Vname(1,idWao)
          Vinfo( 2)=Vname(2,idWao)
          Vinfo( 3)=Vname(3,idWao)
          Vinfo(14)=Vname(4,idWao)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWao,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idWao),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-ocean ice melt/growth rate.
!
        IF (Aout(idWio,ng)) THEN
          Vinfo( 1)=Vname(1,idWio)
          Vinfo( 2)=Vname(2,idWio)
          Vinfo( 3)=Vname(3,idWio)
          Vinfo(14)=Vname(4,idWio)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWio,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idWio),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice melt runoff rate.
!
        IF (Aout(idWro,ng)) THEN
          Vinfo( 1)=Vname(1,idWro)
          Vinfo( 2)=Vname(2,idWro)
          Vinfo( 3)=Vname(3,idWro)
          Vinfo(14)=Vname(4,idWro)
          Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWro,ng),r8)
          status=def_var(ng, iNLM, FIL(10,ng)%ncid,                     &
     &                 FIL(10,ng)%Vid(idWro),                           &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        DO ifile=1,nfile
          CALL netcdf_enddef (ng, iNLM, FIL(ifile,ng)%name,             &
     &                        FIL(ifile,ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF
!
!=======================================================================
!  Open existing filter files, and prepare for appending data.
!=======================================================================
!
      IF (.not.Ldefine) THEN
!
!  Inquire about the contents of filter NetCDF files:  Don't inquire
!  about the dimensions and variables.  Check for consistency.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO

        DO ifile=1,nfile
          ncname=FIL(ifile,ng)%name
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

          IF (Master) WRITE (stdout, 20) TRIM(ncname)
          ifield=0
!
!  Open filter files for read/write.
!
          CALL netcdf_open (ng, iNLM, ncname, 1, FIL(ifile,ng)%ncid)
          IF (exit_flag.ne.NoError) THEN
            IF (Master) WRITE (stdout,50) TRIM(ncname)
            RETURN
          END IF
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
              got_var(idtime)=.TRUE.
              FIL(1,ng)%Vid(idtime) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
              got_var(idFsur)=.TRUE.
              FIL(1,ng)%Vid(idFsur) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
              got_var(idUbar)=.TRUE.
              FIL(1,ng)%Vid(idUbar) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
              got_var(idVbar)=.TRUE.
              FIL(1,ng)%Vid(idVbar) = var_id(i)
# ifdef SOLVE3D
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
              got_var(idUvel)=.TRUE.
              FIL(2,ng)%Vid(idUvel) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
              got_var(idVvel)=.TRUE.
              FIL(3,ng)%Vid(idVvel) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
              got_var(idOvel)=.TRUE.
              FIL(9,ng)%Vid(idOvel) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
              got_var(idWvel)=.TRUE.
              FIL(9,ng)%Vid(idWvel) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
              got_var(idDano)=.TRUE.
              FIL(7,ng)%Vid(idDano) = var_id(i)
#  ifdef LMD_SKPP
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
              got_var(idHsbl)=.TRUE.
              FIL(1,ng)%Vid(idHsbl) = var_id(i)
#  endif
#  ifdef LMD_BKPP
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
              got_var(idHbbl)=.TRUE.
              FIL(1,ng)%Vid(idHbbl) = var_id(i)
#  endif
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
              got_var(idVvis)=.TRUE.
              FIL(6,ng)%Vid(idVvis) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
              got_var(idTdif)=.TRUE.
              FIL(4,ng)%Vid(idTdif) = var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
              got_var(idSdif)=.TRUE.
              FIL(5,ng)%Vid(idSdif) = var_id(i)
#  ifdef ICE_MODEL
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
              got_var(idUice)=.true.
              FIL(10,ng)%Vid(idUice)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
              got_var(idVice)=.true.
              FIL(10,ng)%Vid(idVice)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
              got_var(idAice)=.true.
              FIL(10,ng)%Vid(idAice)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
              got_var(idHice)=.true.
              FIL(10,ng)%Vid(idHice)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
              got_var(idTice)=.true.
              FIL(10,ng)%Vid(idTice)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
              got_var(idHsno)=.true.
              FIL(10,ng)%Vid(idHsno)=var_id(i)
#ifdef MELT_PONDS
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idApond))) THEN
              got_var(idApond)=.true.
              FIL(10,ng)%Vid(idApond)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHpond))) THEN
              got_var(idHpond)=.true.
              FIL(10,ng)%Vid(idHpond)=var_id(i)
#endif
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
              got_var(idAgeice)=.true.
              FIL(10,ng)%Vid(idAgeice)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
              got_var(idIomflx)=.true.
              FIL(10,ng)%Vid(idIomflx)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
              got_var(idTimid)=.true.
              FIL(10,ng)%Vid(idTimid)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
              got_var(idSig11)=.true.
              FIL(10,ng)%Vid(idSig11)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
              got_var(idSig12)=.true.
              FIL(10,ng)%Vid(idSig12)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
              got_var(idSig22)=.true.
              FIL(10,ng)%Vid(idSig22)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
              got_var(idTauiw)=.true.
              FIL(10,ng)%Vid(idTauiw)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
              got_var(idChuiw)=.true.
              FIL(10,ng)%Vid(idChuiw)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
              got_var(idT0mk)=.true.
              FIL(10,ng)%Vid(idT0mk)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
              got_var(idS0mk)=.true.
              FIL(10,ng)%Vid(idS0mk)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
              got_var(idWfr)=.true.
              FIL(10,ng)%Vid(idWfr)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
              got_var(idWai)=.true.
              FIL(10,ng)%Vid(idWai)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
              got_var(idWao)=.true.
              FIL(10,ng)%Vid(idWao)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
              got_var(idWio)=.true.
              FIL(10,ng)%Vid(idWio)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
              got_var(idWro)=.true.
              FIL(10,ng)%Vid(idWro)=var_id(i)
#  endif
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itemp)))) THEN
              got_var(idTsur(itemp))=.TRUE.
              FIL(10,ng)%Vid(idTsur(itemp))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(isalt)))) THEN
              got_var(idTsur(isalt))=.TRUE.
              FIL(10,ng)%Vid(idTsur(isalt))=var_id(i)
#  ifdef BULK_FLUXES
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
              got_var(idLhea)=.TRUE.
              FIL(1,ng)%Vid(idLhea)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
              got_var(idShea)=.TRUE.
              FIL(1,ng)%Vid(idShea)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
              got_var(idLrad)=.TRUE.
              FIL(1,ng)%Vid(idLrad)=var_id(i)
!           ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
!             got_var(idTair)=.TRUE.
!             FIL(1,ng)%Vid(idTair)=var_id(i)
!           ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
!             got_var(idUair)=.TRUE.
!             FIL(1,ng)%Vid(idUair)=var_id(i)
!           ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
!             got_var(idVair)=.TRUE.
!             FIL(1,ng)%Vid(idVair)=var_id(i)
#   ifdef EMINUSP
!           ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
!             got_var(idevap)=.TRUE.
!             FIL(1,ng)%Vid(idevap)=var_id(i)
!           ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
!             got_var(idrain)=.TRUE.
!             FIL(1,ng)%Vid(idrain)=var_id(i)
#   endif
#  endif
#  ifdef SHORTWAVE
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
              got_var(idSrad)=.TRUE.
              FIL(1,ng)%Vid(idSrad)=var_id(i)
#  endif
# endif
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
              got_var(idUsms)=.TRUE.
              FIL(1,ng)%Vid(idUsms)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
              got_var(idVsms)=.TRUE.
              FIL(1,ng)%Vid(idVsms)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
              got_var(idUbms)=.TRUE.
              FIL(1,ng)%Vid(idUbms)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
              got_var(idVbms)=.TRUE.
              FIL(1,ng)%Vid(idVbms)=var_id(i)
            END IF
# ifdef SOLVE3D
            DO itrc=1,NAT
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
                got_var(idTvar(itrc))=.TRUE.
                FIL(8,ng)%Tid(itrc) = var_id(i)
              END IF
            END DO
# endif
          END DO
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Aout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Aout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Aout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Aout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Aout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NAT
          IF (.not.got_var(idTvar(itrc)).and.Aout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))),  &
     &                                  TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano).and.Aout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ICE_MODEL
        IF (.not.got_var(idUice).and.Aout(idUice,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Aout(idVice,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Aout(idAice,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Aout(idHice,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Aout(idTice,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Aout(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef MELT_PONDS
        IF (.not.got_var(idApond).and.Aout(idApond,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idApond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHpond).and.Aout(idHpond,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idHpond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idAgeice).and.Aout(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Aout(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Aout(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Aout(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Aout(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Aout(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Aout(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Aout(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Aout(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Aout(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Aout(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Aout(idWai,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Aout(idWao,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Aout(idWio,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Aout(idWro,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idTsur(itemp)).and.Aout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Aout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Aout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idLhea)), &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Aout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idShea)), &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Aout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idLrad)), &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!       IF (.not.got_var(idTair).and.Aout(idTair,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTair)),          &
!    &                                  TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
!       IF (.not.got_var(idUair).and.Aout(idUair,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUair)),          &
!     &                                 TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
!       IF (.not.got_var(idVair).and.Aout(idVair,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVair)),          &
!     &                                 TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
#   ifdef EMINUSP
!       IF (.not.got_var(idevap).and.Aout(idevap,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idevap)),          &
!    &                                  TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
!       IF (.not.got_var(idrain).and.Aout(idrain,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idrain)),          &
!    &                                  TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
#   endif
#  endif
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Aout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Aout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Aout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Aout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Aout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!       IF (.not.got_var(idUbms).and.Aout(idUbms,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbms)),          &
!    &                                  TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
!       IF (.not.got_var(idVbms).and.Aout(idVbms,ng)) THEN
!         IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbms)),          &
!    &                                  TRIM(ncname)
!         exit_flag=3
!         RETURN
!       END IF
!
!  Set unlimited time record dimension to current value.
!
!     FIL(ifile,ng)%Rindex=RST(ng)%Rindex
!
      END IF

  10  FORMAT (6x,'DEF_FILT - creating filter file:  ',a)
  20  FORMAT (6x,'DEF_FILT - inquiring filter file: ',a)
  30  FORMAT (/,' DEF_FILT - unable to create filter NetCDF file: ',a)
  40  FORMAT (/,' DEF_FILT - unable to find variable: ',a,2x,           &
     &          ' in filter NetCDF file: ',a)
  50  FORMAT (/,' DEF_FILT - unable to open filter NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_filt
#else
      SUBROUTINE def_filt
      RETURN
      END SUBROUTINE def_filt
#endif
