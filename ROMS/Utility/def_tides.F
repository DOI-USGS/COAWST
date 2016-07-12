#include "cppdefs.h"
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      SUBROUTINE def_tides (ng, ldef)
!
!svn $Id: def_tides.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates creates or updates tidal forcing NetCDF file.  !
!  If a file already exists and appropriate,  it will define several   !
!  additional fields needed to detide time-averaged fields.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_tides
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
      logical :: got_var(NV)

      logical :: Ldefine = .FALSE.

      integer, parameter :: Natt = 25

      integer :: i, itrc, j, nvd3, nvd4
      integer :: status, varid

      integer :: DimIDs(32), tharm(2), t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
# endif
      integer :: Vsize(4)

      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_tides.F'
!
!=======================================================================
!  Open existing tidal forcing file and define new variables.
!=======================================================================
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=TIDE(ng)%name
!
      DEFINE : IF (ldef) THEN
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
!  Check if tidal harmonics variables have been already defined.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.Vname(1,idCosW)) THEN
            got_var(idCosW)=.TRUE.
            TIDE(ng)%Vid(idCosW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idSinW)) THEN
            got_var(idSinW)=.TRUE.
            TIDE(ng)%Vid(idSinW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idCos2)) THEN
            got_var(idCos2)=.TRUE.
            TIDE(ng)%Vid(idCos2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idSin2)) THEN
            got_var(idSin2)=.TRUE.
            TIDE(ng)%Vid(idSin2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idSWCW)) THEN
            got_var(idSWCW)=.TRUE.
            TIDE(ng)%Vid(idSWCW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idFsuH)) THEN
            got_var(idFsuH)=.TRUE.
            TIDE(ng)%Vid(idFsuH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idu2dH)) THEN
            got_var(idu2dH)=.TRUE.
            TIDE(ng)%Vid(idu2dH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idv2dH)) THEN
            got_var(idv2dH)=.TRUE.
            TIDE(ng)%Vid(idv2dH)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idu3dH)) THEN
            got_var(idu3dH)=.TRUE.
            TIDE(ng)%Vid(idu3dH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.Vname(1,idv3dH)) THEN
            got_var(idv3dH)=.TRUE.
            TIDE(ng)%Vid(idv3dH)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.Vname(1,idTrcH(itrc))) THEN
              got_var(idTrcH(itrc))=.TRUE.
              TIDE(ng)%Vid(idTrcH(itrc))=var_id(i)
            END IF
          END DO
# endif
        END DO
!
        IF (.not.got_var(idCosW)) Ldefine=.TRUE.
        IF (.not.got_var(idSinW)) Ldefine=.TRUE.
        IF (.not.got_var(idCos2)) Ldefine=.TRUE.
        IF (.not.got_var(idSin2)) Ldefine=.TRUE.
        IF (.not.got_var(idSWCW)) Ldefine=.TRUE.
        IF (.not.got_var(idFsuH).and.Aout(idFsuD,ng)) Ldefine=.TRUE.
        IF (.not.got_var(idu2dH).and.Aout(idu2dD,ng)) Ldefine=.TRUE.
        IF (.not.got_var(idv2dH).and.Aout(idv2dD,ng)) Ldefine=.TRUE.
# ifdef SOLVE3D
        IF (.not.got_var(idu3dH).and.Aout(idu3dD,ng)) Ldefine=.TRUE.
        IF (.not.got_var(idv3dH).and.Aout(idv3dD,ng)) Ldefine=.TRUE.
        DO itrc=1,NAT
          IF (.not.got_var(idTrcH(itrc)).and.                           &
     &        Aout(idTrcD(itrc),ng)) Ldefine=.TRUE.
        END DO
# endif
!
!  Open existing tidal forcing file.
!
        IF (Ldefine) THEN
          CALL netcdf_open (ng, iNLM, ncname, 1, TIDE(ng)%ncid)
          IF (exit_flag.ne.NoError) THEN
            WRITE (stdout,10) TRIM(ncname)
            RETURN
          END IF
        END IF
!
!  Put existing file into define mode so new variables can be added.
!
        IF (Ldefine) THEN
          CALL netcdf_redef (ng, iNLM, ncname, TIDE(ng)%ncid)
          IF (exit_flag.ne.NoError) THEN
            IF (Master) WRITE (stdout,20) TRIM(ncname)
            RETURN
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Inquire and define (if needed) file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0

        IF (Ldefine) THEN
          DO i=1,n_dim
            SELECT CASE (TRIM(ADJUSTL(dim_name(i))))
              CASE ('xi_rho')
                DimIDs( 1)=dim_id(i)
              CASE ('xi_u')
                DimIDs( 2)=dim_id(i)
              CASE ('xi_v')
                DimIDs( 3)=dim_id(i)
              CASE ('eta_rho')
                DimIDs( 5)=dim_id(i)
              CASE ('eta_u')
                DimIDs( 6)=dim_id(i)
              CASE ('eta_v')
                DimIDs( 7)=dim_id(i)
              CASE ('tide_period')
                DimIDs( 8)=dim_id(i)
# ifdef SOLVE3D
              CASE ('s_rho')
                DimIDs( 9)=dim_id(i)
# endif
              CASE ('harmonics')
                DimIDs(12)=dim_id(i)
            END SELECT
          END DO
!
!  Define needed dimensions, if any.
!
# ifdef SOLVE3D
          IF (DimIDs( 9).eq.0) THEN
            status=def_dim(ng, iNLM, TIDE(ng)%ncid, ncname, 's_rho',    &
     &                     N(ng), DimIDs( 9))
            IF (exit_flag.ne.NoError) RETURN
          END IF
# endif
          IF (DimIDs(12).eq.0) THEN
            status=def_dim(ng, iNLM, TIDE(ng)%ncid, ncname, 'harmonics',&
     &                     2*NTC(ng)+1, DimIDs(12))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for tide harmonics variables.
!
        tharm(1)=DimIDs(8)
        tharm(2)=DimIDs(8)
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
# endif
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
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define number of time-accumulated harmonics.
!
        Vinfo( 1)='Hcount'
        Vinfo( 2)='number of time-accumulated tide harmonics'
        status=def_var(ng, iNLM, TIDE(ng)%ncid, varid, nf90_int,        &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model time for accumulated tide harmonic fields.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,1x,a)') 'accumulated harmonics',           &
     &                               TRIM(Vname(2,idtime))
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idtime),   &
     &                 NF_FRST, 1, (/0/), Aval, Vinfo, ncname,          &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define time-accumulated COS(omega(k)*t) harmonics.
!
        IF (.not.got_var(idCosW)) THEN
          Vinfo( 1)=Vname(1,idCosW)
          Vinfo( 2)=Vname(2,idCosW)
          Vinfo( 3)=Vname(3,idCosW)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idCosW), &
     &                   NF_FRST, 1, (/tharm(1)/), Aval, Vinfo, ncname, &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-accumulated SIN(omega(k)*t) harmonics.
!
        IF (.not.got_var(idSinW)) THEN
          Vinfo( 1)=Vname(1,idSinW)
          Vinfo( 2)=Vname(2,idSinW)
          Vinfo( 3)=Vname(3,idSinW)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idSinW), &
     &                   NF_FRST, 1, (/tharm(1)/), Aval, Vinfo, ncname, &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
        IF (.not.got_var(idCos2)) THEN
          Vinfo( 1)=Vname(1,idCos2)
          Vinfo( 2)=Vname(2,idCos2)
          Vinfo( 3)=Vname(3,idCos2)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idCos2), &
     &                   NF_FRST, 2, tharm, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
        IF (.not.got_var(idSin2)) THEN
          Vinfo( 1)=Vname(1,idSin2)
          Vinfo( 2)=Vname(2,idSin2)
          Vinfo( 3)=Vname(3,idSin2)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idSin2), &
     &                   NF_FRST, 2, tharm, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
        IF (.not.got_var(idSWCW)) THEN
          Vinfo( 1)=Vname(1,idSWCW)
          Vinfo( 2)=Vname(2,idSWCW)
          Vinfo( 3)=Vname(3,idSWCW)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idSWCW), &
     &                   NF_FRST, 2, tharm, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define free-surface time-accumulated tide harmonics.
!
        IF (.not.got_var(idFsuH).and.Aout(idFsuD,ng)) THEN
          Vinfo( 1)=Vname(1,idFsuH)
          Vinfo( 2)=Vname(2,idFsuH)
          Vinfo( 3)=Vname(3,idFsuH)
          Vinfo(14)=Vname(4,idFsuH)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsuH,ng),r8)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idFsuH), &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D u-momentum time-accumulated tide harmonics.
!
        IF (.not.got_var(idu2dH).and.Aout(idu2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dH)
          Vinfo( 2)=Vname(2,idu2dH)
          Vinfo( 3)=Vname(3,idu2dH)
          Vinfo(14)=Vname(4,idu2dH)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dH,ng),r8)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idu2dH), &
     &                   NF_FRST, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D v-momentum time-accumulated tide harmonics.
!
        IF (.not.got_var(idv2dH).and.Aout(idv2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dH)
          Vinfo( 2)=Vname(2,idv2dH)
          Vinfo( 3)=Vname(3,idv2dH)
          Vinfo(14)=Vname(4,idv2dH)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dH,ng),r8)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idv2dH), &
     &                   NF_FRST, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D u-momentum time-accumulated tide harmonics.
!
        IF (.not.got_var(idu3dH).and.Aout(idu3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dH)
          Vinfo( 2)=Vname(2,idu3dH)
          Vinfo( 3)=Vname(3,idu3dH)
          Vinfo(14)=Vname(4,idu3dH)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dH,ng),r8)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idu3dH), &
     &                   NF_FRST, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D v-momentum time-accumulated tide harmonics.
!
        IF (.not.got_var(idv3dH).and.Aout(idv3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dH)
          Vinfo( 2)=Vname(2,idv3dH)
          Vinfo( 3)=Vname(3,idv3dH)
          Vinfo(14)=Vname(4,idv3dH)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dH,ng),r8)
          status=def_var(ng, iNLM, TIDE(ng)%ncid, TIDE(ng)%Vid(idv3dH), &
     &                   NF_FRST, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define temperaturea and salinity time-accumulated tide harmonics.
!
        DO itrc=1,NAT
          IF (.not.got_var(idTrcH(itrc)).and.Aout(idTrcD(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTrcH(itrc))
            Vinfo( 2)=Vname(2,idTrcH(itrc))
            Vinfo( 3)=Vname(3,idTrcH(itrc))
            Vinfo(14)=Vname(4,idTrcH(itrc))
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idv3dH,ng),r8)
            status=def_var(ng, iNLM, TIDE(ng)%ncid,                     &
     &                     TIDE(ng)%Vid(idTrcH(itrc)),                  &
     &                     NF_FRST, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, TIDE(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
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
!  Open tidal forcing file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, TIDE(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,10) TRIM(ncname)
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
!  tide variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            TIDE(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCosW))) THEN
            got_var(idCosW)=.TRUE.
            TIDE(ng)%Vid(idCosW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSinW))) THEN
            got_var(idSinW)=.TRUE.
            TIDE(ng)%Vid(idSinW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCos2))) THEN
            got_var(idCos2)=.TRUE.
            TIDE(ng)%Vid(idCos2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSin2))) THEN
            got_var(idSin2)=.TRUE.
            TIDE(ng)%Vid(idSin2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSWCW))) THEN
            got_var(idSWCW)=.TRUE.
            TIDE(ng)%Vid(idSWCW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsuH))) THEN
            got_var(idFsuH)=.TRUE.
            TIDE(ng)%Vid(idFsuH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dH))) THEN
            got_var(idu2dH)=.TRUE.
            TIDE(ng)%Vid(idu2dH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dH))) THEN
            got_var(idv2dH)=.TRUE.
            TIDE(ng)%Vid(idv2dH)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dH))) THEN
            got_var(idu3dH)=.TRUE.
            TIDE(ng)%Vid(idu3dH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dH))) THEN
            got_var(idv3dH)=.TRUE.
            TIDE(ng)%Vid(idv3dH)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTrcH(itrc)))) THEN
              got_var(idTrcH(itrc))=.TRUE.
              TIDE(ng)%Vid(idTrcH(itrc))=var_id(i)
            END IF
          END DO
# endif

        END DO
!
!  Check if tide variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCosW)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idCosW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSinW)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSinW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCos2)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idCos2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSin2)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSin2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWCW)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSWCW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsuH).and.Aout(idFsuD,ng)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idFsuH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dH).and.Aout(idu2dD,ng)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idu2dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dH).and.Aout(idv2dD,ng)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idv2dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idu3dH).and.Aout(idu3dD,ng)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idu3dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dH).and.Aout(idv3dD,ng)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idv3dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NAT
          IF (.not.got_var(idTrcH(itrc)).and.Aout(idTrcD(itrc),ng)) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTrcH(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
      END IF QUERY
!
  10  FORMAT (/,' DEF_TIDES - unable to open tide NetCDF file: ',a)
  20  FORMAT (/,' DEF_TIDES - unable to put in define mode tide',       &
     &        ' NetCDF file: ',a)
  30  FORMAT (/,' DEF_TIDES - unable to find variable: ',a,2x,          &
     &        ' in tide NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_tides
#else
      SUBROUTINE def_tides
      RETURN
      END SUBROUTINE def_tides
#endif
