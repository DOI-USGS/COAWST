#include "cppdefs.h"
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
      SUBROUTINE def_tides (ng, ldef)
!
!svn $Id: def_tides.F 873 2017-10-05 20:27:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates creates least-squares harmomics variables for  !
!  detiding.                                                           !
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
      USE strings_mod, ONLY : FoundError
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
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
      ncname=HAR(ng)%name
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
!  Create a new least-squares detide harmonics file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), HAR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) THEN
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
        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# ifdef SOLVE3D
        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'tide_period',   &
     &                 NTC(ng), DimIDs(13))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, iNLM, HAR(ng)%ncid, ncname, 'harmonics',     &
     &                 2*NTC(ng)+1, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for tide harmonics variables.
!
        tharm(1)=DimIDs(13)
        tharm(2)=DimIDs(13)
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, HAR(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!-----------------------------------------------------------------------
!  Define least-squares detide harmonic variables.
!-----------------------------------------------------------------------
!
!  Define number of time-accumulated harmonics.
!
        Vinfo( 1)='Hcount'
        Vinfo( 2)='number of time-accumulated tide harmonics'
        status=def_var(ng, iNLM, HAR(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define model time for accumulated tide harmonic fields.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'accumulated harmonics ',             &
     &                            TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idtime),     &
     &                 NF_FRST, 1, (/0/), Aval, Vinfo, ncname,          &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define tidal period.
!
        Vinfo( 1)=Vname(1,idTper)
        Vinfo( 2)=Vname(2,idTper)
        Vinfo( 3)=Vname(3,idTper)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idTper),     &
     &                 NF_FRST, 1, (/tharm(1)/), Aval, Vinfo, ncname,   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define time-accumulated COS(omega(k)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idCosW)
        Vinfo( 2)=Vname(2,idCosW)
        Vinfo( 3)=Vname(3,idCosW)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idCosW),     &
     &                 NF_FRST, 1, (/tharm(1)/), Aval, Vinfo, ncname,   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define time-accumulated SIN(omega(k)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idSinW)
        Vinfo( 2)=Vname(2,idSinW)
        Vinfo( 3)=Vname(3,idSinW)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idSinW),     &
     &                 NF_FRST, 1, (/tharm(1)/), Aval, Vinfo, ncname,   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                   __FILE__)) RETURN
!
!  Define time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idCos2)
        Vinfo( 2)=Vname(2,idCos2)
        Vinfo( 3)=Vname(3,idCos2)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idCos2),     &
     &                 NF_FRST, 2, tharm, Aval, Vinfo, ncname,          &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idSin2)
        Vinfo( 2)=Vname(2,idSin2)
        Vinfo( 3)=Vname(3,idSin2)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idSin2),     &
     &                 NF_FRST, 2, tharm, Aval, Vinfo, ncname,          &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
        Vinfo( 1)=Vname(1,idSWCW)
        Vinfo( 2)=Vname(2,idSWCW)
        Vinfo( 3)=Vname(3,idSWCW)
        status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idSWCW),     &
     &                 NF_FRST, 2, tharm, Aval, Vinfo, ncname,          &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define free-surface time-accumulated tide harmonics.
!
        IF (Aout(idFsuD,ng)) THEN
          Vinfo( 1)=Vname(1,idFsuH)
          Vinfo( 2)=Vname(2,idFsuH)
          Vinfo( 3)=Vname(3,idFsuH)
          Vinfo(14)=Vname(4,idFsuH)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsuH,ng),r8)
          status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idFsuH),   &
     &                   NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D u-momentum time-accumulated tide harmonics.
!
        IF (Aout(idu2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dH)
          Vinfo( 2)=Vname(2,idu2dH)
          Vinfo( 3)=Vname(3,idu2dH)
          Vinfo(14)=Vname(4,idu2dH)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dH,ng),r8)
          status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idu2dH),   &
     &                   NF_FRST, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D v-momentum time-accumulated tide harmonics.
!
        IF (Aout(idv2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dH)
          Vinfo( 2)=Vname(2,idv2dH)
          Vinfo( 3)=Vname(3,idv2dH)
          Vinfo(14)=Vname(4,idv2dH)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dH,ng),r8)
          status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idv2dH),   &
     &                   NF_FRST, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D u-momentum time-accumulated tide harmonics.
!
        IF (Aout(idu3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dH)
          Vinfo( 2)=Vname(2,idu3dH)
          Vinfo( 3)=Vname(3,idu3dH)
          Vinfo(14)=Vname(4,idu3dH)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dH,ng),r8)
          status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idu3dH),   &
     &                   NF_FRST, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 3D v-momentum time-accumulated tide harmonics.
!
        IF (.not.got_var(idv3dH).and.Aout(idv3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dH)
          Vinfo( 2)=Vname(2,idv3dH)
          Vinfo( 3)=Vname(3,idv3dH)
          Vinfo(14)=Vname(4,idv3dH)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dH,ng),r8)
          status=def_var(ng, iNLM, HAR(ng)%ncid, HAR(ng)%Vid(idv3dH),   &
     &                   NF_FRST, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define temperaturea and salinity time-accumulated tide harmonics.
!
        DO itrc=1,NAT
          IF (Aout(idTrcD(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTrcH(itrc))
            Vinfo( 2)=Vname(2,idTrcH(itrc))
            Vinfo( 3)=Vname(3,idTrcH(itrc))
            Vinfo(14)=Vname(4,idTrcH(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idv3dH,ng),r8)
            status=def_var(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idTrcH(itrc)),                   &
     &                     NF_FRST, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, HAR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, HAR(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing detide harmonics file, check its contents, and
!  prepare for updating data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=HAR(ng)%name
!
!  Open detide harmonics file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, HAR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) THEN
          WRITE (stdout,40) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname,                        &
     &                         ncid = HAR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname,                          &
     &                       ncid = HAR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  detide harmomics variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            HAR(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTper))) THEN
            got_var(idTper)=.TRUE.
            HAR(ng)%Vid(idTper)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCosW))) THEN
            got_var(idCosW)=.TRUE.
            HAR(ng)%Vid(idCosW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSinW))) THEN
            got_var(idSinW)=.TRUE.
            HAR(ng)%Vid(idSinW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCos2))) THEN
            got_var(idCos2)=.TRUE.
            HAR(ng)%Vid(idCos2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSin2))) THEN
            got_var(idSin2)=.TRUE.
            HAR(ng)%Vid(idSin2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSWCW))) THEN
            got_var(idSWCW)=.TRUE.
            HAR(ng)%Vid(idSWCW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsuH))) THEN
            got_var(idFsuH)=.TRUE.
            HAR(ng)%Vid(idFsuH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dH))) THEN
            got_var(idu2dH)=.TRUE.
            HAR(ng)%Vid(idu2dH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dH))) THEN
            got_var(idv2dH)=.TRUE.
            HAR(ng)%Vid(idv2dH)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dH))) THEN
            got_var(idu3dH)=.TRUE.
            HAR(ng)%Vid(idu3dH)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dH))) THEN
            got_var(idv3dH)=.TRUE.
            HAR(ng)%Vid(idv3dH)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTrcH(itrc)))) THEN
              got_var(idTrcH(itrc))=.TRUE.
              HAR(ng)%Vid(idTrcH(itrc))=var_id(i)
            END IF
          END DO
# endif

        END DO
!
!  Check if detide harmonics variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTper)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTper)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCosW)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idCosW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSinW)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSinW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCos2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idCos2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSin2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSin2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWCW)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSWCW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsuH).and.Aout(idFsuD,ng)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idFsuH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dH).and.Aout(idu2dD,ng)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idu2dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dH).and.Aout(idv2dD,ng)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idv2dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idu3dH).and.Aout(idu3dD,ng)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idu3dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dH).and.Aout(idv3dD,ng)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idv3dH)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NAT
          IF (.not.got_var(idTrcH(itrc)).and.Aout(idTrcD(itrc),ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTrcH(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
      END IF QUERY
!
  10  FORMAT (4x,'DEF_TIDES   - creating  harmonics', t43,              &
     &        ' file, Grid ',i2.2,': ', a)
  20  FORMAT (4x,'DEF_TIDES   - inquiring harmonics', t43,              &
     &        ' file, Grid ',i2.2,': ', a)
  30  FORMAT (/,' DEF_TIDES - unable to create harmonics NetCDF file: ',&
     &        a)
  40  FORMAT (/,' DEF_TIDES - unable to open harmonics NetCDF file: ',a)
  50  FORMAT (/,' DEF_TIDES - unable to find variable: ',a,2x,          &
     &        ' in detide harmonics NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_tides
#else
      SUBROUTINE def_tides
      RETURN
      END SUBROUTINE def_tides
#endif
