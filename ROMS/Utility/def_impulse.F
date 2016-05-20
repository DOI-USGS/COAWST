#include "cppdefs.h"
#ifdef WEAK_CONSTRAINT
      SUBROUTINE def_impulse (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates tangent linear and representer models impulse  !
!  forcing NetCDF file used for weak constraint 4DVAR. It defines its  !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BIOLOGY
      USE mod_biology
# endif
      USE mod_fourdvar
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
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)

      integer, parameter :: Natt = 25

      integer :: i, j, nvd3, nvd4
      integer :: recdim, status, varid

      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_dim

# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
# endif

      real(r8) :: Aval(6)

      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_impulse.F'
!
!=======================================================================
!  Create a new impulse forcing file.
!=======================================================================
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=TLF(ng)%name
!
      DEFINE : IF (LdefTLF(ng)) THEN
        CALL netcdf_create (ng, iTLM, TRIM(ncname), TLF(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
#  endif

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef ECOSIM
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
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
        t2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
#  endif
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered u-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
#  endif
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
#  endif
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        TLF(ng)%Rindex=0
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
        CALL def_info (ng, iTLM, TLF(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define TLM/RPM impulse forcing.
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
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface impulse forcing.
!
        Vinfo( 1)=Vname(1,idZtlf)
        Vinfo( 2)=Vname(2,idZtlf)
        Vinfo( 3)=Vname(3,idZtlf)
        Vinfo(14)=Vname(4,idZtlf)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idZtlf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idZtlf),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D U-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idUbtf)
        Vinfo( 2)=Vname(2,idUbtf)
        Vinfo( 3)=Vname(3,idUbtf)
        Vinfo(14)=Vname(4,idUbtf)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbtf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idUbtf),     &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D V-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idVbtf)
        Vinfo( 2)=Vname(2,idVbtf)
        Vinfo( 3)=Vname(3,idVbtf)
        Vinfo(14)=Vname(4,idVbtf)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbtf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idVbtf),     &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Define 3D U-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idUtlf)
        Vinfo( 2)=Vname(2,idUtlf)
        Vinfo( 3)=Vname(3,idUtlf)
        Vinfo(14)=Vname(4,idUtlf)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUtlf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idUtlf),     &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D V-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idVtlf)
        Vinfo( 2)=Vname(2,idVtlf)
        Vinfo( 3)=Vname(3,idVtlf)
        Vinfo(14)=Vname(4,idVtlf)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVtlf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idVtlf),     &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define tracer type impulse forcing variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTtlf(itrc))
          Vinfo( 2)=Vname(2,idTtlf(itrc))
          Vinfo( 3)=Vname(3,idTtlf(itrc))
          Vinfo(14)=Vname(4,idTtlf(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Tid(itrc),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iTLM, ncname, TLF(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.  Deactive file
!  creation switch.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iTLM, TLF(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
        LdefTLF(ng)=.FALSE.
      END IF DEFINE
!
!=======================================================================
!  Open an existing impulse forcing file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      QUERY : IF (.not.LdefTLF(ng)) THEN
        ncname=TLF(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iTLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iTLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open impulse forcing file for read/write.
!
        CALL netcdf_open (ng, iTLM, ncname, 1, TLF(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,30) TRIM(ncname)
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
!  impulse forcing variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            TLF(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZtlf))) THEN
            got_var(idZtlf)=.TRUE.
            TLF(ng)%Vid(idZtlf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbtf))) THEN
            got_var(idUbtf)=.TRUE.
            TLF(ng)%Vid(idUbtf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbtf))) THEN
            got_var(idVbtf)=.TRUE.
            TLF(ng)%Vid(idVbtf)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUtlf))) THEN
            got_var(idUtlf)=.TRUE.
            TLF(ng)%Vid(idUtlf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVtlf))) THEN
            got_var(idVtlf)=.TRUE.
            TLF(ng)%Vid(idVtlf)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTtlf(itrc)))) THEN
             got_var(idTtlf(itrc))=.TRUE.
             TLF(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if impulse forcing variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZtlf)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbtf)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbtf)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUtlf)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVtlf)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTtlf(itrc))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTtlf(itrc))),  &
     &                                   TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        TLF(ng)%Rindex=rec_size
      END IF QUERY
!
  10  FORMAT (/,' DEF_IMPULSE - unable to create impulse forcing',      &
     &          ' NetCDF file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_IMPULSE - unable to open norm NetCDF file: ',a)
  40  FORMAT (/,' DEF_IMPULSE - unable to find variable: ',a,2x,        &
     &        ' in impulse forcing NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_impulse
#else
      SUBROUTINE def_impulse
      RETURN
      END SUBROUTINE def_impulse
#endif
