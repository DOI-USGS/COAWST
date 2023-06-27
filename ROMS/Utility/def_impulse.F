#include "cppdefs.h"
      MODULE def_impulse_mod
#ifdef WEAK_CONSTRAINT
!
!git $Id$
!svn $Id: def_impulse.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
!
      USE def_dim_mod,  ONLY : def_dim
      USE def_info_mod, ONLY : def_info
      USE def_var_mod,  ONLY : def_var
      USE strings_mod,  ONLY : FoundError
      USE wrt_info_mod, ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_impulse
      PRIVATE :: def_impulse_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_impulse_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_impulse (ng)
!***********************************************************************
!
!  Imported variable declarations.
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
      SELECT CASE (TLF(ng)%IOtype)
        CASE (io_nf90)
          CALL def_impulse_nf90 (ng)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_impulse_pio (ng)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TLF(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_IMPULSE - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_impulse
!
!***********************************************************************
      SUBROUTINE def_impulse_nf90 (ng)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, nvd3, nvd4
      integer :: recdim, status, varid
      integer :: Fcount

      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)

# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_impulse_nf90"
!
      SourceFile=MyFile
!
!=======================================================================
!  Create a new impulse forcing file.
!=======================================================================
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=TLF(ng)%name
!
      IF (Master) THEN
        IF (LdefTLF(ng)) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
      DEFINE : IF (LdefTLF(ng)) THEN
        CALL netcdf_create (ng, iTLM, TRIM(ncname), TLF(ng)%ncid)
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
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, TLF(ng)%ncid, ncname, 'Nbands',        &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iTLM, TLF(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define TLM/RPM impulse forcing.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idtime),     &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        Vinfo(21)=Vname(6,idZtlf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idZtlf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idZtlf),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifndef SOLVE3D
!
!  Define 2D U-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idUbtf)
        Vinfo( 2)=Vname(2,idUbtf)
        Vinfo( 3)=Vname(3,idUbtf)
        Vinfo(14)=Vname(4,idUbtf)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(21)=Vname(6,idUbtf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbtf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idUbtf),     &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 2D V-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idVbtf)
        Vinfo( 2)=Vname(2,idVbtf)
        Vinfo( 3)=Vname(3,idVbtf)
        Vinfo(14)=Vname(4,idVbtf)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(21)=Vname(6,idVbtf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbtf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idVbtf),     &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# else
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
        Vinfo(21)=Vname(6,idUtlf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUtlf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idUtlf),     &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        Vinfo(21)=Vname(6,idVtlf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVtlf,ng),r8)
        status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Vid(idVtlf),     &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTtlf(itrc))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iTLM, TLF(ng)%ncid, TLF(ng)%Tid(itrc),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iTLM, ncname, TLF(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.  Deactive file
!  creation switch.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iTLM, TLF(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!  Open impulse forcing file for read/write.
!
        IF (TLF(ng)%ncid.eq.-1) THEN
          CALL netcdf_open (ng, iTLM, ncname, 1, TLF(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
            WRITE (stdout,50) TRIM(ncname)
            RETURN
          END IF
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iTLM, ncname,                        &
     &                         ncid = TLF(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iTLM, ncname,                          &
     &                       ncid = TLF(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
# ifndef SOLVE3D

          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbtf))) THEN
            got_var(idUbtf)=.TRUE.
            TLF(ng)%Vid(idUbtf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbtf))) THEN
            got_var(idVbtf)=.TRUE.
            TLF(ng)%Vid(idVbtf)=var_id(i)
# else
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
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZtlf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifndef SOLVE3D
        IF (.not.got_var(idUbtf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbtf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# else
        IF (.not.got_var(idUtlf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVtlf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTtlf(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTtlf(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        TLF(ng)%Rindex=rec_size
        Fcount=TLF(ng)%Fcount
        TLF(ng)%Nrec(Fcount)=rec_size
      END IF QUERY
!

  10  FORMAT (2x,'DEF_IMPULSE_NF90 - creating impulse forcing file,',   &
     &        t56,'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_IMPULSE_NF90 - inquiring impulse forcing file,',  &
     &        t56,'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_IMPULSE_NF90 - unable to create impulse forcing', &
     &        ' NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_IMPULSE_NF90 - unable to open norm NetCDF',       &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_IMPULSE_NF90 - unable to find variable: ',a,2x,   &
     &        ' in impulse forcing NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_impulse_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_impulse_pio (ng)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, nvd3, nvd4
      integer :: recdim, status
      integer :: Fcount

      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)

# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_impulse_pio"
!
      SourceFile=MyFile
!
!=======================================================================
!  Create a new impulse forcing file.
!=======================================================================
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=TLF(ng)%name
!
      IF (Master) THEN
        IF (LdefTLF(ng)) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
      DEFINE : IF (LdefTLF(ng)) THEN
        CALL pio_netcdf_create (ng, iTLM, TRIM(ncname), TLF(ng)%pioFile)
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
        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xi_rho',     &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xi_u',       &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xi_v',       &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xi_psi',     &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'eta_rho',    &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'eta_u',      &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'eta_v',      &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'eta_psi',    &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xy_rho',     &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xy_u',       &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xy_v',       &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xyz_rho',    &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xyz_u',      &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xyz_v',      &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xyz_w',      &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 's_rho',      &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 's_w',        &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'tracer',     &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'NST',        &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'Nbed',       &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'xybed',      &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, TLF(ng)%pioFile, ncname, 'Nbands',     &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'Nphy',       &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'Nbac',       &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'Ndom',       &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'Nfec',       &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'boundary',   &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname, 'Nstate',     &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iTLM, TLF(ng)%pioFile, ncname,               &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
        CALL def_info (ng, iTLM, TLF(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define TLM/RPM impulse forcing.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        TLF(ng)%pioVar(idtime)%dkind=PIO_TOUT
        TLF(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, iTLM, TLF(ng)%pioFile,                       &
     &                 TLF(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        Vinfo(21)=Vname(6,idZtlf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idZtlf,ng),r8)
        TLF(ng)%pioVar(idZtlf)%dkind=PIO_FOUT
        TLF(ng)%pioVar(idZtlf)%gtype=r2dvar
!
        status=def_var(ng, iTLM, TLF(ng)%pioFile,                       &
     &                 TLF(ng)%pioVar(idZtlf)%vd,                       &
     &                 PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifndef SOLVE3D
!
!  Define 2D U-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idUbtf)
        Vinfo( 2)=Vname(2,idUbtf)
        Vinfo( 3)=Vname(3,idUbtf)
        Vinfo(14)=Vname(4,idUbtf)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(21)=Vname(6,idUbtf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbtf,ng),r8)
        TLF(ng)%pioVar(idUbtf)%dkind=PIO_FOUT
        TLF(ng)%pioVar(idUbtf)%gtype=u2dvar
!
        status=def_var(ng, iTLM, TLF(ng)%pioFile,                       &
     &                 TLF(ng)%pioVar(idUbtf)%vd,                       &
     &                 PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 2D V-momentum component impulse forcing.
!
        Vinfo( 1)=Vname(1,idVbtf)
        Vinfo( 2)=Vname(2,idVbtf)
        Vinfo( 3)=Vname(3,idVbtf)
        Vinfo(14)=Vname(4,idVbtf)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(21)=Vname(6,idVbtf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbtf,ng),r8)
        TLF(ng)%pioVar(idVbtf)%dkind=PIO_FOUT
        TLF(ng)%pioVar(idVbtf)%gtype=v2dvar
!
        status=def_var(ng, iTLM, TLF(ng)%pioFile,                       &
     &                 TLF(ng)%pioVar(idVbtf)%vd,                       &
     &                 PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# else
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
        Vinfo(21)=Vname(6,idUtlf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUtlf,ng),r8)
        TLF(ng)%pioVar(idUtlf)%dkind=PIO_FOUT
        TLF(ng)%pioVar(idUtlf)%gtype=u3dvar
!
        status=def_var(ng, iTLM, TLF(ng)%pioFile,                       &
     &                 TLF(ng)%pioVar(idUtlf)%vd,                       &
     &                 PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        Vinfo(21)=Vname(6,idVtlf)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVtlf,ng),r8)
        TLF(ng)%pioVar(idVtlf)%dkind=PIO_FOUT
        TLF(ng)%pioVar(idVtlf)%gtype=v3dvar
!
        status=def_var(ng, iTLM, TLF(ng)%pioFile,                       &
     &                 TLF(ng)%pioVar(idVtlf)%vd,                       &
     &                 PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTtlf(itrc))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          TLF(ng)%pioTrc(itrc)%dkind=PIO_FOUT
          TLF(ng)%pioTrc(itrc)%gtype=r3dvar
!
          status=def_var(ng, iTLM, TLF(ng)%pioFile,                     &
     &                   TLF(ng)%pioTrc(itrc)%vd,                       &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, iTLM, ncname, TLF(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.  Deactive file
!  creation switch.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iTLM, TLF(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!  Open impulse forcing file for read/write.
!
        IF (TLF(ng)%pioFile%fh.eq.-1) THEN
          CALL pio_netcdf_open (ng, iTLM, ncname, 1, TLF(ng)%pioFile)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
            WRITE (stdout,50) TRIM(ncname)
            RETURN
          END IF
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, iTLM, ncname,                    &
     &                             pioFile = TLF(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, iTLM, ncname,                      &
     &                           pioFile = TLF(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
            TLF(ng)%pioVar(idtime)%vd=var_desc(i)
            TLF(ng)%pioVar(idtime)%dkind=PIO_TOUT
            TLF(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZtlf))) THEN
            got_var(idZtlf)=.TRUE.
            TLF(ng)%pioVar(idZtlf)%vd=var_desc(i)
            TLF(ng)%pioVar(idZtlf)%dkind=PIO_FOUT
            TLF(ng)%pioVar(idZtlf)%gtype=r2dvar
# ifndef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbtf))) THEN
            got_var(idUbtf)=.TRUE.
            TLF(ng)%pioVar(idUbtf)%vd=var_desc(i)
            TLF(ng)%pioVar(idUbtf)%dkind=PIO_FOUT
            TLF(ng)%pioVar(idUbtf)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbtf))) THEN
            got_var(idVbtf)=.TRUE.
            TLF(ng)%pioVar(idVbtf)%vd=var_desc(i)
            TLF(ng)%pioVar(idVbtf)%dkind=PIO_FOUT
            TLF(ng)%pioVar(idVbtf)%gtype=v2dvar
# else
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUtlf))) THEN
            got_var(idUtlf)=.TRUE.
            TLF(ng)%pioVar(idUtlf)%vd=var_desc(i)
            TLF(ng)%pioVar(idUtlf)%dkind=PIO_FOUT
            TLF(ng)%pioVar(idUtlf)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVtlf))) THEN
            got_var(idVtlf)=.TRUE.
            TLF(ng)%pioVar(idVtlf)%vd=var_desc(i)
            TLF(ng)%pioVar(idVtlf)%dkind=PIO_FOUT
            TLF(ng)%pioVar(idVtlf)%gtype=v3dvar
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTtlf(itrc)))) THEN
             got_var(idTtlf(itrc))=.TRUE.
             TLF(ng)%pioTrc(itrc)%vd=var_desc(i)
             TLF(ng)%pioTrc(itrc)%dkind=PIO_FOUT
             TLF(ng)%pioTrc(itrc)%gtype=r3dvar
            END IF
          END DO
# endif
        END DO
!
!  Check if impulse forcing variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZtlf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifndef SOLVE3D
        IF (.not.got_var(idUbtf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbtf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# else
        IF (.not.got_var(idUtlf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVtlf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTtlf(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTtlf(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        TLF(ng)%Rindex=rec_size
        Fcount=TLF(ng)%Fcount
        TLF(ng)%Nrec(Fcount)=rec_size
      END IF QUERY
!
  10  FORMAT (2x,'DEF_IMPULSE_PIO  - creating impulse forcing file,',   &
     &        t56,'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_IMPULSE_PIO  - inquiring impulse forcing file,',  &
     &        t56,'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_IMPULSE_PIO - unable to create impulse forcing',  &
     &        ' NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_IMPULSE_PIO - unable to open norm NetCDF',        &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_IMPULSE_PIO - unable to find variable: ',a,2x,    &
     &        ' in impulse forcing NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_impulse_pio
# endif
#endif
      END MODULE def_impulse_mod
