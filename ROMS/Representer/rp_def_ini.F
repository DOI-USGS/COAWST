#include "cppdefs.h"
      MODULE rp_def_ini_mod
#if defined TL_IOMS && defined FOUR_DVAR
!
!git $Id$
!svn $Id: rp_def_ini.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module opens existing representer model initial conditions     !
!  file and either define new new variables or inquire its contents.   !
!  It either use the standard NetCDF library or the Parallel-IO (PIO)  !
!  library.                                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE def_dim_mod, ONLY : def_dim
      USE def_var_mod, ONLY : def_var
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: rp_def_ini
      PRIVATE :: rp_def_ini_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: rp_def_ini_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_def_ini (ng)
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
      SELECT CASE (IRP(ng)%IOtype)
        CASE (io_nf90)
          CALL rp_def_ini_nf90 (ng)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL rp_def_ini_pio (ng)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IRP(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' RP_DEF_INI - Illegal output type, io_type = ',i0,       &
     &        /,14x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE rp_def_ini
!
!***********************************************************************
      SUBROUTINE rp_def_ini_nf90 (ng)
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
      logical :: Ldefine = .FALSE.
# ifdef ADJUST_BOUNDARY
      logical :: got_IorJ, got_boundary, got_obc_adjust
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      logical :: got_frc_adjust
# endif
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, status
      integer :: Fcount
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
# endif
      integer :: DimIDs(nDimID)
# ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
#  ifdef SOLVE3D
      integer :: t3dobc(5)
#  endif
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: t4dfrc(4), u4dfrc(4), v4dfrc(4)
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", rp_def_ini_nf90"
!
      SourceFile=MyFile

# if defined ADJUST_BOUNDARY || \
     defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
!
!=======================================================================
!  Open existing representer model initial conditions and define new
!  variables.
!=======================================================================
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=IRP(ng)%name
!
!  Open representer model initialization file for read/write.
!
      CALL netcdf_open (ng, iRPM, ncname, 1, IRP(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, iRPM, ncname,                          &
     &                       ncid = IRP(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, iRPM, ncname,                            &
     &                     ncid = IRP(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Check if surface forcing variables have been already defined.
!
      DO i=1,NV
        got_var(i)=.FALSE.
      END DO
!
      DO i=1,n_var
#  ifdef ADJUST_BOUNDARY
        IF (TRIM(var_name(i)).eq.                                       &
     &      TRIM(Vname(1,idSbry(isFsur)))) THEN
          got_var(idSbry(isFsur))=.TRUE.
          IRP(ng)%Vid(idSbry(isFsur))=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isUbar)))) THEN
          got_var(idSbry(isUbar))=.TRUE.
          IRP(ng)%Vid(idSbry(isUbar))=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isVbar)))) THEN
          got_var(idSbry(isVbar))=.TRUE.
          IRP(ng)%Vid(idSbry(isVbar))=var_id(i)
#   ifdef SOLVE3D
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isUvel)))) THEN
          got_var(idSbry(isUvel))=.TRUE.
          IRP(ng)%Vid(idSbry(isUvel))=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isVvel)))) THEN
          got_var(idSbry(isVvel))=.TRUE.
          IRP(ng)%Vid(idSbry(isVvel))=var_id(i)
#   endif
        END IF
#   ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (TRIM(var_name(i)).eq.                                     &
     &        TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
            got_var(idSbry(isTvar(itrc)))=.TRUE.
            IRP(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
          END IF
        END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
        IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
          got_var(idUsms)=.TRUE.
          IRP(ng)%Vid(idUsms)=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
          got_var(idVsms)=.TRUE.
          IRP(ng)%Vid(idVsms)=var_id(i)
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              IRP(ng)%Vid(idTsur(itrc))=var_id(i)
            END IF
          END IF
        END DO
#  endif
      END DO

#  ifdef ADJUST_BOUNDARY
      IF (.not.got_var(idSbry(isFsur))) Ldefine=.TRUE.
      IF (.not.got_var(idSbry(isUbar))) Ldefine=.TRUE.
      IF (.not.got_var(idSbry(isVbar))) Ldefine=.TRUE.
#   ifdef SOLVE3D
      IF (.not.got_var(idSbry(isUvel))) Ldefine=.TRUE.
      IF (.not.got_var(idSbry(isVvel))) Ldefine=.TRUE.
      DO itrc=1,NT(ng)
        IF (.not.got_var(idSbry(isTvar(itrc)))) Ldefine=.TRUE.
      END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      IF (.not.got_var(idUsms)) Ldefine=.TRUE.
      IF (.not.got_var(idVsms)) Ldefine=.TRUE.
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          IF (.not.got_var(idTsur(itrc))) Ldefine=.TRUE.
        END IF
      END DO
#  endif
!
!  Put existing file into define mode so new variables can be added.
!
      IF (Ldefine) THEN
        CALL netcdf_redef (ng, IRPM, ncname, IRP(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      DEFINE: IF (Ldefine) THEN

#  ifdef ADJUST_BOUNDARY
        got_IorJ=.FALSE.
        got_boundary=.FALSE.
        got_obc_adjust=.FALSE.
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        got_frc_adjust=.FALSE.
#  endif
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
#  ifdef SOLVE3D
            CASE ('s_rho')
              DimIDs( 9)=dim_id(i)
            CASE ('s_w')
              DimIDs(10)=dim_id(i)
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('boundary')
              DimIDs(14)=dim_id(i)
              got_boundary=.TRUE.
            CASE ('IorJ')
              IorJdim=dim_id(i)
              got_IorJ=.TRUE.
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            CASE ('frc_adjust')
              frecdim=dim_id(i)
              got_frc_adjust=.TRUE.
#  endif
#  ifdef ADJUST_BOUNDARY
            CASE ('obc_adjust')
              brecdim=dim_id(i)
              got_obc_adjust=.TRUE.
#  endif
          END SELECT
        END DO

        DimIDs(12)=rec_id
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_boundary) THEN
          status=def_dim(ng, iRPM, IRP(ng)%ncid, ncname, 'boundary',    &
     &                   4, DimIDs(14))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
        IF (.not.got_IorJ) THEN
          status=def_dim(ng, iRPM, IRP(ng)%ncid, ncname, 'IorJ',        &
     &                   IOBOUNDS(ng)%IorJ, IorJdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        IF (.not.got_frc_adjust) THEN
          status=def_dim(ng, iRPM, IRP(ng)%ncid, ncname, 'frc_adjust',  &
     &                   Nfrec(ng), frecdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_obc_adjust) THEN
          status=def_dim(ng, iRPM, IRP(ng)%ncid, ncname, 'obc_adjust',  &
     &                   Nbrec(ng), brecdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#  ifdef ADJUST_BOUNDARY
        t2dobc(1)=IorJdim
        t2dobc(2)=DimIDs(14)
        t2dobc(3)=brecdim
        t2dobc(4)=DimIDs(12)
#   ifdef SOLVE3D
        t3dobc(1)=IorJdim
        t3dobc(2)=DimIDs( 9)
        t3dobc(3)=DimIDs(14)
        t3dobc(4)=brecdim
        t3dobc(5)=DimIDs(12)
#   endif
#  endif
#  ifdef ADJUST_STFLUX
        t4dfrc(1)=DimIDs( 1)
        t4dfrc(2)=DimIDs( 5)
        t4dfrc(3)=frecdim
        t4dfrc(4)=DimIDs(12)
#  endif
#  ifdef ADJUST_WSTRESS
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u4dfrc(1)=DimIDs( 2)
        u4dfrc(2)=DimIDs( 6)
        u4dfrc(3)=frecdim
        u4dfrc(4)=DimIDs(12)
#  endif
#  ifdef ADJUST_WSTRESS
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v4dfrc(1)=DimIDs( 3)
        v4dfrc(2)=DimIDs( 7)
        v4dfrc(3)=frecdim
        v4dfrc(4)=DimIDs(12)
#  endif
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
!  Define additional variables. Notice that these variables have their
!  own fixed time-dimension to allow 4DVAR adjustments at other times
!  in addition to initialization time.
!-----------------------------------------------------------------------

#  ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SOLVE3D
!
!  Define 3D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(ifield),   &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(ifield),   &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables open boundaries.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(ifield), &
     &                     NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        IF (.not.got_var(idUsms)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idUsms),   &
     &                   NF_FOUT, 4, u4dfrc, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (.not.got_var(idVsms)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idVsms),   &
     &                   NF_FOUT, 4, v4dfrc, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Define surface tracer flux.
!
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTsur(itrc)))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, iRPM, IRP(ng)%ncid,                      &
     &                     IRP(ng)%Vid(idTsur(itrc)), NF_FOUT,          &
     &                     4, t4dfrc, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iRPM, ncname, IRP(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
# endif
!
!=======================================================================
!  Open an existing tangent linear initial file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      IF (.not.LdefIRP(ng)) THEN
        ncname=IRP(ng)%name

# if !(defined ADJUST_STFLUX || defined ADJUST_WSTRESS)
!
!  Open representer model initialization file for read/write.
!
        CALL netcdf_open (ng, iRPM, TRIM(ncname), 1, IRP(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iRPM, ncname,                        &
     &                         ncid = IRP(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iRPM, ncname,                          &
     &                         ncid = IRP(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
# endif
!
!  Scan variable list from input NetCDF and activate switches for
!  initialization variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            IRP(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            IRP(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            IRP(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            IRP(ng)%Vid(idVbar)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            IRP(ng)%Vid(idSbry(isFsur))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            IRP(ng)%Vid(idSbry(isUbar))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            IRP(ng)%Vid(idSbry(isVbar))=var_id(i)
# endif
# ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            IRP(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            IRP(ng)%Vid(idVsms)=var_id(i)
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            IRP(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            IRP(ng)%Vid(idVvel)=var_id(i)
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            IRP(ng)%Vid(idSbry(isUvel))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            IRP(ng)%Vid(idSbry(isVvel))=var_id(i)
#  endif
#  if defined BVF_MIXING || defined LMD_MIXING || \
      defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            IRP(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            IRP(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            IRP(ng)%Vid(idSdif)=var_id(i)
#  endif
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              IRP(ng)%Tid(itrc)=var_id(i)
#  ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              IRP(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
#  endif
#  ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              IRP(ng)%Vid(idTsur(itrc))=var_id(i)
#  endif
            END IF
          END DO
# endif
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur)).and.                           &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar)).and.                           &
     &      ANY(Lobc(:,isUbar,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar)).and.                           &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel)).and.                           &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel)).and.                           &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel)).and.                           &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel)).and.                           &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc))).and.                   &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            IF (Master) WRITE (stdout,30)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
#  ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IRP(ng)%Rindex=rec_size
        Fcount=IRP(ng)%Fcount
        IRP(ng)%Nrec(Fcount)=rec_size
      END IF
!
  10  FORMAT (/,' RP_DEF_INI_NF90 - unable to open initial NetCDF',     &
     &        ' file: ',a)
  20  FORMAT (/,' RP_DEF_INI_NF90 - unable to put in define mode',      &
     &        ' initial NetCDF file: ',a)
  30  FORMAT (/,' RP_DEF_INI_NF90 - unable to find variable: ',a,2x,    &
     &        ' in file: ',a)
!
      RETURN
      END SUBROUTINE rp_def_ini_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE rp_def_ini_pio (ng)
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
      logical :: Ldefine = .FALSE.
#  ifdef ADJUST_BOUNDARY
      logical :: got_IorJ, got_boundary, got_obc_adjust
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      logical :: got_frc_adjust
#  endif
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, status
      integer :: Fcount
#  ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
#  endif
      integer :: DimIDs(nDimID)
#  ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
#   ifdef SOLVE3D
      integer :: t3dobc(5)
#   endif
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: t4dfrc(4), u4dfrc(4), v4dfrc(4)
#  endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", rp_def_ini_pio"
!
      SourceFile=MyFile

#  if defined ADJUST_BOUNDARY || \
     defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
!
!=======================================================================
!  Open existing representer model initial conditions and define new
!  variables.
!=======================================================================
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=IRP(ng)%name
!
!  Open representer model initialization file for read/write.
!
      CALL pio_netcdf_open (ng, iRPM, ncname, 1, IRP(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF
!
!  Inquire about the dimensions and check for consistency.
!
      CALL pio_netcdf_check_dim (ng, iRPM, ncname,                      &
     &                           pioFile = IRP(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, iRPM, ncname,                        &
     &                         pioFile = IRP(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Check if surface forcing variables have been already defined.
!
      DO i=1,NV
        got_var(i)=.FALSE.
      END DO
!
      DO i=1,n_var
#   ifdef ADJUST_BOUNDARY
        IF (TRIM(var_name(i)).eq.                                       &
     &      TRIM(Vname(1,idSbry(isFsur)))) THEN
          got_var(idSbry(isFsur))=.TRUE.
          IRP(ng)%Vid(idSbry(isFsur))=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isUbar)))) THEN
          got_var(idSbry(isUbar))=.TRUE.
          IRP(ng)%Vid(idSbry(isUbar))=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isVbar)))) THEN
          got_var(idSbry(isVbar))=.TRUE.
          IRP(ng)%Vid(idSbry(isVbar))=var_id(i)
#    ifdef SOLVE3D
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isUvel)))) THEN
          got_var(idSbry(isUvel))=.TRUE.
          IRP(ng)%Vid(idSbry(isUvel))=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.                                  &
     &           TRIM(Vname(1,idSbry(isVvel)))) THEN
          got_var(idSbry(isVvel))=.TRUE.
          IRP(ng)%Vid(idSbry(isVvel))=var_id(i)
#    endif
        END IF
#    ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (TRIM(var_name(i)).eq.                                     &
     &        TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
            got_var(idSbry(isTvar(itrc)))=.TRUE.
            IRP(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
          END IF
        END DO
#    endif
#   endif
#   ifdef ADJUST_WSTRESS
        IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
          got_var(idUsms)=.TRUE.
          IRP(ng)%Vid(idUsms)=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
          got_var(idVsms)=.TRUE.
          IRP(ng)%Vid(idVsms)=var_id(i)
        END IF
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              IRP(ng)%Vid(idTsur(itrc))=var_id(i)
            END IF
          END IF
        END DO
#   endif
      END DO

#   ifdef ADJUST_BOUNDARY
      IF (.not.got_var(idSbry(isFsur))) Ldefine=.TRUE.
      IF (.not.got_var(idSbry(isUbar))) Ldefine=.TRUE.
      IF (.not.got_var(idSbry(isVbar))) Ldefine=.TRUE.
#    ifdef SOLVE3D
      IF (.not.got_var(idSbry(isUvel))) Ldefine=.TRUE.
      IF (.not.got_var(idSbry(isVvel))) Ldefine=.TRUE.
      DO itrc=1,NT(ng)
        IF (.not.got_var(idSbry(isTvar(itrc)))) Ldefine=.TRUE.
      END DO
#    endif
#   endif
#   ifdef ADJUST_WSTRESS
      IF (.not.got_var(idUsms)) Ldefine=.TRUE.
      IF (.not.got_var(idVsms)) Ldefine=.TRUE.
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          IF (.not.got_var(idTsur(itrc))) Ldefine=.TRUE.
        END IF
      END DO
#   endif
!
!  Put existing file into define mode so new variables can be added.
!
      IF (Ldefine) THEN
        CALL pio_netcdf_redef (ng, IRPM, ncname, IRP(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      DEFINE: IF (Ldefine) THEN

#   ifdef ADJUST_BOUNDARY
        got_IorJ=.FALSE.
        got_boundary=.FALSE.
        got_obc_adjust=.FALSE.
#   endif
#   if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        got_frc_adjust=.FALSE.
#   endif
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
#   ifdef SOLVE3D
            CASE ('s_rho')
              DimIDs( 9)=dim_id(i)
            CASE ('s_w')
              DimIDs(10)=dim_id(i)
#   endif
#   ifdef ADJUST_BOUNDARY
            CASE ('boundary')
              DimIDs(14)=dim_id(i)
              got_boundary=.TRUE.
            CASE ('IorJ')
              IorJdim=dim_id(i)
              got_IorJ=.TRUE.
#   endif
#   if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            CASE ('frc_adjust')
              frecdim=dim_id(i)
              got_frc_adjust=.TRUE.
#   endif
#   ifdef ADJUST_BOUNDARY
            CASE ('obc_adjust')
              brecdim=dim_id(i)
              got_obc_adjust=.TRUE.
#   endif
          END SELECT
        END DO
!
        DimIDs(12)=rec_id
#   ifdef ADJUST_BOUNDARY
        IF (.not.got_boundary) THEN
          status=def_dim(ng, iRPM, IRP(ng)%pioFile, ncname,             &
     &                   'boundary', 4, DimIDs(14))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
        IF (.not.got_IorJ) THEN
          status=def_dim(ng, iRPM, IRP(ng)%pioFile, ncname,             &
     &                   'IorJ', IOBOUNDS(ng)%IorJ, IorJdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        IF (.not.got_frc_adjust) THEN
          status=def_dim(ng, iRPM, IRP(ng)%pioFile, ncname,             &
     &                   'frc_adjust', Nfrec(ng), frecdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef ADJUST_BOUNDARY
        IF (.not.got_obc_adjust) THEN
          status=def_dim(ng, iRPM, IRP(ng)%pioFile, ncname,             &
     &                   'obc_adjust', Nbrec(ng), brecdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#   ifdef ADJUST_BOUNDARY
        t2dobc(1)=IorJdim
        t2dobc(2)=DimIDs(14)
        t2dobc(3)=brecdim
        t2dobc(4)=DimIDs(12)
#    ifdef SOLVE3D
        t3dobc(1)=IorJdim
        t3dobc(2)=DimIDs( 9)
        t3dobc(3)=DimIDs(14)
        t3dobc(4)=brecdim
        t3dobc(5)=DimIDs(12)
#    endif
#   endif
#   ifdef ADJUST_STFLUX
        t4dfrc(1)=DimIDs( 1)
        t4dfrc(2)=DimIDs( 5)
        t4dfrc(3)=frecdim
        t4dfrc(4)=DimIDs(12)
#   endif
#   ifdef ADJUST_WSTRESS
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u4dfrc(1)=DimIDs( 2)
        u4dfrc(2)=DimIDs( 6)
        u4dfrc(3)=frecdim
        u4dfrc(4)=DimIDs(12)
#   endif
#   ifdef ADJUST_WSTRESS
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v4dfrc(1)=DimIDs( 3)
        v4dfrc(2)=DimIDs( 7)
        v4dfrc(3)=frecdim
        v4dfrc(4)=DimIDs(12)
#   endif
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
!  Define additional variables. Notice that these variables have their
!  own fixed time-dimension to allow 4DVAR adjustments at other times
!  in addition to initialization time.
!-----------------------------------------------------------------------

#   ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          IRP(ng)%pioVar(ifield)%dkind=PIO_FOUT
          IRP(ng)%pioVar(ifield)%gtype=r2dobc
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(ifield)%vd, PIO_FOUT,           &
     &                   4, t2dobc, Aval, Vinfo, ncname,                &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          IRP(ng)%pioVar(ifield)%dkind=PIO_FOUT
          IRP(ng)%pioVar(ifield)%gtype=u2dobc
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(ifield)%vd, PIO_FOUT,           &
     &                   4, t2dobc, Aval, Vinfo, ncname,                &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          IRP(ng)%pioVar(ifield)%dkind=PIO_FOUT
          IRP(ng)%pioVar(ifield)%gtype=v2dobc
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(ifield)%vd, PIO_FOUT,           &
     &                   4, t2dobc, Aval, Vinfo, ncname,                &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#    ifdef SOLVE3D
!
!  Define 3D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          IRP(ng)%pioVar(ifield)%dkind=PIO_FOUT
          IRP(ng)%pioVar(ifield)%gtype=u3dobc
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(ifield)%vd, PIO_FOUT,           &
     &                   5, t3dobc, Aval, Vinfo, ncname,                &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          IRP(ng)%pioVar(ifield)%dkind=PIO_FOUT
          IRP(ng)%pioVar(ifield)%gtype=v3dobc
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(ifield)%vd, PIO_FOUT,           &
     &                   5, t3dobc, Aval, Vinfo, ncname,                &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables open boundaries.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            IRP(ng)%pioVar(ifield)%dkind=PIO_FOUT
            IRP(ng)%pioVar(ifield)%gtype=r3dobc
!
            status=def_var(ng, iRPM, IRP(ng)%pioFile,                   &
     &                     IRP(ng)%pioVar(ifield)%vd, PIO_FOUT,         &
     &                     5, t3dobc, Aval, Vinfo, ncname,              &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#    endif
#   endif
#   ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        IF (.not.got_var(idUsms)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          IRP(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          IRP(ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idUsms)%vd, PIO_FOUT,           &
     &                   4, u4dfrc, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (.not.got_var(idVsms)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          IRP(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          IRP(ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idVsms)%vd, PIO_FOUT,           &
     &                   4, v4dfrc, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Define surface tracer flux.
!
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTsur(itrc)))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            IRP(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
            IRP(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
!
            status=def_var(ng, iRPM, IRP(ng)%pioFile,                   &
     &                     IRP(ng)%pioVar(idTsur(itrc))%vd, PIO_FOUT,   &
     &                     4, t4dfrc, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, iRPM, ncname, IRP(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
#  endif
!
!=======================================================================
!  Open an existing tangent linear initial file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      IF (.not.LdefIRP(ng)) THEN
        ncname=IRP(ng)%name

#  if !(defined ADJUST_STFLUX || defined ADJUST_WSTRESS)
!
!  Open representer model initialization file for read/write.
!
        CALL pio_netcdf_open (ng, iRPM, ncname, 1, IRP(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, iRPM, ncname,                    &
     &                             pioFile = IRP(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, iRPM, ncname,                      &
     &                           pioFile = IRP(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
#  endif
!
!  Scan variable list from input NetCDF and activate switches for
!  initialization variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            IRP(ng)%pioVar(idtime)%vd=var_desc(i)
            IRP(ng)%pioVar(idtime)%dkind=PIO_TOUT
            IRP(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            IRP(ng)%pioVar(idFsur)%vd=var_desc(i)
            IRP(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            IRP(ng)%pioVar(idUbar)%vd=var_desc(i)
            IRP(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            IRP(ng)%pioVar(idVbar)%vd=var_desc(i)
            IRP(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idVbar)%gtype=v2dvar
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            IRP(ng)%pioVar(idSbry(isFsur))%vd=var_desc(i)
            IRP(ng)%pioVar(idSbry(isFsur))%dkind=PIO_FOUT
            IRP(ng)%pioVar(idSbry(isFsur))%gtype=r2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            IRP(ng)%pioVar(idSbry(isUbar))%vd=var_desc(i)
            IRP(ng)%pioVar(idSbry(isUbar))%dkind=PIO_FOUT
            IRP(ng)%pioVar(idSbry(isUbar))%gtype=u2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            IRP(ng)%pioVar(idSbry(isVbar))%vd=var_desc(i)
            IRP(ng)%pioVar(idSbry(isVbar))%dkind=PIO_FOUT
            IRP(ng)%pioVar(idSbry(isVbar))%gtype=v2dobc
#  endif
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            IRP(ng)%pioVar(idUsms)%vd=var_desc(i)
            IRP(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            IRP(ng)%pioVar(idVsms)%vd=var_desc(i)
            IRP(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idVsms)%gtype=v2dvar
#  endif
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            IRP(ng)%pioVar(idUvel)%vd=var_desc(i)
            IRP(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            IRP(ng)%pioVar(idVvel)%vd=var_desc(i)
            IRP(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idVvel)%gtype=v3dvar
#   ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            IRP(ng)%pioVar(idSbry(isUvel))%vd=var_desc(i)
            IRP(ng)%pioVar(idSbry(isUvel))%dkind=PIO_FOUT
            IRP(ng)%pioVar(idSbry(isUvel))%gtype=u3dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            IRP(ng)%pioVar(idSbry(isVvel))%vd=var_desc(i)
            IRP(ng)%pioVar(idSbry(isVvel))%dkind=PIO_FOUT
            IRP(ng)%pioVar(idSbry(isVvel))%gtype=v3dobc
#   endif
#   if defined BVF_MIXING || defined LMD_MIXING || \
      defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            IRP(ng)%pioVar(idVvis)%vd=var_desc(i)
            IRP(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idVvis)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            IRP(ng)%pioVar(idTdif)%vd=var_desc(i)
            IRP(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idTdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            IRP(ng)%pioVar(idSdif)%vd=var_desc(i)
            IRP(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            IRP(ng)%pioVar(idSdif)%gtype=w3dvar
#   endif
#  endif
          END IF
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              IRP(ng)%pioTrc(itrc)%vd=var_desc(i)
              IRP(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              IRP(ng)%pioTrc(itrc)%gtype=r3dvar
#   ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              IRP(ng)%pioVar(idSbry(isTvar(itrc)))%vd=var_desc(i)
              IRP(ng)%pioVar(idSbry(isTvar(itrc)))%dkind=PIO_FOUT
              IRP(ng)%pioVar(idSbry(isTvar(itrc)))%gtype=r3dobc
#   endif
#   ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              IRP(ng)%pioVar(idTsur(itrc))%vd=var_desc(i)
              IRP(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              IRP(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
#   endif
            END IF
          END DO
#  endif
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur)).and.                           &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar)).and.                           &
     &      ANY(Lobc(:,isUbar,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar)).and.                           &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel)).and.                           &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel)).and.                           &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel)).and.                           &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel)).and.                           &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc))).and.                   &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            IF (Master) WRITE (stdout,30)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   endif
#   ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   endif
        END DO
#  endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IRP(ng)%Rindex=rec_size
        Fcount=IRP(ng)%Fcount
        IRP(ng)%Nrec(Fcount)=rec_size
      END IF
!
  10  FORMAT (/,' RP_DEF_INI_PIO - unable to open initial NetCDF',      &
     &        ' file: ',a)
  20  FORMAT (/,' RP_DEF_INI_PIO - unable to put in define mode',       &
     &        ' initial NetCDF file: ',a)
  30  FORMAT (/,' RP_DEF_INI_PIO - unable to find variable: ',a,2x,     &
     &        ' in file: ',a)
!
      RETURN
      END SUBROUTINE rp_def_ini_pio
# endif
#endif
      END MODULE rp_def_ini_mod
