#include "cppdefs.h"
      MODULE def_diags_mod
#ifdef DIAGNOSTICS
!
!git $Id$
!svn $Id: def_diags.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module creates output TIME-AVERAGED DIAGNOSTIC file using      !
!  either the standard NetCDF library or the Parallel-IO (PIO)         !
!  library.  It defines its dimensions, attributes, and variables.     !
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
#ifdef INWAVE_MODEL
      USE mod_inwave_params
#endif
!
      USE def_dim_mod,  ONLY : def_dim
      USE def_info_mod, ONLY : def_info
      USE def_var_mod,  ONLY : def_var
      USE strings_mod,  ONLY : FoundError
      USE wrt_info_mod, ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_diags
      PRIVATE :: def_diags_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_diags_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_diags (ng, ldef)
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
      SELECT CASE (DIA(ng)%IOtype)
        CASE (io_nf90)
          CALL def_diags_nf90 (ng, ldef)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_diags_pio (ng, ldef)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) DIA(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_DIAGS - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_diags
!
!***********************************************************************
      SUBROUTINE def_diags_nf90 (ng, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, itrc, ivar, j, nvd3, nvd4, nvd5
      integer :: recdim, status
# if defined WRITE_WATER && defined MASKING
      integer :: xy_pdim, xyz_pdim
# endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
# if defined ECOSIM && defined DIAGNOSTICS_BIO
      integer :: l3dgrd(4), l4dgrd(5)
# endif

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif
# ifdef INWAVE_MODEL
      integer :: r3degrd(4)
      integer :: u3degrd(4)
      integer :: v3degrd(4)
      integer :: t3degrd(4)
# endif
!
      real(r8) :: Aval(6)
!
      character (len= 13)    :: Prefix
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_diags_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=DIA(ng)%name
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
!  Create a new diagnostics NetCDF file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), DIA(ng)%ncid)
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
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xy_psi',        &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

#ifdef INWAVE_MODEL
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname,'ND',             &
                       ND, DimIDs(34))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xyz_psi',       &
     &                 IOBOUNDS(ng)%xy_psi*N(ng), xyz_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
#  ifdef ECOSIM
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Nbands',        &
     &                 NDbands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iNLM, DIA(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
        nvd5=2
#else
        nvd3=3
        nvd4=4
        nvd5=5
#endif
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
# ifdef INWAVE_MODEL
!
!  Define dimension vector for energy wave layer type variables.
!
        r3degrd(1)=DimIDs( 1)
        r3degrd(2)=DimIDs( 5)
        r3degrd(3)=DimIDs(34)
        r3degrd(4)=DimIDs(12)

        u3degrd(1)=DimIDs( 2)
        u3degrd(2)=DimIDs( 6)
        u3degrd(3)=DimIDs(34)
        u3degrd(4)=DimIDs(12)

        v3degrd(1)=DimIDs( 3)
        v3degrd(2)=DimIDs( 7)
        v3degrd(3)=DimIDs(34)
        v3degrd(4)=DimIDs(12)

        t3degrd(1)=DimIDs( 1)
        t3degrd(2)=DimIDs( 5)
        t3degrd(3)=DimIDs(36)
        t3degrd(4)=DimIDs(12)
# endif
# ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#  else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#  endif
#  ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#   if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#   else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#   endif
#  endif
#  if defined ECOSIM && defined DIAGNOSTICS_BIO
!
!  Define dimension vector for spectral light type variables.
!
        l3dgrd(1)=DimIDs( 1)
        l3dgrd(2)=DimIDs( 5)
        l3dgrd(3)=DimIDs(33)
        l3dgrd(4)=DimIDs(12)
!
        l4dgrd(1)=DimIDs( 1)
        l4dgrd(2)=DimIDs( 5)
        l4dgrd(3)=DimIDs( 9)
        l4dgrd(4)=DimIDs(33)
        l4dgrd(5)=DimIDs(12)
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        DIA(ng)%Rindex=0
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
!  Set long name prefix string.
!
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, DIA(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(idtime),     &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define time-averaged free-surface. It is needed for CF compliance for
!  vertical coordinates formula.
!
        Vinfo( 1)=Vname(1,idFsur)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
# if !defined WET_DRY && (defined WRITE_WATER && defined MASKING)
        Vinfo(20)='mask_rho'
# endif
        Vinfo(21)=Vname(6,idFsur)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(idFsur),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef DIAGNOSTICS_UV
!
!  Define 2D momentum diagnostic fields.
!
        DO ivar=1,NDM2d
          ifield=idDu2d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#  endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u2dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
          ifield=idDv2d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#  endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v2dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#  ifdef SOLVE3D
!
!  Define 3D momentum diagnostic fields.
!
        DO ivar=1,NDM3d
          ifield=idDu3d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
          ifield=idDv3d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
# endif
# ifdef DIAGNOSTICS_TS
!
!  Define tracer diagnostic fields.
!
        DO itrc=1,NT(ng)
          DO ivar=1,NDT
            ifield=idDtrc(itrc,ivar)
            IF (Dout(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
#  endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r3dvar,r8)
              status=def_var(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield), NF_FOUT,              &
     &                       nvd4, t3dgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
        END DO
# endif

# ifdef DIAGNOSTICS_BIO
#  if defined BIO_FENNEL || defined HYPOXIA_SRM
!
!  Define 2D biological diagnostic fields.
!
        DO ivar=1,NDbio2d
          ifield=iDbio2(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
#  if defined BIO_FENNEL
!
!  Define 3D biological diagnostic fields.
!
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  elif defined ECOSIM
!
!  Define 3D biological diagnostic fields.
!
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd4, l3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define 4D biological diagnostic fields.
!
        DO ivar=1,NDbio4d
          ifield=iDbio4(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd5, l4dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
# endif

# ifdef DIAGNOSTICS_INWAVE
!
!  Define 2D InWave diagnostic fields.
!
        DO ivar=1,NDInw
          ifield=idIinw(ivar)
!          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(ifield), &
     &                     NF_FOUT, nvd4, r3degrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!          END IF
        END DO
# endif

!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, DIA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, DIA(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing diagnostics file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=DIA(ng)%name
!
!  Open diagnostics file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, DIA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname,                        &
     &                         ncid = DIA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname,                          &
     &                       ncid = DIA(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  diagnostics variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            DIA(ng)%Vid(idtime)=var_id(i)
          END IF
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            DIA(ng)%Vid(idFsur)=var_id(i)
          END IF
# ifdef DIAGNOSTICS_UV
          DO ivar=1,NDM2d
            IF (TRIM(var_name(i)).eq.                                   &
     &          TRIM(Vname(1,idDu2d(ivar)))) THEN
              got_var(idDu2d(ivar))=.TRUE.
              DIA(ng)%Vid(idDu2d(ivar))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idDv2d(ivar)))) THEN
              got_var(idDv2d(ivar))=.TRUE.
              DIA(ng)%Vid(idDv2d(ivar))=var_id(i)
            END IF
          END DO
#  ifdef SOLVE3D
          DO ivar=1,NDM3d
            IF (TRIM(var_name(i)).eq.                                   &
     &          TRIM(Vname(1,idDu3d(ivar)))) THEN
              got_var(idDu3d(ivar))=.TRUE.
              DIA(ng)%Vid(idDu3d(ivar))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idDv3d(ivar)))) THEN
              got_var(idDv3d(ivar))=.TRUE.
              DIA(ng)%Vid(idDv3d(ivar))=var_id(i)
            END IF
          END DO
#  endif
# endif
# ifdef DIAGNOSTICS_TS
          DO itrc=1,NT(ng)
            DO ivar=1,NDT
              ifield=idDtrc(itrc,ivar)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                DIA(ng)%Vid(ifield)=var_id(i)
              END IF
            END DO
          END DO
# endif
# ifdef DIAGNOSTICS_BIO
#  if defined BIO_FENNEL || defined HYPOXIA_SRM
          DO ivar=1,NDbio2d
            ifield=iDbio2(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
#  endif
#  if defined BIO_FENNEL
          DO ivar=1,NDbio3d
            ifield=iDbio3(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
#  elif defined ECOSIM
          DO ivar=1,NDbio3d
            ifield=iDbio3(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
          DO ivar=1,NDbio4d
            ifield=iDbio4(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
#  endif
# endif
# ifdef DIAGNOSTICS_INWAVE
          DO ivar=1,NDInw
            ifield=idIinw(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if diagnostics variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef DIAGNOSTICS_UV
        DO ivar=1,NDM2d
          ifield=idDu2d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng))  THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          ifield=idDv2d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng))  THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  ifdef SOLVE3D
        DO ivar=1,NDM3d
          ifield=idDu3d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          ifield=idDv3d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
# ifdef DIAGNOSTICS_TS
        DO itrc=1,NT(ng)
          DO ivar=1,NDT
            ifield=idDtrc(itrc,ivar)
            IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
              IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_BIO
#  if defined BIO_FENNEL || defined HYPOXIA_SRM
        DO ivar=1,NDbio2d
          ifield=iDbio2(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined BIO_FENNEL
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  elif defined ECOSIM
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO ivar=1,NDbio4d
          ifield=iDbio4(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
# ifdef DIAGNOSTICS_INWAVE
        DO ivar=1,NDInw
          ifield=idIinw(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nDIA(ng)) THEN
          IF (ndefDIA(ng).gt.0) THEN
            DIA(ng)%Rindex=((ntstart(ng)-1)-                            &
     &                      ndefDIA(ng)*((ntstart(ng)-1)/ndefDIA(ng)))/ &
     &                     nDIA(ng)
          ELSE
            DIA(ng)%Rindex=(ntstart(ng)-1)/nDIA(ng)
          END IF
        ELSE
          DIA(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial averaged time.
!
      IF (ntsDIA(ng).eq.1) THEN
        DIAtime(ng)=time(ng)-0.5_r8*REAL(nDIA(ng),r8)*dt(ng)
      ELSE
        DIAtime(ng)=time(ng)+REAL(ntsDIA(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nDIA(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (2x,'DEF_DIAGS_NF90   - creating diagnostics file,',t56,   &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_DIAGS_NF90   - inquiring diagnostics file,',t56,  &
     &        'Grid ',i2.2,': ', a)
  30  FORMAT (/,' DEF_DIAGS_NF90 - unable to create diagnostics NetCDF',&
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_DIAGS_NF90 - unable to open diagnostics NetCDF',  &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_DIAGS_NF90 - unable to find variable: ',a,2x,     &
     &        ' in diagnostics NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_diags_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_diags_pio (ng, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, itrc, ivar, j, nvd3, nvd4, nvd5
      integer :: recdim, status
#  if defined WRITE_WATER && defined MASKING
      integer :: xy_pdim, xyz_pdim
#  endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
#  if defined ECOSIM && defined DIAGNOSTICS_BIO
      integer :: l3dgrd(4), l4dgrd(5)
#  endif

#  ifdef SOLVE3D
#   ifdef SEDIMENT
      integer :: b3dgrd(4)
#   endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#  endif
!
      real(r8) :: Aval(6)
!
      character (len= 13)    :: Prefix
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_diags_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=DIA(ng)%name
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
!  Create a new diagnostics NetCDF file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, iNLM, TRIM(ncname), DIA(ng)%pioFile)
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
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xi_rho',     &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xi_u',       &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xi_v',       &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xi_psi',     &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'eta_rho',    &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'eta_u',      &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'eta_v',      &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'eta_psi',    &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xy_psi',     &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xy_rho',     &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xy_u',       &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xy_v',       &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xyz_psi',    &
     &                 IOBOUNDS(ng)%xy_psi*N(ng), xyz_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xyz_rho',    &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xyz_u',      &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xyz_v',      &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xyz_w',      &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 's_rho',      &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 's_w',        &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'tracer',     &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SEDIMENT
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'NST',        &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Nbed',       &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'xybed',      &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif
#   ifdef ECOSIM
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Nbands',     &
     &                 NDbands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Nphy',       &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Nbac',       &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Ndom',       &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Nfec',       &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'boundary',   &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname, 'Nstate',     &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iNLM, DIA(ng)%pioFile, ncname,               &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#  if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
        nvd5=2
#  else
        nvd3=3
        nvd4=4
        nvd5=5
#  endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#  if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
#   endif
#  else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
#   endif
#  endif
!
!  Define dimension vectors for staggered u-momentum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
#   endif
#  else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
#   endif
#  endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
#   endif
#  else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
#   endif
#  endif
#  ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
#   if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#   else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#   endif
#   ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#    if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#    else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#    endif
#   endif
#   if defined ECOSIM && defined DIAGNOSTICS_BIO
!
!  Define dimension vector for spectral light type variables.
!
        l3dgrd(1)=DimIDs( 1)
        l3dgrd(2)=DimIDs( 5)
        l3dgrd(3)=DimIDs(33)
        l3dgrd(4)=DimIDs(12)
!
        l4dgrd(1)=DimIDs( 1)
        l4dgrd(2)=DimIDs( 5)
        l4dgrd(3)=DimIDs( 9)
        l4dgrd(4)=DimIDs(33)
        l4dgrd(5)=DimIDs(12)
#   endif
#  endif
!
!  Initialize unlimited time record dimension.
!
        DIA(ng)%Rindex=0
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
!  Set long name prefix string.
!
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, DIA(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        DIA(ng)%pioVar(idtime)%dkind=PIO_TOUT
        DIA(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, iNLM, DIA(ng)%pioFile,                       &
     &                 DIA(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define time-averaged free-surface. It is needed for CF compliance for
!  vertical coordinates formula.
!
        Vinfo( 1)=Vname(1,idFsur)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
#  if !defined WET_DRY && (defined WRITE_WATER && defined MASKING)
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(21)=Vname(6,idFsur)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        AVG(ng)%pioVar(idFsur)%dkind=PIO_FOUT
        AVG(ng)%pioVar(idFsur)%gtype=r2dvar
!
        status=def_var(ng, iNLM, DIA(ng)%pioFile,                       &
     &                 DIA(ng)%pioVar(idFsur)%vd,                       &
     &                 PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef DIAGNOSTICS_UV
!
!  Define 2D momentum diagnostic fields.
!
        DO ivar=1,NDM2d
          ifield=idDu2d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u2dvar,r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=u2dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          ifield=idDv2d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v2dvar,r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=v2dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#   ifdef SOLVE3D
!
!  Define 3D momentum diagnostic fields.
!
        DO ivar=1,NDM3d
          ifield=idDu3d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=u3dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          ifield=idDv3d(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=v3dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
#  endif
#  ifdef DIAGNOSTICS_TS
!
!  Define tracer diagnostic fields.
!
        DO itrc=1,NT(ng)
          DO ivar=1,NDT
            ifield=idDtrc(itrc,ivar)
            IF (Dout(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
#   endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r3dvar,r8)
              DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
              DIA(ng)%pioVar(ifield)%gtype=r3dvar
!
              status=def_var(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield)%vd,                 &
     &                       PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo,       &
     &                       ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
        END DO
#  endif

#  ifdef DIAGNOSTICS_BIO
#   if defined BIO_FENNEL || defined HYPOXIA_SRM
!
!  Define 2D biological diagnostic fields.
!
        DO ivar=1,NDbio2d
          ifield=iDbio2(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=r2dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
#   if defined BIO_FENNEL
!
!  Define 3D biological diagnostic fields.
!
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=r3dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#   elif defined ECOSIM
!
!  Define 3D biological diagnostic fields.
!
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=l3dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd4, l3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define 4D biological diagnostic fields.
!
        DO ivar=1,NDbio4d
          ifield=iDbio4(ivar)
          IF (Dout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
            DIA(ng)%pioVar(ifield)%gtype=l4dvar
!
            status=def_var(ng, iNLM, DIA(ng)%pioFile,                   &
     &                     DIA(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd5, l4dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, iNLM, ncname, DIA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, DIA(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing diagnostics file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=DIA(ng)%name
!
!  Open diagnostics file for read/write.
!
        CALL pio_netcdf_open (ng, iNLM, ncname, 1, DIA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, iNLM, ncname,                    &
     &                             pioFile = DIA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, iNLM, ncname,                      &
     &                           pioFile = DIA(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  diagnostics variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            DIA(ng)%pioVar(idtime)%vd=var_desc(i)
            DIA(ng)%pioVar(idtime)%dkind=PIO_TOUT
            DIA(ng)%pioVar(idtime)%gtype=0
          END IF
#  ifdef DIAGNOSTICS_UV
          DO ivar=1,NDM2d
            IF (TRIM(var_name(i)).eq.                                   &
     &          TRIM(Vname(1,idDu2d(ivar)))) THEN
              got_var(idDu2d(ivar))=.TRUE.
              DIA(ng)%pioVar(idDu2d(ivar))%vd=var_desc(i)
              DIA(ng)%pioVar(idDu2d(ivar))%dkind=PIO_FOUT
              DIA(ng)%pioVar(idDu2d(ivar))%gtype=u2dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idDv2d(ivar)))) THEN
              got_var(idDv2d(ivar))=.TRUE.
              DIA(ng)%pioVar(idDv2d(ivar))%vd=var_desc(i)
              DIA(ng)%pioVar(idDv2d(ivar))%dkind=PIO_FOUT
              DIA(ng)%pioVar(idDv2d(ivar))%gtype=v2dvar
            END IF
          END DO
#   ifdef SOLVE3D
          DO ivar=1,NDM3d
            IF (TRIM(var_name(i)).eq.                                   &
     &          TRIM(Vname(1,idDu3d(ivar)))) THEN
              got_var(idDu3d(ivar))=.TRUE.
              DIA(ng)%pioVar(idDu3d(ivar))%vd=var_desc(i)
              DIA(ng)%pioVar(idDu3d(ivar))%dkind=PIO_FOUT
              DIA(ng)%pioVar(idDu3d(ivar))%gtype=u3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idDv3d(ivar)))) THEN
              got_var(idDv3d(ivar))=.TRUE.
              DIA(ng)%pioVar(idDv3d(ivar))%vd=var_desc(i)
              DIA(ng)%pioVar(idDv3d(ivar))%dkind=PIO_FOUT
              DIA(ng)%pioVar(idDv3d(ivar))%gtype=v3dvar
            END IF
          END DO
#   endif
#  endif
#  ifdef DIAGNOSTICS_TS
          DO itrc=1,NT(ng)
            DO ivar=1,NDT
              ifield=idDtrc(itrc,ivar)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                DIA(ng)%pioVar(ifield)%vd=var_desc(i)
                DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
                DIA(ng)%pioVar(ifield)%gtype=r3dvar
              END IF
            END DO
          END DO
#  endif
#  ifdef DIAGNOSTICS_BIO
#   if defined BIO_FENNEL || defined HYPOXIA_SRM
          DO ivar=1,NDbio2d
            ifield=iDbio2(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%pioVar(ifield)%vd=var_desc(i)
              DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
              DIA(ng)%pioVar(ifield)%gtype=r2dvar
            END IF
          END DO
#   endif
#   if defined BIO_FENNEL
          DO ivar=1,NDbio3d
            ifield=iDbio3(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%pioVar(ifield)%vd=var_desc(i)
              DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
              DIA(ng)%pioVar(ifield)%gtype=r3dvar
            END IF
          END DO
#   elif defined ECOSIM
          DO ivar=1,NDbio3d
            ifield=iDbio3(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%pioVar(ifield)%vd=var_desc(i)
              DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
              DIA(ng)%pioVar(ifield)%gtype=l3dvar
            END IF
          END DO
          DO ivar=1,NDbio4d
            ifield=iDbio4(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              DIA(ng)%pioVar(ifield)%vd=var_desc(i)
              DIA(ng)%pioVar(ifield)%dkind=PIO_FOUT
              DIA(ng)%pioVar(ifield)%gtype=l4dvar
            END IF
          END DO
#   endif
#  endif
        END DO
!
!  Check if diagnostics variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef DIAGNOSTICS_UV
        DO ivar=1,NDM2d
          ifield=idDu2d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng))  THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          ifield=idDv2d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng))  THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   ifdef SOLVE3D
        DO ivar=1,NDM3d
          ifield=idDu3d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          ifield=idDv3d(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif
#  endif
#  ifdef DIAGNOSTICS_TS
        DO itrc=1,NT(ng)
          DO ivar=1,NDT
            ifield=idDtrc(itrc,ivar)
            IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
              IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END DO
        END DO
#  endif
#  ifdef DIAGNOSTICS_BIO
#   if defined BIO_FENNEL || defined HYPOXIA_SRM
        DO ivar=1,NDbio2d
          ifield=iDbio2(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif
#   if defined BIO_FENNEL
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   elif defined ECOSIM
        DO ivar=1,NDbio3d
          ifield=iDbio3(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO ivar=1,NDbio4d
          ifield=iDbio4(ivar)
          IF (.not.got_var(ifield).and.Dout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif
#  endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nDIA(ng)) THEN
          IF (ndefDIA(ng).gt.0) THEN
            DIA(ng)%Rindex=((ntstart(ng)-1)-                            &
     &                      ndefDIA(ng)*((ntstart(ng)-1)/ndefDIA(ng)))/ &
     &                     nDIA(ng)
          ELSE
            DIA(ng)%Rindex=(ntstart(ng)-1)/nDIA(ng)
          END IF
        ELSE
          DIA(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial averaged time.
!
      IF (ntsDIA(ng).eq.1) THEN
        DIAtime(ng)=time(ng)-0.5_r8*REAL(nDIA(ng),r8)*dt(ng)
      ELSE
        DIAtime(ng)=time(ng)+REAL(ntsDIA(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nDIA(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (2x,'DEF_DIAGS_PIO    - creating diagnostics file,',t56,   &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_DIAGS_PIO    - inquiring diagnostics file,',t56,  &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_DIAGS_PIO - unable to create diagnostics NetCDF', &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_DIAGS_PIO - unable to open diagnostics NetCDF',   &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_DIAGS_PIO - unable to find variable: ',a,2x,      &
     &        ' in diagnostics NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_diags_pio
# endif
#endif
      END MODULE def_diags_mod
