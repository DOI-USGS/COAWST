#include "cppdefs.h"
      MODULE def_norm_mod
#ifdef FOUR_DVAR
!
!git $Id$
!svn $Id: def_norm.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates error covariance normalization factors file    !
!  using either the standard NetCDF library or the Parallel-IO (PIO)   !
!  library.  It defines its dimensions, attributes, and variables.     !
!                                                                      !
!  The NetCDF files are  used for variational data assimilation. Four  !
!  different files can be generates for initial, model, opem boundary, !
!  and surface forcing error covariances.                              !
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
      PUBLIC  :: def_norm
      PRIVATE :: def_norm_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_norm_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_norm (ng, model, ifile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ifile
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
      SELECT CASE (NRM(ifile,ng)%IOtype)
        CASE (io_nf90)
          CALL def_norm_nf90 (ng, model, ifile)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_norm_pio (ng, model, ifile)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) NRM(ifile,ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_NORM - Illegal output file type, io_type = ',i0,    &
     &        /,12x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_norm
!
!***********************************************************************
      SUBROUTINE def_norm_nf90 (ng, model, ifile)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ifile
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, nvd3, nvd4
      integer :: recdim, status, varid
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim, ifield
# endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef ADJUST_BOUNDARY
      integer :: t2dobc(3)
# endif
# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
#  ifdef ADJUST_BOUNDARY
      integer :: t3dobc(4)
#  endif
# endif
!
      real(r8) :: Aval(6)
!
      character (len=60 )    :: Text
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_norm_nf90"
!
      SourceFile=MyFile
!
!=======================================================================
!  Create a new background covariace normalization file.
!=======================================================================
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=NRM(ifile,ng)%name
!
      DEFINE : IF (LdefNRM(ifile,ng)) THEN
        CALL netcdf_create (ng, iTLM, TRIM(ncname), NRM(ifile,ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xi_rho',  &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xi_u',    &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xi_v',    &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xi_psi',  &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'eta_rho', &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'eta_u',   &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'eta_v',   &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'eta_psi', &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
        IF (ifile.eq.3) THEN
          status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'IorJ',  &
     &                   IOBOUNDS(ng)%IorJ, IorJdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xy_rho',  &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xy_u',    &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xy_v',    &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xyz_rho', &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xyz_u',   &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xyz_v',   &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xyz_w',   &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 's_rho',   &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 's_w',     &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'tracer',  &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'NST',     &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'Nbed',    &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'xybed',   &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, NRM(ifile,ng)%ncid, ncname, 'Nbands',  &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'Nphy',    &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'Nbac',    &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'Ndom',    &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'Nfec',    &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'boundary',&
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname, 'Nstate',  &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iTLM, NRM(ifile,ng)%ncid, ncname,            &
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
# ifdef ADJUST_BOUNDARY
        IF (ifile.eq.3) THEN
          t2dobc(1)=IorJdim
          t2dobc(2)=DimIDs(14)
          t2dobc(3)=DimIDs(12)
#  ifdef SOLVE3D
          t3dobc(1)=IorJdim
          t3dobc(2)=DimIDs( 9)
          t3dobc(3)=DimIDs(14)
          t3dobc(4)=DimIDs(12)
#  endif
        END IF
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
        NRM(ifile,ng)%Rindex=0
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
        CALL def_info (ng, iTLM, NRM(ifile,ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define background covariance normalization variables.
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
        status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                    &
     &                 NRM(ifile,ng)%Vid(idtime),                       &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Initial conditions or model error covariance normalization factors.
!-----------------------------------------------------------------------
!
        IF ((ifile.eq.1).or.(ifile.eq.2)) THEN
          IF (ifile.eq.1) THEN
            Text='initial conditions error covariance normalization'
          ELSE IF (ifile.eq.2) THEN
            Text='model error covariance normalization'
          END IF
!
!  Define free-surface normalization factor.
!
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idFsur)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idFsur),                     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 2D U-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idUbar)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idUbar),                     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 2D V-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idVbar)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idVbar),                     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
!
!  Define 3D U-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idUvel)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idUvel),                     &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 3D V-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idVvel)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idVvel),                     &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define tracer type normalization factors.
!
          DO itrc=1,NT(ng)
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),20) TRIM(Vname(2,idTvar(itrc))), TRIM(Text)
            Vinfo( 3)='nondimensional'
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),30) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                &
     &                     NRM(ifile,ng)%Vid(idTvar(itrc)),             &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Boundary conditions error covariance normalization factors.
!-----------------------------------------------------------------------
!
        ELSE IF (ifile.eq.3) THEN
          Text='error covariance normalization'
!
!  Define free-surface open boundaries.
!
          IF (ANY(Lobc(:,isFsur,ng))) THEN
            ifield=idSbry(isFsur)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                &
     &                     NRM(ifile,ng)%Vid(ifield),                   &
     &                     NF_FOUT, 3, t2dobc, Aval, Vinfo, ncname,     &
     &                   SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define 2D U-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isUbar,ng))) THEN
            ifield=idSbry(isUbar)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                &
     &                     NRM(ifile,ng)%Vid(ifield),                   &
     &                     NF_FOUT, 3, t2dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define 2D V-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isVbar,ng))) THEN
            ifield=idSbry(isVbar)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                &
     &                     NRM(ifile,ng)%Vid(ifield),                   &
     &                     NF_FOUT, 3, t2dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#  ifdef SOLVE3D
!
!  Define 3D U-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isUvel,ng))) THEN
            ifield=idSbry(isUvel)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                &
     &                     NRM(ifile,ng)%Vid(ifield),                   &
     &                     NF_FOUT, 4, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define 3D V-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isVvel,ng))) THEN
            ifield=idSbry(isVvel)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                &
     &                     NRM(ifile,ng)%Vid(ifield),                   &
     &                     NF_FOUT, 4, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define tracer type variables open boundaries.
!
          DO itrc=1,NT(ng)
            IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
              ifield=idSbry(isTvar(itrc))
              Vinfo( 1)=Vname(1,ifield)
              WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
              DO i=1,NST
                IF (itrc.eq.idsed(i)) THEN
                  WRITE (Vinfo(19),30) 1000.0_r8*Sd50(i,ng)
                END IF
              END DO
#   endif
              Vinfo(21)=Vname(6,ifield)
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,              &
     &                       NRM(ifile,ng)%Vid(ifield),                 &
     &                       NF_FOUT, 4, t3dobc, Aval, Vinfo, ncname,   &
     &                       SetFillVal = .FALSE.)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Surface forcing error covariance normalization factors.
!-----------------------------------------------------------------------
!
        ELSE IF (ifile.eq.4) THEN
          Text='error covariance normalization'

#  ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress normalization factors.
!
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idUsms)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idUsms),                     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define surface V-momentum stress normalization factors.
!
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idVsms)), TRIM(Text)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,                  &
     &                   NRM(ifile,ng)%Vid(idVsms),                     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Define surface tracer fluxes.
!
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              Vinfo( 1)=Vname(1,idTsur(itrc))
              WRITE (Vinfo( 2),20) TRIM(Vname(2,idTsur(itrc))),         &
     &                             TRIM(Text)
              Vinfo( 3)='nondimensional'
              Vinfo(14)=Vname(4,idTsur(itrc))
              Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
#   endif
              Vinfo(21)=Vname(6,idTsur(itrc))
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r2dvar,r8)
              status=def_var(ng, iTLM, NRM(ifile,ng)%ncid,              &
     &                       NRM(ifile,ng)%Vid(idTsur(itrc)), NF_FOUT,  &
     &                       nvd3, t2dgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#  endif
# endif
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, NRM(ifile,ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.  Deactive file
!  creation switch.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, NRM(ifile,ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        LdefNRM(ifile,ng)=.FALSE.
      END IF DEFINE
!
!=======================================================================
!  Open an existing normalization file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      QUERY : IF (.not.LdefNRM(ifile,ng)) THEN
        ncname=NRM(ifile,ng)%name
!
!  Open normalization file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, NRM(ifile,ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,40) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = NRM(ifile,ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = NRM(ifile,ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  normalization variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            NRM(ifile,ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            NRM(ifile,ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            NRM(ifile,ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            NRM(ifile,ng)%Vid(idVbar)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            NRM(ifile,ng)%Vid(idSbry(isFsur))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            NRM(ifile,ng)%Vid(idSbry(isUbar))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            NRM(ifile,ng)%Vid(idSbry(isVbar))=var_id(i)
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            NRM(ifile,ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            NRM(ifile,ng)%Vid(idVvel)=var_id(i)
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            NRM(ifile,ng)%Vid(idSbry(isUvel))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            NRM(ifile,ng)%Vid(idSbry(isVvel))=var_id(i)
#  endif
# endif
# ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            NRM(ifile,ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            NRM(ifile,ng)%Vid(idVsms)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              NRM(ifile,ng)%Vid(idTvar(itrc))=var_id(i)
#  ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              NRM(ifile,ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
#  endif
#  ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              NRM(ifile,ng)%Vid(idTsur(itrc))=var_id(i)
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
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSvar(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSvar(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms).and.(ifile.eq.4)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.(ifile.eq.4)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.(ifile.le.2)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc))).and.(ifile.eq.3)) THEN
            IF (Master) WRITE (stdout,50)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
#  ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.(ifile.eq.4).and.          &
     &        Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
# endif
!
!  Set unlimited time record dimension to zero to allow other programs
!  to process and write normalization factors for different variable.
!
        NRM(ifile,ng)%Rindex=0
      END IF QUERY
!
  10  FORMAT (/,' DEF_NORM_NF90 - unable to create norm NetCDF',        &
     &        ' file: ',a)
  20  FORMAT (a,', ',a)
  30  FORMAT (1pe11.4,1x,'millimeter')
  40  FORMAT (/,' DEF_NORM_NF90 - unable to open norm NetCDF file: ',a)
  50  FORMAT (/,' DEF_NORM_NF90 - unable to find variable: ',a,2x,      &
     &        ' in norm NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_norm_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_norm_pio (ng, model, ifile)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ifile
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, nvd3, nvd4
      integer :: recdim, status
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim, ifield
# endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef ADJUST_BOUNDARY
      integer :: t2dobc(3)
# endif
# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
#  ifdef ADJUST_BOUNDARY
      integer :: t3dobc(4)
#  endif
# endif
!
      real(r8) :: Aval(6)
!
      character (len=60 )    :: Text
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_norm_pio"
!
      SourceFile=MyFile
!
!=======================================================================
!  Create a new background covariace normalization file.
!=======================================================================
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=NRM(ifile,ng)%name
!
      DEFINE : IF (LdefNRM(ifile,ng)) THEN
        CALL pio_netcdf_create (ng, iTLM, TRIM(ncname),                 &
     &                          NRM(ifile,ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xi_rho', IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xi_u',IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xi_v', IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xi_psi', IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'eta_rho', IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'eta_u', IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'eta_v', IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'eta_psi', IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
        IF (ifile.eq.3) THEN
          status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,       &
     &                   'IorJ', IOBOUNDS(ng)%IorJ, IorJdim)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xy_rho', IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xy_u', IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,
     &                 'xy_v', IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xyz_rho', IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xyz_u', IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xyz_v', IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xyz_w', IOBOUNDS(ng)%xy_rho*(N(ng)+1),          &
     &                 DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 's_rho', N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 's_w', N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'tracer', NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'NST', NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'Nbed', Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'xybed', IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'Nbands', NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'Nphy', Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'Nbac', Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,
     &                 'Ndom', Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'Nfec', Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'boundary', 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 'Nstate', NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iTLM, NRM(ifile,ng)%pioFile, ncname,         &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
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
# ifdef ADJUST_BOUNDARY
        IF (ifile.eq.3) THEN
          t2dobc(1)=IorJdim
          t2dobc(2)=DimIDs(14)
          t2dobc(3)=DimIDs(12)
#  ifdef SOLVE3D
          t3dobc(1)=IorJdim
          t3dobc(2)=DimIDs( 9)
          t3dobc(3)=DimIDs(14)
          t3dobc(4)=DimIDs(12)
#  endif
        END IF
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
        NRM(ifile,ng)%Rindex=0
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
        CALL def_info (ng, iTLM, NRM(ifile,ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define background covariance normalization variables.
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
        NRM(ifile,ng)%pioVar(idtime)%dkind=PIO_TOUT
        NRM(ifile,ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,                 &
     &                 NRM(ifile,ng)%pioVar(idtime)%vd,                 &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Initial conditions or model error covariance normalization factors.
!-----------------------------------------------------------------------
!
        IF ((ifile.eq.1).or.(ifile.eq.2)) THEN
          IF (ifile.eq.1) THEN
            Text='initial conditions error covariance normalization'
          ELSE IF (ifile.eq.2) THEN
            Text='model error covariance normalization'
          END IF
!
!  Define free-surface normalization factor.
!
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idFsur)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          NRM(ifile,ng)%pioVar(idFsur)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idFsur)%gtype=r2dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idFsur)%vd,               &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 2D U-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idUbar)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          NRM(ifile,ng)%pioVar(idUbar)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idUbar)%gtype=u2dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idUbar)%vd,               &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 2D V-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idVbar)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          NRM(ifile,ng)%pioVar(idVbar)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idVbar)%gtype=v2dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idVbar)%vd,               &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
!
!  Define 3D U-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idUvel)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          NRM(ifile,ng)%pioVar(idUvel)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idUvel)%gtype=u3dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idUvel)%vd,               &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define 3D V-momentum normalization factor.
!
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idVvel)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          NRM(ifile,ng)%pioVar(idVvel)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idVvel)%gtype=v3dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idVvel)%vd,               &
     &                 PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define tracer type normalization factors.
!
          DO itrc=1,NT(ng)
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),20) TRIM(Vname(2,idTvar(itrc))), TRIM(Text)
            Vinfo( 3)='nondimensional'
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),30) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            NRM(ifile,ng)%pioTrc(itrc)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioTrc(itrc)%gtype=r3dvar
!
            status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,             &
     &                     NRM(ifile,ng)%pioTrc(itrc)%vd,               &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Boundary conditions error covariance normalization factors.
!-----------------------------------------------------------------------
!
        ELSE IF (ifile.eq.3) THEN
          Text='error covariance normalization'
!
!  Define free-surface open boundaries.
!
          IF (ANY(Lobc(:,isFsur,ng))) THEN
            ifield=idSbry(isFsur)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=r2dobc
!
            status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,             &
     &                     NRM(ifile,ng)%pioVar(ifield)%vd,             &
     &                     PIO_FOUT, 3, t2dobc, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define 2D U-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isUbar,ng))) THEN
            ifield=idSbry(isUbar)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=u2dobc
!
            status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,             &
     &                     NRM(ifile,ng)%pioVar(ifield)%vd,             &
     &                     PIO_FOUT, 3, t2dobc, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define 2D V-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isVbar,ng))) THEN
            ifield=idSbry(isVbar)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=v2dobc
!
            status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,             &
     &                     NRM(ifile,ng)%pioVar(ifield)%vd,             &
     &                     PIO_FOUT, 3, t2dobc, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#  ifdef SOLVE3D
!
!  Define 3D U-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isUvel,ng))) THEN
            ifield=idSbry(isUvel)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=u3dobc
!
            status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,             &
     &                     NRM(ifile,ng)%pioVar(ifield)%vd,             &
     &                     PIO_FOUT, 4, t3dobc, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define 3D V-momentum component open boundaries.
!
          IF (ANY(Lobc(:,isVvel,ng))) THEN
            ifield=idSbry(isVvel)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=v3dobc
!
            status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,             &
     &                     NRM(ifile,ng)%pioVar(ifield)%vd,             &
     &                     PIO_FOUT, 4, t3dobc, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define tracer type variables open boundaries.
!
          DO itrc=1,NT(ng)
            IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
              ifield=idSbry(isTvar(itrc))
              Vinfo( 1)=Vname(1,ifield)
              WRITE (Vinfo( 2),20) TRIM(Vname(2,ifield)), TRIM(Text)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
              DO i=1,NST
                IF (itrc.eq.idsed(i)) THEN
                  WRITE (Vinfo(19),30) 1000.0_r8*Sd50(i,ng)
                END IF
              END DO
#   endif
              Vinfo(21)=Vname(6,ifield)
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
              NRM(ifile,ng)%pioVar(ifield)%gtype=r3dobc
!
              status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,           &
     &                       NRM(ifile,ng)%pioVar(ifield)%vd,           &
     &                       PIO_FOUT, 4, t3dobc, Aval, Vinfo, ncname,  &
     &                       SetFillVal = .FALSE.)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Surface forcing error covariance normalization factors.
!-----------------------------------------------------------------------
!
        ELSE IF (ifile.eq.4) THEN
          Text='error covariance normalization'

#  ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress normalization factors.
!
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idUsms)), TRIM(Text)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          NRM(ifile,ng)%pioVar(idUsms)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idUsms)%vd,               &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define surface V-momentum stress normalization factors.
!
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),20) TRIM(Vname(2,idVsms)), TRIM(Text)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          NRM(ifile,ng)%pioVar(idVsms)%dkind=PIO_FOUT
          NRM(ifile,ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,               &
     &                   NRM(ifile,ng)%pioVar(idVsms)%vd,               &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Define surface tracer fluxes.
!
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              Vinfo( 1)=Vname(1,idTsur(itrc))
              WRITE (Vinfo( 2),20) TRIM(Vname(2,idTsur(itrc))),         &
     &                             TRIM(Text)
              Vinfo( 3)='nondimensional'
              Vinfo(14)=Vname(4,idTsur(itrc))
              Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
#   endif
              Vinfo(21)=Vname(6,idTsur(itrc))
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r2dvar,r8)
              NRM(ifile,ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              NRM(ifile,ng)%pioVar(idTsur(itrc))%gtype=r2dvar
!
              status=def_var(ng, iTLM, NRM(ifile,ng)%pioFile,           &
     &                       NRM(ifile,ng)%pioVar(idTsur(itrc))%vd,     &
     &                       PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo,       &
     &                       ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#  endif
# endif
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname,                      &
     &                          NRM(ifile,ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.  Deactive file
!  creation switch.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, NRM(ifile,ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        LdefNRM(ifile,ng)=.FALSE.
      END IF DEFINE
!
!=======================================================================
!  Open an existing normalization file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      QUERY : IF (.not.LdefNRM(ifile,ng)) THEN
        ncname=NRM(ifile,ng)%name
!
!  Open normalization file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1,                     &
     &                        NRM(ifile,ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,40) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = NRM(ifile,ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = NRM(ifile,ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  normalization variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            NRM(ifile,ng)%pioVar(idtime)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idtime)%dkind=PIO_TOUT
            NRM(ifile,ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            NRM(ifile,ng)%pioVar(idFsur)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idFsur)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            NRM(ifile,ng)%pioVar(idUbar)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idUbar)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            NRM(ifile,ng)%pioVar(idVbar)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idVbar)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idVbar)%gtype=v2dvar
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            ifield=idSbry(isFsur)
            got_var(ifield)=.TRUE.
            NRM(ifile,ng)%pioVar(ifield)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=r2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            ifield=idSbry(isUbar)
            got_var(ifield)=.TRUE.
            NRM(ifile,ng)%pioVar(ifield)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=u2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            ifield=idSbry(isVbar)
            got_var(ifield)=.TRUE.
            NRM(ifile,ng)%pioVar(ifield)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=v2dobc
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            NRM(ifile,ng)%pioVar(idUvel)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idUvel)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            NRM(ifile,ng)%pioVar(idVvel)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idVvel)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idVvel)%gtype=v3dvar
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            ifield=idSbry(isUvel)
            got_var(ifield)=.TRUE.
            NRM(ifile,ng)%pioVar(ifield)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=u3dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            ifield=idSbry(isVvel)
            got_var(ifield)=.TRUE.
            NRM(ifile,ng)%pioVar(ifield)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(ifield)%gtype=v3dobc
#  endif
# endif
# ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            NRM(ifile,ng)%pioVar(idUsms)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idUsms)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            NRM(ifile,ng)%pioVar(idVsms)%vd=var_desc(i)
            NRM(ifile,ng)%pioVar(idVsms)%dkind=PIO_FOUT
            NRM(ifile,ng)%pioVar(idVsms)%gtype=v2dvar
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              NRM(ifile,ng)%pioTrc(itrc)%vd=var_desc(i)
              NRM(ifile,ng)%pioTrc(itrc)%dkind=PIO_FOUT
              NRM(ifile,ng)%pioTrc(itrc)%gtype=r3dvar
#  ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              ifield=idSbry(isTvar(itrc))
              got_var(ifield)=.TRUE.
              NRM(ifile,ng)%pioVar(ifield)%vd=var_desc(i)
              NRM(ifile,ng)%pioVar(ifield)%dkind=PIO_FOUT
              NRM(ifile,ng)%pioVar(ifield)%gtype=r3dobc
#  endif
#  ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              NRM(ifile,ng)%pioVar(idTsur(itrc))%vd=var_desc(i)
              NRM(ifile,ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              NRM(ifile,ng)%pioVar(idTsur(itrc))%gtype=r2dvar
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
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.(ifile.le.2)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSvar(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel)).and.(ifile.eq.3)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idSvar(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms).and.(ifile.eq.4)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.(ifile.eq.4)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.(ifile.le.2)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc))).and.(ifile.eq.3)) THEN
            IF (Master) WRITE (stdout,50)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
#  ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.(ifile.eq.4).and.          &
     &        Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
# endif
!
!  Set unlimited time record dimension to zero to allow other programs
!  to process and write normalization factors for different variable.
!
        NRM(ifile,ng)%Rindex=0
      END IF QUERY
!
  10  FORMAT (/,' DEF_NORM_PIO - unable to create norm NetCDF',         &
     &        ' file: ',a)
  20  FORMAT (a,', ',a)
  30  FORMAT (1pe11.4,1x,'millimeter')
  40  FORMAT (/,' DEF_NORM_PIO - unable to open norm NetCDF file: ',a)
  50  FORMAT (/,' DEF_NORM_PIO - unable to find variable: ',a,2x,       &
     &        ' in norm NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_norm_pio
# endif
#endif
      END MODULE def_norm_mod
