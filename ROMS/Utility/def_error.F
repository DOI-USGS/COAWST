#include "cppdefs.h"
      MODULE def_error_mod
#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
!
!git $Id$
!svn $Id: def_error.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module creates full posterior error covariance (diagonal)      !
!  matrix for weak constraint 4DVar data assimilation using either     !
!  the standard NetCDF library or the Parallel-IO (PIO) library. It    !
!  defines its dimensions, attributes, and variables.                  !
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
# ifdef SEDIMENT
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
      PUBLIC  :: def_error
      PRIVATE :: def_error_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_error_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_error (ng)
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
      SELECT CASE (ERR(ng)%IOtype)
        CASE (io_nf90)
          CALL def_error_nf90 (ng)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_error_pio (ng)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) ERR(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_ERROR - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_error
!
!***********************************************************************
      SUBROUTINE def_error_nf90 (ng)
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
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, nrec, nvd, nvd3, nvd4
      integer :: recdim, status, varid
      integer :: NinnerDim
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
# endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
# endif

# ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#  ifdef ADJUST_BOUNDARY
      integer :: t3dobc(5)
#  endif
#  ifdef ADJUST_STFLUX
      integer :: t3dfrc(4)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      integer :: u3dfrc(4), v3dfrc(4)
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_error_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=ERR(ng)%name
!
      IF (Master) THEN
        IF (LdefERR(ng)) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new posterior error covariance matrix file.
!=======================================================================
!
      DEFINE : IF (LdefERR(ng)) THEN
        CALL netcdf_create (ng, iTLM, TRIM(ncname), ERR(ng)%ncid)
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
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'IorJ',          &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'N',             &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, ERR(ng)%ncid, ncname, 'Nbands',        &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'Ninner',        &
     &                 Ninner, NinnerDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'frc_adjust',    &
     &                 Nfrec(ng), DimIDs(30))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname, 'obc_adjust',    &
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iTLM, ERR(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        frecdim=DimIDs(30)
# endif
# ifdef ADJUST_BOUNDARY
        brecdim=DimIDs(31)
# endif
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
#  ifdef ADJUST_STFLUX
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
        t3dfrc(4)=DimIDs(12)
#  endif
# endif
# ifdef ADJUST_BOUNDARY
        t2dobc(1)=IorJdim
        t2dobc(2)=DimIDs(14)
        t2dobc(3)=brecdim
        t2dobc(4)=DimIDs(12)
#  ifdef SOLVE3D
        t3dobc(1)=IorJdim
        t3dobc(2)=DimIDs( 9)
        t3dobc(3)=DimIDs(14)
        t3dobc(4)=brecdim
        t3dobc(5)=DimIDs(12)
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
#  ifdef ADJUST_WSTRESS
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
        u3dfrc(4)=DimIDs(12)
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
#  ifdef ADJUST_WSTRESS
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
        v3dfrc(4)=DimIDs(12)
#  endif
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
# endif
!
!  Initialize unlimited time record dimension.
!
        ERR(ng)%Rindex=0
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
        CALL def_info (ng, iTLM, ERR(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define inner-loop Lanczos vectors tridiagonal matrix.
!
        Vinfo( 1)='zLanczos_coef'
        Vinfo( 2)='inner-loop Lanczos vector tridiagonal matrix'
        status=def_var(ng, iTLM, ERR(ng)%ncid, varid, NF_TYPE,          &
     &                 2, (/NinnerDim,NinnerDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define inverse of inner-loop Lanczos vectors tridiagonal matrix.
!
        Vinfo( 1)='zLanczos_inv'
        Vinfo( 2)='inverse inner-loop Lanczos vector tridiagonal matrix'
        status=def_var(ng, iTLM, ERR(ng)%ncid, varid, NF_TYPE,          &
     &                 2, (/NinnerDim,NinnerDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define inner-loop Lanczos vectors tridiagonal matrix inversion error.
!  We need to get an identity matrix withing roundoff.
!
        Vinfo( 1)='zLanczos_err'
        Vinfo( 2)=                                                      &
          'inner-loop Lanczos vector tridiagonal matrix inversion error'
        status=def_var(ng, iTLM, ERR(ng)%ncid, varid, NF_TYPE,          &
     &                 2, (/NinnerDim,NinnerDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idtime),     &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface error variance.
!
        Vinfo( 1)=Vname(1,idFsur)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idFsur))
        Vinfo( 3)='meter2'
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(21)=Vname(6,idFsur)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idFsur),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries error variance.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D U-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idUbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbar))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(21)=Vname(6,idUbar)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idUbar),     &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define 2D U-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D V-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idVbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbar))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(21)=Vname(6,idVbar)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idVbar),     &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define 2D V-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# ifdef SOLVE3D
!
!  Define 3D U-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idUvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUvel))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(21)=Vname(6,idUvel)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idUvel),     &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define 3D U-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(ifield),   &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 3D V-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idVvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvel))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(21)=Vname(6,idVvel)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idVvel),     &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define 3D V-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(ifield),   &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define tracer type variables error variance.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTvar(itrc)))
          IF (itrc.eq.itemp) THEN
            Vinfo( 3)='Celsius2'
          ELSE IF (itrc.eq.isalt) THEN
            Vinfo( 3)='nondimensional'
          ELSE
            WRITE (Vinfo( 3),50) TRIM(Vname(3,idTvar(itrc)))
          END IF
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),60) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTvar(itrc))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Tid(itrc),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

#  ifdef ADJUST_BOUNDARY
!
!  Define tracer type variables open boundaries error variance.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
            IF (itrc.eq.itemp) THEN
              Vinfo( 3)='Celsius2'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo( 3)='nondimensional'
            ELSE
              WRITE (Vinfo( 3),50) TRIM(Vname(3,idTvar(itrc)))
            END IF
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),60) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#   endif
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(ifield), &
     &                     NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Define surface tracer fluxes error variance.
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
            WRITE (Vinfo( 3),50) TRIM(Vname(3,idTsur(itrc)))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, iTLM, ERR(ng)%ncid,                      &
     &                     ERR(ng)%Vid(idTsur(itrc)),                   &
     &                     NF_FOUT, nvd4, t3dfrc, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress error variance.
!
        Vinfo( 1)=Vname(1,idUsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
        Vinfo( 3)='meter4 second-4'
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(21)=Vname(6,idUsms)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u2dvar,r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idUsms),     &
     &                 NF_FOUT, nvd4, u3dfrc, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define surface V-momentum stress error variance.
!
        Vinfo( 1)=Vname(1,idVsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)='meter4 second4'
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(21)=Vname(6,idVsms)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v2dvar,r8)
        status=def_var(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idVsms),     &
     &                 NF_FOUT, nvd4, v3dfrc, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iTLM, ncname, ERR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iTLM, ERR(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing posterior error covariance matrix, check its
!  contents, and prepare for appending data.
!=======================================================================
!
      QUERY: IF (.not.LdefERR(ng)) THEN
        ncname=ERR(ng)%name
!
!  Open posterior error covariance file for read/write.
!
        CALL netcdf_open (ng, iTLM, ncname, 1, ERR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,70) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iTLM, ncname,                        &
     &                         ncid = ERR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iTLM, ncname,                          &
     &                       ncid = ERR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  posterior error covariance matrix variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            ERR(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            ERR(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            ERR(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            ERR(ng)%Vid(idVbar)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            ERR(ng)%Vid(idSbry(isFsur))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            ERR(ng)%Vid(idSbry(isUbar))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            ERR(ng)%Vid(idSbry(isVbar))=var_id(i)
# endif
# ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            ERR(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            ERR(ng)%Vid(idVsms)=var_id(i)
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            ERR(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            ERR(ng)%Vid(idVvel)=var_id(i)
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            ERR(ng)%Vid(idSbry(isUvel))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            ERR(ng)%Vid(idSbry(isVvel))=var_id(i)
#  endif
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              ERR(ng)%Tid(itrc)=var_id(i)
#  ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              ERR(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
#  endif
#  ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              ERR(ng)%Vid(idTsur(itrc))=var_id(i)
#  endif
            END IF
          END DO
# endif
        END DO
!
!  Check if posterior error covariance matrix variables are available
!  in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,80) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,80)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
#  ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,80) TRIM(Vname(1,idTsur(itrc))),  &
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
        ERR(ng)%Rindex=rec_size
      END IF QUERY
!
  10  FORMAT (2x,'DEF_ERROR_NF90   - creating error file,',t56,         &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_ERROR_NF90   - inquiring error file,',t56,        &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_ERROR_NF90 - unable to create 4DVar error NetCDF' &
     &          ' file:',1x,a)
# if defined POSTERIOR_ERROR_I
  40  FORMAT (a,', initial posterior error variance')
# elif defined POSTERIOR_ERROR_F
  40  FORMAT (a,', final posterior error variance')
# endif
  50  FORMAT ('(',a,')^2')
  60  FORMAT (1pe11.4,1x,'millimeter')
  70  FORMAT (/,' DEF_ERROR_NF90 - unable to open error NetCDF',        &
     &        ' file: ',a)
  80  FORMAT (/,' DEF_ERROR_NF90 - unable to find variable: ',a,2x,     &
     &        ' in 4DVar error NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_error_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_error_pio (ng)
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
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, nrec, nvd, nvd3, nvd4
      integer :: recdim, status
      integer :: NinnerDim
#  ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
#  endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
#  ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
#  endif
#  ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#   ifdef ADJUST_BOUNDARY
      integer :: t3dobc(5)
#   endif
#   ifdef ADJUST_STFLUX
      integer :: t3dfrc(4)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      integer :: u3dfrc(4), v3dfrc(4)
#  endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_error_pio"
!
      TYPE (Var_desc_t) :: varDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=ERR(ng)%name
!
      IF (Master) THEN
        IF (LdefERR(ng)) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new posterior error covariance matrix file.
!=======================================================================
!
      DEFINE : IF (LdefERR(ng)) THEN
        CALL pio_netcdf_create (ng, iTLM, TRIM(ncname), ERR(ng)%pioFile)
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
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xi_rho',     &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xi_u',       &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xi_v',       &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xi_psi',     &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'eta_rho',    &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'eta_u',      &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'eta_v',      &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'eta_psi',    &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'IorJ',       &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xy_rho',     &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xy_u',       &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xy_v',       &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xyz_rho',    &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xyz_u',      &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xyz_v',      &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xyz_w',      &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'N',          &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 's_rho',      &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 's_w',        &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'tracer',     &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SEDIMENT
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'NST',        &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Nbed',       &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'xybed',      &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif

#   ifdef ECOSIM
        status=def_dim(ng, iNLM, ERR(ng)%pioFile, ncname, 'Nbands',     &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Nphy',       &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Nbac',       &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Ndom',       &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Nfec',       &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'boundary',   &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FOUR_DVAR
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Nstate',     &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'Ninner',     &
     &                 Ninner, NinnerDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'frc_adjust', &
     &                 Nfrec(ng), DimIDs(30))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname, 'obc_adjust', &
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, iTLM, ERR(ng)%pioFile, ncname,               &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        frecdim=DimIDs(30)
#  endif
#  ifdef ADJUST_BOUNDARY
        brecdim=DimIDs(31)
#  endif
!
!  Set number of dimensions for output variables.
!
#  if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#  else
        nvd3=3
        nvd4=4
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
#   ifdef ADJUST_STFLUX
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
        t3dfrc(4)=DimIDs(12)
#   endif
#  endif
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
#   ifdef ADJUST_WSTRESS
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
        u3dfrc(4)=DimIDs(12)
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
#   ifdef ADJUST_WSTRESS
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
        v3dfrc(4)=DimIDs(12)
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
#  endif
!
!  Initialize unlimited time record dimension.
!
        ERR(ng)%Rindex=0
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
        CALL def_info (ng, iTLM, ERR(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define inner-loop Lanczos vectors tridiagonal matrix.
!
        Vinfo( 1)='zLanczos_coef'
        Vinfo( 2)='inner-loop Lanczos vector tridiagonal matrix'
        status=def_var(ng, iTLM, ERR(ng)%pioFile, varDesc, PIO_TYPE,    &
     &                 2, (/NinnerDim,NinnerDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define inverse of inner-loop Lanczos vectors tridiagonal matrix.
!
        Vinfo( 1)='zLanczos_inv'
        Vinfo( 2)='inverse inner-loop Lanczos vector tridiagonal matrix'
        status=def_var(ng, iTLM, ERR(ng)%pioFile, varDesc, PIO_TYPE,    &
     &                 2, (/NinnerDim,NinnerDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define inner-loop Lanczos vectors tridiagonal matrix inversion error.
!  We need to get an identity matrix withing roundoff.
!
        Vinfo( 1)='zLanczos_err'
        Vinfo( 2)='inner-loop Lanczos vector tridiagonal matrix '//     &
     &            'inversion error'
        status=def_var(ng, iTLM, ERR(ng)%pioFile, varDesc, PIO_TYPE,    &
     &                 2, (/NinnerDim,NinnerDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        ERR(ng)%pioVar(idtime)%dkind=PIO_TOUT
        ERR(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface error variance.
!
        Vinfo( 1)=Vname(1,idFsur)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idFsur))
        Vinfo( 3)='meter2'
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(21)=Vname(6,idFsur)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        ERR(ng)%pioVar(idFsur)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idFsur)%gtype=r2dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idFsur)%vd,                       &
     &                 PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries error variance.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          ERR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          ERR(ng)%pioVar(ifield)%gtype=r2dobc
!
          status=def_var(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 4, t2dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D U-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idUbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbar))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(21)=Vname(6,idUbar)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        ERR(ng)%pioVar(idUbar)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idUbar)%gtype=u2dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idUbar)%vd,                       &
     &                 PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define 2D U-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          ERR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          ERR(ng)%pioVar(ifield)%gtype=u2dobc
!
          status=def_var(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 4, t2dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D V-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idVbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbar))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(21)=Vname(6,idVbar)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        ERR(ng)%pioVar(idVbar)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idVbar)%gtype=v2dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idVbar)%vd,                       &
     &                 PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define 2D V-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          ERR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          ERR(ng)%pioVar(ifield)%gtype=v2dobc
!
          status=def_var(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 4, t2dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef SOLVE3D
!
!  Define 3D U-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idUvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUvel))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#   endif
        Vinfo(21)=Vname(6,idUvel)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        ERR(ng)%pioVar(idUvel)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idUvel)%gtype=u3dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idUvel)%vd,                       &
     &                 PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef ADJUST_BOUNDARY
!
!  Define 3D U-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          ERR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          ERR(ng)%pioVar(ifield)%gtype=u3dobc
!
          status=def_var(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 5, t3dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define 3D V-momentum component error variance.
!
        Vinfo( 1)=Vname(1,idVvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvel))
        Vinfo( 3)='meter2 second-2'
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#   endif
        Vinfo(21)=Vname(6,idVvel)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        ERR(ng)%pioVar(idVvel)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idVvel)%gtype=v3dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idVvel)%vd,                       &
     &                 PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef ADJUST_BOUNDARY
!
!  Define 3D V-momentum component open boundaries error variance.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='meter2 second-2'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,ifield)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          ERR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          ERR(ng)%pioVar(ifield)%gtype=v3dobc
!
          status=def_var(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 5, t3dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define tracer type variables error variance.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTvar(itrc)))
          IF (itrc.eq.itemp) THEN
            Vinfo( 3)='Celsius2'
          ELSE IF (itrc.eq.isalt) THEN
            Vinfo( 3)='nondimensional'
          ELSE
            WRITE (Vinfo( 3),50) TRIM(Vname(3,idTvar(itrc)))
          END IF
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),60) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#   endif
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTvar(itrc))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          ERR(ng)%pioTrc(itrc)%dkind=PIO_FOUT
          ERR(ng)%pioTrc(itrc)%gtype=r3dvar
!
          status=def_var(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioTrc(itrc)%vd,                       &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

#   ifdef ADJUST_BOUNDARY
!
!  Define tracer type variables open boundaries error variance.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
            IF (itrc.eq.itemp) THEN
              Vinfo( 3)='Celsius2'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo( 3)='nondimensional'
            ELSE
              WRITE (Vinfo( 3),50) TRIM(Vname(3,idTvar(itrc)))
            END IF
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),60) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#    endif
            Vinfo(21)=Vname(6,ifield)
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            ERR(ng)%pioVar(ifield)%dkind=PIO_FOUT
            ERR(ng)%pioVar(ifield)%gtype=r3dobc
!
            status=def_var(ng, iTLM, ERR(ng)%pioFile,                   &
     &                     ERR(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, 5, t3dobc, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
#   ifdef ADJUST_STFLUX
!
!  Define surface tracer fluxes error variance.
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
            WRITE (Vinfo( 3),50) TRIM(Vname(3,idTsur(itrc)))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            ERR(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
            ERR(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
!
            status=def_var(ng, iTLM, ERR(ng)%pioFile,                   &
     &                     ERR(ng)%pioVar(idTsur(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, t3dfrc, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress error variance.
!
        Vinfo( 1)=Vname(1,idUsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
        Vinfo( 3)='meter4 second-4'
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#   endif
        Vinfo(21)=Vname(6,idUsms)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u2dvar,r8)
        ERR(ng)%pioVar(idUsms)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idUsms)%gtype=u2dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idUsms)%vd,                       &
     &                 PIO_FOUT, nvd4, u3dfrc, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define surface V-momentum stress error variance.
!
        Vinfo( 1)=Vname(1,idVsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)='meter4 second4'
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#   endif
        Vinfo(21)=Vname(6,idVsms)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v2dvar,r8)
        ERR(ng)%pioVar(idVsms)%dkind=PIO_FOUT
        ERR(ng)%pioVar(idVsms)%gtype=v2dvar
!
        status=def_var(ng, iTLM, ERR(ng)%pioFile,                       &
     &                 ERR(ng)%pioVar(idVsms)%vd,                       &
     &                 PIO_FOUT, nvd4, v3dfrc, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, iTLM, ncname, ERR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iTLM, ERR(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing posterior error covariance matrix, check its
!  contents, and prepare for appending data.
!=======================================================================
!
      QUERY: IF (.not.LdefERR(ng)) THEN
        ncname=ERR(ng)%name
!
!  Open posterior error covariance file for read/write.
!
        CALL pio_netcdf_open (ng, iTLM, ncname, 1, ERR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,70) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, iTLM, ncname,                    &
     &                             pioFile = ERR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, iTLM, ncname,                      &
     &                           pioFile = ERR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  posterior error covariance matrix variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            ERR(ng)%pioVar(idtime)%vd=var_desc(i)
            ERR(ng)%pioVar(idtime)%dkind=PIO_TOUT
            ERR(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            ERR(ng)%pioVar(idFsur)%vd=var_desc(i)
            ERR(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            ERR(ng)%pioVar(idUbar)%vd=var_desc(i)
            ERR(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            ERR(ng)%pioVar(idVbar)%vd=var_desc(i)
            ERR(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idVbar)%gtype=v2dvar
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            ERR(ng)%pioVar(idSbry(isFsur))%vd=var_desc(i)
            ERR(ng)%pioVar(idSbry(isFsur))%dkind=PIO_FOUT
            ERR(ng)%pioVar(idSbry(isFsur))%gtype=r2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            ERR(ng)%pioVar(idSbry(isUbar))%vd=var_desc(i)
            ERR(ng)%pioVar(idSbry(isUbar))%dkind=PIO_FOUT
            ERR(ng)%pioVar(idSbry(isUbar))%gtype=u2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            ERR(ng)%pioVar(idSbry(isVbar))%vd=var_desc(i)
            ERR(ng)%pioVar(idSbry(isVbar))%dkind=PIO_FOUT
            ERR(ng)%pioVar(idSbry(isVbar))%gtype=v2dobc
#  endif
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            ERR(ng)%pioVar(idUsms)%vd=var_desc(i)
            ERR(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            ERR(ng)%pioVar(idVsms)%vd=var_desc(i)
            ERR(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idVsms)%gtype=v2dvar
#  endif
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            ERR(ng)%pioVar(idUvel)%vd=var_desc(i)
            ERR(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            ERR(ng)%pioVar(idVvel)%vd=var_desc(i)
            ERR(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            ERR(ng)%pioVar(idVvel)%gtype=v3dvar
#   ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            ERR(ng)%pioVar(idSbry(isUvel))%vd=var_desc(i)
            ERR(ng)%pioVar(idSbry(isUvel))%dkind=PIO_FOUT
            ERR(ng)%pioVar(idSbry(isUvel))%gtype=u3dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            ERR(ng)%pioVar(idSbry(isVvel))%vd=var_desc(i)
            ERR(ng)%pioVar(idSbry(isVvel))%dkind=PIO_FOUT
            ERR(ng)%pioVar(idSbry(isVvel))%gtype=v3dobc
#   endif
#  endif
          END IF
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              ERR(ng)%pioTrc(itrc)%vd=var_desc(i)
              ERR(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              ERR(ng)%pioTrc(itrc)%gtype=r3dvar
#   ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              ERR(ng)%pioVar(idSbry(isTvar(itrc)))%vd=var_desc(i)
              ERR(ng)%pioVar(idSbry(isTvar(itrc)))%dkind=PIO_FOUT
              ERR(ng)%pioVar(idSbry(isTvar(itrc)))%gtype=r3dobc
#   endif
#   ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              ERR(ng)%pioVar(idTsur(itrc))%vd=var_desc(i)
              ERR(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              ERR(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
#   endif
            END IF
          END DO
#  endif
        END DO
!
!  Check if posterior error covariance matrix variables are available
!  in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel))) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,80) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,80)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   endif
#   ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,80) TRIM(Vname(1,idTsur(itrc))),  &
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
        ERR(ng)%Rindex=rec_size
      END IF QUERY
!
  10  FORMAT (2x,'DEF_ERROR_PIO    - creating error file,',t56,         &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_ERROR_PIO    - inquiring error file,',t56,        &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_ERROR_PIO - unable to create 4DVar error NetCDF'  &
     &          ' file:',1x,a)
#  if defined POSTERIOR_ERROR_I
  40  FORMAT (a,', initial posterior error variance')
#  elif defined POSTERIOR_ERROR_F
  40  FORMAT (a,', final posterior error variance')
#  endif
  50  FORMAT ('(',a,')^2')
  60  FORMAT (1pe11.4,1x,'millimeter')
  70  FORMAT (/,' DEF_ERROR_PIO - unable to open error NetCDF',         &
     &        ' file: ',a)
  80  FORMAT (/,' DEF_ERROR_PIO - unable to find variable: ',a,2x,      &
     &        ' in 4DVar error NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_error_pio
# endif
#endif
      END MODULE def_error_mod
