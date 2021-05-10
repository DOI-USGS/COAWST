#include "cppdefs.h"
      SUBROUTINE def_state (ng, model, label, Lcreate, S, Lclose)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates an ocean state NetCDF file according to input  !
!  derived (TYPE T_IO) structure S. It defines its dimensions,         !
!  attributes, and variables.                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     label        Identification label (string)                       !
!     ldef         Switch to define a new file (T) or restart (F)      !
!     S            File derived type structure (TYPE T_IO)             !
!     Lclose       Switch to close (.TRUE.) NetCDF after created.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#ifdef SEDIMENT
      USE mod_sediment
#endif
!
      USE def_var_mod, ONLY : def_var
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lcreate, Lclose
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: label
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, nrec, nvd, nvd3, nvd4
      integer :: recdim, status, varid
#ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
#endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
#ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
#endif

      integer :: def_dim

#ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# ifdef ADJUST_BOUNDARY
      integer :: t3dobc(5)
# endif
# ifdef ADJUST_STFLUX
      integer :: t3dfrc(4)
# endif
#endif
#ifdef ADJUST_WSTRESS
      integer :: u3dfrc(4), v3dfrc(4)
#endif
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      IF (Master) THEN
        IF (Lcreate) THEN
          WRITE (stdout,10) TRIM(label), ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(label), ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new ocean state NetCDF file.
!=======================================================================
!
      DEFINE : IF (Lcreate) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), S(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(label), TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xi_rho',         &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xi_u',           &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xi_v',           &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xi_psi',         &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'eta_rho',        &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'eta_u',          &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'eta_v',          &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'eta_psi',        &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'IorJ',           &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xy_rho',         &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xy_u',           &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xy_v',           &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xyz_rho',        &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xyz_u',          &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xyz_v',          &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xyz_w',          &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'N',              &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 's_rho',          &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 's_w',            &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'tracer',         &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SEDIMENT
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'NST',            &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'Nbed',           &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'xybed',          &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef ECOSIM
        status=def_dim(ng, iNLM, S(ng)%ncid, ncname, 'Nbands',          &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'Nphy',           &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'Nbac',           &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'Ndom',           &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'Nfec',           &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

        status=def_dim(ng, model, S(ng)%ncid, ncname, 'boundary',       &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef FOUR_DVAR
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'Nstate',         &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'frc_adjust',     &
     &                 Nfrec(ng), DimIDs(30))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, S(ng)%ncid, ncname, 'obc_adjust',     &
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

        status=def_dim(ng, model, S(ng)%ncid, ncname,                   &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        frecdim=DimIDs(30)
#endif
#ifdef ADJUST_BOUNDARY
        brecdim=DimIDs(31)
#endif
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#else
        nvd3=3
        nvd4=4
#endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
# endif
#else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
# endif
# ifdef ADJUST_STFLUX
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
        t3dfrc(4)=DimIDs(12)
# endif
#endif
#ifdef ADJUST_BOUNDARY
        t2dobc(1)=IorJdim
        t2dobc(2)=DimIDs(14)
        t2dobc(3)=brecdim
        t2dobc(4)=DimIDs(12)
# ifdef SOLVE3D
        t3dobc(1)=IorJdim
        t3dobc(2)=DimIDs( 9)
        t3dobc(3)=DimIDs(14)
        t3dobc(4)=brecdim
        t3dobc(5)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered u-momentum type variables.
!
#if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
# endif
#else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
# endif
# ifdef ADJUST_WSTRESS
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
        u3dfrc(4)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
#if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
# endif
#else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
# endif
# ifdef ADJUST_WSTRESS
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
        v3dfrc(4)=DimIDs(12)
# endif
#endif
#ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        S(ng)%Rindex=0
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
        CALL def_info (ng, model, S(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idtime),        &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo,ncname,      &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idFsur),        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(ifield),      &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
!
!  Define 2D U-momentum component.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbar),        &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef ADJUST_BOUNDARY
!
!  Define 2D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(ifield),      &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
!
!  Define 2D V-momentum component.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbar),        &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef ADJUST_BOUNDARY
!
!  Define 2D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(ifield),      &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
#ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        Vinfo( 1)=Vname(1,idUvel)
        Vinfo( 2)=Vname(2,idUvel)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUvel),        &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define 3D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(ifield),      &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 3D V-momentum component.
!
        Vinfo( 1)=Vname(1,idVvel)
        Vinfo( 2)=Vname(2,idVvel)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVvel),        &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define 3D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(ifield),      &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
# endif
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Tid(itrc),        &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

# ifdef ADJUST_BOUNDARY
!
!  Define tracer type variables open boundaries.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)='nondimensional'
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(ifield),    &
     &                     NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
# ifdef ADJUST_STFLUX
!
!  Define surface tracer fluxes.
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)='nondimensional'
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idTsur(itrc)),                     &
     &                     NF_FOUT, nvd4, t3dfrc, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
#endif
#ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        Vinfo( 1)=Vname(1,idUsms)
        Vinfo( 2)=Vname(2,idUsms)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u2dvar,r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUsms),        &
     &                 NF_FOUT, nvd4, u3dfrc, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define surface V-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v2dvar,r8)
        status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVsms),        &
     &                 NF_FOUT, nvd4, v3dfrc, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, S(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, S(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        IF (Lclose) THEN
          CALL netcdf_close (ng, model, S(ng)%ncid, ncname, .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

      END IF DEFINE
!
!=======================================================================
!  Open an existing ocean state NetCDF file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      QUERY: IF (.not.Lcreate) THEN
        ncname=S(ng)%name
!
!  Open ocean state file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, S(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(label), TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = S(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = S(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Hessian eigenvectors variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            S(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            S(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            S(ng)%Vid(idVbar)=var_id(i)
#ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            S(ng)%Vid(idSbry(isFsur))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            S(ng)%Vid(idSbry(isUbar))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            S(ng)%Vid(idSbry(isVbar))=var_id(i)
#endif
#ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            S(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            S(ng)%Vid(idVsms)=var_id(i)
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            S(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            S(ng)%Vid(idVvel)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            S(ng)%Vid(idSbry(isUvel))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            S(ng)%Vid(idSbry(isVvel))=var_id(i)
# endif
#endif
          END IF
#ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              S(ng)%Tid(itrc)=var_id(i)
# ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              S(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
# endif
# ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              S(ng)%Vid(idTsur(itrc))=var_id(i)
# endif
            END IF
          END DO
#endif
        END DO
!
!  Check if ocean state variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur))) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar))) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar))) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel))) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(label), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel))) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(label), TRIM(ncname)
          RETURN
        END IF
# endif
#endif
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(label), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
# ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,60)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(label), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
# endif
# ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(label), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
# endif
        END DO
#endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        S(ng)%Rindex=rec_size
      END IF QUERY
!
  10  FORMAT (6x,'DEF_STATE - creating  ',a,t43,' file, Grid ',i2.2,    &
     &        ': ', a)
  20  FORMAT (6x,'DEF_STATE - inquiring ',a,t43,' file, Grid ',i2.2,    &
     &        ': ', a)
  30  FORMAT (/,' DEF_STATE - unable to create ',a,' NetCDF file:',1x,a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_STATE - unable to open ',a,' NetCDF file: ',a)
  60  FORMAT (/,' DEF_STATE - unable to find variable: ',a,2x,' in ',a, &
     &        ' NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_state
