#include "cppdefs.h"
#if defined IS4DVAR              || \
    (defined WEAK_CONSTRAINT     && \
     (defined POSTERIOR_EOFS     || defined POSTERIOR_ERROR_I || \
      defined POSTERIOR_ERROR_F))
      SUBROUTINE def_hessian (ng)
!
!svn $Id: def_hessian.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates the Hessian eigenvectors NetCDF file used for  !
!  the pre-conditioning of the 4DVar conjugate gradient algorithm. It  !
!  defines its dimensions, attributes, and variables.                  !
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
      USE mod_netcdf
      USE mod_scalars
# ifdef SEDIMENT
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
      logical :: got_var(NV)

      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, nrec, nvd, nvd3, nvd4
      integer :: recdim, status, varid
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
# endif
# if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
      integer :: NpostDim
# endif
      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
# endif
      integer :: Vsize(4)

      integer :: def_dim

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

      real(r8) :: Aval(6)

      character (len=80) :: fname, ncname
      character (len=120) :: Vinfo(Natt)
!
      SourceFile='def_hessian.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=HSSname(ng)
!
      IF (Master) THEN
        IF (LdefHSS(ng)) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new Hessian eigenvectors file.
!=======================================================================
!
      DEFINE : IF (LdefHSS(ng)) THEN
        CALL netcdf_create (ng, iADM, TRIM(ncname), ncHSSid(ng))
        IF (exit_flag.ne.NoError) THEN
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
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xi_rho',         &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xi_u',           &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xi_v',           &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xi_psi',         &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'eta_rho',        &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'eta_u',          &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'eta_v',          &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'eta_psi',        &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'IorJ',           &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xy_rho',         &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xy_u',           &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xy_v',           &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xyz_rho',        &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xyz_u',          &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xyz_v',          &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xyz_w',          &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
#  endif

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'N',              &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 's_rho',          &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 's_w',            &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'tracer',         &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'NST',            &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Nbed',           &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'xybed',          &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Nphy',           &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Nbac',           &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Ndom',           &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Nfec',           &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'boundary',       &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Nstate',         &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'Nposterior',     &
     &                 NpostI+1, NpostDim)
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'frc_adjust',     &
     &                 Nfrec(ng), DimIDs(30))
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iADM, ncHSSid(ng), ncname, 'obc_adjust',     &
     &                 Nbrec(ng), DimIDs(31))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iADM, ncHSSid(ng), ncname,                   &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

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
!  Define dimension vectors for staggered u-momemtum type variables.
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
!  Define dimension vectors for staggered v-momemtum type variables.
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
!  Define dimension vector for staggered w-momemtum type variables.
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
        tHSSindx(ng)=0
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
        CALL def_info (ng, iADM, ncHSSid(ng), ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='number of converged Ritz eigenvalues'
        status=def_var(ng, iADM, ncHSSid(ng), varid, nf90_int,          &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Ritz eigenvalues.
!
        Vinfo( 1)='Ritz'
        Vinfo( 2)='Ritz eigenvalues'
        status=def_var(ng, iADM, ncHSSid(ng), varid, NF_TYPE,           &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define accuracy of Ritz eigenvalues .
!
        Vinfo( 1)='Ritz_error'
        Vinfo( 2)='accuracy of Ritz eigenvalues'
        status=def_var(ng, iADM, ncHSSid(ng), varid, NF_TYPE,           &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

# ifdef POSTERIOR_EOFS
!
!  Define posterior analysis error covariance matrix trace.
!
        Vinfo( 1)='ae_trace'
        Vinfo( 2)='posterior analysis error covariance matrix trace'
        status=def_var(ng, iADM, ncHSSid(ng), varid, NF_TYPE,           &
     &                 1, (/NpostDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
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
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idtime,ng),        &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo,ncname,      &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idFsur))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idFsur,ng),        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iADM, ncHSSid(ng), hssVid(ifield,ng),      &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define 2D U-momentum component.
!
        Vinfo( 1)=Vname(1,idUbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbar))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idUbar,ng),        &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define 2D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iADM, ncHSSid(ng), hssVid(ifield,ng),      &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define 2D V-momentum component.
!
        Vinfo( 1)=Vname(1,idVbar)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbar))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idVbar,ng),        &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Define 2D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iADM, ncHSSid(ng), hssVid(ifield,ng),      &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        Vinfo( 1)=Vname(1,idUvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUvel))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idUvel,ng),        &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define 3D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iADM, ncHSSid(ng), hssVid(ifield,ng),      &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
!
!  Define 3D V-momentum component.
!
        Vinfo( 1)=Vname(1,idVvel)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvel))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idVvel,ng),        &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

#  ifdef ADJUST_BOUNDARY
!
!  Define 3D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iADM, ncHSSid(ng), hssVid(ifield,ng),      &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTvar(itrc)))
          Vinfo( 3)='nondimensional'
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iADM, ncHSSid(ng), hssTid(itrc,ng),        &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO

#  ifdef ADJUST_BOUNDARY
!
!  Define tracer type variables open boundaries.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),40) TRIM(Vname(2,ifield))
            Vinfo( 3)='nondimensional'
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#   endif
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iADM, ncHSSid(ng), hssVid(ifield,ng),    &
     &                     NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Define surface tracer fluxes.
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
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
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, iADM, ncHSSid(ng),                       &
     &                     hssVid(idTsur(itrc),ng),                     &
     &                     NF_FOUT, nvd4, t3dfrc, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        Vinfo( 1)=Vname(1,idUsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u2dvar,r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idUsms,ng),        &
     &                 NF_FOUT, nvd4, u3dfrc, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface V-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)='nondimensional'
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v2dvar,r8)
        status=def_var(ng, iADM, ncHSSid(ng), hssVid(idVsms,ng),        &
     &                 NF_FOUT, nvd4, v3dfrc, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iADM, ncname, ncHSSid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iADM, ncHSSid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing Hessian eigenvectors file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      QUERY: IF (.not.LdefHSS(ng)) THEN
        ncname=HSSname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iADM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iADM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open Hessian eigenvectors file for read/write.
!
        CALL netcdf_open (ng, iADM, ncname, 1, ncHSSid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,60) TRIM(ncname)
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
!  Hessian eigenvectors variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            hssVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            hssVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            hssVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            hssVid(idVbar,ng)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            hssVid(idSbry(isFsur),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            hssVid(idSbry(isUbar),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            hssVid(idSbry(isVbar),ng)=var_id(i)
# endif
# ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            hssVid(idUsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            hssVid(idVsms,ng)=var_id(i)
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            hssVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            hssVid(idVvel,ng)=var_id(i)
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            hssVid(idSbry(isUvel),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            hssVid(idSbry(isVvel),ng)=var_id(i)
#  endif
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              hssTid(itrc,ng)=var_id(i)
#  ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              hssVid(idSbry(isTvar(itrc)),ng)=var_id(i)
#  endif
#  ifdef ADJUST_STFLUX
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              hssVid(idTsur(itrc),ng)=var_id(i)
#  endif
            END IF
          END DO
# endif
        END DO
!
!  Check if Hessian eigenvectors variables are available in input
!  NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,70)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
#  ifdef ADJUST_STFLUX
          IF (.not.got_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
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
        tHSSindx(ng)=rec_size
      END IF QUERY
!
  10  FORMAT (3x,'DEF_HESSIAN  - creating Hessian  file: ',a)
  20  FORMAT (3x,'DEF_HESSIAN  - inquiring Hessian file: ',a)
  30  FORMAT (/,' DEF_HESSIAN - unable to create Hessian NetCDF file:', &
     &        1x,a)
  40  FORMAT (a,', Hessian eigenvectors')
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_HESSIAN - unable to open Hessian NetCDF file: ',a)
  70  FORMAT (/,' DEF_HESSIAN - unable to find variable: ',a,2x,        &
     &        ' in Hessian NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_hessian
#else
      SUBROUTINE def_hessian
      RETURN
      END SUBROUTINE def_hessian
#endif
