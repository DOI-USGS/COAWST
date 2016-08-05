#include "cppdefs.h"
      SUBROUTINE def_his (ng, ldef)
!
!svn $Id: def_his.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates history NetCDF file, it defines its            !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
#ifdef VEGETATION
      USE mod_vegetation
      USE mod_vegarr
#endif
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
#ifdef INWAVE_MODEL
      USE mod_inwave_params
#endif
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

      integer, parameter :: Natt = 25

      integer :: i, j, ifield, itrc, nvd3, nvd4, varid
      integer :: recdim, status
#ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
#endif
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
#ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
#endif

#if defined VEGETATION && (defined VEG_DRAG || defined VEG_BIOMASS)
      integer :: DimIDs(35)
#elif defined ESTUARYBGC 
      integer :: DimIDs(34)
#elif defined INWAVE_MODEL 
      integer :: DimIDs(33)
#else
      integer :: DimIDs(32)
#endif
      integer :: Vsize(4)

      integer :: def_dim

#ifdef INWAVE_MODEL
      integer :: r3degrd(4)
      integer :: u3degrd(4)
      integer :: v3degrd(4)
      integer :: t3degrd(4)
#endif
#ifdef SOLVE3D
# ifdef SEDIMENT
      integer :: b3dgrd(4)
# endif
# if defined ESTUARYBGC
      integer :: s3dgrd(4)
# endif
#if defined VEGETATION 
# if defined VEG_DRAG || defined VEG_BIOMASS
      integer :: v3pgrd(4)
# endif
#endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# ifdef ADJUST_BOUNDARY
      integer :: t3dobc(5)
# endif
#endif
#ifdef WET_DRY
      integer :: sp2dgrd(3)
#endif
      real(r8) :: Aval(6)

      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_his.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=HIS(ng)%name
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
!  Create a new history file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), HIS(ng)%ncid)
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
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

#ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'IorJ',          &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (exit_flag.ne.NoError) RETURN
#endif

#if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef INWAVE_MODEL
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname,'ND',             &
                       ND(ng), DimIDs(33))
#endif

#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'N',             &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

# ifdef SEDIMENT
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# if defined ESTUARYBGC && defined SPECTRAL_LIGHT
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'NBAND',         &
     &                 NBAND, DimIDs(34))
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef ECOSIM
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
# endif
#endif

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

#ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
#endif

#if defined VEGETATION
# if defined VEG_DRAG || defined VEG_BIOMASS
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'NVEG',          &
     &                 NVEG, DimIDs(35))
        IF (exit_flag.ne.NoError) RETURN
# endif 
#endif

#ifdef ADJUST_BOUNDARY
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'obc_adjust',    &
     &                 Nbrec(ng), DimIDs(31))
        IF (exit_flag.ne.NoError) RETURN
#endif

        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
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
#ifdef WET_DRY
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        sp2dgrd(1)=DimIDs( 4)
        sp2dgrd(2)=DimIDs( 8)
        sp2dgrd(3)=DimIDs(12)
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
#endif
!
!  Define dimensions for the plant array variables
!
#ifdef VEGETATION 
# if defined VEG_DRAG || defined VEG_BIOMASS
        v3pgrd(1)=DimIDs( 1)
        v3pgrd(2)=DimIDs( 5)
        v3pgrd(3)=DimIDs(35)
        v3pgrd(4)=DimIDs(12)
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
# ifdef INWAVE_MODEL
!
!  Define dimension vector for energy wave layer type variables.
!
        r3degrd(1)=DimIDs( 1)
        r3degrd(2)=DimIDs( 5)
        r3degrd(3)=DimIDs(33)
        r3degrd(4)=DimIDs(12)

        u3degrd(1)=DimIDs( 2)
        u3degrd(2)=DimIDs( 6)
        u3degrd(3)=DimIDs(33)
        u3degrd(4)=DimIDs(12)

        v3degrd(1)=DimIDs( 3)
        v3degrd(2)=DimIDs( 7)
        v3degrd(3)=DimIDs(33)
        v3degrd(4)=DimIDs(12)

        t3degrd(1)=DimIDs( 1)
        t3degrd(2)=DimIDs( 5)
        t3degrd(3)=DimIDs(33)
        t3degrd(4)=DimIDs(12)
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
# ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#  if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#  else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#  endif
# endif
# if defined ESTUARYBGC 
!
!  Define dimension vector for spectral light type variables.
!
        s3dgrd(1)=DimIDs( 1)
        s3dgrd(2)=DimIDs( 5)
        s3dgrd(3)=DimIDs(34)
        s3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        HIS(ng)%Rindex=0
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
        CALL def_info (ng, iNLM, HIS(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
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
        status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

#ifdef SO_SEMI
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='real Ritz eigenvalues'
        status=def_var(ng, iNLM, HIS(ng)%ncid, varid, NF_TYPE,          &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ng, iNLM, HIS(ng)%ncid, varid, NF_TYPE,          &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#endif
#if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
        IF (Hout(idbath,ng)) THEN
          Vinfo( 1)=Vname(1,idbath)
          Vinfo( 2)=Vname(2,idbath)
          Vinfo( 3)=Vname(3,idbath)
          Vinfo(14)=Vname(4,idbath)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idbath),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef WET_DRY
!
!  Define wet/dry mask on PSI-points.
!
        Vinfo( 1)=Vname(1,idPwet)
        Vinfo( 2)=Vname(2,idPwet)
        Vinfo( 3)=Vname(3,idPwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idPwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idPwet,ng),r8)
        status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPwet),     &
     &                 NF_FOUT, nvd3, sp2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on RHO-points.
!
        Vinfo( 1)=Vname(1,idRwet)
        Vinfo( 2)=Vname(2,idRwet)
        Vinfo( 3)=Vname(3,idRwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idRwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRwet,ng),r8)
        status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRwet),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on U-points.
!
        Vinfo( 1)=Vname(1,idUwet)
        Vinfo( 2)=Vname(2,idUwet)
        Vinfo( 3)=Vname(3,idUwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idUwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUwet,ng),r8)
        status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUwet),     &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on V-points.
!
        Vinfo( 1)=Vname(1,idVwet)
        Vinfo( 2)=Vname(2,idVwet)
        Vinfo( 3)=Vname(3,idVwet)
        Vinfo(14)=Vname(4,idVwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVwet,ng),r8)
        status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVwet),     &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SOLVE3D
!
!  Define time-varying depth of RHO-points.
!
        IF (Hout(idpthR,ng)) THEN
          Vinfo( 1)=Vname(1,idpthR)
          WRITE (Vinfo( 2),40) Vname(2,idpthR)
          Vinfo( 3)=Vname(3,idpthR)
          Vinfo(14)=Vname(4,idpthR)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthR,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idpthR),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Hout(idpthU,ng)) THEN
          Vinfo( 1)=Vname(1,idpthU)
          WRITE (Vinfo( 2),40) Vname(2,idpthU)
          Vinfo( 3)=Vname(3,idpthU)
          Vinfo(14)=Vname(4,idpthU)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthU,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idpthU),   &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Hout(idpthV,ng)) THEN
          Vinfo( 1)=Vname(1,idpthV)
          WRITE (Vinfo( 2),40) Vname(2,idpthV)
          Vinfo( 3)=Vname(3,idpthV)
          Vinfo(14)=Vname(4,idpthV)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthV,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idpthV),   &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Hout(idpthW,ng)) THEN
          Vinfo( 1)=Vname(1,idpthW)
          WRITE (Vinfo( 2),40) Vname(2,idpthW)
          Vinfo( 3)=Vname(3,idpthW)
          Vinfo(14)=Vname(4,idpthW)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthW,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idpthW),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idFsur),   &
#ifdef WET_DRY
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
#else
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
#endif
          IF (exit_flag.ne.NoError) RETURN
#if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRzet)
          Vinfo( 2)=Vname(2,idRzet)
          Vinfo( 3)=Vname(3,idRzet)
          Vinfo(14)=Vname(4,idRzet)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRzet),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#endif
        END IF
#ifdef ADJUST_BOUNDARY
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
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Define 2D U-momentum component.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbar),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRu2d)
          Vinfo( 2)=Vname(2,idRu2d)
          Vinfo( 3)=Vname(3,idRu2d)
          Vinfo(14)=Vname(4,idRu2d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRu2d),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRuct)
          Vinfo( 2)=Vname(2,idRuct)
          Vinfo( 3)=Vname(3,idRuct)
          Vinfo(14)=Vname(4,idRuct)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRuct),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#  endif
          Vinfo( 1)=Vname(1,idUfx1)
          Vinfo( 2)=Vname(2,idUfx1)
          Vinfo( 3)=Vname(3,idUfx1)
          Vinfo(14)=Vname(4,idUfx1)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUfx1),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo( 1)=Vname(1,idUfx2)
          Vinfo( 2)=Vname(2,idUfx2)
          Vinfo( 3)=Vname(3,idUfx2)
          Vinfo(14)=Vname(4,idUfx2)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUfx2),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
#endif
        END IF
#ifdef ADJUST_BOUNDARY
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
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Define 2D V-momentum component.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbar),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRv2d)
          Vinfo( 2)=Vname(2,idRv2d)
          Vinfo( 3)=Vname(3,idRv2d)
          Vinfo(14)=Vname(4,idRv2d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRv2d),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRvct)
          Vinfo( 2)=Vname(2,idRvct)
          Vinfo( 3)=Vname(3,idRvct)
          Vinfo(14)=Vname(4,idRvct)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRvct),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#  endif
          Vinfo( 1)=Vname(1,idVfx1)
          Vinfo( 2)=Vname(2,idVfx1)
          Vinfo( 3)=Vname(3,idVfx1)
          Vinfo(14)=Vname(4,idVfx1)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVfx1),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo( 1)=Vname(1,idVfx2)
          Vinfo( 2)=Vname(2,idVfx2)
          Vinfo( 3)=Vname(3,idVfx2)
          Vinfo(14)=Vname(4,idVfx2)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVfx2),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
#endif
        END IF
#ifdef ADJUST_BOUNDARY
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
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(ifield),   &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Hout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)='barotropic_eastward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu2dE),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Hout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)='barotropic_northward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv2dN),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUvel),   &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRu3d)
          Vinfo( 2)=Vname(2,idRu3d)
          Vinfo( 3)=Vname(3,idRu3d)
          Vinfo(14)=Vname(4,idRu3d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRu3d),   &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
        END IF
# ifdef ADJUST_BOUNDARY
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
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(ifield),   &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define 3D V-momentum component.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvel),   &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRv3d)
          Vinfo( 2)=Vname(2,idRv3d)
          Vinfo( 3)=Vname(3,idRv3d)
          Vinfo(14)=Vname(4,idRv3d)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRv3d),   &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
        END IF
# ifdef ADJUST_BOUNDARY
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
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(ifield),   &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Hout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)='eastward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu3dE),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Hout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)='northward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv3dN),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)='upward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWvel),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idOvel),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Tid(itrc),   &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
          END IF
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
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(ifield), &
     &                     NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idDano),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# ifdef NEMURO_SED1
!
!  Define PON in sediment.
!
        IF (Hout(idPONsed,ng)) THEN
          Vinfo( 1)=Vname(1,idPONsed)
          Vinfo( 2)=Vname(2,idPONsed)
          Vinfo( 3)=Vname(3,idPONsed)
          Vinfo(14)=Vname(4,idPONsed)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPONsed,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPONsed), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define OPAL in sediment.
!
        IF (Hout(idOPALsed,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALsed)
          Vinfo( 2)=Vname(2,idOPALsed)
          Vinfo( 3)=Vname(3,idOPALsed)
          Vinfo(14)=Vname(4,idOPALsed)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOPALsed,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idOPALsed),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define DENIT in sediment.
!
        IF (Hout(idDENITsed,ng)) THEN
          Vinfo( 1)=Vname(1,idDENITsed)
          Vinfo( 2)=Vname(2,idDENITsed)
          Vinfo( 3)=Vname(3,idDENITsed)
          Vinfo(14)=Vname(4,idDENITsed)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDENITsed,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid,                        &
     &                   HIS(ng)%Vid(idDENITsed),                       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define buried PON in sediment.
!
        IF (Hout(idPONbur,ng)) THEN
          Vinfo( 1)=Vname(1,idPONbur)
          Vinfo( 2)=Vname(2,idPONbur)
          Vinfo( 3)=Vname(3,idPONbur)
          Vinfo(14)=Vname(4,idPONbur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPONbur,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPONbur), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define buried OPAL in sediment.
!
        IF (Hout(idOPALbur,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALbur)
          Vinfo( 2)=Vname(2,idOPALbur)
          Vinfo( 3)=Vname(3,idOPALbur)
          Vinfo(14)=Vname(4,idOPALbur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOPALbur,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idOPALbur),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef PRIMARY_PROD
!
!  Define Net primary productivity.
!
        IF (Hout(idNPP,ng)) THEN
          Vinfo( 1)=Vname(1,idNPP)
          Vinfo( 2)=Vname(2,idNPP)
          Vinfo( 3)=Vname(3,idNPP)
          Vinfo(14)=Vname(4,idNPP)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idNPP,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idNPP),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idHsbl),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idHbbl),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Define out KPP nonlocal transport.
!
        DO itrc=1,NAT
          IF (Hout(idGhat(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idGhat(itrc))
            Vinfo( 2)=Vname(2,idGhat(itrc))
            Vinfo( 3)=Vname(3,idGhat(itrc))
            Vinfo(14)=Vname(4,idGhat(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idGhat(itrc),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idGhat(itrc)), NF_FOUT,          &
     &                     nvd4, w3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
!
!  Define sea surface salinity flux correction
!
        IF (Hout(idSSSf,ng)) THEN
          Vinfo( 1)=Vname(1,idSSSf)
          Vinfo( 2)=Vname(2,idSSSf)
          Vinfo( 3)=Vname(3,idSSSf)
          Vinfo(14)=Vname(4,idSSSf)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSSSf,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSSSf),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvis),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTdif),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSdif),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idMtke),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
#  ifdef FORWARD_WRITE
          Vinfo( 1)=Vname(1,idVmKK)
          Vinfo( 2)=Vname(2,idVmKK)
          Vinfo( 3)=Vname(3,idVmKK)
          Vinfo(14)=Vname(4,idVmKK)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVmKK),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#  endif
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idMtls),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
#  ifdef FORWARD_WRITE
          Vinfo( 1)=Vname(1,idVmLS)
          Vinfo( 2)=Vname(2,idVmLS)
          Vinfo( 3)=Vname(3,idVmLS)
          Vinfo(14)=Vname(4,idVmLS)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVmLS),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#  endif
#  if defined FORWARD_WRITE && defined GLS_MIXING
          Vinfo( 1)=Vname(1,idVmKP)
          Vinfo( 2)=Vname(2,idVmKP)
          Vinfo( 3)=Vname(3,idVmKP)
          Vinfo(14)=Vname(4,idVmKP)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKP,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVmKP),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
#  endif
        END IF
# endif
# if defined RUNOFF
!
!  Define surface air pressure.
!
        IF (Hout(idRunoff,ng)) THEN
          Vinfo( 1)=Vname(1,idRunoff)
          Vinfo( 2)=Vname(2,idRunoff)
          Vinfo( 3)=Vname(3,idRunoff)
          Vinfo(14)=Vname(4,idRunoff)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRunoff,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idRunoff), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
  || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
!
!  Define surface air pressure.
!
        IF (Hout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM2OCN_FLUXES
!
!  Define surface u-wind component.
!
        IF (Hout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface v-wind component.
!
        IF (Hout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Eastward surface winds at RHO-points.
!
        IF (Hout(idUairE,ng)) THEN
          Vinfo( 1)=Vname(1,idUairE)
          Vinfo( 2)=Vname(2,idUairE)
          Vinfo( 3)=Vname(3,idUairE)
          Vinfo(14)=Vname(4,idUairE)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUairE,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUairE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Northward surface winds at RHO-points.
!
        IF (Hout(idVairN,ng)) THEN
          Vinfo( 1)=Vname(1,idVairN)
          Vinfo( 2)=Vname(2,idVairN)
          Vinfo( 3)=Vname(3,idVairN)
          Vinfo(14)=Vname(4,idVairN)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVairN,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVairN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
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
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idTsur(itrc)), NF_FOUT,          &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLhea),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idShea),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLrad),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef EMINUSP
!
!  Define E-P flux (as computed by bulk_flux.F).
!
        IF (Hout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idEmPf),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idevap),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idrain),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
# ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSrad),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
#endif
!
# if defined ESTUARYBGC && defined SPECTRAL_LIGHT
!  Define Photosynthetically Available Radiation (PAR)
!
        IF (Hout(idPARo,ng)) THEN
          Vinfo( 1)=Vname(1,idPARo)
          Vinfo( 2)=Vname(2,idPARo)
          Vinfo( 3)=Vname(3,idPARo)
          Vinfo(14)=Vname(4,idPARo)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPARo,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPARo),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!  Define Spectral Photosynthetically Available Radiation (PARs)
!
        IF (Hout(idPARs,ng)) THEN
          Vinfo( 1)=Vname(1,idPARs)
          Vinfo( 2)=Vname(2,idPARs)
          Vinfo( 3)=Vname(3,idPARs)
          Vinfo(14)=Vname(4,idPARs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPARs,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPARs),   &
     &                   NF_FOUT, nvd4, s3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  Define Spectral attenuation (SpKd)
!
        IF (Hout(idSpKd,ng)) THEN
          Vinfo( 1)=Vname(1,idSpKd)
          Vinfo( 2)=Vname(2,idSpKd)
          Vinfo( 3)=Vname(3,idSpKd)
          Vinfo(14)=Vname(4,idSpKd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSpKd,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSpKd),   &
     &                   NF_FOUT, nvd4, s3dgrd, Aval, Vinfo, ncname)
        END IF
# endif
!
# if defined ESTUARYBGC && defined SAV_BIOMASS
!
!  Define Dissolved Inorganic Nitrogen in water column
!
        IF (Hout(iddinw,ng)) THEN 
          Vinfo( 1)=Vname(1,iddinw)
          Vinfo( 2)=Vname(2,iddinw)
          Vinfo( 3)=Vname(3,iddinw)
          Vinfo(14)=Vname(4,iddinw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,iddinw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(iddinw),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  Define Dissolved Inorganic Nitrogen in sediment column
!
        IF (Hout(iddins,ng)) THEN 
          Vinfo( 1)=Vname(1,iddins)
          Vinfo( 2)=Vname(2,iddins)
          Vinfo( 3)=Vname(3,iddins)
          Vinfo(14)=Vname(4,iddins)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,iddins,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(iddins),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  Define Dissolved Oxygen in water column
!
        IF (Hout(iddowc,ng)) THEN 
          Vinfo( 1)=Vname(1,iddowc)
          Vinfo( 2)=Vname(2,iddowc)
          Vinfo( 3)=Vname(3,iddowc)
          Vinfo(14)=Vname(4,iddowc)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,iddowc,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(iddowc),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  Define Dissolved Inorganic Nitrogen in water column due to SAV 
!
        IF (Hout(idwsvl,ng)) THEN 
          Vinfo( 1)=Vname(1,idwsvl)
          Vinfo( 2)=Vname(2,idwsvl)
          Vinfo( 3)=Vname(3,idwsvl)
          Vinfo(14)=Vname(4,idwsvl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idwsvl,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idwsvl),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!
!  Above ground biomass
!
        IF (Hout(idsagb,ng)) THEN
          Vinfo( 1)=Vname(1,idsagb)
          Vinfo( 2)=Vname(2,idsagb)
          Vinfo( 3)=Vname(3,idsagb)
          Vinfo(14)=Vname(4,idsagb)
          Vinfo(16)=Vname(1,idsagb)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsagb,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsagb),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        END IF
!
!  Above ground biomass
!
        IF (Hout(idsbgb,ng)) THEN
          Vinfo( 1)=Vname(1,idsbgb)
          Vinfo( 2)=Vname(2,idsbgb)
          Vinfo( 3)=Vname(3,idsbgb)
          Vinfo(14)=Vname(4,idsbgb)
          Vinfo(16)=Vname(1,idsbgb)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsbgb,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsbgb),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        END IF
!
!  PP
!
        IF (Hout(idsvpp,ng)) THEN
          Vinfo( 1)=Vname(1,idsvpp)
          Vinfo( 2)=Vname(2,idsvpp)
          Vinfo( 3)=Vname(3,idsvpp)
          Vinfo(14)=Vname(4,idsvpp)
          Vinfo(16)=Vname(1,idsvpp)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsvpp,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsvpp),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  AGM
!
        IF (Hout(idsvam,ng)) THEN
          Vinfo( 1)=Vname(1,idsvam)
          Vinfo( 2)=Vname(2,idsvam)
          Vinfo( 3)=Vname(3,idsvam)
          Vinfo(14)=Vname(4,idsvam)
          Vinfo(16)=Vname(1,idsvam)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsvam,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsvam),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  AGAR
!
        IF (Hout(idsgar,ng)) THEN
          Vinfo( 1)=Vname(1,idsgar)
          Vinfo( 2)=Vname(2,idsgar)
          Vinfo( 3)=Vname(3,idsgar)
          Vinfo(14)=Vname(4,idsgar)
          Vinfo(16)=Vname(1,idsgar)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsgar,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsgar),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  AGBR
!
        IF (Hout(idsvbr,ng)) THEN
          Vinfo( 1)=Vname(1,idsvbr)
          Vinfo( 2)=Vname(2,idsvbr)
          Vinfo( 3)=Vname(3,idsvbr)
          Vinfo(14)=Vname(4,idsvbr)
          Vinfo(16)=Vname(1,idsvbr)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsvbr,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsvbr),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
!  SEARS
!
        IF (Hout(idsvrs,ng)) THEN
          Vinfo( 1)=Vname(1,idsvrs)
          Vinfo( 2)=Vname(2,idsvrs)
          Vinfo( 3)=Vname(3,idsvrs)
          Vinfo(14)=Vname(4,idsvrs)
          Vinfo(16)=Vname(1,idsvrs)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsvrs,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsvrs),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
! AGBG
! 
        IF (Hout(idsvbg,ng)) THEN
          Vinfo( 1)=Vname(1,idsvbg)
          Vinfo( 2)=Vname(2,idsvbg)
          Vinfo( 3)=Vname(3,idsvbg)
          Vinfo(14)=Vname(4,idsvbg)
          Vinfo(16)=Vname(1,idsvbg)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsvbg,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsvbg),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
! BGAG
! 
        IF (Hout(idsvag,ng)) THEN
          Vinfo( 1)=Vname(1,idsvag)
          Vinfo( 2)=Vname(2,idsvag)
          Vinfo( 3)=Vname(3,idsvag)
          Vinfo(14)=Vname(4,idsvag)
          Vinfo(16)=Vname(1,idsvag)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsvag,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsvag),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
! BGR 
! 
        IF (Hout(idsbgr,ng)) THEN
          Vinfo( 1)=Vname(1,idsbgr)
          Vinfo( 2)=Vname(2,idsbgr)
          Vinfo( 3)=Vname(3,idsbgr)
          Vinfo(14)=Vname(4,idsbgr)
          Vinfo(16)=Vname(1,idsbgr)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsbgr,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsbgr),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
! BGM
! 
        IF (Hout(idsbgm,ng)) THEN
          Vinfo( 1)=Vname(1,idsbgm)
          Vinfo( 2)=Vname(2,idsbgm)
          Vinfo( 3)=Vname(3,idsbgm)
          Vinfo(14)=Vname(4,idsbgm)
          Vinfo(16)=Vname(1,idsbgm)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idsbgm,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idsbgm),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        END IF
!
# endif
!
#if defined VEGETATION 
#include "vegetation_def_his.h"
#endif
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUsms),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVsms),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbms),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbms),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        
#ifdef INWAVE_MODEL
!
!  Define 3D Action density wave energy.
!
        IF (Hout(idACen,ng)) THEN
          Vinfo( 1)=Vname(1,idACen)
          Vinfo( 2)=Vname(2,idACen)
          Vinfo( 3)=Vname(3,idACen)
          Vinfo(14)=Vname(4,idACen)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idACen),   &
     &                   NF_FOUT, nvd4, r3degrd, Aval, Vinfo, ncname)
        END IF
!
!  Define 3D wave celerity xi-direction.
!
        IF (Hout(idACcx,ng)) THEN
          Vinfo( 1)=Vname(1,idACcx)
          Vinfo( 2)=Vname(2,idACcx)
          Vinfo( 3)=Vname(3,idACcx)
          Vinfo(14)=Vname(4,idACcx)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idACcx),   &
     &                   NF_FOUT, nvd4, u3degrd, Aval, Vinfo, ncname)
        END IF
!
!  Define 3D wave celerity eta-direction.
!
        IF (Hout(idACcy,ng)) THEN
          Vinfo( 1)=Vname(1,idACcy)
          Vinfo( 2)=Vname(2,idACcy)
          Vinfo( 3)=Vname(3,idACcy)
          Vinfo(14)=Vname(4,idACcy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idACcy),   &
     &                   NF_FOUT, nvd4,v3degrd, Aval, Vinfo, ncname)
        END IF
!
!  Define 3D wave celerity theta-direction.
!
        IF (Hout(idACct,ng)) THEN
          Vinfo( 1)=Vname(1,idACct)
          Vinfo( 2)=Vname(2,idACct)
          Vinfo( 3)=Vname(3,idACct)
          Vinfo(14)=Vname(4,idACct)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idACct),   &
     &                   NF_FOUT, nvd4, r3degrd, Aval, Vinfo, ncname)
        END IF
!
!  Define Peak period.
!
        IF (Hout(idACtp,ng)) THEN
          Vinfo( 1)=Vname(1,idACtp)
          Vinfo( 2)=Vname(2,idACtp)
          Vinfo( 3)=Vname(3,idACtp)
          Vinfo(14)=Vname(4,idACtp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idACtp),   &
     &                   NF_FOUT, nvd4, t3degrd, Aval, Vinfo, ncname)
        END IF
#endif
#ifdef SOLVE3D
# ifdef UV_KIRBY
!
!  Define coupling U-current given by Kirby and Chen (1989).
!
        IF (Hout(idUwav,ng)) THEN
          Vinfo( 1)=Vname(1,idUwav)
          Vinfo( 2)=Vname(2,idUwav)
          Vinfo( 3)=Vname(3,idUwav)
          Vinfo(14)=Vname(4,idUwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUwav,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUwav),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        END IF
!
!  Define coupling V-current given by Kirby and Chen (1989).
!
        IF (Hout(idVwav,ng)) THEN
          Vinfo( 1)=Vname(1,idVwav)
          Vinfo( 2)=Vname(2,idVwav)
          Vinfo( 3)=Vname(3,idVwav)
          Vinfo(14)=Vname(4,idVwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVwav,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVWav),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        END IF
# endif
#endif
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Define bottom U-current stress.
!
        IF (Hout(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          Vinfo( 2)=Vname(2,idUbrs)
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbrs,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbrs),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-current stress.
!
        IF (Hout(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          Vinfo( 2)=Vname(2,idVbrs)
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbrs,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbrs),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (Hout(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          Vinfo( 2)=Vname(2,idUbws)
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbws,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbws),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (Hout(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          Vinfo( 2)=Vname(2,idVbws)
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbws,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbws),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (Hout(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          Vinfo( 2)=Vname(2,idUbcs)
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbcs,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbcs),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (Hout(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          Vinfo( 2)=Vname(2,idVbcs)
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbcs,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbcs),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define maximum wave and current bottom stress magnitude.
!
        IF (Hout(idUVwc,ng)) THEN
          Vinfo( 1)=Vname(1,idUVwc)
          Vinfo( 2)=Vname(2,idUVwc)
          Vinfo( 3)=Vname(3,idUVwc)
          Vinfo(14)=Vname(4,idUVwc)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUVwc,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUVwc),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (Hout(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          Vinfo( 2)=Vname(2,idUbot)
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbot,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbot),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (Hout(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          Vinfo( 2)=Vname(2,idVbot)
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbot,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbot),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (Hout(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          Vinfo( 2)=Vname(2,idUbur)
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbur,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbur),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (Hout(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          Vinfo( 2)=Vname(2,idVbvr)
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbvr,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbvr),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Define Bedload U-direction.
!
        DO i=1,NST
          IF (Hout(idUbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idUbld(i))
            Vinfo( 2)=Vname(2,idUbld(i))
            Vinfo( 3)=Vname(3,idUbld(i))
            Vinfo(14)=Vname(4,idUbld(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idUbld(i)), NF_FOUT,             &
     &                     nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Define Bedload V-direction.
!
          IF (Hout(idVbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idVbld(i))
            Vinfo( 2)=Vname(2,idVbld(i))
            Vinfo( 3)=Vname(3,idVbld(i))
            Vinfo(14)=Vname(4,idVbld(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idVbld(i)), NF_FOUT,             &
     &                     nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          IF (Hout(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            Vinfo( 2)=Vname(2,idfrac(i))
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idfrac(i),ng))
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idfrac(i)), NF_FOUT,             &
     &                     nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define sediment mass of each size class in each bed layer.
!
        DO i=1,NST
          IF (Hout(idBmas(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBmas(i))
            Vinfo( 2)=Vname(2,idBmas(i))
            Vinfo( 3)=Vname(3,idBmas(i))
            Vinfo(14)=Vname(4,idBmas(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idBmas(i),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idBmas(i)), NF_FOUT,             &
     &                     nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (Hout(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            Vinfo( 2)=Vname(2,idSbed(i))
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idSbed(i),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idSbed(i)), NF_FOUT,             &
     &                     nvd4, b3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Define exposed sediment layer properties.
!
        DO i=1,MBOTP
          IF (Hout(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            Vinfo( 2)=Vname(2,idBott(i))
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idBott(i)), NF_FOUT,             &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
#endif
#ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        IF (Hout(idUice,ng)) THEN
          Vinfo( 1)=Vname(1,idUice)
          Vinfo( 2)=Vname(2,idUice)
          Vinfo( 3)=Vname(3,idUice)
          Vinfo(14)=Vname(4,idUice)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUice,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUice),   &
     &                     NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Hout(idVice,ng)) THEN
          Vinfo( 1)=Vname(1,idVice)
          Vinfo( 2)=Vname(2,idVice)
          Vinfo( 3)=Vname(3,idVice)
          Vinfo(14)=Vname(4,idVice)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVice,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVice),   &
     &                     NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Eastward ice momentum component at RHO-points.
!
        IF (Hout(idUiceE,ng)) THEN
          Vinfo( 1)=Vname(1,idUiceE)
          Vinfo( 2)=Vname(2,idUiceE)
          Vinfo( 3)=Vname(3,idUiceE)
          Vinfo(14)=Vname(4,idUiceE)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUiceE,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUiceE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Northward ice momentum component at RHO-points.
!
        IF (Hout(idViceN,ng)) THEN
          Vinfo( 1)=Vname(1,idViceN)
          Vinfo( 2)=Vname(2,idViceN)
          Vinfo( 3)=Vname(3,idViceN)
          Vinfo(14)=Vname(4,idViceN)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idViceN,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idViceN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice concentration.
!
        IF (Hout(idAice,ng)) THEN
          Vinfo( 1)=Vname(1,idAice)
          Vinfo( 2)=Vname(2,idAice)
          Vinfo( 3)=Vname(3,idAice)
          Vinfo(14)=Vname(4,idAice)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAice,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idAice),   &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice average thickness.
!
        IF (Hout(idHice,ng)) THEN
          Vinfo( 1)=Vname(1,idHice)
          Vinfo( 2)=Vname(2,idHice)
          Vinfo( 3)=Vname(3,idHice)
          Vinfo(14)=Vname(4,idHice)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHice,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idHice),   &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Define ice/snow surface temperature.
!
        IF (Hout(idTice,ng)) THEN
          Vinfo( 1)=Vname(1,idTice)
          Vinfo( 2)=Vname(2,idTice)
          Vinfo( 3)=Vname(3,idTice)
          Vinfo(14)=Vname(4,idTice)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTice,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTice),   &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define snow thickness.
!
        IF (Hout(idHsno,ng)) THEN
          Vinfo( 1)=Vname(1,idHsno)
          Vinfo( 2)=Vname(2,idHsno)
          Vinfo( 3)=Vname(3,idHsno)
          Vinfo(14)=Vname(4,idHsno)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsno,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idHsno),   &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#ifdef MELT_PONDS
!
!  Define surface water fraction (on ice).
!
        IF (Hout(idHpond,ng)) THEN
          Vinfo( 1)=Vname(1,idApond)
          Vinfo( 2)=Vname(2,idApond)
          Vinfo( 3)=Vname(3,idApond)
          Vinfo(14)=Vname(4,idApond)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idApond,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idApond),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface water thickness (on ice).
!
        IF (Hout(idHpond,ng)) THEN
          Vinfo( 1)=Vname(1,idHpond)
          Vinfo( 2)=Vname(2,idHpond)
          Vinfo( 3)=Vname(3,idHpond)
          Vinfo(14)=Vname(4,idHpond)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHpond,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idHpond),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Define ice age.
!
        IF (Hout(idAgeice,ng)) THEN
          Vinfo( 1)=Vname(1,idAgeice)
          Vinfo( 2)=Vname(2,idAgeice)
          Vinfo( 3)=Vname(3,idAgeice)
          Vinfo(14)=Vname(4,idAgeice)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAgeice,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idAgeice), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-ocean mass flux
!
        IF (Hout(idIomflx,ng)) THEN
          Vinfo( 1)=Vname(1,idIomflx)
          Vinfo( 2)=Vname(2,idIomflx)
          Vinfo( 3)=Vname(3,idIomflx)
          Vinfo(14)=Vname(4,idIomflx)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIomflx,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idIomflx), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice internal temperature.
!
        IF (Hout(idTimid,ng)) THEN
          Vinfo( 1)=Vname(1,idTimid)
          Vinfo( 2)=Vname(2,idTimid)
          Vinfo( 3)=Vname(3,idTimid)
          Vinfo(14)=Vname(4,idTimid)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTimid,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTimid),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress component 11.
!
        IF (Hout(idSig11,ng)) THEN
          Vinfo( 1)=Vname(1,idSig11)
          Vinfo( 2)=Vname(2,idSig11)
          Vinfo( 3)=Vname(3,idSig11)
          Vinfo(14)=Vname(4,idSig11)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig11,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSig11),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress component 22.
!
        IF (Hout(idSig22,ng)) THEN
          Vinfo( 1)=Vname(1,idSig22)
          Vinfo( 2)=Vname(2,idSig22)
          Vinfo( 3)=Vname(3,idSig22)
          Vinfo(14)=Vname(4,idSig22)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig22,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSig22),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress component 12.
!
        IF (Hout(idSig12,ng)) THEN
          Vinfo( 1)=Vname(1,idSig12)
          Vinfo( 2)=Vname(2,idSig12)
          Vinfo( 3)=Vname(3,idSig12)
          Vinfo(14)=Vname(4,idSig12)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig12,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSig12),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water friction velocity.
!
        IF (Hout(idTauiw,ng)) THEN
          Vinfo( 1)=Vname(1,idTauiw)
          Vinfo( 2)=Vname(2,idTauiw)
          Vinfo( 3)=Vname(3,idTauiw)
          Vinfo(14)=Vname(4,idTauiw)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTauiw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTauiw),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water momentum transfer coefficient.
!
        IF (Hout(idChuiw,ng)) THEN
          Vinfo( 1)=Vname(1,idChuiw)
          Vinfo( 2)=Vname(2,idChuiw)
          Vinfo( 3)=Vname(3,idChuiw)
          Vinfo(14)=Vname(4,idChuiw)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idChuiw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idChuiw),  &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define salinity of molecular sub-layer under ice.
!
        IF (Hout(idS0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idS0mk)
          Vinfo( 2)=Vname(2,idS0mk)
          Vinfo( 3)=Vname(3,idS0mk)
          Vinfo(14)=Vname(4,idS0mk)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idS0mk,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idS0mk), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define temperature of molecular sub-layer under ice.
!
        IF (Hout(idT0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idT0mk)
          Vinfo( 2)=Vname(2,idT0mk)
          Vinfo( 3)=Vname(3,idT0mk)
          Vinfo(14)=Vname(4,idT0mk)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idT0mk,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idT0mk), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice freeze wfr
!
        IF (Hout(idWfr,ng)) THEN
          Vinfo( 1)=Vname(1,idWfr)
          Vinfo( 2)=Vname(2,idWfr)
          Vinfo( 3)=Vname(3,idWfr)
          Vinfo(14)=Vname(4,idWfr)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWfr,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWfr),    &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice freeze/melt wai
!
        IF (Hout(idWai,ng)) THEN
          Vinfo( 1)=Vname(1,idWai)
          Vinfo( 2)=Vname(2,idWai)
          Vinfo( 3)=Vname(3,idWai)
          Vinfo(14)=Vname(4,idWai)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWai,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWai),    &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice freeze/melt wao
!
        IF (Hout(idWao,ng)) THEN
          Vinfo( 1)=Vname(1,idWao)
          Vinfo( 2)=Vname(2,idWao)
          Vinfo( 3)=Vname(3,idWao)
          Vinfo(14)=Vname(4,idWao)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWao,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWao),    &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice freeze/melt wio
!
        IF (Hout(idWio,ng)) THEN
          Vinfo( 1)=Vname(1,idWio)
          Vinfo( 2)=Vname(2,idWio)
          Vinfo( 3)=Vname(3,idWio)
          Vinfo(14)=Vname(4,idWio)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWio,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWio),    &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice runoff wro
!
        IF (Hout(idWro,ng)) THEN
          Vinfo( 1)=Vname(1,idWro)
          Vinfo( 2)=Vname(2,idWro)
          Vinfo( 3)=Vname(3,idWro)
          Vinfo(14)=Vname(4,idWro)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWro,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWro),    &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice divergence wdiv
!
        IF (Hout(idWdiv,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiv)
          Vinfo( 2)=Vname(2,idWdiv)
          Vinfo( 3)=Vname(3,idWdiv)
          Vinfo(14)=Vname(4,idWdiv)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiv,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdiv),   &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef WEC_MELLOR
!
!  Define 2D radiation stress, Sxx-component.
!
        IF (Hout(idW2xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xx)
          Vinfo( 2)=Vname(2,idW2xx)
          Vinfo( 3)=Vname(3,idW2xx)
          Vinfo(14)=Vname(4,idW2xx)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW2xx),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress, Sxy-component.
!
        IF (Hout(idW2xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xy)
          Vinfo( 2)=Vname(2,idW2xy)
          Vinfo( 3)=Vname(3,idW2xy)
          Vinfo(14)=Vname(4,idW2xy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW2xy),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress, Syy-component.
!
        IF (Hout(idW2yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2yy)
          Vinfo( 2)=Vname(2,idW2yy)
          Vinfo( 3)=Vname(3,idW2yy)
          Vinfo(14)=Vname(4,idW2yy)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW2yy),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!

# ifdef SOLVE3D
!
!  Define 3D radiation stress, Sxx-component.
!
        IF (Hout(idW3xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xx)
          Vinfo( 2)=Vname(2,idW3xx)
          Vinfo( 3)=Vname(3,idW3xx)
          Vinfo(14)=Vname(4,idW3xx)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3xx),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Sxy-component.
!
        IF (Hout(idW3xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xy)
          Vinfo( 2)=Vname(2,idW3xy)
          Vinfo( 3)=Vname(3,idW3xy)
          Vinfo(14)=Vname(4,idW3xy)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3xy),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Syy-component.
!
        IF (Hout(idW3yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3yy)
          Vinfo( 2)=Vname(2,idW3yy)
          Vinfo( 3)=Vname(3,idW3yy)
          Vinfo(14)=Vname(4,idW3yy)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3yy),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Szx-component.
!
        IF (Hout(idW3zx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zx)
          Vinfo( 2)=Vname(2,idW3zx)
          Vinfo( 3)=Vname(3,idW3zx)
          Vinfo(14)=Vname(4,idW3zx)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3zx),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress, Szy-component.
!
        IF (Hout(idW3zy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zy)
          Vinfo( 2)=Vname(2,idW3zy)
          Vinfo( 3)=Vname(3,idW3zy)
          Vinfo(14)=Vname(4,idW3zy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3zy),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
#endif
#ifdef WEC
!
!  Define 2D total WEC u-stress.
!
        IF (Hout(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          Vinfo( 2)=Vname(2,idU2rs)
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idU2rs),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D total WEC v-stress.
!
        IF (Hout(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          Vinfo( 2)=Vname(2,idV2rs)
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idV2rs),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D u-Stokes drift velocity.
!
        IF (Hout(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          Vinfo( 2)=Vname(2,idU2Sd)
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idU2Sd),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D v-Stokes drift velocity.
!
        IF (Hout(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          Vinfo( 2)=Vname(2,idV2Sd)
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idV2Sd),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# ifdef SOLVE3D
!
!  Define 3D total WEC u-stress.
!
        IF (Hout(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          Vinfo( 2)=Vname(2,idU3rs)
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idU3rs),   &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D total WEC v-stress.
!
        IF (Hout(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          Vinfo( 2)=Vname(2,idV3rs)
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idV3rs),   &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D u-Stokes velocity.
!
        IF (Hout(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          Vinfo( 2)=Vname(2,idU3Sd)
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idU3Sd),   &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D v-Stokes velocity.
!
        IF (Hout(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          Vinfo( 2)=Vname(2,idV3Sd)
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idV3Sd),   &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D omega-Stokes velocity.
!
        IF (Hout(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          Vinfo( 2)=Vname(2,idW3Sd)
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3Sd,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3Sd),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D w-Stokes velocity test
!
        IF (Hout(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          Vinfo( 2)=Vname(2,idW3St)
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3St,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idW3St),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# endif
#endif
#ifdef WAVES_HEIGHT
!
!  Define wind-induced significat wave height.
!
        IF (Hout(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          Vinfo( 2)=Vname(2,idWamp)
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWamp,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWamp),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef WAVES_LENGTH
!
!  Define wind-induced mean wavelength.
!
        IF (Hout(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          Vinfo( 2)=Vname(2,idWlen)
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlen,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWlen),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef WAVES_LENGTHP
!
!  Define wind-induced peak wavelength.
!
        IF (Hout(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          Vinfo( 2)=Vname(2,idWlep)
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlep,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWlep),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef WAVES_DIR
!
!  Define wind-induced wave direction.
!
        IF (Hout(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          Vinfo( 2)=Vname(2,idWdir)
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdir,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdir),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef WAVES_TOP_PERIOD
!
!  Define wind-induced surface wave period.
!
        IF (Hout(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          Vinfo( 2)=Vname(2,idWptp)
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWptp,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWptp),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#endif
#ifdef WAVES_BOT_PERIOD
!
!  Define wind-induced bottom wave period.
!
        IF (Hout(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          Vinfo( 2)=Vname(2,idWpbt)
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWpbt,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWpbt),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#endif
#if defined BBL_MODEL || defined WAVES_OCEAN
!
!  Define wind-induced bottom orbital velocity.
!
        IF (Hout(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          Vinfo( 2)=Vname(2,idWorb)
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWorb,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWorb),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
!
!  Define wave dissipation due to bottom friction.
!
        IF (Hout(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          Vinfo( 2)=Vname(2,idWdif)
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdif,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdif),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WDISS_INWAVE
!
!  Define wave dissipation due to breaking.
!
        IF (Hout(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdib)
          Vinfo( 2)=Vname(2,idWdib)
          Vinfo( 3)=Vname(3,idWdib)
          Vinfo(14)=Vname(4,idWdib)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdib,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdib),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define wave dissipation due to white capping.
!
        IF (Hout(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          Vinfo( 2)=Vname(2,idWdiw)
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdiw),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined ROLLER_SVENDSEN
!
!  Define percent wave breaking.
!
        IF (Hout(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          Vinfo( 2)=Vname(2,idWbrk)
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbrk,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWbrk),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined WAVES_DSPR
!
!  Define waves directional spreading.
!
        IF (Hout(idWvds,ng)) THEN
          Vinfo( 1)=Vname(1,idWvds)
          Vinfo( 2)=Vname(2,idWvds)
          Vinfo( 3)=Vname(3,idWvds)
          Vinfo(14)=Vname(4,idWvds)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvds,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWvds),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define waves directional spreading.
!
        IF (Hout(idWvqp,ng)) THEN
          Vinfo( 1)=Vname(1,idWvqp)
          Vinfo( 2)=Vname(2,idWvqp)
          Vinfo( 3)=Vname(3,idWvqp)
          Vinfo(14)=Vname(4,idWvqp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvqp,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWvqp),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined WEC_ROLLER
!
!  Define Roller dissipation.
!
        IF (Hout(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          Vinfo( 2)=Vname(2,idWdis)
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdis,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWdis),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef ROLLER_RENIERS
!
!  Define roller action density.
!
        IF (Hout(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          Vinfo( 2)=Vname(2,idWrol)
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWrol,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWrol),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined WEC_VF
!
!  Define WEC quasi-static sea level adjustment.
!
        IF (Hout(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          Vinfo( 2)=Vname(2,idWztw)
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWztw,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWztw),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define WEC quasi-static pressure.
!
        IF (Hout(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          Vinfo( 2)=Vname(2,idWqsp)
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWqsp,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWqsp),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define WEC Bernoulli head.
!
        IF (Hout(idWbeh,ng)) THEN
          Vinfo( 1)=Vname(1,idWbeh)
          Vinfo( 2)=Vname(2,idWbeh)
          Vinfo( 3)=Vname(3,idWbeh)
          Vinfo(14)=Vname(4,idWbeh)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbeh,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWbeh),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, HIS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, HIS(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing history file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=HIS(ng)%name
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
!  Open history file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, HIS(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
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
!  history variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            HIS(ng)%Vid(idtime)=var_id(i)
#if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            HIS(ng)%Vid(idbath)=var_id(i)
#endif
#if defined WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            HIS(ng)%Vid(idPwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            HIS(ng)%Vid(idRwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            HIS(ng)%Vid(idUwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            HIS(ng)%Vid(idVwet)=var_id(i)
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            HIS(ng)%Vid(idpthR)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            HIS(ng)%Vid(idpthU)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            HIS(ng)%Vid(idpthV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            HIS(ng)%Vid(idpthW)=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            HIS(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            HIS(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            HIS(ng)%Vid(idVbar)=var_id(i)
#ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            HIS(ng)%Vid(idSbry(isFsur))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            HIS(ng)%Vid(idSbry(isUbar))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            HIS(ng)%Vid(idSbry(isVbar))=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            HIS(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            HIS(ng)%Vid(idv2dN)=var_id(i)
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            HIS(ng)%Vid(idRzet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            HIS(ng)%Vid(idRu2d)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            HIS(ng)%Vid(idRv2d)=var_id(i)
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRuct))) THEN
            got_var(idRuct)=.TRUE.
            HIS(ng)%Vid(idRuct)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRvct))) THEN
            got_var(idRvct)=.TRUE.
            HIS(ng)%Vid(idRvct)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            HIS(ng)%Vid(idUfx1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            HIS(ng)%Vid(idUfx2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            HIS(ng)%Vid(idVfx1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            HIS(ng)%Vid(idVfx2)=var_id(i)
#  ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            HIS(ng)%Vid(idRu3d)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            HIS(ng)%Vid(idRv3d)=var_id(i)
#  endif
# endif
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            HIS(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            HIS(ng)%Vid(idVvel)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            HIS(ng)%Vid(idSbry(isUvel))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            HIS(ng)%Vid(idSbry(isVvel))=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            HIS(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            HIS(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            HIS(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            HIS(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            HIS(ng)%Vid(idDano)=var_id(i)
# ifdef NEMURO_SED1
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONsed))) THEN
            got_var(idPONsed)=.TRUE.
            HIS(ng)%Vid(idPONsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALsed))) THEN
            got_var(idOPALsed)=.TRUE.
            HIS(ng)%Vid(idOPALsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDENITsed))) THEN
            got_var(idDENITsed)=.TRUE.
            HIS(ng)%Vid(idDENITsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONbur))) THEN
            got_var(idPONbur)=.TRUE.
            HIS(ng)%Vid(idPONbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALbur))) THEN
            got_var(idOPALbur)=.TRUE.
            HIS(ng)%Vid(idOPALbur)=var_id(i)
# endif
# ifdef PRIMARY_PROD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNPP))) THEN
            got_var(idNPP)=.TRUE.
            HIS(ng)%Vid(idNPP)=var_id(i)
# endif
# ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            HIS(ng)%Vid(idHsbl)=var_id(i)
# endif
# ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            HIS(ng)%Vid(idHbbl)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSSSf))) THEN
            got_var(idSSSf)=.TRUE.
            HIS(ng)%Vid(idSSSf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            HIS(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            HIS(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            HIS(ng)%Vid(idSdif)=var_id(i)
# if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            HIS(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            HIS(ng)%Vid(idMtls)=var_id(i)
# endif
# if defined RUNOFF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRunoff))) THEN
            got_var(idRunoff)=.TRUE.
            HIS(ng)%Vid(idRunoff)=var_id(i)
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
  || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            HIS(ng)%Vid(idPair)=var_id(i)
# endif
# if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            HIS(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            HIS(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUairE))) THEN
            got_var(idUairE)=.TRUE.
            HIS(ng)%Vid(idUairE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVairN))) THEN
            got_var(idVairN)=.TRUE.
            HIS(ng)%Vid(idVairN)=var_id(i)
# endif
# ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            HIS(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            HIS(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            HIS(ng)%Vid(idLrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            HIS(ng)%Vid(idTair)=var_id(i)
#  ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            HIS(ng)%Vid(idEmPf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            HIS(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            HIS(ng)%Vid(idrain)=var_id(i)
#  endif
# endif
# ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            HIS(ng)%Vid(idSrad)=var_id(i)
# endif
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            HIS(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            HIS(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            HIS(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            HIS(ng)%Vid(idVbms)=var_id(i)
#ifdef SOLVE3D
# ifdef UV_KIRBY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwav))) THEN
            got_var(idUwav)=.TRUE.
            HIS(ng)%Vid(idUwav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwav))) THEN
            got_var(idVwav)=.TRUE.
            HIS(ng)%Vid(idVwav)=var_id(i)
# endif
# ifdef BBL_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            HIS(ng)%Vid(idUbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            HIS(ng)%Vid(idVbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            HIS(ng)%Vid(idUbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            HIS(ng)%Vid(idVbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            HIS(ng)%Vid(idUbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            HIS(ng)%Vid(idVbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            HIS(ng)%Vid(idUVwc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            HIS(ng)%Vid(idUbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            HIS(ng)%Vid(idVbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            HIS(ng)%Vid(idUbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            HIS(ng)%Vid(idVbvr)=var_id(i)
# endif
#endif
#ifdef ICE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.true.
            HIS(ng)%Vid(idUice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.true.
            HIS(ng)%Vid(idVice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiceE))) THEN
            got_var(idUiceE)=.true.
            HIS(ng)%Vid(idUiceE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViceN))) THEN
            got_var(idViceN)=.true.
            HIS(ng)%Vid(idViceN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.true.
            HIS(ng)%Vid(idAice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.true.
            HIS(ng)%Vid(idHice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.true.
            HIS(ng)%Vid(idTice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.true.
            HIS(ng)%Vid(idHsno)=var_id(i)
# if defined ICE_BIO
         ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIcePhL))) THEN
            got_var(idIcePhL)=.true.
            HIS(ng)%Vid(idIcePhL)=var_id(i)
         ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIceNO3))) THEN
            got_var(idIceNO3)=.true.
            HIS(ng)%Vid(idIceNO3)=var_id(i)
         ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIceNH4))) THEN
            got_var(idIceNH4)=.true.
            HIS(ng)%Vid(idIceNH4)=var_id(i)
         ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIceLog))) THEN
            got_var(idIceLog)=.true.
            HIS(ng)%Vid(idIceLog)=var_id(i)
# endif
# ifdef MELT_PONDS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idApond))) THEN
            got_var(idApond)=.true.
            HIS(ng)%Vid(idApond)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHpond))) THEN
            got_var(idHpond)=.true.
            HIS(ng)%Vid(idHpond)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.true.
            HIS(ng)%Vid(idAgeice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
            got_var(idIomflx)=.true.
            HIS(ng)%Vid(idIomflx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.true.
            HIS(ng)%Vid(idTimid)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.true.
            HIS(ng)%Vid(idSig11)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.true.
            HIS(ng)%Vid(idSig12)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.true.
            HIS(ng)%Vid(idSig22)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.true.
            HIS(ng)%Vid(idTauiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.true.
            HIS(ng)%Vid(idChuiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.true.
            HIS(ng)%Vid(idT0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.true.
            HIS(ng)%Vid(idS0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
            got_var(idWfr)=.true.
            HIS(ng)%Vid(idWfr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
            got_var(idWai)=.true.
            HIS(ng)%Vid(idWai)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
            got_var(idWao)=.true.
            HIS(ng)%Vid(idWao)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
            got_var(idWio)=.true.
            HIS(ng)%Vid(idWio)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
            got_var(idWro)=.true.
            HIS(ng)%Vid(idWro)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiv))) THEN
            got_var(idWdiv)=.true.
            HIS(ng)%Vid(idWdiv)=var_id(i)
#endif
#if defined WEC_MELLOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xx))) THEN
            got_var(idW2xx)=.TRUE.
            HIS(ng)%Vid(idW2xx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xy))) THEN
            got_var(idW2xy)=.TRUE.
            HIS(ng)%Vid(idW2xy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2yy))) THEN
            got_var(idW2yy)=.TRUE.
            HIS(ng)%Vid(idW2yy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            HIS(ng)%Vid(idU2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            HIS(ng)%Vid(idV2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            HIS(ng)%Vid(idU2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            HIS(ng)%Vid(idV2Sd)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xx))) THEN
            got_var(idW3xx)=.TRUE.
            HIS(ng)%Vid(idW3xx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xy))) THEN
            got_var(idW3xy)=.TRUE.
            HIS(ng)%Vid(idW3xy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3yy))) THEN
            got_var(idW3yy)=.TRUE.
            HIS(ng)%Vid(idW3yy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zx))) THEN
            got_var(idW3zx)=.TRUE.
            HIS(ng)%Vid(idW3zx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zy))) THEN
            got_var(idW3zy)=.TRUE.
            HIS(ng)%Vid(idW3zy)=var_id(i)
# endif
#endif
#ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            HIS(ng)%Vid(idU2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            HIS(ng)%Vid(idV2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            HIS(ng)%Vid(idU2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            HIS(ng)%Vid(idV2Sd)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            HIS(ng)%Vid(idU3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            HIS(ng)%Vid(idV3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            HIS(ng)%Vid(idU3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            HIS(ng)%Vid(idV3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            HIS(ng)%Vid(idW3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            HIS(ng)%Vid(idW3St)=var_id(i)
# endif
#endif
#ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            HIS(ng)%Vid(idWamp)=var_id(i)
#endif
#ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            HIS(ng)%Vid(idWlen)=var_id(i)
#endif
#ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            HIS(ng)%Vid(idWlep)=var_id(i)
#endif
#ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            HIS(ng)%Vid(idWdir)=var_id(i)
#endif
#ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            HIS(ng)%Vid(idWptp)=var_id(i)
#endif
#ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            HIS(ng)%Vid(idWpbt)=var_id(i)
#endif
#if defined BBL_MODEL || defined WAVES_OCEAN
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            HIS(ng)%Vid(idWorb)=var_id(i)
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            HIS(ng)%Vid(idWdif)=var_id(i)
#endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            HIS(ng)%Vid(idWdib)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            HIS(ng)%Vid(idWdiw)=var_id(i)
#endif
#ifdef ROLLER_SVENDSEN
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            HIS(ng)%Vid(idWbrk)=var_id(i)
#endif
#ifdef WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            HIS(ng)%Vid(idWdis)=var_id(i)
#endif
#ifdef ROLLER_RENIERS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            HIS(ng)%Vid(idWrol)=var_id(i)
#endif
#ifdef WAVES_DSPR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvds))) THEN
            got_var(idWvds)=.TRUE.
            HIS(ng)%Vid(idWvds)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvqp))) THEN
            got_var(idWvqp)=.TRUE.
            HIS(ng)%Vid(idWvqp)=var_id(i)
#endif
#ifdef WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            HIS(ng)%Vid(idWztw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            HIS(ng)%Vid(idWqsp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            HIS(ng)%Vid(idWbeh)=var_id(i)
#endif

#ifdef INWAVE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACen))) THEN
            got_var(idACen)=.TRUE.
            HIS(ng)%Vid(idACen)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACcx))) THEN
            got_var(idACcx)=.TRUE.
            HIS(ng)%Vid(idACcx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACcy))) THEN
            got_var(idACcy)=.TRUE.
            HIS(ng)%Vid(idACcy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACct))) THEN
            got_var(idACct)=.TRUE.
            HIS(ng)%Vid(idACct)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACtp))) THEN
            got_var(idACtp)=.TRUE.
            HIS(ng)%Vid(idACtp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            HIS(ng)%Vid(idWdib)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            HIS(ng)%Vid(idWdiw)=var_id(i)
#endif

          END IF
!
#if defined VEGETATION 
# if defined VEG_DRAG || defined VEG_BIOMASS
          DO itrc=1,NVEGP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idvprp(itrc)))) THEN
              got_var(idvprp(itrc))=.TRUE.
              HIS(ng)%Vid(itrc)=var_id(i)
            END IF
          END DO
# endif 
#endif
! 
#ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              HIS(ng)%Tid(itrc)=var_id(i)
# ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              HIS(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
# endif
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              HIS(ng)%Vid(idTsur(itrc))=var_id(i)
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idGhat(itrc)))) THEN
              got_var(idGhat(itrc))=.TRUE.
              HIS(ng)%Vid(idGhat(itrc))=var_id(i)
# endif
            END IF
          END DO
# ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.                                   &
     &               TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              HIS(ng)%Vid(idfrac(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              HIS(ng)%Vid(idBmas(itrc))=var_id(i)
#  ifdef BEDLOAD
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              HIS(ng)%Vid(idUbld(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              HIS(ng)%Vid(idVbld(itrc))=var_id(i)
#  endif
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              HIS(ng)%Vid(idSbed(itrc))=var_id(i)
            END IF
          END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              HIS(ng)%Vid(idBott(itrc))=var_id(i)
            END IF
          END DO
# endif
#endif
        END DO
!
!  Check if history variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath).and.Hout(idbath,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idbath)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined WET_DRY
        IF (.not.got_var(idPwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        IF (.not.got_var(idpthR).and.Hout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Hout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Hout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Hout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef ADJUST_BOUNDARY
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
#endif
        IF (.not.got_var(idu2dE).and.Hout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Hout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
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
# endif
        IF (.not.got_var(idu3dE).and.Hout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Hout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef NEMURO_SED1
        IF (.not.got_var(idPONsed).and.Hout(idPONsed,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPONsed)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALsed).and.Hout(idOPALsed,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOPALsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDENITsed).and.Hout(idDENITsed,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDENITsed)),      &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPONbur).and.Hout(idPONbur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPONbur)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALbur).and.Hout(idOPALbur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOPALbur)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef PRIMARY_PROD
        IF (.not.got_var(idNPP).and.Hout(idNPP,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idNPP)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idSSSf).and.Hout(idSSSf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSSSf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SALINITY
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined RUNOFF
        IF (.not.got_var(idRunoff).and.Hout(idRunoff,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRunoff)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
  || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
        IF (.not.got_var(idPair).and.Hout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM2OCN_FLUXES
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUairE).and.Hout(idUairE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUairE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVairN).and.Hout(idVairN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVairN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef EMINUSP
        IF (.not.got_var(idEmPf).and.Hout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
# ifdef UV_KIRBY
        IF (.not.got_var(idUwav).and.Hout(idUwav,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwav).and.Hout(idVwav,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.Hout(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.Hout(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.Hout(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.Hout(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.Hout(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVwc).and.Hout(idUVwc,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUVwc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.Hout(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.Hout(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.Hout(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.Hout(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.Hout(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef ICE_MODEL
        IF (.not.got_var(idUice).and.Hout(idUice,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Hout(idVice,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUiceE).and.Hout(idUiceE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUiceE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViceN).and.Hout(idViceN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idViceN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Hout(idAice,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Hout(idHice,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Hout(idTice,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Hout(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined ICE_BIO
        IF (.not.got_var(idIcePhL).and.Hout(idIcePhL,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idIcePhL)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF            
# endif
# ifdef MELT_PONDS
        IF (.not.got_var(idApond).and.Hout(idApond,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idApond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHpond).and.Hout(idHpond,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHpond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idAgeice).and.Hout(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Hout(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Hout(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Hout(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Hout(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Hout(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Hout(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Hout(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Hout(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Hout(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Hout(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Hout(idWai,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Hout(idWao,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Hout(idWio,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Hout(idWro,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiv).and.Hout(idWdiv,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWdiv)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined WEC_MELLOR
        IF (.not.got_var(idW2xx).and.Hout(idW2xx,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW2xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2xy).and.Hout(idW2xy,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW2xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2yy).and.Hout(idW2yy,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW2yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idW3xx).and.Hout(idW3xx,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3xy).and.Hout(idW3xy,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3yy).and.Hout(idW3yy,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zx).and.Hout(idW3zx,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3zx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zy).and.Hout(idW3zy,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3zy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef WEC
        IF (.not.got_var(idU2Sd).and.Hout(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Hout(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Hout(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Hout(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idU3Sd).and.Hout(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Hout(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Hout(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Hout(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.Hout(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.Hout(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.Hout(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.Hout(idWlen,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWlen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.Hout(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.Hout(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.Hout(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.Hout(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined BBL_MODEL || defined WAVES_OCEAN
        IF (.not.got_var(idWorb).and.Hout(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (.not.got_var(idWdif).and.Hout(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (.not.got_var(idWdib).and.Hout(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.Hout(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef ROLLER_SVENDSEN
        IF (.not.got_var(idWbrk).and.Hout(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WEC_ROLLER
        IF (.not.got_var(idWdis).and.Hout(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef ROLLER_RENIERS
        IF (.not.got_var(idWrol).and.Hout(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WAVES_DSPR
        IF (.not.got_var(idWvds).and.Hout(idWvds,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvds)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvqp).and.Hout(idWvqp,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvqp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef WEC_VF
        IF (.not.got_var(idWztw).and.Hout(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.Hout(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.Hout(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWbeh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
# ifdef ADJUST_BOUNDARY
          ELSE IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,70)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
# endif
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
          IF (.not.got_var(idGhat(itrc)).and.Hout(idGhat(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idGhat(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
# endif
        END DO
# ifdef VEGETATION 
#  if defined VEG_DRAG || defined VEG_BIOMASS 
        DO i=1,NVEGP 
          IF (.not.got_var(idvprp(i)).and.Hout(idvprp(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idvprp(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO 
#  endif
# endif 
# ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i)).and.Hout(idfrac(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          IF(.not.got_var(idBmas(i)).and.Hout(idBmas(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef BEDLOAD
          IF (.not.got_var(idUbld(i)).and.Hout(idUbld(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.Hout(idVbld(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i)).and.Hout(idSbed(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO i=1,MBOTP
          IF (.not.got_var(idBott(i)).and.Hout(idBott(i),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefHIS(ng).gt.0) THEN
          HIS(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefHIS(ng)*((ntstart(ng)-1)/ndefHIS(ng)))/   &
     &                   nHIS(ng)
        ELSE
          HIS(ng)%Rindex=(ntstart(ng)-1)/nHIS(ng)
        END IF
        HIS(ng)%Rindex=MIN(HIS(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (6x,'DEF_HIS   - creating history file, Grid ',i2.2,': ',  &
     &        a)
  20  FORMAT (6x,'DEF_HIS   - inquiring history file, Grid ',i2.2,': ', &
     &        a)
  30  FORMAT (/,' DEF_HIS - unable to create history NetCDF file: ',a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_HIS - unable to open history NetCDF file: ',a)
  70  FORMAT (/,' DEF_HIS - unable to find variable: ',a,2x,            &
     &        ' in history NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_his
