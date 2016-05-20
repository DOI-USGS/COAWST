#include "cppdefs.h"
      SUBROUTINE def_rst (ng)
!
!svn $Id: def_rst.F 788 2016-05-05 16:14:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates restart NetCDF file, it defines its            !
!  dimensions, attributes, and variables.                              !
!                                                                      !
#if defined PERFECT_RESTART && defined WRITE_WATER && defined MASKING
!  Currently, perfect restart is not compatible with writing           !
!  only water points.                                                  !
!                                                                      !
#endif
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

      integer :: i, j, nvd3, nvd4, nvd5
      integer :: recdim, status, varid

#if defined VEGETATION && (defined VEG_DRAG || defined VEG_BIOMASS)
      integer :: DimIDs(35)
#else
      integer :: DimIDs(32)
#endif
      integer :: r2dgrd(4), ru2dgrd(4), rv2dgrd(4)
      integer :: sp2dgrd(3), sr2dgrd(3), su2dgrd(3), sv2dgrd(3)
      integer :: sr3dgrd(4), su3dgrd(4), sv3dgrd(4)
      integer :: t2dgrd(4), u2dgrd(4), v2dgrd(4)
      integer :: Vsize(4)
#if defined VEGETATION 
# if defined VEG_DRAG || defined VEG_BIOMASS
      integer :: v3pgrd(4)
# endif
#endif
      integer :: def_dim

#ifdef SOLVE3D
      integer :: itrc

      integer :: k3dgrd(5), t3dgrd(5)
      integer :: r3dgrd(4), ru3dgrd(5), rv3dgrd(5)
      integer :: u3dgrd(5), v3dgrd(5), w3dgrd(4)
#endif
      real(r8) :: Aval(6)

      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_rst.F'
!
!=======================================================================
!  Create a new restart NetCDF file.
!=======================================================================
!
!  Activate creation of restart NetCDF file.  Create a new restart
!  file if during a restart run, the restart filename "RST(ng)%name"
!  is different than the initial filename "INI(ng)%name".
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=RST(ng)%name
      Ldefine=.FALSE.
      IF (((nrrec(ng).eq.0).and.(iic(ng).eq.ntstart(ng))).or.           &
     &    ((nrrec(ng).ne.0).and.                                        &
     &     (TRIM(ncname).ne.TRIM(INI(ng)%name)))) THEN
        Ldefine=.TRUE.
      END IF
!
      DEFINE : IF (Ldefine) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), RST(ng)%ncid)
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
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef SOLVE3D
# if !defined PERFECT_RESTART && \
     (defined WRITE_WATER && defined MASKING)
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'N',             &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

# ifdef SEDIMENT
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#  if !defined PERFECT_RESTART && \
      (defined WRITE_WATER && defined MASKING)
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
# if defined VEGETATION
#  if defined VEG_DRAG || defined VEG_BIOMASS
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'NVEG',          &
     &                 NVEG, DimIDs(35))
        IF (exit_flag.ne.NoError) RETURN
#  endif 
# endif
!
# if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef ECOSIM
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
# endif
#endif

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

#ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef PERFECT_RESTART
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'two',           &
     &                 2, DimIDs(30))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'three',         &
     &                 3, DimIDs(31))
        IF (exit_flag.ne.NoError) RETURN
#endif

        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
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
#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
        sr2dgrd(1)=DimIDs(17)
        sr2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
        r3dgrd(1)=DimIDs(20)
        r3dgrd(2)=DimIDs(12)
# endif
#else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        sr2dgrd(1)=DimIDs( 1)
        sr2dgrd(2)=DimIDs( 5)
        sr2dgrd(3)=DimIDs(12)
# ifdef PERFECT_RESTART
        t2dgrd(3)=DimIDs(31)
        t2dgrd(4)=DimIDs(12)
# else
        t2dgrd(3)=DimIDs(12)
# endif
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        r3dgrd(1)=DimIDs( 1)
        r3dgrd(2)=DimIDs( 5)
        r3dgrd(3)=DimIDs( 9)
#  ifdef PERFECT_RESTART
        t3dgrd(4)=DimIDs(30)
        t3dgrd(5)=DimIDs(12)
#  else
        t3dgrd(4)=DimIDs(12)
#  endif
        r3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        sp2dgrd(1)=DimIDs( 4)
        sp2dgrd(2)=DimIDs( 8)
        sp2dgrd(3)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momentum type variables.
!
#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
# endif
#else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
# ifdef PERFECT_RESTART
        u2dgrd(3)=DimIDs(31)
        u2dgrd(4)=DimIDs(12)
# else
        u2dgrd(3)=DimIDs(12)
# endif
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
#  ifdef PERFECT_RESTART
        u3dgrd(4)=DimIDs(30)
        u3dgrd(5)=DimIDs(12)
#  else
        u3dgrd(4)=DimIDs(12)
#  endif
# endif
#endif
!
!  Define dimensions for the plant array variables
!
#if !defined PERFECT_RESTART && defined VEGETATION     
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
#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
# endif
#else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
# ifdef PERFECT_RESTART
        v2dgrd(3)=DimIDs(31)
        v2dgrd(4)=DimIDs(12)
# else
        v2dgrd(3)=DimIDs(12)
# endif
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
#  ifdef PERFECT_RESTART
        v3dgrd(4)=DimIDs(30)
        v3dgrd(5)=DimIDs(12)
#  else
        v3dgrd(4)=DimIDs(12)
#  endif
# endif
#endif
#ifdef PERFECT_RESTART
!
!  Define dimension vectors for RHS free-surface equation.
!
# if !defined PERFECT_RESTART && \
     (defined WRITE_WATER && defined MASKING)
        r2dgrd(1)=DimIDs(17)
        r2dgrd(2)=DimIDs(12)
# else
        r2dgrd(1)=DimIDs( 1)
        r2dgrd(2)=DimIDs( 5)
        r2dgrd(3)=DimIDs(30)
        r2dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for RHS u-momentum equation.
!
# if !defined PERFECT_RESTART && \
     (defined WRITE_WATER && defined MASKING)
        ru2dgrd(1)=DimIDs(18)
        ru2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        ru2dgrd(1)=DimIDs(21)
        ru2dgrd(2)=DimIDs(12)
#  endif
# else
        ru2dgrd(1)=DimIDs( 2)
        ru2dgrd(2)=DimIDs( 6)
        ru2dgrd(3)=DimIDs(30)
        ru2dgrd(4)=DimIDs(12)
#  ifdef SOLVE3D
        ru3dgrd(1)=DimIDs( 2)
        ru3dgrd(2)=DimIDs( 6)
        ru3dgrd(3)=DimIDs(10)
        ru3dgrd(4)=DimIDs(30)
        ru3dgrd(5)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for RHS v-momentum equation.
!
# if !defined PERFECT_RESTART && \
     (defined WRITE_WATER && defined MASKING)
        rv2dgrd(1)=DimIDs(19)
        rv2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        rv3dgrd(1)=DimIDs(22)
        rv3dgrd(2)=DimIDs(12)
#  endif
# else
        rv2dgrd(1)=DimIDs( 3)
        rv2dgrd(2)=DimIDs( 7)
        rv2dgrd(3)=DimIDs(30)
        rv2dgrd(4)=DimIDs(12)
#  ifdef SOLVE3D
        rv3dgrd(1)=DimIDs( 3)
        rv3dgrd(2)=DimIDs( 7)
        rv3dgrd(3)=DimIDs(10)
        rv3dgrd(4)=DimIDs(30)
        rv3dgrd(5)=DimIDs(12)
#  endif
# endif
#endif
#ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
# if !defined PERFECT_RESTART && \
     (defined WRITE_WATER && defined MASKING)
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#  ifdef PERFECT_RESTART
        k3dgrd(1)=DimIDs(23)
        k3dgrd(2)=DimIDs(12)
#  endif
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#  ifdef PERFECT_RESTART
        k3dgrd(1)=DimIDs( 1)
        k3dgrd(2)=DimIDs( 5)
        k3dgrd(3)=DimIDs(10)
        k3dgrd(4)=DimIDs(30)
        k3dgrd(5)=DimIDs(12)
#  endif
# endif
#endif
!
!  Define dimension vector for sediment, radiation stress variables.
!
#if !defined PERFECT_RESTART && \
    (defined WRITE_WATER && defined MASKING)
        su2dgrd(1)=DimIDs(18)
        su2dgrd(2)=DimIDs(12)
        sv2dgrd(1)=DimIDs(19)
        sv2dgrd(2)=DimIDs(12)
#else
        su2dgrd(1)=DimIDs( 2)
        su2dgrd(2)=DimIDs( 6)
        su2dgrd(3)=DimIDs(12)
        sv2dgrd(1)=DimIDs( 3)
        sv2dgrd(2)=DimIDs( 7)
        sv2dgrd(3)=DimIDs(12)
#endif
#ifdef SOLVE3D
# if !defined PERFECT_RESTART && \
     (defined WRITE_WATER && defined MASKING)
        sr3dgrd(1)=DimIDs(24)
        sr3dgrd(2)=DimIDs(12)
        su3dgrd(1)=DimIDs(21)
        su3dgrd(2)=DimIDs(12)
        sv3dgrd(1)=DimIDs(22)
        sv3dgrd(2)=DimIDs(12)
# else
        sr3dgrd(1)=DimIDs( 1)
        sr3dgrd(2)=DimIDs( 5)
        sr3dgrd(3)=DimIDs(16)
        sr3dgrd(4)=DimIDs(12)
        su3dgrd(1)=DimIDs( 2)
        su3dgrd(2)=DimIDs( 6)
        su3dgrd(3)=DimIDs( 9)
        su3dgrd(4)=DimIDs(12)
        sv3dgrd(1)=DimIDs( 3)
        sv3dgrd(2)=DimIDs( 7)
        sv3dgrd(3)=DimIDs( 9)
        sv3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        RST(ng)%Rindex=0
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
        CALL def_info (ng, iNLM, RST(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------

#ifdef PERFECT_RESTART
!
!  Define time-stepping indices.
!
# ifdef SOLVE3D
        Vinfo( 1)='nstp'
        Vinfo( 2)='3D equations time level index, nstp'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='nrhs'
        Vinfo( 2)='3D equations time level index, nrhs'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='nnew'
        Vinfo( 2)='3D equations time level index, nnew'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
        Vinfo( 1)='kstp'
        Vinfo( 2)='3D equations time level index, kstp'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='krhs'
        Vinfo( 2)='3D equations time level index, krhs'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='knew'
        Vinfo( 2)='3D equations time level index, knew'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# ifdef ICE_MODEL
        Vinfo( 1)='linew'
        Vinfo( 2)='ice equations time level index, linew'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='liold'
        Vinfo( 2)='ice equations time level index, liold'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='liunw'
        Vinfo( 2)='ice equations time level index, liunw'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='liuol'
        Vinfo( 2)='ice equations time level index, liuol'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='lienw'
        Vinfo( 2)='ice equations time level index, lienw'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='lieol'
        Vinfo( 2)='ice equations time level index, lieol'
        status=def_var(ng, iNLM, RST(ng)%ncid, varid, nf90_int,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
#endif
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
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

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
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idPwet),     &
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
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRwet),     &
     &                 NF_FOUT, nvd3, sr2dgrd, Aval, Vinfo, ncname,     &
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
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUwet),     &
     &                 NF_FOUT, nvd3, su2dgrd, Aval, Vinfo, ncname,     &
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
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVwet),     &
     &                 NF_FOUT, nvd3, sv2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#endif
!
#if defined VEGETATION 
#include "vegetation_def_rst.h"
#endif
!
#if defined SEDIMENT && defined SED_MORPH
!
!  Define  time-varying bathymetry
!
        Vinfo( 1)=Vname(1,idbath)
        Vinfo( 2)=Vname(2,idbath)
        Vinfo( 3)=Vname(3,idbath)
        Vinfo(14)=Vname(4,idbath)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idbath),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
# if !defined WET_DRY && defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# else
        Vinfo(20)='mask_rho'
# endif
#endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
#ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),     &
# ifdef WET_DRY
     &                 NF_FRST, nvd4, t2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
# else
     &                 NF_FRST, nvd4, t2dgrd, Aval, Vinfo, ncname)
# endif
#else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),     &
# ifdef WET_DRY
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
# else
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
#endif
        IF (exit_flag.ne.NoError) RETURN

#ifdef PERFECT_RESTART
!
!  Define RHS of free-surface equation.
!
        Vinfo( 1)=Vname(1,idRzet)
        Vinfo( 2)=Vname(2,idRzet)
        Vinfo( 3)=Vname(3,idRzet)
        Vinfo(14)=Vname(4,idRzet)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_rho'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRzet),     &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)=Vname(3,idUbar)
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
# if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# else
        Vinfo(20)='mask_u'
# endif
#endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
#ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),     &
     &                 NF_FRST, nvd4, u2dgrd, Aval, Vinfo, ncname)
#else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),     &
     &                 NF_FRST, nvd3, u2dgrd, Aval, Vinfo, ncname)
#endif
        IF (exit_flag.ne.NoError) RETURN

#ifdef PERFECT_RESTART
!
!  Define RHS of 2D momentum equation in the XI-direction.
!
        Vinfo( 1)=Vname(1,idRu2d)
        Vinfo( 2)=Vname(2,idRu2d)
        Vinfo( 3)=Vname(3,idRu2d)
        Vinfo(14)=Vname(4,idRu2d)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_u'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRu2d,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRu2d),     &
     &                 NF_FRST, nvd4, ru2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)=Vname(3,idVbar)
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
# if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# else
        Vinfo(20)='mask_v'
# endif
#endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
#ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),     &
     &                 NF_FRST, nvd4, v2dgrd, Aval, Vinfo, ncname)
#else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),     &
     &                 NF_FRST, nvd3, v2dgrd, Aval, Vinfo, ncname)
#endif
        IF (exit_flag.ne.NoError) RETURN

#ifdef PERFECT_RESTART
!
!  Define RHS of 2D momentum equation in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idRv2d)
        Vinfo( 2)=Vname(2,idRv2d)
        Vinfo( 3)=Vname(3,idRv2d)
        Vinfo(14)=Vname(4,idRv2d)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_v'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRv2d,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRv2d),     &
     &                 NF_FRST, nvd4, rv2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUvel)
        Vinfo( 2)=Vname(2,idUvel)
        Vinfo( 3)=Vname(3,idUvel)
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_u'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUvel),     &
     &                 NF_FRST, nvd5, u3dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUvel),     &
     &                 NF_FRST, nvd4, u3dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN

# ifdef PERFECT_RESTART
!
!  Define RHS of 3D momentum equation in the XI-direction.
!  Although this variable is a U-points, a negative value is used
!  here to set "s_w" in the "coordinate" attribute.  The k=0 index
!  is used during coupling in step2d.
!
        Vinfo( 1)=Vname(1,idRu3d)
        Vinfo( 2)=Vname(2,idRu3d)
        Vinfo( 3)=Vname(3,idRu3d)
        Vinfo(14)=Vname(4,idRu3d)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_u'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(-u3dvar,r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRu3d),     &
     &                 NF_FRST, nvd5, ru3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define 3D momentum component in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVvel)
        Vinfo( 2)=Vname(2,idVvel)
        Vinfo( 3)=Vname(3,idVvel)
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_v'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvel),     &
     &                 NF_FRST, nvd5, v3dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvel),     &
     &                 NF_FRST, nvd4, v3dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN

# ifdef PERFECT_RESTART
!
!  Define RHS of 3D momentum equation in the ETA-direction.
!  Although this variable is a V-points, a negative value is used
!  here to set "s_w" in the "coordinate" attribute.  The k=0 index
!  is used during coupling in step2d.
!
        Vinfo( 1)=Vname(1,idRv3d)
        Vinfo( 2)=Vname(2,idRv3d)
        Vinfo( 3)=Vname(3,idRv3d)
        Vinfo(14)=Vname(4,idRv3d)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_v'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(-v3dvar,r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRv3d),     &
     &                 NF_FRST, nvd5, rv3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
          Vinfo(20)='mask_rho'
#  endif
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
# ifdef PERFECT_RESTART
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Tid(itrc),     &
     &                   NF_FRST, nvd5, t3dgrd, Aval, Vinfo, ncname)
# else
          status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Tid(itrc),     &
     &                   NF_FRST, nvd4, t3dgrd, Aval, Vinfo, ncname)
# endif
          IF (exit_flag.ne.NoError) RETURN
        END DO

# if defined PERFECT_RESTART && defined ICE_MODEL
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
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
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idTsur(itrc)), NF_FRST,            &
     &                   nvd3, sr2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Define surface U-momentum stress.
!
        Vinfo( 1)=Vname(1,idUsms)
        Vinfo( 2)=Vname(2,idUsms)
        Vinfo( 3)=Vname(3,idUsms)
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUsms),     &
     &                   NF_FRST, nvd3, su2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface V-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)=Vname(3,idVsms)
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVsms),     &
     &                   NF_FRST, nvd3, sv2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined BENTHIC
!------------------------------------------
!  Define benthic tracer type variables.
!------------------------------------------   
        DO itrc=1,NBeT(ng)
          rstBid(itrc,ng) = var_id(i)
          IF (Hout(idBeTvar(itrc),ng)) THEN

            Vinfo( 1)=Vname(1,idBeTvar(itrc))
            Vinfo( 2)=Vname(2,idBeTvar(itrc))
            Vinfo(3)=Vname(3,idBeTvar(itrc))
            Vinfo(14)=Vname(4,idBeTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng,iNLM,RST(ng)%ncid,RST(ng)%Bid(itrc),NF_FOUT, &
     &                     nvd3,t2dgrd,Aval,Vinfo,ncname)

!Will need to switch this if have more than one depth level for benthos 
!         status=def_var(ncRSTid(ng),hisBid(itrc,ng),NF_FOUT, &
!    &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif

!------------------------------------------
!  Define ice bio tracer type variables.
!------------------------------------------  
# if defined ICE_BIO
#  ifdef CLIM_ICE_1D
        DO itrc=1,NIceT(ng)
          IF (Hout(idIceBeTvar(itrc),ng)) THEN

            Vinfo( 1)=Vname(1,idIceBeTvar(itrc))
            Vinfo( 2)=Vname(2,idIceBeTvar(itrc))
            Vinfo(3)=Vname(3,idIceBeTvar(itrc))
            Vinfo(14)=Vname(4,idIceBeTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)  
            status=def_var(ng,iNLM,RST(ng)%ncid,RST(ng)%IceBid(itrc),      &
     &                   NF_FOUT, nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
# endif     
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_rho'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idDano),     &
     &                 NF_FRST, nvd4, r3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef NEMURO_SED1
!
!  Define PON in sediment.
!
        Vinfo( 1)=Vname(1,idPONsed)
        Vinfo( 2)=Vname(2,idPONsed)
        Vinfo( 3)=Vname(3,idPONsed)
        Vinfo(14)=Vname(4,idPONsed)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idPONsed,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idPONsed),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define OPAL in sediment.
!
        Vinfo( 1)=Vname(1,idOPALsed)
        Vinfo( 2)=Vname(2,idOPALsed)
        Vinfo( 3)=Vname(3,idOPALsed)
        Vinfo(14)=Vname(4,idOPALsed)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idOPALsed,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idOPALsed),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define DENIT in sediment.
!
        Vinfo( 1)=Vname(1,idDENITsed)
        Vinfo( 2)=Vname(2,idDENITsed)
        Vinfo( 3)=Vname(3,idDENITsed)
        Vinfo(14)=Vname(4,idDENITsed)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idDENITsed,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idDENITsed),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define buried PON in sediment.
!
        Vinfo( 1)=Vname(1,idPONbur)
        Vinfo( 2)=Vname(2,idPONbur)
        Vinfo( 3)=Vname(3,idPONbur)
        Vinfo(14)=Vname(4,idPONbur)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idPONbur,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idPONbur),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define buried OPAL in sediment.
!
        Vinfo( 1)=Vname(1,idOPALbur)
        Vinfo( 2)=Vname(2,idOPALbur)
        Vinfo( 3)=Vname(3,idOPALbur)
        Vinfo(14)=Vname(4,idOPALbur)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idOPALbur,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idOPALbur),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        Vinfo( 1)=Vname(1,idHsbl)
        Vinfo( 2)=Vname(2,idHsbl)
        Vinfo( 3)=Vname(3,idHsbl)
        Vinfo(14)=Vname(4,idHsbl)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHsbl),     &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        Vinfo( 1)=Vname(1,idHbbl)
        Vinfo( 2)=Vname(2,idHbbl)
        Vinfo( 3)=Vname(3,idHbbl)
        Vinfo(14)=Vname(4,idHbbl)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHbbl),     &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined PERFECT_RESTART && defined LMD_NONLOCAL
!
!  Define out KPP nonlocal transport.
!
        DO itrc=1,NAT
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
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                     RST(ng)%Vid(idGhat(itrc)), NF_FRST,          &
     &                     nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
# if defined BVF_MIXING || defined LMD_MIXING || \
     defined GLS_MIXING || defined MY25_MIXING
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        Vinfo( 2)=Vname(2,idVvis)
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvis),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        Vinfo( 2)=Vname(2,idTdif)
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTdif),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)

        IF (exit_flag.ne.NoError) RETURN

#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        Vinfo( 2)=Vname(2,idSdif)
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSdif),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
!
!  Define turbulent kinetic energy.
!
        Vinfo( 1)=Vname(1,idMtke)
        Vinfo( 2)=Vname(2,idMtke)
        Vinfo( 3)=Vname(3,idMtke)
        Vinfo(14)=Vname(4,idMtke)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idMtke),     &
     &                 NF_FRST, nvd5, k3dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)

        IF (exit_flag.ne.NoError) RETURN
!
!  Define turbulent kinetic energy time length scale.
!
        Vinfo( 1)=Vname(1,idMtls)
        Vinfo( 2)=Vname(2,idMtls)
        Vinfo( 3)=Vname(3,idMtls)
        Vinfo(14)=Vname(4,idMtls)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idMtls),     &
     &                 NF_FRST, nvd5, k3dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical mixing turbulent length scale.
!
        Vinfo( 1)=Vname(1,idVmLS)
        Vinfo( 2)=Vname(2,idVmLS)
        Vinfo( 3)=Vname(3,idVmLS)
        Vinfo(14)=Vname(4,idVmLS)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVmLS),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define turbulent kinetic energy vertical diffusion coefficient.
!
        Vinfo( 1)=Vname(1,idVmKK)
        Vinfo( 2)=Vname(2,idVmKK)
        Vinfo( 3)=Vname(3,idVmKK)
        Vinfo(14)=Vname(4,idVmKK)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVmKK),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

#  ifdef GLS_MIXING
!
!  Define turbulent length scale vertical diffusion coefficient.
!
        Vinfo( 1)=Vname(1,idVmKP)
        Vinfo( 2)=Vname(2,idVmKP)
        Vinfo( 3)=Vname(3,idVmKP)
        Vinfo(14)=Vname(4,idVmKP)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#    else
        Vinfo(20)='mask_rho'
#    endif
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVmKP,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVmKP),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Define Bedload U-direction.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idUbld(i))
          Vinfo( 2)=Vname(2,idUbld(i))
          Vinfo( 3)=Vname(3,idUbld(i))
          Vinfo(14)=Vname(4,idUbld(i))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#    else
          Vinfo(20)='mask_u'
#    endif
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idUbld(i)), NF_FRST,               &
     &                   nvd3, su2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
!  Define Bedload V-direction.
!
          Vinfo( 1)=Vname(1,idVbld(i))
          Vinfo( 2)=Vname(2,idVbld(i))
          Vinfo( 3)=Vname(3,idVbld(i))
          Vinfo(14)=Vname(4,idVbld(i))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#    else
          Vinfo(20)='mask_v'
#    endif
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idVbld(i)), NF_FRST,               &
     &                   nvd3, sv2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#  endif
!
!  Define sediment fraction of each size class in each bed layer.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idfrac(i))
          Vinfo( 2)=Vname(2,idfrac(i))
          Vinfo( 3)=Vname(3,idfrac(i))
          Vinfo(14)=Vname(4,idfrac(i))
          Vinfo(16)=Vname(1,idtime)
          WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
          Vinfo(20)='mask_rho'
#   endif
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idfrac(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idfrac(i)), NF_FRST,               &
     &                   nvd4, sr3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Define sediment mass of each size class in each bed layer.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idBmas(i))
          Vinfo( 2)=Vname(2,idBmas(i))
          Vinfo( 3)=Vname(3,idBmas(i))
          Vinfo(14)=Vname(4,idBmas(i))
          Vinfo(16)=Vname(1,idtime)
          WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
          Vinfo(20)='mask_rho'
#   endif
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idBmas(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idBmas(i)), NF_FRST,               &
     &                   nvd4, sr3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Define sediment properties in each bed layer.
!
        DO i=1,MBEDP
          Vinfo( 1)=Vname(1,idSbed(i))
          Vinfo( 2)=Vname(2,idSbed(i))
          Vinfo( 3)=Vname(3,idSbed(i))
          Vinfo(14)=Vname(4,idSbed(i))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#   else
          Vinfo(20)='mask_rho'
#   endif
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSbed(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idSbed(i)), NF_FRST,               &
     &                   nvd4, sr3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  define exposed sediment layer properties. Notice that only the
!  first four properties (mean grain diameter, mean grain density,
!  mean settling velocity, mean critical erosion stress,
!  ripple length and ripple height) are written.
!
        DO i=1,6
          Vinfo( 1)=Vname(1,idBott(i))
          Vinfo( 2)=Vname(2,idBott(i))
          Vinfo( 3)=Vname(3,idBott(i))
          Vinfo(14)=Vname(4,idBott(i))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#   else
          Vinfo(20)='mask_rho'
#   endif
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idBott(i),ng),r8)
          status=def_var(ng, iNLM, RST(ng)%ncid,                        &
     &                   RST(ng)%Vid(idBott(i)), NF_FRST,               &
     &                   nvd3, sr2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
#endif
#ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUice)
        Vinfo( 2)=Vname(2,idUice)
        Vinfo( 3)=Vname(3,idUice)
        Vinfo(14)=Vname(4,idUice)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUice,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUice),     &
     &                 NF_FRST, nvd4, ru2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUice),     &
     &                 NF_FRST, nvd3, u2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D ice momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVice)
        Vinfo( 2)=Vname(2,idVice)
        Vinfo( 3)=Vname(3,idVice)
        Vinfo(14)=Vname(4,idVice)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVice,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVice),     &
     &                 NF_FRST, nvd4, rv2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVice),     &
     &                 NF_FRST, nvd3, v2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice concentration.
!
        Vinfo( 1)=Vname(1,idAice)
        Vinfo( 2)=Vname(2,idAice)
        Vinfo( 3)=Vname(3,idAice)
        Vinfo(14)=Vname(4,idAice)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idAice,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idAice),     &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idAice),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice average thickness.
!
        Vinfo( 1)=Vname(1,idHice)
        Vinfo( 2)=Vname(2,idHice)
        Vinfo( 3)=Vname(3,idHice)
        Vinfo(14)=Vname(4,idHice)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHice,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHice),     &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHice),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice/snow surface temperature.
!
        Vinfo( 1)=Vname(1,idTice)
        Vinfo( 2)=Vname(2,idTice)
        Vinfo( 3)=Vname(3,idTice)
        Vinfo(14)=Vname(4,idTice)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTice,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTice),     &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define snow thickness.
!
        Vinfo( 1)=Vname(1,idHsno)
        Vinfo( 2)=Vname(2,idHsno)
        Vinfo( 3)=Vname(3,idHsno)
        Vinfo(14)=Vname(4,idHsno)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHsno,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHsno),     &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHsno),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
# ifdef MELT_PONDS
!
!  Define surface water fraction (on ice).
!
        Vinfo( 1)=Vname(1,idApond)
        Vinfo( 2)=Vname(2,idApond)
        Vinfo( 3)=Vname(3,idApond)
        Vinfo(14)=Vname(4,idApond)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idApond,ng),r8)
#  ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idApond),    &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
#  else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idApond),    &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
#  endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface water thickness (on ice).
!
        Vinfo( 1)=Vname(1,idHpond)
        Vinfo( 2)=Vname(2,idHpond)
        Vinfo( 3)=Vname(3,idHpond)
        Vinfo(14)=Vname(4,idHpond)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHpond,ng),r8)
#  ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHpond),    &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
#  else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHpond),    &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
#  endif
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define ice age.
!
        Vinfo( 1)=Vname(1,idAgeice)
        Vinfo( 2)=Vname(2,idAgeice)
        Vinfo( 3)=Vname(3,idAgeice)
        Vinfo(14)=Vname(4,idAgeice)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idAgeice,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idAgeice),   &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idAgeice),   &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice internal temperature.
!
        Vinfo( 1)=Vname(1,idTimid)
        Vinfo( 2)=Vname(2,idTimid)
        Vinfo( 3)=Vname(3,idTimid)
        Vinfo(14)=Vname(4,idTimid)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTimid,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTimid),    &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTimid),    &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define internal ice stress component 11.
!
        Vinfo( 1)=Vname(1,idSig11)
        Vinfo( 2)=Vname(2,idSig11)
        Vinfo( 3)=Vname(3,idSig11)
        Vinfo(14)=Vname(4,idSig11)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSig11,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig11),    &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig11),    &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define internal ice stress component 22.
!
        Vinfo( 1)=Vname(1,idSig22)
        Vinfo( 2)=Vname(2,idSig22)
        Vinfo( 3)=Vname(3,idSig22)
        Vinfo(14)=Vname(4,idSig22)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSig22,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig22),    &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig22),    &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define internal ice stress component 12.
!
        Vinfo( 1)=Vname(1,idSig12)
        Vinfo( 2)=Vname(2,idSig12)
        Vinfo( 3)=Vname(3,idSig12)
        Vinfo(14)=Vname(4,idSig12)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSig12,ng),r8)
# ifdef PERFECT_RESTART
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig12),    &
     &                 NF_FRST, nvd4, r2dgrd, Aval, Vinfo, ncname)
# else
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig12),    &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-water friction velocity.
!
        Vinfo( 1)=Vname(1,idTauiw)
        Vinfo( 2)=Vname(2,idTauiw)
        Vinfo( 3)=Vname(3,idTauiw)
        Vinfo(14)=Vname(4,idTauiw)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTauiw,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTauiw),    &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-water momentum transfer coefficient.
!
        Vinfo( 1)=Vname(1,idChuiw)
        Vinfo( 2)=Vname(2,idChuiw)
        Vinfo( 3)=Vname(3,idChuiw)
        Vinfo(14)=Vname(4,idChuiw)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idChuiw,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idChuiw),    &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define salinity of molecular sub-layer under ice.
!
        Vinfo( 1)=Vname(1,idS0mk)
        Vinfo( 2)=Vname(2,idS0mk)
        Vinfo( 3)=Vname(3,idS0mk)
        Vinfo(14)=Vname(4,idS0mk)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idS0mk,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idS0mk),     &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define temperature of molecular sub-layer under ice.
!
        Vinfo( 1)=Vname(1,idT0mk)
        Vinfo( 2)=Vname(2,idT0mk)
        Vinfo( 3)=Vname(3,idT0mk)
        Vinfo(14)=Vname(4,idT0mk)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idT0mk,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idT0mk),     &
     &                 NF_FRST, nvd3, sr2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef NCEP_FLUXES
!
!  Define NCEP gustiness squared.
!
        Vinfo( 1)=Vname(1,idWg2d)
        Vinfo( 2)=Vname(2,idWg2d)
        Vinfo( 3)=Vname(3,idWg2d)
        Vinfo(14)=Vname(4,idWg2d)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWg2d),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define NCEP air-ice/ocean momentum transfer coefficient.
!
        Vinfo( 1)=Vname(1,idCdd)
        Vinfo( 2)=Vname(2,idCdd)
        Vinfo( 3)=Vname(3,idCdd)
        Vinfo(14)=Vname(4,idCdd)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idCdd),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define NCEP air-ice/ocean specific heat transfer coefficient.
!
        Vinfo( 1)=Vname(1,idChd)
        Vinfo( 2)=Vname(2,idChd)
        Vinfo( 3)=Vname(3,idChd)
        Vinfo(14)=Vname(4,idChd)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idChd),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define NCEP air-ice/ocean latent heat transfer coefficient.
!
        Vinfo( 1)=Vname(1,idCed)
        Vinfo( 2)=Vname(2,idCed)
        Vinfo( 3)=Vname(3,idCed)
        Vinfo(14)=Vname(4,idCed)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idCed),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model gustiness squared.
!
        Vinfo( 1)=Vname(1,idWg2m)
        Vinfo( 2)=Vname(2,idWg2m)
        Vinfo( 3)=Vname(3,idWg2m)
        Vinfo(14)=Vname(4,idWg2m)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWg2m),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model air-ice/ocean momentum transfer coefficient.
!
        Vinfo( 1)=Vname(1,idCdm)
        Vinfo( 2)=Vname(2,idCdm)
        Vinfo( 3)=Vname(3,idCdm)
        Vinfo(14)=Vname(4,idCdm)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idCdm),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model air-ice/ocean specific heat transfer coefficient.
!
        Vinfo( 1)=Vname(1,idChm)
        Vinfo( 2)=Vname(2,idChm)
        Vinfo( 3)=Vname(3,idChm)
        Vinfo(14)=Vname(4,idChm)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idChm),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model air-ice/ocean latent heat transfer coefficient.
!
        Vinfo( 1)=Vname(1,idCem)
        Vinfo( 2)=Vname(2,idCem)
        Vinfo( 3)=Vname(3,idCem)
        Vinfo(14)=Vname(4,idCem)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idCem),      &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

!
!  Define NCEP near-surface air density.
!
        Vinfo( 1)=Vname(1,idRhoa)
        Vinfo( 2)=Vname(2,idRhoa)
        Vinfo( 3)=Vname(3,idRhoa)
        Vinfo(14)=Vname(4,idRhoa)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRhoa),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
#endif
#ifdef WEC
!
!  Define 2D u-stokes velocity.
!
        Vinfo( 1)=Vname(1,idU2Sd)
        Vinfo( 2)=Vname(2,idU2Sd)
        Vinfo( 3)=Vname(3,idU2Sd)
        Vinfo(14)=Vname(4,idU2Sd)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_u'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idU2Sd),     &
     &                 NF_FRST, nvd3, su2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D v-stokes velocity.
!
        Vinfo( 1)=Vname(1,idV2Sd)
        Vinfo( 2)=Vname(2,idV2Sd)
        Vinfo( 3)=Vname(3,idV2Sd)
        Vinfo(14)=Vname(4,idV2Sd)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  else
        Vinfo(20)='mask_v'
#  endif
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idV2Sd),     &
     &                 NF_FRST, nvd3, sv2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Define 3D u-stokes velocity.
!
        Vinfo( 1)=Vname(1,idU3Sd)
        Vinfo( 2)=Vname(2,idU3Sd)
        Vinfo( 3)=Vname(3,idU3Sd)
        Vinfo(14)=Vname(4,idU3Sd)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_u'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idU3Sd),     &
     &                 NF_FRST, nvd4, su3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D v-stokes velocity.
!
        Vinfo( 1)=Vname(1,idV3Sd)
        Vinfo( 2)=Vname(2,idV3Sd)
        Vinfo( 3)=Vname(3,idV3Sd)
        Vinfo(14)=Vname(4,idV3Sd)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_v'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idV3Sd),     &
     &                 NF_FRST, nvd4, sv3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D omega-stokes velocity.
!
        Vinfo( 1)=Vname(1,idW3Sd)
        Vinfo( 2)=Vname(2,idW3Sd)
        Vinfo( 3)=Vname(3,idW3Sd)
        Vinfo(14)=Vname(4,idW3Sd)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idW3Sd,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idW3Sd),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

!
!  Define 3D w-stokes velocity test
!
        Vinfo( 1)=Vname(1,idW3St)
        Vinfo( 2)=Vname(2,idW3St)
        Vinfo( 3)=Vname(3,idW3St)
        Vinfo(14)=Vname(4,idW3St)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if defined PERFECT_RESTART
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#   else
        Vinfo(20)='mask_rho'
#   endif
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idW3St,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idW3St),     &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# endif
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, RST(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, RST(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing restart file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.Ldefine) THEN
        ncname=RST(ng)%name
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
!  Open restart file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, RST(ng)%ncid)
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
!  restart variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            RST(ng)%Vid(idtime)=var_id(i)
#if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
            got_var(idbath)=.TRUE.
            RST(ng)%Vid(idbath)=var_id(i)
#endif
#if defined WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            RST(ng)%Vid(idRwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            RST(ng)%Vid(idUwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            RST(ng)%Vid(idVwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            RST(ng)%Vid(idPwet)=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            RST(ng)%Vid(idFsur)=var_id(i)
#ifdef PERFECT_RESTART
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            RST(ng)%Vid(idRzet)=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            RST(ng)%Vid(idUbar)=var_id(i)
#ifdef PERFECT_RESTART
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            RST(ng)%Vid(idRu2d)=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            RST(ng)%Vid(idVbar)=var_id(i)
#ifdef PERFECT_RESTART
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            RST(ng)%Vid(idRv2d)=var_id(i)
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            RST(ng)%Vid(idUvel)=var_id(i)
# ifdef PERFECT_RESTART
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            RST(ng)%Vid(idRu3d)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            RST(ng)%Vid(idVvel)=var_id(i)
# ifdef PERFECT_RESTART
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            RST(ng)%Vid(idRv3d)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            RST(ng)%Vid(idDano)=var_id(i)
# ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            RST(ng)%Vid(idHsbl)=var_id(i)
# endif
# ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            RST(ng)%Vid(idHbbl)=var_id(i)
# endif
# if defined PERFECT_RESTART && defined LMD_NONLOCAL
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idGhat(itemp)))) THEN
            got_var(idGhat(itemp))=.TRUE.
            RST(ng)%Vid(idGhat(itemp))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idGhat(isalt)))) THEN
            got_var(idGhat(isalt))=.TRUE.
            RST(ng)%Vid(idGhat(isalt))=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            RST(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            RST(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            RST(ng)%Vid(idSdif)=var_id(i)
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            RST(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            RST(ng)%Vid(idMtls)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVmLS))) THEN
            got_var(idVmLS)=.TRUE.
            RST(ng)%Vid(idVmLS)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVmKK))) THEN
            got_var(idVmKK)=.TRUE.
            RST(ng)%Vid(idVmKK)=var_id(i)
#  ifdef GLS_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVmKP))) THEN
            got_var(idVmKP)=.TRUE.
            RST(ng)%Vid(idVmKP)=var_id(i)
#  endif
# endif
#endif
#if defined WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            RST(ng)%Vid(idU2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            RST(ng)%Vid(idV2Sd)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            RST(ng)%Vid(idU3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            RST(ng)%Vid(idV3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            RST(ng)%Vid(idW3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            RST(ng)%Vid(idW3St)=var_id(i)
# endif
#endif
#if defined NEMURO_SED1
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONsed))) THEN
            got_var(idPONsed)=.TRUE.
            RST(ng)%Vid(idPONsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALsed))) THEN
            got_var(idOPALsed)=.TRUE.
            RST(ng)%Vid(idOPALsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDENITsed))) THEN
            got_var(idDENITsed)=.TRUE.
            RST(ng)%Vid(idDENITsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONbur))) THEN
            got_var(idPONbur)=.TRUE.
            RST(ng)%Vid(idPONbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALbur))) THEN
            got_var(idOPALbur)=.TRUE.
            RST(ng)%Vid(idOPALbur)=var_id(i)
#endif
#if defined PERFECT_RESTART && defined ICE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            RST(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            RST(ng)%Vid(idVsms)=var_id(i)
#endif
          END IF
#ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              RST(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
# if defined PERFECT_RESTART && defined ICE_MODEL
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              RST(ng)%Vid(idTsur(itrc))=var_id(i)
            END IF
          END DO
# endif
# ifdef BENTHIC
          DO itrc=1,NBeT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBeTvar(itrc)))) THEN
              got_var(idBeTvar(itrc))=.TRUE.
              rstBid(itrc,ng)=var_id(i)
            END IF
          END DO
# endif
# ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.                                   &
     &          TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.TRUE.
              RST(ng)%Vid(idfrac(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idBmas(itrc)))) THEN
              got_var(idBmas(itrc))=.TRUE.
              RST(ng)%Vid(idBmas(itrc))=var_id(i)
#  ifdef BEDLOAD
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              RST(ng)%Vid(idUbld(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              RST(ng)%Vid(idVbld(itrc))=var_id(i)
#  endif
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.TRUE.
              RST(ng)%Vid(idSbed(itrc))=var_id(i)
            END IF
          END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,6
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              got_var(idBott(itrc))=.TRUE.
              RST(ng)%Vid(idBott(itrc))=var_id(i)
            END IF
          END DO
# endif
#endif
#ifdef ICE_MODEL
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.TRUE.
            RST(ng)%Vid(idUice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.TRUE.
            RST(ng)%Vid(idVice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.TRUE.
            RST(ng)%Vid(idAice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.TRUE.
            RST(ng)%Vid(idHice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.TRUE.
            RST(ng)%Vid(idTice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.TRUE.
            RST(ng)%Vid(idHsno)=var_id(i)
# ifdef MELT_PONDS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idApond))) THEN
            got_var(idApond)=.TRUE.
            RST(ng)%Vid(idApond)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHpond))) THEN
            got_var(idHpond)=.TRUE.
            RST(ng)%Vid(idHpond)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.TRUE.
            RST(ng)%Vid(idAgeice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.TRUE.
            RST(ng)%Vid(idTimid)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.TRUE.
            RST(ng)%Vid(idSig11)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.TRUE.
            RST(ng)%Vid(idSig22)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.TRUE.
            RST(ng)%Vid(idSig12)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.TRUE.
            RST(ng)%Vid(idTauiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.TRUE.
            RST(ng)%Vid(idChuiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.TRUE.
            RST(ng)%Vid(idS0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.TRUE.
            RST(ng)%Vid(idT0mk)=var_id(i)
          END IF
#endif
#ifdef NCEP_FLUXES
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWg2d))) THEN
            got_var(idWg2d)=.TRUE.
            RST(ng)%Vid(idWg2d)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCdd))) THEN
            got_var(idCdd)=.TRUE.
            RST(ng)%Vid(idCdd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChd))) THEN
            got_var(idChd)=.TRUE.
            RST(ng)%Vid(idChd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCed))) THEN
            got_var(idCed)=.TRUE.
            RST(ng)%Vid(idCed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWg2m))) THEN
            got_var(idWg2m)=.TRUE.
            RST(ng)%Vid(idWg2m)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCdm))) THEN
            got_var(idCdm)=.TRUE.
            RST(ng)%Vid(idCdm)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChm))) THEN
            got_var(idChm)=.TRUE.
            RST(ng)%Vid(idChm)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCem))) THEN
            got_var(idCem)=.TRUE.
            RST(ng)%Vid(idCem)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRhoa))) THEN
            got_var(idRhoa)=.TRUE.
            RST(ng)%Vid(idRhoa)=var_id(i)
          END IF
#endif
        END DO
#if defined VEGETATION 
# if defined VEG_DRAG || defined VEG_BIOMASS 
        DO i=1,NVEGP
            IF (TRIM(var_name(i)).eq.                                   &
     &          TRIM(Vname(1,idvprp(i)))) THEN
              got_var(idvprp(i))=.TRUE.
              RST(ng)%Vid(idvprp(i))=var_id(i)
          END IF
        END DO
#  endif
#endif 
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#if defined SEDIMENT && defined SED_MORPH
        IF (.not.got_var(idbath)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idbath)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined WET_DRY
        IF (.not.got_var(idRwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idPwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef PERFECT_RESTART
        IF (.not.got_var(idRzet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRzet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef PERFECT_RESTART
        IF (.not.got_var(idRu2d)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRu2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef PERFECT_RESTART
        IF (.not.got_var(idRv2d)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRv2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef PERFECT_RESTART
        IF (.not.got_var(idRu3d)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRu3d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef PERFECT_RESTART
        IF (.not.got_var(idRv3d)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRv3d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined PERFECT_RESTART && \
     (defined GLS_MIXING     || defined MY25_MIXING)
        IF (.not.got_var(idMtke)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVmLS)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVmLS)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVmKK)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVmKK)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef GLS_MIXING
        IF (.not.got_var(idVmKP)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVmKP)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
#endif
#if defined WEC
        IF (.not.got_var(idU2Sd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idU3Sd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef NEMURO_SED1
        IF (.not.got_var(idPONsed)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idPONsed)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALsed)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idOPALsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDENITsed)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idDENITsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPONbur)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idPONbur)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALbur)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idOPALbur)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# ifdef VEGETATION 
#  if defined VEG_DRAG || defined VEG_BIOMASS 
        DO i=1,NVEGP
          IF (.not.got_var(idvprp(i)).and.Hout(idvprp(i),ng)) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idvprp(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif 
# ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idfrac(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idBmas(i))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idBmas(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef BEDLOAD
          IF (.not.got_var(idUbld(i))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSbed(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO i=1,6
          IF (.not.got_var(idBott(i))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idBott(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
#ifdef ICE_MODEL
        IF (.not.got_var(idUice)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef MELT_PONDS
        IF (.not.got_var(idApond)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idApond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHpond)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idHpond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idAgeice)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef NCEP_FLUXES
        IF (.not.got_var(idWg2d)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWg2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
        RETURN
        END IF
        IF (.not.got_var(idCdd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idCdd)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChd)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idChd)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCed)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idCed)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWg2m)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idWg2m)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCdm)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idCdm)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChm)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idChm)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idCem)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idCem)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRhoa)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRhoa)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  Set unlimited time record dimension to current value.
!
        IF (LcycleRST(ng)) THEN
          RST(ng)%Rindex=0
        ELSE
          RST(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
  10  FORMAT (/,' DEF_RST - unable to create restart NetCDF file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_RST - unable to open restart NetCDF file: ',a)
  40  FORMAT (/,' DEF_RST - unable to find variable: ',a,2x,            &
     &        ' in restart NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_rst
