#include "cppdefs.h"
      SUBROUTINE wrt_info (ng, model, ncid, ncname)
!
!svn $Id: wrt_info.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     ncid     NetCDF file ID (integer)                                !
!     ncname   NetCDF file name (string)                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
      USE mod_grid
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#ifdef SEDIMENT
      USE mod_sediment
#endif
#if defined SOLVE3D && defined TS_PSOURCE
      USE mod_sources
#endif

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
#ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
#endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.

      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status, varid

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
      integer :: ifield = 0

      real(r8) :: scale
#ifdef SOLVE3D
# ifdef TS_DIF4
      real(r8), dimension(NT(ng)) :: diff
# endif
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
#endif
#ifdef STATIONS
      real(r8), dimension(Nstation(ng)) :: Zpos, wrk
#endif

      character (len=1) :: char1
      character (len=1), dimension(NT(ng)) :: char1_trc
!
      SourceFile='wrt_info.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname, ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Time stepping parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes',                &
     &                      ntimes(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndtfast',               &
     &                      ndtfast(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dt',                    &
     &                      dt(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dtfast',                &
     &                      dtfast(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dstart',                &
     &                      dstart, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#if defined NETCDF4 && defined DEFLATE
      CALL netcdf_put_ivar (ng, model, ncname, 'shuffle',               &
     &                      shuffle, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'deflate',               &
     &                      deflate, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'deflate_level',         &
     &                      deflate_level, (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

      CALL netcdf_put_ivar (ng, model, ncname, 'nHIS',                  &
     &                      nHIS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefHIS',               &
     &                      ndefHIS(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nRST',                  &
     &                      nRST(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef AVERAGES
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsAVG',                &
     &                      ntsAVG(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nAVG',                  &
     &                      nAVG(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefAVG',               &
     &                      ndefAVG(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef ADJOINT
      CALL netcdf_put_ivar (ng, model, ncname, 'nADJ',                  &
     &                      nADJ(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefADJ',               &
     &                      ndefADJ(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef TANGENT
      CALL netcdf_put_ivar (ng, model, ncname, 'nTLM',                  &
     &                      nTLM(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefTLM',               &
     &                      ndefTLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef ADJUST_BOUNDARY
      CALL netcdf_put_ivar (ng, model, ncname, 'nOBC',                  &
     &                      nOBC(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      CALL netcdf_put_ivar (ng, model, ncname, 'nSFF',                  &
     &                      nSFF(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef PROPAGATOR
      WRITE (char1,'(l1)') LrstGST
      CALL netcdf_put_svar (ng, model, ncname, 'LrstGST',               &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'MaxIterGST',            &
     &                      MaxIterGST, (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nGST',                  &
     &                      nGST, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN


      CALL netcdf_put_fvar (ng, model, ncname, 'Ritz_tol',              &
     &                      Ritz_tol, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef DIAGNOSTICS
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsDIA',                &
     &                      ntsDIA(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nDIA',                  &
     &                      nDIA(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefDIA',               &
     &                      ndefDIA(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef STATIONS
      CALL netcdf_put_ivar (ng, model, ncname, 'nSTA',                  &
     &                      nSTA(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef FOUR_DVAR
      CALL netcdf_put_ivar (ng, model, ncname, 'Nouter',                &
     &                      Nouter, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Ninner',                &
     &                      Ninner, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Falpha',                &
     &                      Falpha, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Fbeta',                 &
     &                      Fbeta, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Fgamma',                &
     &                      Fgamma, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_tnu2',               &
     &                      nl_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_tnu2',               &
     &                      ad_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_tnu2',               &
     &                      tl_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

#endif

#if defined SOLVE3D && defined TS_DIF4
      DO itrc=1,NT(ng)
        diff(itrc)=nl_tnu4(itrc,ng)*nl_tnu4(itrc,ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_tnu4',               &
     &                      diff, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      DO itrc=1,NT(ng)
        diff(itrc)=ad_tnu4(itrc,ng)*ad_tnu4(itrc,ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_tnu4',               &
     &                      diff, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      DO itrc=1,NT(ng)
        diff(itrc)=tl_tnu4(itrc,ng)*tl_tnu4(itrc,ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_tnu4',               &
     &                      diff, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef UV_VIS2
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_visc2',              &
     &                      nl_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_visc2',              &
     &                      ad_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_visc2',              &
     &                      tl_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef UV_VIS4
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_visc4',              &
     &                      nl_visc4(ng)**2, (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_visc4',              &
     &                      ad_visc4(ng)**2, (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_visc4',              &
     &                      tl_visc4(ng)**2, (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
      CALL netcdf_put_fvar (ng, model, ncname, 'tkenu2',                &
     &                      tkenu2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef TKE_DIF4
      CALL netcdf_put_fvar (ng, model, ncname, 'tkenu4',                &
     &                      tkenu4(ng)**2, (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Akt_bak',               &
     &                      Akt_bak(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Akv_bak',               &
     &                      Akv_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# if defined MY25_MIXING || defined GLS_MIXING
      CALL netcdf_put_fvar (ng, model, ncname, 'Akk_bak',               &
     &                      Akk_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Akp_bak',               &
     &                      Akp_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef FORWARD_MIXING
!
!  Basic state vertical mixing scale used in adjoint-based applications.
!
#  ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_Akt_fac',            &
     &                      ad_Akt_fac(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_Akt_fac',            &
     &                      tl_Akt_fac(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_Akv_fac',            &
     &                      ad_Akv_fac(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_Akv_fac',            &
     &                      tl_Akv_fac(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
#endif
!
!  Drag coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg',                  &
     &                      rdrg(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg2',                 &
     &                      rdrg2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'Zob',                   &
     &                      Zob(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Zos',                   &
     &                      Zos(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_p',                 &
     &                      gls_p(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_m',                 &
     &                      gls_m(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_n',                 &
     &                      gls_n(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_cmu0',              &
     &                      gls_cmu0(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c1',                &
     &                      gls_c1(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c2',                &
     &                      gls_c2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3m',               &
     &                      gls_c3m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3p',               &
     &                      gls_c3p(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigk',              &
     &                      gls_sigk(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigp',              &
     &                      gls_sigp(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Kmin',              &
     &                      gls_Kmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Pmin',              &
     &                      gls_Pmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Charnok_alpha',         &
     &                      charnok_alpha(ng), (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Zos_hsig_alpha',        &
     &                      zos_hsig_alpha(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'sz_alpha',              &
     &                      sz_alpha(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'CrgBan_cw',             &
     &                      crgban_cw(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Nudging inverse time scales used in various tasks.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudg',                 &
     &                      Znudg(ng)/sec2day, (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudg',                &
     &                      M2nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'M3nudg',                &
     &                      M3nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tnudg',                 &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef NUDGING
!
!  Nudging inverse time scales used in data assimilation.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudass',               &
     &                      Znudass(ng)/sec2day, (/0/), (/0/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudass',              &
     &                      M2nudass(ng)/sec2day, (/0/), (/0/),         &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'M3nudass',              &
     &                      M3nudass(ng)/sec2day, (/0/), (/0/),         &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudass(itrc,ng)/sec2day
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tnudass',               &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef NUDGING_COFF
!
!  Open boundary nudging, inverse time scales.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_in',              &
     &                      FSobc_in(ng,:), (/1/), (/4/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_out',             &
     &                      FSobc_out(ng,:), (/1/), (/4/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_in',              &
     &                      M2obc_in(ng,:), (/1/), (/4/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_out',             &
     &                      M2obc_out(ng,:), (/1/), (/4/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      DO ibry=1,4
        DO itrc=1,NT(ng)
          Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
        END DO
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_in',               &
     &                      Tobc, (/1,1/), (/NT(ng),4/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      DO ibry=1,4
        DO itrc=1,NT(ng)
          Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
        END DO
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_out',              &
     &                      Tobc, (/1,1/), (/NT(ng),4/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_in',              &
     &                      M3obc_in(ng,:), (/1/), (/4/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_out',             &
     &                      M3obc_out(ng,:), (/1/), (/4/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif
!
!  Equation of State parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rho0',                  &
     &                      rho0, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#if defined SOLVE3D && defined PROPAGATOR
      CALL netcdf_put_fvar (ng, model, ncname, 'bvf_bak',               &
     &                      bvf_bak, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined SOLVE3D     && \
   (!defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR)
      CALL netcdf_put_fvar (ng, model, ncname, 'R0',                    &
     &                      R0(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Tcoef',                 &
     &                      Tcoef(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Scoef',                 &
     &                      Scoef(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SOLVE3D
# ifdef BODYFORCE
!
!  Body force parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'levsfrc',               &
     &                      levsfrc(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'levbfrc',               &
     &                      levbfrc(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif
!
!  Slipperiness parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gamma2',                &
     &                      gamma2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#if defined SOLVE3D && defined TS_PSOURCE
!
!  Logical switches indicating which tracer variables are processed
!  during point Sources/Sinks.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',            &
     &                      LtracerSrc(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef FOUR_DVAR
!
!  4DVAR assimilation parameters.
!
# ifdef ADJUST_STFLUX
      CALL netcdf_put_lvar (ng, model, ncname, 'Lstflux',               &
     &                      Lstflux(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef ADJUST_BOUNDARY
      CALL netcdf_put_lvar (ng, model, ncname, 'Lobc',                  &
     &                      Lobc(:,:,ng), (/1,1/), (/4,NstateVar(ng)/), &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifndef IS4DVAR_SENSITIVITY
      CALL netcdf_put_lvar (ng, model, ncname, 'LhessianEV',            &
     &                      LhessianEV, (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
      CALL netcdf_put_lvar (ng, model, ncname, 'LhotStart',             &
     &                      LhotStart, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

      CALL netcdf_put_lvar (ng, model, ncname, 'Lprecond',              &
     &                      Lprecond, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Lritz',                 &
     &                      Lritz, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
      IF (Lprecond.and.(NritzEV.gt.0)) THEN
        CALL netcdf_put_ivar (ng, model, ncname, 'NritzEV',             &
     &                        NritzEV, (/0/), (/0/),                    &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
# endif
# if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
        CALL netcdf_put_ivar (ng, model, ncname, 'NpostI',              &
     &                        NpostI, (/0/), (/0/),                     &
     &                        ncid = ncid)
# endif
# ifndef IS4DVAR_SENSITIVITY
      CALL netcdf_put_fvar (ng, model, ncname, 'GradErr',               &
     &                      GradErr, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HevecErr',              &
     &                      HevecErr, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

      CALL netcdf_put_ivar (ng, model, ncname, 'Nmethod',               &
     &                      Nmethod(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Rscheme',               &
     &                      Rscheme(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Nrandom',               &
     &                      Nrandom, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Hgamma',                &
     &                      Hgamma(1), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef WEAK_CONSTRAINT
      CALL netcdf_put_fvar (ng, model, ncname, 'HgammaM',               &
     &                      Hgamma(2), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef ADJUST_BOUNDARY
      CALL netcdf_put_fvar (ng, model, ncname, 'HgammaB',               &
     &                      Hgamma(3), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef ADJUST_STFLUX
      CALL netcdf_put_fvar (ng, model, ncname, 'HgammaF',               &
     &                      Hgamma(4), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'Vgamma',                &
     &                      Vgamma(1), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
      CALL netcdf_put_fvar (ng, model, ncname, 'VgammaM',               &
     &                      Vgamma(2), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ADJUST_BOUNDARY
      CALL netcdf_put_fvar (ng, model, ncname, 'VgammaB',               &
     &                      Vgamma(3), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

      CALL netcdf_put_fvar (ng, model, ncname, 'Hdecay',                &
     &                      Hdecay(1,:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'Vdecay',                &
     &                      Vdecay(1,:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WEAK_CONSTRAINT
      CALL netcdf_put_fvar (ng, model, ncname, 'HdecayM',               &
     &                      Hdecay(2,:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'VdecayM',               &
     &                      Vdecay(2,:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef ADJUST_BOUNDARY
      CALL netcdf_put_fvar (ng, model, ncname, 'HdecayB',               &
     &                      HdecayB(:,:,ng),                            &
     &                      (/1,1/), (/NstateVar(ng),4/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'VdecayB',               &
     &                      VdecayB(:,:,ng),                            &
     &                      (/1,1/), (/NstateVar(ng),4/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef RPM_RELAXATION
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_M2diff',             &
     &                      tl_M2diff(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_M3diff',             &
     &                      tl_M3diff(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'tl_Tdiff',              &
     &                      tl_Tdiff(:,ng), (/1/), (/NT(ng)/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef BALANCE_OPERATOR
#  ifdef ZETA_ELLIPTIC
      CALL netcdf_put_ivar (ng, model, ncname, 'Nbico',                 &
     &                      Nbico(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
#  endif

      CALL netcdf_put_lvar (ng, model, ncname, 'Lbalance',              &
     &                      balance, (/1/), (/NstateVar(ng)/),          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'LNM_flag',              &
     &                      LNM_flag, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'LNM_depth',             &
     &                      LNM_depth(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dTdz_min',              &
     &                      dTdz_min(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ml_depth',              &
     &                      ml_depth(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
!
!  Adjoint sensitivity parameters.
!
      WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isFsur)
      CALL netcdf_put_svar (ng, model, ncname, 'Lzeta',                 &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isUbar)
      CALL netcdf_put_svar (ng, model, ncname, 'Lubar',                 &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isVbar)
      CALL netcdf_put_svar (ng, model, ncname, 'Lvbar',                 &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isUvel)
      CALL netcdf_put_svar (ng, model, ncname, 'Luvel',                 &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isVvel)
      CALL netcdf_put_svar (ng, model, ncname, 'Lvvel',                 &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      DO i=1,NT(ng)
        WRITE (char1_trc(i),'(l1)') SCALARS(ng)%Lstate(isTvar(i))
      END DO
      CALL netcdf_put_svar (ng, model, ncname, 'Ltracer',               &
     &                      char1_trc, (/1/), (/NT(ng)/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'KstrS',                 &
     &                      KstrS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'KendS',                 &
     &                      KendS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SO_SEMI
#  ifndef SO_SEMI_WHITE
      CALL netcdf_put_fvar (ng, model, ncname, 'SO_decay',              &
     &                      SO_decay(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

      CALL netcdf_put_fvar (ng, model, ncname, 'SO_trace',              &
     &                      TRnorm(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      WRITE (char1,'(l1)') SCALARS(ng)%SOstate(isUstr)
      CALL netcdf_put_svar (ng, model, ncname, 'Lso_Ustr',              &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      WRITE (char1,'(l1)') SCALARS(ng)%SOstate(isVstr)
      CALL netcdf_put_svar (ng, model, ncname, 'Lso_Vstr',              &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      DO i=1,NT(ng)
        WRITE (char1_trc(i),'(l1)') SCALARS(ng)%SOstate(isTsur(i))
      END DO
      CALL netcdf_put_svar (ng, model, ncname, 'Lso_tracer',            &
     &                      char1_trc, (/1/), (/NT(ng)/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif

      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_Ustr',           &
     &                      SO_sdev(isUstr,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN


      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_Vstr',           &
     &                      SO_sdev(isUstr,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      DO itrc=1,NT(ng)
        nudg(itrc)=SO_sdev(isTsur(itrc),ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_stflx',          &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
#endif

#if defined BIOLOGY && defined SOLVE3D
# if defined BIO_FENNEL
#  include <fennel_wrt.h>
# elif defined ECOSIM
#  include <ecosim_wrt.h>
# elif defined NEMURO
#  include <nemuro_wrt.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_wrt.h>
# elif defined NPZD_IRON
#  include <npzd_iron_wrt.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_wrt.h>
# endif
#endif

#ifdef SEDIMENT
# include <sediment_wrt.h>
#endif
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch.
!
      WRITE (char1,'(l1)') spherical
      CALL netcdf_put_svar (ng, model, ncname, 'spherical',             &
     &                      char1, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Domain Length.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'xl',                    &
     &                      xl(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'el',                    &
     &                      el(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'Vtransform',            &
     &                      Vtransform(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Vstretching',           &
     &                      Vstretching(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'theta_s',               &
     &                      theta_s(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'theta_b',               &
     &                      theta_b(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Tcline',                &
     &                      Tcline(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'hc',                    &
     &                      hc(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL netcdf_put_fvar (ng, model, ncname, 's_rho',                 &
     &                      SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 's_w',                   &
     &                      SCALARS(ng)%sc_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_r',                  &
     &                      SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_w',                  &
     &                      SCALARS(ng)%Cs_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'user',                &
     &                        user, (/1/), (/Nuser/),                   &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF

#ifdef STATIONS
!
!  Stations positions.
!
      IF (ncid.eq.ncSTAid(ng)) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'Ipos',                &
     &                        SCALARS(ng)%SposX(:), (/1/),              &
     &                        (/Nstation(ng)/), ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'Jpos',                &
     &                        SCALARS(ng)%SposY(:), (/1/),              &
     &                        (/Nstation(ng)/), ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif

#ifdef NO_WRITE_GRID
      GRID_VARS : IF (ncid.eq.ncSTAid(ng)) THEN
#else
      GRID_VARS : IF (ncid.ne.ncFLTid(ng)) THEN
#endif
#if !(defined SED_MORPH && defined SEDIMENT)
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_r8
          IF (ncid.ne.ncSTAid(ng)) THEN
            IF (find_string(var_name, n_var, 'h', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'h', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
# ifdef STATIONS
          ELSE
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%h,                      &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            CALL netcdf_put_fvar (ng, model, ncname, 'h',               &
     &                            wrk, (/1/), (/Nstation(ng)/),         &
     &                            ncid = ncid)
            IF (exit_flag.ne.NoError) RETURN
# endif
          END IF
        END IF
#endif
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'f', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#endif
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'f', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'pm', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#endif
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'pm', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'pn', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#endif
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'pn', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.ncSTAid(ng)) THEN
              IF (find_string(var_name, n_var, 'lon_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
#ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%lonr,                 &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, 'lon_rho',       &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (exit_flag.ne.NoError) RETURN
#endif
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.ncSTAid(ng)) THEN
              IF (find_string(var_name, n_var, 'lat_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
#ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%latr,                 &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, 'lat_rho',       &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (exit_flag.ne.NoError) RETURN
#endif
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.ncSTAid(ng)) THEN
              IF (find_string(var_name, n_var, 'x_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
#ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%xr,                     &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, 'x_rho',         &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (exit_flag.ne.NoError) RETURN
#endif
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.ncSTAid(ng)) THEN
              IF (find_string(var_name, n_var, 'y_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
#ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%yr,                   &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, 'y_rho',         &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (exit_flag.ne.NoError) RETURN
# endif
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % pmask,                    &
#endif
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_p', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_p', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % pmask,                    &
#endif
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_p', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.ncSTAid(ng)) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_p', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF

#ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_r8
          IF (ncid.ne.ncSTAid(ng)) THEN
            IF (find_string(var_name, n_var, 'angle', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % angler,                     &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'angle', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'angle', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
# ifdef STATIONS
          ELSE
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%angler,                 &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            CALL netcdf_put_fvar (ng, model, ncname, 'angle',           &
     &                            wrk, (/1/), (/Nstation(ng)/),         &
     &                            ncid = ncid)
            IF (exit_flag.ne.NoError) RETURN
# endif
          END IF
        END IF
#endif

#ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_rho', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_u', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_v', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_psi', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_psi', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
!
!  Adjoint sensitivity spatial scope mask at RHO-, U-, and V-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'scope_rho', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % Rscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'scope_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'scope_rho', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'scope_u', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % umask,                      &
# endif
     &                           GRID(ng) % Uscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'scope_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'scope_u', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'scope_v', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % Vscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'scope_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'scope_v', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
#endif
      END IF GRID_VARS
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (exit_flag.ne.NoError) RETURN

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' WRT_INFO - error while writing variable: ',a,/,       &
     &        12x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO - error while inquiring ID for variable: ',  &
     &        a,/,12x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO - unable to synchronize to disk file: ',     &
     &        /,12x,a)

      RETURN
      END SUBROUTINE wrt_info
