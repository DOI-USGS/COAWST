#include "cppdefs.h"
      MODULE wrt_info_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#if defined FLOATS && defined FLOAT_BIOLOGY
      USE mod_behavior
#endif
#ifdef BIOLOGY
      USE mod_biology
# ifdef ECOSIM
      USE mod_eclight
# endif
#endif
#ifdef GRID_EXTRACT
      USE mod_extract
#endif
      USE mod_grid
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#ifdef INWAVE_MODEL
      USE mod_inwave_vars
      USE mod_inwave_bound
#endif
#ifdef SEDIMENT
      USE mod_sediment
#endif
#ifdef VEGETATION 
      USE mod_vegetation
#endif
#ifdef PROPAGATOR
      USE mod_storage
#endif
      USE mod_sources
!
#ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
#endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      INTERFACE wrt_info
        MODULE PROCEDURE wrt_info_nf90
#if defined PIO_LIB && defined DISTRIBUTE
        MODULE PROCEDURE wrt_info_pio
#endif
      END INTERFACE wrt_info
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_info_nf90 (ng, model, ncid, ncname)
!***********************************************************************
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file using the standard NetCDF-3 or NetCDF-4 library.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     ncname       NetCDF filename (string)                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod,  ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status, varid

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
      integer :: ifield = 0
!
      real(dp) :: scale
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
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_info_nf90"
!
      SourceFile=MyFile
!
      IF (ncid.ne.XTR(ng)%ncid) THEN
        LBi=LBOUND(GRID(ng)%h,DIM=1)
        UBi=UBOUND(GRID(ng)%h,DIM=1)
        LBj=LBOUND(GRID(ng)%h,DIM=2)
        UBj=UBOUND(GRID(ng)%h,DIM=2)
#ifdef GRID_EXTRACT
      ELSE
        LBi=LBOUND(EXTRACT(ng)%h,DIM=1)
        UBi=UBOUND(EXTRACT(ng)%h,DIM=1)
        LBj=LBOUND(EXTRACT(ng)%h,DIM=2)
        UBj=UBOUND(EXTRACT(ng)%h,DIM=2)
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time stepping parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes',                &
     &                      ntimes(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndtfast',               &
     &                      ndtfast(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dt',                    &
     &                      dt(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dtfast',                &
     &                      dtfast(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dstart',                &
     &                      dstart, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef RBL4DVAR_FCT_SENSITIVITY
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes_ana',            &
     &                      ntimes_ana(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes_fct',            &
     &                      ntimes_fct(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined HDF5 && defined DEFLATE
      CALL netcdf_put_ivar (ng, model, ncname, 'shuffle',               &
     &                      shuffle, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'deflate',               &
     &                      deflate, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'deflate_level',         &
     &                      deflate_level, (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

      CALL netcdf_put_ivar (ng, model, ncname, 'nHIS',                  &
     &                      nHIS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefHIS',               &
     &                      ndefHIS(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef GRID_EXTRACT
      CALL netcdf_put_ivar (ng, model, ncname, 'nXTR',                  &
     &                      nXTR(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefXTR',               &
     &                      ndefXTR(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ExtractFlag',           &
     &                      ExtractFlag(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

      CALL netcdf_put_ivar (ng, model, ncname, 'nRST',                  &
     &                      nRST(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsAVG',                &
     &                      ntsAVG(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nAVG',                  &
     &                      nAVG(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefAVG',               &
     &                      ndefAVG(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef ADJOINT
      CALL netcdf_put_ivar (ng, model, ncname, 'nADJ',                  &
     &                      nADJ(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefADJ',               &
     &                      ndefADJ(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef TANGENT
      CALL netcdf_put_ivar (ng, model, ncname, 'nTLM',                  &
     &                      nTLM(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefTLM',               &
     &                      ndefTLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef ADJUST_BOUNDARY
      CALL netcdf_put_ivar (ng, model, ncname, 'nOBC',                  &
     &                      nOBC(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      CALL netcdf_put_ivar (ng, model, ncname, 'nSFF',                  &
     &                      nSFF(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef PROPAGATOR
      CALL netcdf_put_lvar (ng, model, ncname, 'LmultiGST',             &
     &                      LmultiGST, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'LrstGST',               &
     &                      LrstGST, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'MaxIterGST',            &
     &                      MaxIterGST, (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nGST',                  &
     &                      nGST, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'NEV',                   &
     &                      NEV, (/0/), (/0/),                          &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'NCV',                   &
     &                      NCV, (/0/), (/0/),                          &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Ritz_tol',              &
     &                      Ritz_tol, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef DIAGNOSTICS
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsDIA',                &
     &                      ntsDIA(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'nDIA',                  &
     &                      nDIA(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'ndefDIA',               &
     &                      ndefDIA(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef STATIONS
      CALL netcdf_put_ivar (ng, model, ncname, 'nSTA',                  &
     &                      nSTA(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef FOUR_DVAR
      CALL netcdf_put_ivar (ng, model, ncname, 'Nouter',                &
     &                      Nouter, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Ninner',                &
     &                      Ninner, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Falpha',                &
     &                      Falpha, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Fbeta',                 &
     &                      Fbeta, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Fgamma',                &
     &                      Fgamma, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_tnu2',               &
     &                      nl_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_tnu2',               &
     &                      ad_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_tnu2',               &
     &                      tl_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

#endif

#if defined SOLVE3D && defined TS_DIF4
      DO itrc=1,NT(ng)
        diff(itrc)=nl_tnu4(itrc,ng)*nl_tnu4(itrc,ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_tnu4',               &
     &                      diff, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJOINT
      DO itrc=1,NT(ng)
        diff(itrc)=ad_tnu4(itrc,ng)*ad_tnu4(itrc,ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_tnu4',               &
     &                      diff, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      DO itrc=1,NT(ng)
        diff(itrc)=tl_tnu4(itrc,ng)*tl_tnu4(itrc,ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_tnu4',               &
     &                      diff, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#ifdef UV_VIS2
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_visc2',              &
     &                      nl_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_visc2',              &
     &                      ad_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_visc2',              &
     &                      tl_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#ifdef UV_VIS4
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_visc4',              &
     &                      nl_visc4(ng)**2, (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_visc4',              &
     &                      ad_visc4(ng)**2, (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_visc4',              &
     &                      tl_visc4(ng)**2, (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
      CALL netcdf_put_fvar (ng, model, ncname, 'tkenu2',                &
     &                      tkenu2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef TKE_DIF4
      CALL netcdf_put_fvar (ng, model, ncname, 'tkenu4',                &
     &                      tkenu4(ng)**2, (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#if defined UV_VIS2 || defined UV_VIS4
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSponge',             &
     &                      LuvSponge(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if (defined TS_DIF2 || defined TS_DIF4) && defined SOLVE3D
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSponge',         &
     &                      LtracerSponge(:,ng), (/1/), (/NT(ng)/),     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Akt_bak',               &
     &                      Akt_bak(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Akv_bak',               &
     &                      Akv_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined MY25_MIXING || defined GLS_MIXING
      CALL netcdf_put_fvar (ng, model, ncname, 'Akk_bak',               &
     &                      Akk_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Akp_bak',               &
     &                      Akp_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef FORWARD_MIXING
!
!  Basic state vertical mixing scale used in adjoint-based applications.
!
#  ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_Akt_fac',            &
     &                      ad_Akt_fac(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_Akt_fac',            &
     &                      tl_Akt_fac(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ADJOINT
      CALL netcdf_put_fvar (ng, model, ncname, 'ad_Akv_fac',            &
     &                      ad_Akv_fac(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_Akv_fac',            &
     &                      tl_Akv_fac(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
#endif
!
!  Drag coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg',                  &
     &                      rdrg(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg2',                 &
     &                      rdrg2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'Zob',                   &
     &                      Zob(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Zos',                   &
     &                      Zos(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_p',                 &
     &                      gls_p(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_m',                 &
     &                      gls_m(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_n',                 &
     &                      gls_n(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_cmu0',              &
     &                      gls_cmu0(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c1',                &
     &                      gls_c1(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c2',                &
     &                      gls_c2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3m',               &
     &                      gls_c3m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3p',               &
     &                      gls_c3p(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigk',              &
     &                      gls_sigk(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigp',              &
     &                      gls_sigp(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Kmin',              &
     &                      gls_Kmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Pmin',              &
     &                      gls_Pmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Charnok_alpha',         &
     &                      charnok_alpha(ng), (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Zos_hsig_alpha',        &
     &                      zos_hsig_alpha(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'sz_alpha',              &
     &                      sz_alpha(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'CrgBan_cw',             &
     &                      crgban_cw(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef WEC
      CALL netcdf_put_fvar (ng, model, ncname, 'wec_alpha',             &
     &                      wec_alpha(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Nudging inverse time scales used in various tasks.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudg',                 &
     &                      Znudg(ng)/sec2day, (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudg',                &
     &                      M2nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'M3nudg',                &
     &                      M3nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tnudg',                 &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifndef DEBUGGING
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_in',            &
     &                        FSobc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_out',           &
     &                        FSobc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_in',            &
     &                        M2obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_out',           &
     &                        M2obc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_in',             &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_out',            &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_in',            &
     &                        M3obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_out',           &
     &                        M3obc_out(ng,:), (/1/), (/4/),            &
     &                      ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
      END IF
#endif
!
!  Equation of State parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rho0',                  &
     &                      rho0, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#if defined SOLVE3D && defined PROPAGATOR
      CALL netcdf_put_fvar (ng, model, ncname, 'bvf_bak',               &
     &                      bvf_bak, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined SOLVE3D     && \
   (!defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR)
      CALL netcdf_put_fvar (ng, model, ncname, 'R0',                    &
     &                      R0(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Tcoef',                 &
     &                      Tcoef(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Scoef',                 &
     &                      Scoef(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef SOLVE3D
# ifdef BODYFORCE
!
!  Body force parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'levsfrc',               &
     &                      levsfrc(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'levbfrc',               &
     &                      levbfrc(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif
!
!  Slipperiness parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gamma2',                &
     &                      gamma2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSrc',                &
     &                      LuvSrc(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'LwSrc',                 &
     &                      LwSrc(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef SOLVE3D
!
!  Logical switches to activate tracer point Sources/Sinks.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',            &
     &                      LtracerSrc(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Logical switches to process climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LsshCLM',               &
     &                      LsshCLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',                &
     &                      Lm2CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef SOLVE3D
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm3CLM',                &
     &                      Lm3CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerCLM',            &
     &                      LtracerCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Logical switches for nudging climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',           &
     &                      LnudgeM2CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef SOLVE3D
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM3CLM',           &
     &                      LnudgeM3CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeTCLM',            &
     &                      LnudgeTCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef FOUR_DVAR
!
!  4DVAR assimilation parameters.
!
# ifdef ADJUST_STFLUX
      CALL netcdf_put_lvar (ng, model, ncname, 'Lstflux',               &
     &                      Lstflux(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef ADJUST_BOUNDARY
      CALL netcdf_put_lvar (ng, model, ncname, 'Lobc',                  &
     &                      Lobc(:,:,ng), (/1,1/), (/4,NstateVar(ng)/), &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifndef I4DVAR_ANA_SENSITIVITY
      CALL netcdf_put_lvar (ng, model, ncname, 'LhessianEV',            &
     &                      LhessianEV, (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef WEAK_CONSTRAINT
      CALL netcdf_put_lvar (ng, model, ncname, 'LhotStart',             &
     &                      LhotStart, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

      CALL netcdf_put_lvar (ng, model, ncname, 'Lprecond',              &
     &                      Lprecond, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Lritz',                 &
     &                      Lritz, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef WEAK_CONSTRAINT
      IF (Lprecond.and.(NritzEV.gt.0)) THEN
        CALL netcdf_put_ivar (ng, model, ncname, 'NritzEV',             &
     &                        NritzEV, (/0/), (/0/),                    &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
# endif
# if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
        CALL netcdf_put_ivar (ng, model, ncname, 'NpostI',              &
     &                        NpostI, (/0/), (/0/),                     &
     &                        ncid = ncid)
# endif
# if defined ARRAY_MODES              || \
     defined I4DVAR_ANA_SENSITIVITY   || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined R4DVAR_ANA_SENSITIVITY
        CALL netcdf_put_ivar (ng, model, ncname, 'Nimpact',             &
     &                        Nimpact, (/0/), (/0/),                    &
     &                        ncid = ncid)
# endif
# ifndef I4DVAR_ANA_SENSITIVITY
      CALL netcdf_put_fvar (ng, model, ncname, 'GradErr',               &
     &                      GradErr, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HevecErr',              &
     &                      HevecErr, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      CALL netcdf_put_ivar (ng, model, ncname, 'Nmethod',               &
     &                      Nmethod(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Rscheme',               &
     &                      Rscheme(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Nrandom',               &
     &                      Nrandom, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Hgamma',                &
     &                      Hgamma(1), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef WEAK_CONSTRAINT
      CALL netcdf_put_fvar (ng, model, ncname, 'HgammaM',               &
     &                      Hgamma(2), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
      CALL netcdf_put_fvar (ng, model, ncname, 'HgammaB',               &
     &                      Hgamma(3), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_STFLUX
      CALL netcdf_put_fvar (ng, model, ncname, 'HgammaF',               &
     &                      Hgamma(4), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'Vgamma',                &
     &                      Vgamma(1), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef WEAK_CONSTRAINT
      CALL netcdf_put_fvar (ng, model, ncname, 'VgammaM',               &
     &                      Vgamma(2), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ADJUST_BOUNDARY
      CALL netcdf_put_fvar (ng, model, ncname, 'VgammaB',               &
     &                      Vgamma(3), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef MULTI_SCALE_B
#  ifndef NONUNIFORM_SCALES
      CALL netcdf_put_fvar (ng, model, ncname, 'HdecayX',               &
     &                      HdecayX(1,:,:,ng),                          &
     &                      (/1,1/), (/NstateVar(ng),Nscale(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HdecayY',               &
     &                      HdecayY(1,:,:,ng),                          &
     &                      (/1,1/), (/NstateVar(ng),Nscale(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# else
      CALL netcdf_put_fvar (ng, model, ncname, 'Hdecay',                &
     &                      Hdecay(1,:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'Vdecay',                &
     &                      Vdecay(1,:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      IF (NSA.eq.2) THEN
# ifdef MULTI_SCALE_B
#  ifndef NONUNIFORM_SCALES
        CALL netcdf_put_fvar (ng, model, ncname, 'HdecayMX',            &
     &                        HdecayX(2,:,:,ng),                        &
     &                        (/1,1/), (/NstateVar(ng),Nscale(ng)/),    &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'HdecayMY',            &
     &                        HdecayY(2,:,:,ng),                        &
     &                        (/1,1/), (/NstateVar(ng),Nscale(ng)/),    &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# else
        CALL netcdf_put_fvar (ng, model, ncname, 'HdecayM',             &
     &                        Hdecay(2,:,ng), (/1/), (/NstateVar(ng)/), &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
        CALL netcdf_put_fvar (ng, model, ncname, 'VdecayM',             &
     &                        Vdecay(2,:,ng), (/1/), (/NstateVar(ng)/), &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
      END IF

# ifdef ADJUST_BOUNDARY
#  ifdef MULTI_SCALE_B
      CALL netcdf_put_fvar (ng, model, ncname, 'HdecayB',               &
     &                      HdecayB(:,:,:,ng),                          &
     &                      (/1,1,1/), (/NstateVar(ng),Nscale(ng),4/),  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  else
      CALL netcdf_put_fvar (ng, model, ncname, 'HdecayB',               &
     &                      HdecayB(:,:,ng),                            &
     &                      (/1,1/), (/NstateVar(ng),4/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'VdecayB',               &
     &                      VdecayB(:,:,ng),                            &
     &                      (/1,1/), (/NstateVar(ng),4/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef RPM_RELAXATION
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_M2diff',             &
     &                      tl_M2diff(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'tl_M3diff',             &
     &                      tl_M3diff(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'tl_Tdiff',              &
     &                      tl_Tdiff(:,ng), (/1/), (/NT(ng)/),          &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'LNM_flag',              &
     &                      LNM_flag, (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'LNM_depth',             &
     &                      LNM_depth(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'dTdz_min',              &
     &                      dTdz_min(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'ml_depth',              &
     &                      ml_depth(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef STD_MODEL
!
!  Background error covariance standard deviation modeling parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Sigma_max',             &
     &                      Sigma_max(:,ng), (/1/), (/NstateVar(ng)/),  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Sigma_ml',              &
     &                      Sigma_ml(:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Sigma_do',              &
     &                      Sigma_do(:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Sigma_dz',              &
     &                      Sigma_dz(:,ng), (/1/), (/NstateVar(ng)/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifndef COMPUTE_MLD
      CALL netcdf_put_fvar (ng, model, ncname, 'mld_uniform',           &
     &                      mld_uniform(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
#endif

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint sensitivity parameters.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'Lzeta',                 &
     &                      SCALARS(ng)%Lstate(isFsur),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Lubar',                 &
     &                      SCALARS(ng)%Lstate(isUbar),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Lvbar',                 &
     &                      SCALARS(ng)%Lstate(isVbar),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      CALL netcdf_put_lvar (ng, model, ncname, 'Luvel',                 &
     &                      SCALARS(ng)%Lstate(isUvel),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Lvvel',                 &
     &                      SCALARS(ng)%Lstate(isVvel),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Ltracer',               &
     &                      SCALARS(ng)%Lstate(isTvar(:)),              &
     &                      (/1/), (/NT(ng)/),                          &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'KstrS',                 &
     &                      KstrS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'KendS',                 &
     &                      KendS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#if defined FORCING_SV || defined SO_SEMI || defined STOCHASTIC_OPT
!
!  Singular Forcing Vectors or Stochastic Optimals state switches.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'Fzeta',                 &
     &                      SCALARS(ng)%Fstate(isFsur),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifndef SOLVE3D
      CALL netcdf_put_lvar (ng, model, ncname, 'Fubar',                 &
     &                      SCALARS(ng)%Fstate(isUbar),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Fvbar',                 &
     &                      SCALARS(ng)%Fstate(isVbar),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# else

      CALL netcdf_put_lvar (ng, model, ncname, 'Fuvel',                 &
     &                      SCALARS(ng)%Fstate(isUvel),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Fvvel',                 &
     &                      SCALARS(ng)%Fstate(isVvel),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Ftracer',               &
     &                      SCALARS(ng)%Fstate(isTvar(:)),              &
     &                      (/1/), (/NT(ng)/),                          &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      CALL netcdf_put_lvar (ng, model, ncname, 'Fsustr',                &
     &                      SCALARS(ng)%Fstate(isUstr),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_lvar (ng, model, ncname, 'Fsvstr',                &
     &                      SCALARS(ng)%Fstate(isVstr),                 &
     &                      (/0/), (/0/),                               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      CALL netcdf_put_lvar (ng, model, ncname, 'Fstflx',                &
     &                      SCALARS(ng)%Fstate(isTsur(:)),              &
     &                      (/1/), (/NT(ng)/),                          &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#ifdef SO_SEMI
# ifndef SO_SEMI_WHITE
      CALL netcdf_put_fvar (ng, model, ncname, 'SO_decay',              &
     &                      SO_decay(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      CALL netcdf_put_fvar (ng, model, ncname, 'SO_trace',              &
     &                      TRnorm(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_zeta',           &
     &                      SO_sdev(isFsur,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifndef SOLVE3D
      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_ubar',           &
     &                      SO_sdev(isUbar,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_vbar',           &
     &                      SO_sdev(isUbar,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# else

      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_uvel',           &
     &                      SO_sdev(isUvel,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_vvel',           &
     &                      SO_sdev(isVvel,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      DO itrc=1,NT(ng)
        nudg(itrc)=SO_sdev(isTvar(itrc),ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_tracer',         &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_sustr',          &
     &                      SO_sdev(isUstr,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN


      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_svstr',          &
     &                      SO_sdev(isVstr,ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      DO itrc=1,NT(ng)
        nudg(itrc)=SO_sdev(isTsur(itrc),ng)
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'SOsdev_stflx',          &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

#if defined BIOLOGY && defined SOLVE3D
# if defined BIO_FENNEL
#  include <fennel_wrt.h>
# elif defined ECB
#  include <ecb_wrt.h>
# elif defined ESTUARYBGC
#  include <estuarybgc_wrt.h>
# elif defined ECOSIM
#  include <ecosim_wrt.h>
# elif defined HYPOXIA_SRM
#  include <hypoxia_srm_wrt.h>
# elif defined NEMURO
#  include <nemuro_wrt.h>
# elif defined BIO_UMAINE
#  include <umaine_wrt.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_wrt.h>
# elif defined NPZD_IRON
#  include <npzd_iron_wrt.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_wrt.h>
# elif defined RED_TIDE
#  include <red_tide_wrt.h>
# endif
#endif

#if defined FLOATS && defined FLOAT_BIOLOGY
# if defined FLOAT_OYSTER
#  include <oyster_floats_wrt.h>
# endif
#endif

#ifdef SEDIMENT
# include <sediment_wrt.h>
#endif

#ifdef INWAVE_MODEL
# include <inwave_wrt.h>
#endif

!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'spherical',             &
     &                      spherical, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Domain Length.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'xl',                    &
     &                      xl(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'el',                    &
     &                      el(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'Vtransform',            &
     &                      Vtransform(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'Vstretching',           &
     &                      Vstretching(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'theta_s',               &
     &                      theta_s(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'theta_b',               &
     &                      theta_b(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Tcline',                &
     &                      Tcline(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'hc',                    &
     &                      hc(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids. The value
!  is arbitrary but is set to unity so it can be used as logical during
!  post-processing.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'grid',                  &
     &                      (/1/), (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL netcdf_put_fvar (ng, model, ncname, 's_rho',                 &
     &                      SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 's_w',                   &
     &                      SCALARS(ng)%sc_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_r',                  &
     &                      SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_w',                  &
     &                      SCALARS(ng)%Cs_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Depths of horizontal slices.
!
      IF (Nslice.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'z_slice',             &
     &                        Zslice, (/1/), (/Nslice/),                &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'user',                &
     &                        user, (/1/), (/Nuser/),                   &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#ifdef STATIONS
!
!  Stations positions.
!
      IF (ncid.eq.STA(ng)%ncid) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'Ipos',                &
     &                        SCALARS(ng)%SposX(:), (/1/),              &
     &                        (/Nstation(ng)/), ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, model, ncname, 'Jpos',                &
     &                        SCALARS(ng)%SposY(:), (/1/),              &
     &                        (/Nstation(ng)/), ncid = ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Write out grid tiled variables.
!-----------------------------------------------------------------------
!
#ifdef NO_WRITE_GRID
      GRID_VARS : IF (ncid.eq.STA(ng)%ncid) THEN
#else
      GRID_VARS : IF (ncid.ne.FLT(ng)%ncid) THEN
#endif
#if !(defined SED_MORPH && defined SEDIMENT)
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idtopo)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idtopo,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtopo)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtopo)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idtopo)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idtopo,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
#  endif
     &                           EXTRACT(ng) % h,                       &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtopo)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtopo)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# endif
# ifdef STATIONS
          ELSE IF (ncid.eq.STA(ng)%ncid) THEN
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%h,                      &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            CALL netcdf_put_fvar (ng, model, ncname, Vname(1,idtopo),   &
     &                            wrk, (/1/), (/Nstation(ng)/),         &
     &                            ncid = ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
          END IF
        END IF
#endif
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idfcor)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idfcor,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#endif
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idfcor)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfcor)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idfcor)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idfcor,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           EXTRACT(ng) % f,                       &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idfcor)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfcor)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#endif
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idpmdx)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idpmdx,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#endif
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpmdx)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpmdx)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idpmdx)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idpmdx,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           EXTRACT(ng) % pm,                      &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpmdx)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpmdx)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#endif
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idpndy)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idpndy,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#endif
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpndy)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpndy)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idpndy)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idpndy,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
# endif
     &                           EXTRACT(ng) % pn,                      &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpndy)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpndy)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#endif
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % rmask,                 &
# endif
     &                             EXTRACT(ng) % lonr,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
#ifdef STATIONS
            ELSE IF (ncid.eq.STA(ng)%ncid) THEN
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%lonr,                 &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, Vname(1,idLonR), &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % rmask,                 &
# endif
     &                             EXTRACT(ng) % latr,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
#ifdef STATIONS
            ELSE IF (ncid.eq.STA(ng)%ncid) THEN
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%latr,                 &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, Vname(1,idLatR), &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#endif
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % rmask,                 &
# endif
     &                             EXTRACT(ng) % xr,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
#ifdef STATIONS
            ELSE IF (ncid.eq.STA(ng)%ncid) THEN
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%xr,                   &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, Vname(1,idXgrR), &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % rmask,                    &
#endif
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % rmask,                 &
# endif
     &                             EXTRACT(ng) % yr,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
#ifdef STATIONS
            ELSE IF (ncid.eq.STA(ng)%ncid) THEN
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%yr,                   &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL netcdf_put_fvar (ng, model, ncname, Vname(1,idYgrR), &
     &                              wrk, (/1/), (/Nstation(ng)/),       &
     &                              ncid = ncid)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % umask,                 &
# endif
     &                             EXTRACT(ng) % lonu,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % umask,                 &
# endif
     &                             EXTRACT(ng) % latu,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % umask,                 &
# endif
     &                             EXTRACT(ng) % xu,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % umask,                    &
#endif
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % umask,                 &
# endif
     &                             EXTRACT(ng) % yu,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % vmask,                 &
# endif
     &                             EXTRACT(ng) % lonv,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % vmask,                 &
# endif
     &                             EXTRACT(ng) % latv,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % vmask,                 &
# endif
     &                             EXTRACT(ng) % xv,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % vmask,                    &
#endif
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % vmask,                 &
# endif
     &                             EXTRACT(ng) % yv,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % pmask,                    &
#endif
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % pmask,                 &
# endif
     &                             EXTRACT(ng) % lonp,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % pmask,                    &
#endif
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % pmask,                 &
# endif
     &                             EXTRACT(ng) % latp,                  &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % pmask,                    &
#endif
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % pmask,                 &
# endif
     &                             EXTRACT(ng) % xp,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
#ifdef MASKING
     &                             GRID(ng) % pmask,                    &
#endif
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#ifdef GRID_EXTRACT
            ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             EXTRACT(ng) % pmask,                 &
# endif
     &                             EXTRACT(ng) % yp,                    &
     &                             SetFillVal = .FALSE.,                &
     &                             ExtractField = -1)
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
#endif
            END IF
          END IF
        END IF

#ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idangR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idangR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % angler,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idangR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idangR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idangR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idangR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
#  endif
     &                           EXTRACT(ng) % angler,                  &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idangR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idangR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# endif
# ifdef STATIONS
          ELSE IF (ncid.eq.STA(ng)%ncid) THEN
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%angler,                 &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            CALL netcdf_put_fvar (ng, model, ncname, Vname(1,idangR),   &
     &                            wrk, (/1/), (/Nstation(ng)/),         &
     &                            ncid = ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
          END IF
        END IF
#endif

#ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           EXTRACT(ng) % rmask,                   &
     &                           EXTRACT(ng) % rmask,                   &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# endif
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskU)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskU,               &
     &                           varid, 0, u2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskU)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskU,               &
     &                           varid, 0, u2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           EXTRACT(ng) % umask,                   &
     &                           EXTRACT(ng) % umask,                   &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# endif
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskV)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskV,               &
     &                           varid, 0, v2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskV)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskV,               &
     &                           varid, 0, v2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           EXTRACT(ng) % vmask,                   &
     &                           EXTRACT(ng) % vmask,                   &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# endif
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskP)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskP,               &
     &                           varid, 0, p2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskP)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskP)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskP,               &
     &                           varid, 0, p2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           EXTRACT(ng) % pmask,                   &
     &                           EXTRACT(ng) % pmask,                   &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskP)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
# endif
          END IF
        END IF
#endif

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint sensitivity spatial scope mask at RHO-, U-, and V-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idscoR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idscoR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % Rscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
# ifdef GRID_EXTRACT
        ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idscoR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idscoR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
#  endif
     &                           EXTRACT(ng) % Rscope,                  &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
# endif
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idscoU)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idscoU,               &
     &                           varid, 0, u2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % umask,                      &
# endif
     &                           GRID(ng) % Uscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
# ifdef GRID_EXTRACT
        ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idscoU)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idscoU,               &
     &                           varid, 0, u2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           EXTRACT(ng) % umask,                   &
#  endif
     &                           EXTRACT(ng) % Uscope,                  &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
# endif
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idscoV)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idscoV,               &
     &                           varid, 0, v2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % vmask,                      &
# endif
     &                           GRID(ng) % Vscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
# ifdef GRID_EXTRACT
        ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idscoV)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idscoV,               &
     &                           varid, 0, v2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           EXTRACT(ng) % vmask,                   &
#  endif
     &                           EXTRACT(ng) % Vscope,                  &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
# endif
        END IF
#endif
#ifdef UV_DRAG_GRID
!
!  Spatially bottom friction parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
# if defined UV_LOGDRAG || defined BBL_MODEL
            IF (find_string(var_name, n_var, TRIM(Vname(1,idZoBL)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idZoBL,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#  endif
     &                           GRID(ng) % ZoBot,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idZoBL)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idZoBL)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idZoBL)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idZoBL,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
#   endif
     &                           EXTRACT(ng) % ZoBot,                   &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idZoBL)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idZoBL)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  endif
# endif
          END IF
# ifdef UV_LDRAG
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idragL)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idragL,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#  endif
     &                           GRID(ng) % rdrag,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragL)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragL)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idragL)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idragL,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
#   endif
     &                           EXTRACT(ng) % rdrag,                   &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragL)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragL)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  endif
          END IF
# endif
# ifdef UV_QDRAG
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idragQ)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idragQ,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#  endif
     &                           GRID(ng) % rdrag2,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragQ)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragQ)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  ifdef GRID_EXTRACT
          ELSE IF (ncid.eq.XTR(ng)%ncid) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idragQ)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idragQ,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           EXTRACT(ng) % rmask,                   &
#   endif
     &                           EXTRACT(ng) % rdrag2,                  &
     &                           SetFillVal = .FALSE.,                  &
     &                           ExtractField = -1)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragQ)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragQ)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  endif
          END IF
# endif
        END IF
#endif
      END IF GRID_VARS

#if !(defined OUT_NETCDF4 && defined DELAYED_SYNC_NF90)
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
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
  10  FORMAT (/,' WRT_INFO_NF90 - error while writing variable: ',a,/,  &
     &        17x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO_NF90 - error while inquiring ID for',        &
     &        ' variable: ',a,/,17x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO_NF90 - unable to synchronize to disk file:', &
     &        /,17x,a)
!
      RETURN
      END SUBROUTINE wrt_info_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_info_pio (ng, model, pioFile, ncname)
!***********************************************************************
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file using the standard NetCDF-3 or NetCDF-4 library.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     pioFile      PIO file descriptor structure, TYPE(File_desc_t)    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     ncname       PIO filename (string)                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status
      integer :: ival
      integer :: FileH, MY_FOUT

      integer :: ifield = 0
!
      real(dp) :: scale
# ifdef SOLVE3D
#  ifdef TS_DIF4
      real(r8), dimension(NT(ng)) :: diff
#  endif
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
# endif
# ifdef STATIONS
      real(r8), dimension(Nstation(ng)) :: Zpos, wrk
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_info_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
      TYPE (My_VarDesc) :: pioVar
!
      SourceFile=MyFile
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
!  Get NetCDF file handler from descriptor.
!
      FileH=ABS(pioFile%fh)
!
!  Time stepping parameters.
!
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntimes',            &
     &                          ntimes(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndtfast',           &
     &                          ndtfast(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dt',                &
     &                          dt(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dtfast',            &
     &                          dtfast(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dstart',            &
     &                          dstart, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef RBL4DVAR_FCT_SENSITIVITY
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntimes_ana',        &
     &                          ntimes_ana(ng), (/0/), (/0/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntimes_fct',        &
     &                          ntimes_fct(ng), (/0/), (/0/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined HDF5 && defined DEFLATE
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'shuffle',           &
     &                          shuffle, (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'deflate',           &
     &                          deflate, (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'deflate_level',     &
     &                          deflate_level, (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nHIS',              &
     &                          nHIS(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefHIS',           &
     &                          ndefHIS(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef GRID_EXTRACT
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nXTR',              &
     &                          nXTR(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefXTR',           &
     &                          ndefXTR(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ExtractFlag',       &
     &                          ExtractFlag(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nRST',              &
     &                          nRST(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntsAVG',            &
     &                          ntsAVG(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nAVG',              &
     &                          nAVG(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefAVG',           &
     &                          ndefAVG(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJOINT
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nADJ',              &
     &                          nADJ(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefADJ',           &
     &                          ndefADJ(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef TANGENT
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nTLM',              &
     &                          nTLM(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefTLM',           &
     &                          ndefTLM(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nOBC',              &
     &                          nOBC(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nSFF',              &
     &                          nSFF(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef PROPAGATOR
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LmultiGST',         &
     &                          LmultiGST, (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LrstGST',           &
     &                          LrstGST, (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'MaxIterGST',        &
     &                          MaxIterGST, (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nGST',              &
     &                          nGST, (/0/), (/0/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'NEV',               &
     &                          NEV, (/0/), (/0/),                      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'NCV',               &
     &                          NCV, (/0/), (/0/),                      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Ritz_tol',          &
     &                          Ritz_tol, (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef DIAGNOSTICS
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntsDIA',            &
     &                          ntsDIA(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nDIA',              &
     &                          nDIA(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefDIA',           &
     &                          ndefDIA(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef STATIONS
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nSTA',              &
     &                          nSTA(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef FOUR_DVAR
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Nouter',            &
     &                          Nouter, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Ninner',            &
     &                          Ninner, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Falpha',            &
     &                          Falpha, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Fbeta',             &
     &                          Fbeta, (/0/), (/0/),                    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Fgamma',            &
     &                          Fgamma, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Horizontal mixing coefficients.
!
# if defined SOLVE3D && defined TS_DIF2
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'nl_tnu2',           &
     &                          nl_tnu2(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJOINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ad_tnu2',           &
     &                          ad_tnu2(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_tnu2',           &
     &                          tl_tnu2(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

# endif

# if defined SOLVE3D && defined TS_DIF4
      DO itrc=1,NT(ng)
        diff(itrc)=nl_tnu4(itrc,ng)*nl_tnu4(itrc,ng)
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'nl_tnu4',           &
     &                          diff, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJOINT
      DO itrc=1,NT(ng)
        diff(itrc)=ad_tnu4(itrc,ng)*ad_tnu4(itrc,ng)
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ad_tnu4',           &
     &                          diff, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      DO itrc=1,NT(ng)
        diff(itrc)=tl_tnu4(itrc,ng)*tl_tnu4(itrc,ng)
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_tnu4',           &
     &                          diff, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef UV_VIS2
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'nl_visc2',          &
     &                          nl_visc2(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJOINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ad_visc2',          &
     &                          ad_visc2(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_visc2',          &
     &                          tl_visc2(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef UV_VIS4
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'nl_visc4',          &
     &                          nl_visc4(ng)**2, (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJOINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ad_visc4',          &
     &                          ad_visc4(ng)**2, (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_visc4',          &
     &                          tl_visc4(ng)**2, (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
#  ifdef TKE_DIF2
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tkenu2',            &
     &                          tkenu2(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef TKE_DIF4
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tkenu4',            &
     &                          tkenu4(ng)**2, (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined UV_VIS2 || defined UV_VIS4
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LuvSponge',         &
     &                          LuvSponge(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if (defined TS_DIF2 || defined TS_DIF4) && defined SOLVE3D
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LtracerSponge',     &
     &                          LtracerSponge(:,ng), (/1/), (/NT(ng)/), &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akt_bak',           &
     &                          Akt_bak(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akv_bak',           &
     &                          Akv_bak(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined MY25_MIXING || defined GLS_MIXING
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akk_bak',           &
     &                          Akk_bak(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akp_bak',           &
     &                          Akp_bak(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef FORWARD_MIXING
!
!  Basic state vertical mixing scale used in adjoint-based applications.
!
#   ifdef ADJOINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ad_Akt_fac',        &
     &                          ad_Akt_fac(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   if defined TANGENT || defined TL_IOMS
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_Akt_fac',        &
     &                          tl_Akt_fac(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef ADJOINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ad_Akv_fac',        &
     &                          ad_Akv_fac(ng), (/0/), (/0/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   if defined TANGENT || defined TL_IOMS
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_Akv_fac',        &
     &                          tl_Akv_fac(ng), (/0/), (/0/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
# endif
!
!  Drag coefficients.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rdrg',              &
     &                          rdrg(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rdrg2',             &
     &                          rdrg2(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Zob',               &
     &                          Zob(ng), (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Zos',               &
     &                          Zos(ng), (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_p',             &
     &                          gls_p(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_m',             &
     &                          gls_m(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_n',             &
     &                          gls_n(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_cmu0',          &
     &                          gls_cmu0(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c1',            &
     &                          gls_c1(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c2',            &
     &                          gls_c2(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c3m',           &
     &                          gls_c3m(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c3p',           &
     &                          gls_c3p(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_sigk',          &
     &                          gls_sigk(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_sigp',          &
     &                          gls_sigp(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_Kmin',          &
     &                          gls_Kmin(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_Pmin',          &
     &                          gls_Pmin(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Charnok_alpha',     &
     &                          charnok_alpha(ng), (/0/), (/0/),        &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Zos_hsig_alpha',    &
     &                          zos_hsig_alpha(ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'sz_alpha',          &
     &                          sz_alpha(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'CrgBan_cw',         &
     &                          crgban_cw(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef WEC
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'wec_alpha',         &
     &                          wec_alpha(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Nudging inverse time scales used in various tasks.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Znudg',             &
     &                          Znudg(ng)/sec2day, (/0/), (/0/),        &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'M2nudg',            &
     &                          M2nudg(ng)/sec2day, (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'M3nudg',            &
     &                          M3nudg(ng)/sec2day, (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tnudg',             &
     &                          nudg, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifndef DEBUGGING
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'FSobc_in',        &
     &                            FSobc_in(ng,:), (/1/), (/4/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'FSobc_out',       &
     &                            FSobc_out(ng,:), (/1/), (/4/),        &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M2obc_in',        &
     &                            M2obc_in(ng,:), (/1/), (/4/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M2obc_out',       &
     &                            M2obc_out(ng,:), (/1/), (/4/),        &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tobc_in',         &
     &                            Tobc, (/1,1/), (/NT(ng),4/),          &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tobc_out',        &
     &                            Tobc, (/1,1/), (/NT(ng),4/),          &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M3obc_in',        &
     &                            M3obc_in(ng,:), (/1/), (/4/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M3obc_out',       &
     &                            M3obc_out(ng,:), (/1/), (/4/),        &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
      END IF
# endif
!
!  Equation of State parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rho0',              &
     &                          rho0, (/0/), (/0/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined SOLVE3D && defined PROPAGATOR
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'bvf_bak',           &
     &                          bvf_bak, (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SOLVE3D     && \
   (!defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR)
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'R0',                &
     &                          R0(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tcoef',             &
     &                          Tcoef(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Scoef',             &
     &                          Scoef(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef SOLVE3D
#  ifdef BODYFORCE
!
!  Body force parameters.
!
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'levsfrc',           &
     &                          levsfrc(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'levbfrc',           &
     &                          levbfrc(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
!
!  Slipperiness parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gamma2',            &
     &                          gamma2(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LuvSrc',            &
     &                          LuvSrc(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LwSrc',             &
     &                          LwSrc(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
!
!  Logical switches to activate tracer point Sources/Sinks.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',        &
     &                          LtracerSrc(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Logical switches to process climatology fields.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LsshCLM',           &
     &                          LsshCLM(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',            &
     &                          Lm2CLM(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lm3CLM',            &
     &                          Lm3CLM(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LtracerCLM',        &
     &                          LtracerCLM(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Logical switches for nudging climatology fields.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',       &
     &                          LnudgeM2CLM(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LnudgeM3CLM',       &
     &                          LnudgeM3CLM(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LnudgeTCLM',        &
     &                          LnudgeTCLM(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef FOUR_DVAR
!
!  4DVAR assimilation parameters.
!
#  ifdef ADJUST_STFLUX
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lstflux',           &
     &                          Lstflux(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef ADJUST_BOUNDARY
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lobc',              &
     &                          Lobc(:,:,ng),                           &
     &                          (/1,1/), (/4,NstateVar(ng)/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifndef I4DVAR_ANA_SENSITIVITY
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LhessianEV',        &
     &                          LhessianEV, (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef WEAK_CONSTRAINT
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LhotStart',         &
     &                          LhotStart, (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lprecond',          &
     &                          Lprecond, (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lritz',             &
     &                          Lritz, (/0/), (/0/),                    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef WEAK_CONSTRAINT
      IF (Lprecond.and.(NritzEV.gt.0)) THEN
        CALL pio_netcdf_put_ivar (ng, model, ncname, 'NritzEV',         &
     &                            NritzEV, (/0/), (/0/),                &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#  endif
#  if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
        CALL pio_netcdf_put_ivar (ng, model, ncname, 'NpostI',          &
     &                            NpostI, (/0/), (/0/),                 &
     &                            pioFile = pioFile)
#  endif
#  if defined ARRAY_MODES              || \
     defined I4DVAR_ANA_SENSITIVITY   || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined R4DVAR_ANA_SENSITIVITY
        CALL pio_netcdf_put_ivar (ng, model, ncname, 'Nimpact',         &
     &                            Nimpact, (/0/), (/0/),                &
     &                            pioFile = pioFile)
#  endif
#  ifndef I4DVAR_ANA_SENSITIVITY
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'GradErr',           &
     &                          GradErr, (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HevecErr',          &
     &                          HevecErr, (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Nmethod',           &
     &                          Nmethod(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Rscheme',           &
     &                          Rscheme(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Nrandom',           &
     &                          Nrandom, (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Hgamma',            &
     &                          Hgamma(1), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef WEAK_CONSTRAINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HgammaM',           &
     &                          Hgamma(2), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ADJUST_BOUNDARY
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HgammaB',           &
     &                          Hgamma(3), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ADJUST_STFLUX
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HgammaF',           &
     &                          Hgamma(4), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Vgamma',            &
     &                          Vgamma(1), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef WEAK_CONSTRAINT
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'VgammaM',           &
     &                          Vgamma(2), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef ADJUST_BOUNDARY
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'VgammaB',           &
     &                          Vgamma(3), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef MULTI_SCALE_B
#   ifndef NONUNIFORM_SCALES
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayX',           &
     &                          HdecayX(1,:,:,ng),                      &
     &                          (/1,1/), (/NstateVar(ng),Nscale(ng)/),  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayY',           &
     &                          HdecayY(1,:,:,ng),                      &
     &                          (/1,1/), (/NstateVar(ng),Nscale(ng)/),  &
     &                          pioFile = pioFile)
#   endif
#  else
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Hdecay',            &
     &                          Hdecay(1,:,ng),                         &
     &                          (/1/), (/NstateVar(ng)/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Vdecay',            &
     &                          Vdecay(1,:,ng),                         &
     &                          (/1/), (/NstateVar(ng)/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

      IF (NSA.eq.2) THEN
#  ifdef MULTI_SCALE_B
#   ifndef NONUNIFORM_SCALES
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayMX',        &
     &                            HdecayX(2,:,:,ng),                    &
     &                            (/1,1/),                              &
     &                            (/NstateVar(ng),Nscale(ng)/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayMY',        &
     &                            HdecayY(2,:,:,ng),                    &
     &                            (/1,1/),                              &
     &                            (/NstateVar(ng),Nscale(ng)/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  else
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayM',         &
     &                            Hdecay(2,:,ng),                       &
     &                            (/1/), (/NstateVar(ng)/),             &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'VdecayM',         &
     &                            Vdecay(2,:,ng),                       &
     &                            (/1/), (/NstateVar(ng)/),             &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
#   ifdef MULTI_SCALE_B
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayB',           &
     &                          HdecayB(:,:,:,ng),                      &
     &                          (/1,1,1/),                              &
     &                          (/NstateVar(ng),Nscale(ng),4/),         &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   else
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'HdecayB',           &
     &                          HdecayB(:,:,ng),                        &
     &                          (/1,1/), (/NstateVar(ng),4/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

#   ifdef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'VdecayB',           &
     &                          VdecayB(:,:,ng),                        &
     &                          (/1,1/), (/NstateVar(ng),4/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef RPM_RELAXATION
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_M2diff',         &
     &                          tl_M2diff(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_M3diff',         &
     &                          tl_M3diff(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'tl_Tdiff',          &
     &                          tl_Tdiff(:,ng), (/1/), (/NT(ng)/),      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef BALANCE_OPERATOR
#   ifdef ZETA_ELLIPTIC
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Nbico',             &
     &                          Nbico(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
#   endif

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lbalance',          &
     &                          balance, (/1/), (/NstateVar(ng)/),      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'LNM_flag',          &
     &                          LNM_flag, (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'LNM_depth',         &
     &                          LNM_depth(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dTdz_min',          &
     &                          dTdz_min(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'ml_depth',          &
     &                          ml_depth(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef STD_MODEL
!
!  Background error covariance standard deviation modeling parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Sigma_max',         &
     &                          Sigma_max(:,ng),                        &
     &                          (/1/), (/NstateVar(ng)/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Sigma_ml',          &
     &                          Sigma_ml(:,ng),                         &
     &                          (/1/), (/NstateVar(ng)/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Sigma_do',          &
     &                          Sigma_do(:,ng),                         &
     &                          (/1/), (/NstateVar(ng)/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Sigma_dz',          &
     &                          Sigma_dz(:,ng),                         &
     &                          (/1/), (/NstateVar(ng)/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifndef COMPUTE_MLD
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'mld_uniform',       &
     &                          mld_uniform(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
# endif

# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint sensitivity parameters.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lzeta',             &
     &                          SCALARS(ng)%Lstate(isFsur),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lubar',             &
     &                          SCALARS(ng)%Lstate(isUbar),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lvbar',             &
     &                          SCALARS(ng)%Lstate(isVbar),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Luvel',             &
     &                          SCALARS(ng)%Lstate(isUvel),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lvvel',             &
     &                          SCALARS(ng)%Lstate(isVvel),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Ltracer',           &
     &                          SCALARS(ng)%Lstate(isTvar(:)),          &
     &                          (/1/), (/NT(ng)/),                      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'KstrS',             &
     &                          KstrS(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'KendS',             &
     &                          KendS(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined FORCING_SV || defined SO_SEMI || defined STOCHASTIC_OPT
!
!  Singular Forcing Vectors or Stochastic Optimals state switches.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fzeta',             &
     &                          SCALARS(ng)%Fstate(isFsur),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifndef SOLVE3D
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fubar',             &
     &                          SCALARS(ng)%Fstate(isUbar),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fvbar',             &
     &                          SCALARS(ng)%Fstate(isVbar),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  else

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fuvel',             &
     &                          SCALARS(ng)%Fstate(isUvel),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fvvel',             &
     &                          SCALARS(ng)%Fstate(isVvel),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Ftracer',           &
     &                          SCALARS(ng)%Fstate(isTvar(:)),          &
     &                          (/1/), (/NT(ng)/),                      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fsustr',            &
     &                          SCALARS(ng)%Fstate(isUstr),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fsvstr',            &
     &                          SCALARS(ng)%Fstate(isVstr),             &
     &                          (/0/), (/0/),                           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Fstflx',            &
     &                          SCALARS(ng)%Fstate(isTsur(:)),          &
     &                          (/1/), (/NT(ng)/),                      &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef SO_SEMI
#  ifndef SO_SEMI_WHITE
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SO_decay',          &
     &                          SO_decay(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SO_trace',          &
     &                          TRnorm(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_zeta',       &
     &                          SO_sdev(isFsur,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifndef SOLVE3D
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_ubar',       &
     &                          SO_sdev(isUbar,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_vbar',       &
     &                          SO_sdev(isUbar,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  else

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_uvel',       &
     &                          SO_sdev(isUvel,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_vvel',       &
     &                          SO_sdev(isVvel,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      DO itrc=1,NT(ng)
        nudg(itrc)=SO_sdev(isTvar(itrc),ng)
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_tracer',     &
     &                          nudg, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_sustr',      &
     &                          SO_sdev(isUstr,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN


      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_svstr',      &
     &                          SO_sdev(isVstr,ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
      DO itrc=1,NT(ng)
        nudg(itrc)=SO_sdev(isTsur(itrc),ng)
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'SOsdev_stflx',      &
     &                          nudg, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined BIOLOGY && defined SOLVE3D
#  if defined BIO_FENNEL
#   include <fennel_wrt_pio.h>
#  elif defined ECB
#   include <ecb_wrt_pio.h>
#  elif defined ECOSIM
#   include <ecosim_wrt_pio.h>
#  elif defined HYPOXIA_SRM
#   include <hypoxia_srm_wrt_pio.h>
#  elif defined NEMURO
#   include <nemuro_wrt_pio.h>
#  elif defined NPZD_FRANKS
#   include <npzd_Franks_wrt_pio.h>
#  elif defined NPZD_IRON
#   include <npzd_iron_wrt_pio.h>
#  elif defined NPZD_POWELL
#   include <npzd_Powell_wrt_pio.h>
#  elif defined RED_TIDE
#   include <red_tide_wrt_pio.h>
#  endif
# endif

# if defined FLOATS && defined FLOAT_BIOLOGY
#  if defined FLOAT_OYSTER
#   include <oyster_floats_wrt_pio.h>
#  endif
# endif

# ifdef SEDIMENT
#  include <sediment_wrt_pio.h>
# endif
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'spherical',         &
     &                          spherical, (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Domain Length.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'xl',                &
     &                          xl(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'el',                &
     &                          el(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
!
!  S-coordinate parameters.
!
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Vtransform',        &
     &                          Vtransform(ng), (/0/), (/0/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Vstretching',       &
     &                          Vstretching(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'theta_s',           &
     &                          theta_s(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'theta_b',           &
     &                          theta_b(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tcline',            &
     &                          Tcline(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'hc',                &
     &                          hc(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids. The value
!  is arbitrary but is set to unity so it can be used as logical during
!  post-processing.
!
      ival=1
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'grid',              &
     &                          ival, (/0/), (/0/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 's_rho',             &
     &                          SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 's_w',               &
     &                          SCALARS(ng)%sc_w(0:),                   &
     &                          (/1/), (/N(ng)+1/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Cs_r',              &
     &                          SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Cs_w',              &
     &                          SCALARS(ng)%Cs_w(0:),                   &
     &                          (/1/), (/N(ng)+1/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Depths of horizontal slices.
!
      IF (Nslice.gt.0) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'z_slice',         &
     &                            Zslice, (/1/), (/Nslice/),            &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'user',            &
     &                            user, (/1/), (/Nuser/),               &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef STATIONS
!
!  Stations positions.
!
      IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'Ipos',            &
     &                            SCALARS(ng)%SposX(:), (/1/),          &
     &                            (/Nstation(ng)/), pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL pio_netcdf_put_fvar (ng, model, ncname, 'Jpos',            &
     &                            SCALARS(ng)%SposY(:), (/1/),          &
     &                            (/Nstation(ng)/), pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Write out grid tiled variables.
!-----------------------------------------------------------------------
!
# ifdef NO_WRITE_GRID
      GRID_VARS : IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
# else
      GRID_VARS : IF (FileH.ne.ABS(FLT(ng)%pioFile%fh)) THEN
# endif
# if !(defined SED_MORPH && defined SEDIMENT)
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idtopo)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idtopo,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#  endif
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtopo)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtopo)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  ifdef STATIONS
          ELSE
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%h,                      &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            CALL pio_netcdf_put_fvar (ng, model, ncname,                &
     &                                Vname(1,idtopo),                  &
     &                                wrk, (/1/), (/Nstation(ng)/),     &
     &                                pioFile = pioFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
          END IF
        END IF
# endif
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idfcor)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idfcor,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idfcor)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfcor)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idpmdx)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idpmdx,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpmdx)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpmdx)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idpndy)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idpndy,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpndy)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpndy)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % rmask,                    &
# endif
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
# ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%lonr,                 &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,idLonR),                &
     &                                  wrk, (/1/), (/Nstation(ng)/),   &
     &                                  pioFile = pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % rmask,                    &
# endif
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
# ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%latr,                 &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  Vname(1,idLatR),                &
     &                                  wrk, (/1/), (/Nstation(ng)/),   &
     &                                  pioFile = pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % rmask,                    &
# endif
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
# ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%xr,                     &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  TRIM(Vname(1,idXgrR)),          &
     &                                  wrk, (/1/), (/Nstation(ng)/),   &
     &                                  pioFile = pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
# endif
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % rmask,                    &
# endif
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
# ifdef STATIONS
            ELSE
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, GRID(ng)%yr,                   &
     &                            Nstation(ng), SCALARS(ng)%SposX,      &
     &                            SCALARS(ng)%SposY, wrk)
              CALL pio_netcdf_put_fvar (ng, model, ncname,              &
     &                                  TRIM(Vname(1,idYgrR)),          &
     &                                  wrk, (/1/), (/Nstation(ng)/),   &
     &                                  pioFile = pioFile)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  endif
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % umask,                    &
# endif
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % umask,                    &
# endif
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % umask,                    &
# endif
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % umask,                    &
# endif
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % vmask,                    &
# endif
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % vmask,                    &
# endif
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrV,              &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % vmask,                    &
# endif
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % vmask,                    &
# endif
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % pmask,                    &
# endif
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % pmask,                    &
# endif
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % pmask,                    &
# endif
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_psi', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
# ifdef MASKING
     &                             GRID(ng) % pmask,                    &
# endif
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF

# ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idangR)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idangR,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#  endif
     &                           GRID(ng) % angler,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idangR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idangR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  ifdef STATIONS
          ELSE
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%angler,                 &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            CALL pio_netcdf_put_fvar (ng, model, ncname, 'angle',       &
     &                                wrk, (/1/), (/Nstation(ng)/),     &
     &                                pioFile = pioFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
          END IF
        END IF
# endif

# ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskR)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskR,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskU)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=u2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskU,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskV)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=v2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskV,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskP)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=p2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskP,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskP)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
# endif

# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
!
!  Adjoint sensitivity spatial scope mask at RHO-, U-, and V-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idscoR)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idscoR,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#  endif
     &                           GRID(ng) % Rscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idscoU)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=u2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idscoU,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % umask,                      &
#  endif
     &                           GRID(ng) % Uscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idscoV)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=v2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idscoV,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % vmask,                      &
#  endif
     &                           GRID(ng) % Vscope,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idscoV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idscoV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
# endif
# ifdef UV_DRAG_GRID
!
!  Spatially bottom friction parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
#  if defined UV_LOGDRAG || defined BBL_MODEL
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              Vname(1,idZoBL), pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idZoBL,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#   endif
     &                           GRID(ng) % ZoBot,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idZoBL)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idZoBL)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  endif
#  ifdef UV_LDRAG
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              Vname(1,idragL), pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idragL,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#   endif
     &                           GRID(ng) % rdrag,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragL)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragL)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  endif
#  ifdef UV_QDRAG
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              Vname(1,idragQ), pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idragQ,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#   endif
     &                           GRID(ng) % rdrag2,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idragQ)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idragQ)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
#  endif
          END IF
        END IF
# endif
      END IF GRID_VARS

# ifndef DELAYED_SYNC_PIO
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, ncname, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
  10  FORMAT (/,' WRT_INFO_PIO - error while writing variable: ',a,/,   &
     &        16x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO_PIO - error while inquiring ID for',         &
     &        ' variable: ',a,/,16x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO_PIO - unable to synchronize to disk file:',  &
     &        /,16x,a)
!
      RETURN
      END SUBROUTINE wrt_info_pio
#endif
      END MODULE wrt_info_mod
