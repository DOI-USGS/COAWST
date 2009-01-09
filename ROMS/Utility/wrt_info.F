#include "cppdefs.h"
      FUNCTION wrt_info (ng, model, ncid, wrtThread, ncname)
!
!svn $Id: wrt_info.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number.                                     !
!     model    Calling model identifier.                               !
!     ncid     NetCDF file ID.                                         !
!     ncname   NetCDF file name.                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     wrt_info   Error flag (integer).                                 !
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
!
#ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
!
#endif
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: wrtThread

      integer, intent(in) :: ng, model, ncid

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.

      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, varid

      integer, dimension(2) :: start, total

      integer :: nf_fwrite2d
      integer :: wrt_info

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
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
      character (len=1), dimension(NT(ng)) :: char1_trc
#endif
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
      IF (wrtThread) THEN
!
!  Time stepping parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'ntimes', varid)
        wrt_info=nf90_put_var(ncid, varid, ntimes(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ntimes', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ndtfast', varid)
        wrt_info=nf90_put_var(ncid, varid, ndtfast(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ndtfast', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'dt', varid)
        wrt_info=nf90_put_var(ncid, varid, dt(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'dt', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'dtfast', varid)
        wrt_info=nf90_put_var(ncid, varid, dtfast(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'dtfast', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'dstart', varid)
        wrt_info=nf90_put_var(ncid, varid, dstart)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'dstart', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'nHIS', varid)
        wrt_info=nf90_put_var(ncid, varid, nHIS(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nHIS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ndefHIS', varid)
        wrt_info=nf90_put_var(ncid, varid, ndefHIS(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ndefHIS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'nRST', varid)
        wrt_info=nf90_put_var(ncid, varid, nRST(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nRST', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#ifdef AVERAGES
        wrt_info=nf90_inq_varid(ncid, 'ntsAVG', varid)
        wrt_info=nf90_put_var(ncid, varid, ntsAVG(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ntsAVG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'nAVG', varid)
        wrt_info=nf90_put_var(ncid, varid, nAVG(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nAVG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ndefAVG', varid)
        wrt_info=nf90_put_var(ncid, varid, ndefAVG(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ndefAVG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef ADJOINT
        wrt_info=nf90_inq_varid(ncid, 'nADJ', varid)
        wrt_info=nf90_put_var(ncid, varid, nADJ(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nADJ', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ndefADJ', varid)
        wrt_info=nf90_put_var(ncid, varid, ndefADJ(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ndefADJ', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef TANGENT
        wrt_info=nf90_inq_varid(ncid, 'nTLM', varid)
        wrt_info=nf90_put_var(ncid, varid, nTLM(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nTLM', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ndefTLM', varid)
        wrt_info=nf90_put_var(ncid, varid, ndefTLM(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ndefTLM', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        wrt_info=nf90_inq_varid(ncid, 'nSFF', varid)
        wrt_info=nf90_put_var(ncid, varid, nSFF(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nSFF', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef PROPAGATOR
        wrt_info=nf90_inq_varid(ncid, 'LrstGST', varid)
        WRITE (char1,'(l1)') LrstGST
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'LrstGST', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'MaxIterGST', varid)
        wrt_info=nf90_put_var(ncid, varid, MaxIterGST)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'MaxIterGST', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'nGST', varid)
        wrt_info=nf90_put_var(ncid, varid, nGST)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nGST', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Ritz_tol', varid)
        wrt_info=nf90_put_var(ncid, varid, Ritz_tol)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Ritz_tol', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef DIAGNOSTICS
        wrt_info=nf90_inq_varid(ncid, 'ntsDIA', varid)
        wrt_info=nf90_put_var(ncid, varid, ntsDIA(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ntsDIA', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'nDIA', varid)
        wrt_info=nf90_put_var(ncid, varid, nDIA(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nDIA', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ndefDIA', varid)
        wrt_info=nf90_put_var(ncid, varid, ndefDIA(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ndefDIA', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef STATIONS
        wrt_info=nf90_inq_varid(ncid, 'nSTA', varid)
        wrt_info=nf90_put_var(ncid, varid, nSTA(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'nSTA', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef FOUR_DVAR
        wrt_info=nf90_inq_varid(ncid, 'Nouter', varid)
        wrt_info=nf90_put_var(ncid, varid, Nouter)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Nouter', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Ninner', varid)
        wrt_info=nf90_put_var(ncid, varid, Ninner)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Ninner', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
        wrt_info=nf90_inq_varid(ncid, 'Falpha', varid)
        wrt_info=nf90_put_var(ncid, varid, Falpha)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Falpha', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Fbeta', varid)
        wrt_info=nf90_put_var(ncid, varid, Fbeta)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Fbeta', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Fgamma', varid)
        wrt_info=nf90_put_var(ncid, varid, Fgamma)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Fgamma', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
        wrt_info=nf90_inq_varid(ncid, 'tnu2', varid)
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, tnu2(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tnu2', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#if defined SOLVE3D && defined TS_DIF4
        wrt_info=nf90_inq_varid(ncid, 'tnu4', varid)
        DO itrc=1,NT(ng)
          diff(itrc)=tnu4(itrc,ng)*tnu4(itrc,ng)
        END DO
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, diff, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tnu4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef UV_VIS2
        wrt_info=nf90_inq_varid(ncid, 'visc2', varid)
        wrt_info=nf90_put_var(ncid, varid, visc2(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'visc2', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#ifdef UV_VIS4
        wrt_info=nf90_inq_varid(ncid, 'visc4', varid)
        wrt_info=nf90_put_var(ncid, varid, visc4(ng)**2)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'visc4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
        wrt_info=nf90_inq_varid(ncid, 'tkenu2', varid)
        wrt_info=nf90_put_var(ncid, varid, tkenu2(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tkenu2', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif

# ifdef TKE_DIF4
        wrt_info=nf90_inq_varid(ncid, 'tkenu4', varid)
        wrt_info=nf90_put_var(ncid, varid, tkenu4(ng)**2)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tkenu4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
#endif
#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
        wrt_info=nf90_inq_varid(ncid, 'Akt_bak', varid)
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, Akt_bak(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Akt_bak', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Akv_bak', varid)
        wrt_info=nf90_put_var(ncid, varid, Akv_bak(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Akv_bak', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# if defined MY25_MIXING || defined GLS_MIXING
        wrt_info=nf90_inq_varid(ncid, 'Akk_bak', varid)
        wrt_info=nf90_put_var(ncid, varid, Akk_bak(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Akk_bak', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Akp_bak', varid)
        wrt_info=nf90_put_var(ncid, varid, Akp_bak(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Akp_bak', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
#endif
!
!  Drag coefficients.
!
        wrt_info=nf90_inq_varid(ncid, 'rdrg', varid)
        wrt_info=nf90_put_var(ncid, varid, rdrg(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'rdrg', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'rdrg2', varid)
        wrt_info=nf90_put_var(ncid, varid, rdrg2(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'rdrg2', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'Zob', varid)
        wrt_info=nf90_put_var(ncid, varid, Zob(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Zob', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Zos', varid)
        wrt_info=nf90_put_var(ncid, varid, Zos(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Zos', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif
#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'gls_p', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_p(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_p', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_m', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_m(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_m', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_n', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_n(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_n', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_cmu0', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_cmu0(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_cmu0', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_c1', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_c1(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_c1', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_c2', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_c2(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_c2', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_c3m', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_c3m(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_c3m', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_c3p', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_c3p(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_c3p', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_sigk', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_sigk(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_sigk', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_sigp', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_sigp(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_sigp', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_Kmin', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_Kmin(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_Kmin', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'gls_Pmin', varid)
        wrt_info=nf90_put_var(ncid, varid, gls_Pmin(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gls_Pmin', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif
!
!  Nudging inverse time scales used in various tasks.
!
        wrt_info=nf90_inq_varid(ncid, 'Znudg', varid)
        wrt_info=nf90_put_var(ncid, varid, Znudg(ng)/sec2day)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Znudg', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'M2nudg', varid)
        wrt_info=nf90_put_var(ncid, varid, M2nudg(ng)/sec2day)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M2nudg', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'M3nudg', varid)
        wrt_info=nf90_put_var(ncid, varid, M3nudg(ng)/sec2day)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M3nudg', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Tnudg', varid)
        DO itrc=1,NT(ng)
          nudg(itrc)=Tnudg(itrc,ng)/sec2day
        END DO
        start(1)=1
        total(1)=NT(ng) 
        wrt_info=nf90_put_var(ncid, varid, nudg, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Tnudg', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif
#ifdef NUDGING
!
!  Nudging inverse time scales used in data assimilation.
!
        wrt_info=nf90_inq_varid(ncid, 'Znudass', varid)
        wrt_info=nf90_put_var(ncid, varid, Znudass(ng)/sec2day)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Znudass', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'M2nudass', varid)
        wrt_info=nf90_put_var(ncid, varid, M2nudass(ng)/sec2day)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M2nudass', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'M3nudass', varid)
        wrt_info=nf90_put_var(ncid, varid, M3nudass(ng)/sec2day)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M3nudass', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Tnudass', varid)
        DO itrc=1,NT(ng)
          nudg(itrc)=Tnudass(itrc,ng)/sec2day
        END DO
        start(1)=1
        total(1)=NT(ng) 
        wrt_info=nf90_put_var(ncid, varid, nudg, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Tnudass', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
#endif
#ifdef NUDGING_COFF
!
!  Open boundary nudging, inverse time scales.
!
        wrt_info=nf90_inq_varid(ncid, 'FSobc_in', varid)
        start(1)=1
        total(1)=4
        wrt_info=nf90_put_var(ncid, varid, FSobc_in(ng,:),              &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'FSobc_in', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'FSobc_out', varid)
        start(1)=1
        total(1)=4
        wrt_info=nf90_put_var(ncid, varid, FSobc_out(ng,:),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'FSobc_out', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'M2obc_in', varid)
        start(1)=1
        total(1)=4
        wrt_info=nf90_put_var(ncid, varid, M2obc_in(ng,:),              &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M2obc_in', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'M2obc_out', varid)
        start(1)=1
        total(1)=4
        wrt_info=nf90_put_var(ncid, varid, M2obc_out(ng,:),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M2obc_out', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# ifdef SOLVE3D
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        wrt_info=nf90_inq_varid(ncid, 'Tobc_in', varid)
        start(1)=1
        total(1)=NT(ng)
        start(2)=1
        total(2)=4
        wrt_info=nf90_put_var(ncid, varid, Tobc, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Tobc_in', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        wrt_info=nf90_inq_varid(ncid, 'Tobc_out', varid)
        start(1)=1
        total(1)=NT(ng)
        start(2)=1
        total(2)=4
        wrt_info=nf90_put_var(ncid, varid, Tobc, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Tobc_out', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'M3obc_in', varid)
        start(1)=1
        total(1)=4
        wrt_info=nf90_put_var(ncid, varid, M3obc_in(ng,:),              &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M3obc_in', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'M3obc_out', varid)
        start(1)=1
        total(1)=4
        wrt_info=nf90_put_var(ncid, varid, M3obc_out(ng,:),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'M3obc_out', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
#endif
!
!  Equation of State parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'rho0', varid)
        wrt_info=nf90_put_var(ncid, varid, rho0)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'rho0', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#if defined SOLVE3D && (defined FOUR_DVAR || defined PROPAGATOR)
        wrt_info=nf90_inq_varid(ncid, 'bvf_bak', varid)
        wrt_info=nf90_put_var(ncid, varid, bvf_bak)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'bvf_bak', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif

#if defined SOLVE3D     && \
   (!defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR)
        wrt_info=nf90_inq_varid(ncid, 'R0', varid)
        wrt_info=nf90_put_var(ncid, varid, R0(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'R0', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Tcoef', varid)
        wrt_info=nf90_put_var(ncid, varid, Tcoef(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Tcoef', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Scoef', varid)
        wrt_info=nf90_put_var(ncid, varid, Scoef(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Scoef', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif
#ifdef SOLVE3D
# ifdef BODYFORCE
!
!  Body force parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'levsfrc', varid)
        wrt_info=nf90_put_var(ncid, varid, levsfrc(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'levsfrc', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'levbfrc', varid)
        wrt_info=nf90_put_var(ncid, varid, levbfrc(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'levbfrc', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
#endif
!
!  Slipperiness parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'gamma2', varid)
        wrt_info=nf90_put_var(ncid, varid, gamma2(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'gamma2', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#ifdef FOUR_DVAR
!
!  Incremental 4DVAR parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'GradErr', varid)
        wrt_info=nf90_put_var(ncid, varid, Hgamma)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'GradErr', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HevecErr', varid)
        wrt_info=nf90_put_var(ncid, varid, Hgamma)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HevecErr', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Nmethod', varid)
        wrt_info=nf90_put_var(ncid, varid, Nmethod(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Nmethod', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Rscheme', varid)
        wrt_info=nf90_put_var(ncid, varid, Rscheme(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Rscheme', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Nrandom', varid)
        wrt_info=nf90_put_var(ncid, varid, Nrandom)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Nrandom', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Hgamma', varid)
        wrt_info=nf90_put_var(ncid, varid, Hgamma)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Hgamma', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'Vgamma', varid)
        wrt_info=nf90_put_var(ncid, varid, Vgamma)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Vgamma', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif

        wrt_info=nf90_inq_varid(ncid, 'Hdecay', varid)
        start(1)=1
        total(1)=NstateVar(ng)
        wrt_info=nf90_put_var(ncid, varid, Hdecay(:,ng),                &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Hdecay', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'Vdecay', varid)
        start(1)=1
        total(1)=NstateVar(ng)
        wrt_info=nf90_put_var(ncid, varid, Vdecay(:,ng),                &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Vdecay', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif

# ifdef RPM_RELAXATION
        wrt_info=nf90_inq_varid(ncid, 'tl_M2diff', varid)
        wrt_info=nf90_put_var(ncid, varid, tl_M2diff(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tl_M2diff', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#  ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'tl_M3diff', varid)
        wrt_info=nf90_put_var(ncid, varid, tl_M3diff(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tl_M3diff', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'tl_Tdiff', varid)
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, tl_Tdiff(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'tl_Tdiff', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#  endif
# endif
#endif

#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
!
!  Adjoint sensitivity parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'Lzeta', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isFsur)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Lstate(isFsur)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Lubar', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isUbar)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Lstate(isUbar)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Lvbar', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isVbar)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Lstate(isVbar)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'Luvel', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isUvel)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Lstate(isUvel)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Lvvel', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%Lstate(isVvel)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Lstate(isVvel)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Ltracer', varid)
        DO i=1,NT(ng)
          WRITE (char1_trc(i),'(l1)') SCALARS(ng)%Lstate(isTvar(i))
        END DO
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, char1_trc, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Lstate(isTvar)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'KstrS', varid)
        wrt_info=nf90_put_var(ncid, varid, KstrS(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'KstrS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'KendS', varid)
        wrt_info=nf90_put_var(ncid, varid, KendS(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'KendS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
# ifdef SO_SEMI
#  ifndef SO_SEMI_WHITE
        wrt_info=nf90_inq_varid(ncid, 'SO_decay', varid)
        wrt_info=nf90_put_var(ncid, varid, SO_decay(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SO_decay', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#  endif

        wrt_info=nf90_inq_varid(ncid, 'SO_trace', varid)
        wrt_info=nf90_put_var(ncid, varid, TRnorm(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SO_trace', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Lso_Ustr', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%SOstate(isUstr)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SOstate(isUstr)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Lso_Vstr', varid)
        WRITE (char1,'(l1)') SCALARS(ng)%SOstate(isVstr)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SOstate(isVstr)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#  ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'Lso_tracer', varid)
        DO i=1,NT(ng)
          WRITE (char1_trc(i),'(l1)') SCALARS(ng)%SOstate(isTsur(i))
        END DO
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, char1_trc, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SOstate(isTsur)', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#  endif

        wrt_info=nf90_inq_varid(ncid, 'SOsdev_Ustr', varid)
        wrt_info=nf90_put_var(ncid, varid, SO_sdev(isUstr,ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SOsdev_Ustr', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'SOsdev_Vstr', varid)
        wrt_info=nf90_put_var(ncid, varid, SO_sdev(isUstr,ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SOsdev_Vstr', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

#  ifdef SOLVE3D
        wrt_info=nf90_inq_varid(ncid, 'SOsdev_stflx', varid)
        DO itrc=1,NT(ng)
          nudg(itrc)=SO_sdev(isTsur(itrc),ng)
        END DO
        start(1)=1
        total(1)=NT(ng)
        wrt_info=nf90_put_var(ncid, varid, nudg, start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SOsdev_stflx', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#  endif
# endif
#endif

#if defined BIOLOGY && defined SOLVE3D
# ifdef BIO_FASHAM
!
!  Write out Fasham type biological model parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'BioIter', varid)
        wrt_info=nf90_put_var(ncid, varid, BioIter(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BioIter', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'AttSW', varid)
        wrt_info=nf90_put_var(ncid, varid, AttSW(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'AttSW', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'AttChl', varid)
        wrt_info=nf90_put_var(ncid, varid, AttChl(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'AttChl', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PARfrac', varid)
        wrt_info=nf90_put_var(ncid, varid, PARfrac(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PARfrac', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Vp0', varid)
        wrt_info=nf90_put_var(ncid, varid, Vp0(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Vp0', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'I_thNH4', varid)
        wrt_info=nf90_put_var(ncid, varid, I_thNH4(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'I_thNH4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'D_p5NH4', varid)
        wrt_info=nf90_put_var(ncid, varid, D_p5NH4(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'D_p5NH4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'NitriR', varid)
        wrt_info=nf90_put_var(ncid, varid, NitriR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'NitriR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_NO3', varid)
        wrt_info=nf90_put_var(ncid, varid, K_NO3(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_NO3', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_NH4', varid)
        wrt_info=nf90_put_var(ncid, varid, K_NH4(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_NH4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_Phy', varid)
        wrt_info=nf90_put_var(ncid, varid, K_Phy(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_Phy', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Chl2C_m', varid)
        wrt_info=nf90_put_var(ncid, varid, Chl2C_m(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Chl2C_m', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ChlMin', varid)
        wrt_info=nf90_put_var(ncid, varid, ChlMin(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ChlMin', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyCN', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyCN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyCN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyIP', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyIP(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyIP', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyIS', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyIS(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyIS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyMin', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyMin(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyMin', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyMR', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyMR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyMR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooAE_N', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooAE_N(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooAE_N', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooBM', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooBM(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooBM', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooCN', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooCN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooCN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooER', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooER(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooER', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooGR', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooGR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooGR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooMin', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooMin(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooMin', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooMR', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooMR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooMR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'LDeRRN', varid)
        wrt_info=nf90_put_var(ncid, varid, LDeRRN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'LDeRRN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'LDeRRC', varid)
        wrt_info=nf90_put_var(ncid, varid, LDeRRC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'LDeRRC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'CoagR', varid)
        wrt_info=nf90_put_var(ncid, varid, CoagR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'CoagR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'SDeRRN', varid)
        wrt_info=nf90_put_var(ncid, varid, SDeRRN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SDeRRN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'SDeRRC', varid)
        wrt_info=nf90_put_var(ncid, varid, SDeRRC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'SDeRRC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'wPhy', varid)
        wrt_info=nf90_put_var(ncid, varid, wPhy(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'wPhy', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'wLDet', varid)
        wrt_info=nf90_put_var(ncid, varid, wLDet(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'wLDet', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'wSDet', varid)
        wrt_info=nf90_put_var(ncid, varid, wSDet(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'wSDet', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'pCO2air', varid)
        wrt_info=nf90_put_var(ncid, varid, pCO2air(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'pCO2air', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# elif defined NPZD_FRANKS
!
!  Write out NPZD (Franks et al., 1986) biological model parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'BioIter', varid)
        wrt_info=nf90_put_var(ncid, varid, BioIter(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BioIter', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_ext', varid)
        wrt_info=nf90_put_var(ncid, varid, K_ext(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_ext', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_NO3', varid)
        wrt_info=nf90_put_var(ncid, varid, K_NO3(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_NO3', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_Phy', varid)
        wrt_info=nf90_put_var(ncid, varid, K_Phy(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_Phy', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Vm_NO3', varid)
        wrt_info=nf90_put_var(ncid, varid, Vm_NO3(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Vm_NO3', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyMR', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyMR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyMR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooGR', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooGR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooGR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooMR', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooMR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooMR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooMD', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooMD(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooMD', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooGA', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooGA(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooGA', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooEC', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooEC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooEC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'DetRR', varid)
        wrt_info=nf90_put_var(ncid, varid, DetRR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'DetRR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'wDet', varid)
        wrt_info=nf90_put_var(ncid, varid, wDet(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'wDet', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# elif defined NPZD_POWELL
!
!  Write out NPZD (Powell et al., 2006) biological model parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'BioIter', varid)
        wrt_info=nf90_put_var(ncid, varid, BioIter(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BioIter', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'AttSW', varid)
        wrt_info=nf90_put_var(ncid, varid, AttSW(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'AttSW', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'AttPhy', varid)
        wrt_info=nf90_put_var(ncid, varid, AttPhy(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'AttPhy', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyIS', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyIS(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyIS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Vm_NO3', varid)
        wrt_info=nf90_put_var(ncid, varid, Vm_NO3(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Vm_NO3', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyMRD', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyMRD(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyMRD', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyMRN', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyMRN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyMRN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'K_NO3', varid)
        wrt_info=nf90_put_var(ncid, varid, K_NO3(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'K_NO3', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Ivlev', varid)
        wrt_info=nf90_put_var(ncid, varid, Ivlev(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Ivlev', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooGR', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooGR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooGR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooEED', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooEED(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooEED', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooEEN', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooEEN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooEEN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooMRD', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooMRD(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooMRD', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ZooMRN', varid)
        wrt_info=nf90_put_var(ncid, varid, ZooMRN(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ZooMRN', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'DetRR', varid)
        wrt_info=nf90_put_var(ncid, varid, DetRR(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'DetRR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'wPhy', varid)
        wrt_info=nf90_put_var(ncid, varid, wPhy(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'wPhy', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'wDet', varid)
        wrt_info=nf90_put_var(ncid, varid, wDet(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'wDet', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# elif defined ECOSIM
!
!  Write out EcoSim bio-optical model parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'BioIter', varid)
        wrt_info=nf90_put_var(ncid, varid, BioIter(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BioIter', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RtUVR_flag', varid)
        WRITE (char1,'(l1)') RtUVR_flag(ng)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RtUVR_flag', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'NFIX_flag', varid)
        WRITE (char1,'(l1)') NFIX_flag(ng)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'NFIX_flag', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Regen_flag', varid)
        WRITE (char1,'(l1)') Regen_flag(ng)
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Regen_flag', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        start(1)=1
        total(1)=Nphy

        wrt_info=nf90_inq_varid(ncid, 'HsNO3', varid)
        wrt_info=nf90_put_var(ncid, varid, HsNO3(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsNO3', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsNH4', varid)
        wrt_info=nf90_put_var(ncid, varid, HsNH4(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsNH4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsSiO', varid)
        wrt_info=nf90_put_var(ncid, varid, HsSiO(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsSiO', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsPO4', varid)
        wrt_info=nf90_put_var(ncid, varid, HsPO4(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsPO4', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsFe', varid)
        wrt_info=nf90_put_var(ncid, varid, HsFe(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsFe', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'GtALG_max', varid)
        wrt_info=nf90_put_var(ncid, varid, GtALG_max(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'GtALG_max', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyTbase', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyTbase(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyTbase', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'PhyTfac', varid)
        wrt_info=nf90_put_var(ncid, varid, PhyTfac(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'PhyTfac', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'BET_', varid)
        wrt_info=nf90_put_var(ncid, varid, BET_(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BET_', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'maxC2nALG', varid)
        wrt_info=nf90_put_var(ncid, varid, maxC2nALG(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'maxC2nALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'minC2nALG', varid)
        wrt_info=nf90_put_var(ncid, varid, minC2nALG(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'minC2nALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2nALGminABS', varid)
        wrt_info=nf90_put_var(ncid, varid, C2nALGminABS(:,ng),          &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2nALGminABS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'maxC2SiALG', varid)
        wrt_info=nf90_put_var(ncid, varid, maxC2SiALG(:,ng),            &
     &                        start, total)
          IF (wrt_info.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'maxC2SiALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'minC2SiALG', varid)
        wrt_info=nf90_put_var(ncid, varid, minC2SiALG(:,ng),            &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'minC2SiALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2SiALGminABS', varid)
        wrt_info=nf90_put_var(ncid, varid, C2SiALGminABS(:,ng),         &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2SiALGminABS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'maxC2pALG', varid)
        wrt_info=nf90_put_var(ncid, varid, maxC2pALG(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'maxC2pALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'minC2pALG', varid)
        wrt_info=nf90_put_var(ncid, varid, minC2pALG(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'minC2pALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2pALGminABS', varid)
        wrt_info=nf90_put_var(ncid, varid, C2pALGminABS(:,ng),          &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2pALGminABS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'maxC2FeALG', varid)
        wrt_info=nf90_put_var(ncid, varid, maxC2FeALG(:,ng),            &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'maxC2FeALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'minC2FeALG', varid)
        wrt_info=nf90_put_var(ncid, varid, minC2FeALG(:,ng),            &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'minC2FeALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2FeALGminABS', varid)
        wrt_info=nf90_put_var(ncid, varid, C2FeALGminABS(:,ng),         &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2FeALGminABS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'qu_yld', varid)
        wrt_info=nf90_put_var(ncid, varid, qu_yld(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'qu_yld', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'E0_comp', varid)
        wrt_info=nf90_put_var(ncid, varid, E0_comp(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'E0_comp', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'E0_inhib', varid)
        wrt_info=nf90_put_var(ncid, varid, E0_inhib(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'E0_inhib', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'inhib_fac', varid)
        wrt_info=nf90_put_var(ncid, varid, inhib_fac(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'inhib_fac', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2CHL_max', varid)
        wrt_info=nf90_put_var(ncid, varid, C2CHL_max(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2CHL_max', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxC2Cl', varid)
        wrt_info=nf90_put_var(ncid, varid, mxC2Cl(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxC2Cl', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_C2Cl', varid)
        wrt_info=nf90_put_var(ncid, varid, b_C2Cl(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_C2Cl', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxC2Cn', varid)
        wrt_info=nf90_put_var(ncid, varid, mxC2Cn(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxC2Cn', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_C2Cn', varid)
        wrt_info=nf90_put_var(ncid, varid, b_C2Cn(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_C2Cn', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxPacEff', varid)
        wrt_info=nf90_put_var(ncid, varid, mxPacEff(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxPacEff', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_PacEff', varid)
        wrt_info=nf90_put_var(ncid, varid, b_PacEff(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_PacEff', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxChlB', varid)
        wrt_info=nf90_put_var(ncid, varid, mxChlB(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxChlB', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_ChlB', varid)
        wrt_info=nf90_put_var(ncid, varid, b_ChlB(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_ChlB', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxChlC', varid)
        wrt_info=nf90_put_var(ncid, varid, mxChlC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxChlC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_ChlC', varid)
        wrt_info=nf90_put_var(ncid, varid, b_ChlC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_ChlC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxPSC', varid)
        wrt_info=nf90_put_var(ncid, varid, mxPSC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxPSC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_PSC', varid)
        wrt_info=nf90_put_var(ncid, varid, b_PSC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_PSC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxPPC', varid)
        wrt_info=nf90_put_var(ncid, varid, mxPPC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxPPC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_PPC', varid)
        wrt_info=nf90_put_var(ncid, varid, b_PPC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_PPC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxLPUb', varid)
        wrt_info=nf90_put_var(ncid, varid, mxLPUb(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxLPUb', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_LPUb', varid)
        wrt_info=nf90_put_var(ncid, varid, b_LPUb(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_LPUb', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'mxHPUb', varid)
        wrt_info=nf90_put_var(ncid, varid, mxHPUb(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'mxHPUb', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'b_HPUb', varid)
        wrt_info=nf90_put_var(ncid, varid, b_HPUb(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'b_HPUb', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'FecDOC', varid)
        wrt_info=nf90_put_var(ncid, varid, FecDOC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'FecDOC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        start(2)=1
        total(2)=Nfec

        wrt_info=nf90_inq_varid(ncid, 'FecPEL', varid)
        wrt_info=nf90_put_var(ncid, varid, FecPEL(:,:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'FecPEL', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'FecCYC', varid)
        wrt_info=nf90_put_var(ncid, varid, FecCYC(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'FecCYC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ExALG', varid)
        wrt_info=nf90_put_var(ncid, varid, ExALG(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ExALG', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'WS', varid)
        wrt_info=nf90_put_var(ncid, varid, WS(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'WS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsGRZ', varid)
        wrt_info=nf90_put_var(ncid, varid, HsGRZ(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsGRZ', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'MinRefuge', varid)
        wrt_info=nf90_put_var(ncid, varid, MinRefuge(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'MinRefuge', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RefugeDep', varid)
        wrt_info=nf90_put_var(ncid, varid, RefugeDep(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RefugeDep', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Norm_Vol', varid)
        wrt_info=nf90_put_var(ncid, varid, Norm_Vol(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Norm_Vol', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Norm_Surf', varid)
        wrt_info=nf90_put_var(ncid, varid, Norm_Surf(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Norm_Surf', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsDOP', varid)
        wrt_info=nf90_put_var(ncid, varid, HsDOP(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsDOP', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2pALKPHOS', varid)
        wrt_info=nf90_put_var(ncid, varid, C2pALKPHOS(:,ng),            &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2pALKPHOS', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsDON', varid)
        wrt_info=nf90_put_var(ncid, varid, HsDON(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsDON', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2nNupDON', varid)
        wrt_info=nf90_put_var(ncid, varid, C2nNupDON(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2nNupDON', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        start(1)=1
        total(1)=Nbac

        wrt_info=nf90_inq_varid(ncid, 'HsDOC_ba', varid)
        wrt_info=nf90_put_var(ncid, varid, HsDOC_ba(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsDOC_ba', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'GtBAC_max', varid)
        wrt_info=nf90_put_var(ncid, varid, GtBAC_max(:,ng),             &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'GtBAC_max', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'BacTbase', varid)
        wrt_info=nf90_put_var(ncid, varid, BacTbase(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BacTbase', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'BacTfac', varid)
        wrt_info=nf90_put_var(ncid, varid, BacTfac(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BacTfac', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2nBAC', varid)
        wrt_info=nf90_put_var(ncid, varid, C2nBAC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2nBAC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2pBAC', varid)
        wrt_info=nf90_put_var(ncid, varid, C2pBAC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2pBAC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'C2FeBAC', varid)
        wrt_info=nf90_put_var(ncid, varid, C2FeBAC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'C2FeBAC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'BacDOC', varid)
        wrt_info=nf90_put_var(ncid, varid, BacDOC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BacDOC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'BacPEL', varid)
        wrt_info=nf90_put_var(ncid, varid, BacPEL(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BacPEL', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'BacCYC', varid)
        wrt_info=nf90_put_var(ncid, varid, BacCYC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'BacCYC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ExBAC_c', varid)
        wrt_info=nf90_put_var(ncid, varid, ExBAC_c(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ExBAC_c', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'ExBacC2N', varid)
        wrt_info=nf90_put_var(ncid, varid, ExBacC2N(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'ExBacC2N', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Bac_Ceff', varid)
        wrt_info=nf90_put_var(ncid, varid, Bac_Ceff(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Bac_Ceff', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RtNIT', varid)
        wrt_info=nf90_put_var(ncid, varid, RtNIT(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RtNIT', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'HsNIT', varid)
        wrt_info=nf90_put_var(ncid, varid, HsNIT(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'HsNIT', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        start(1)=1
        total(1)=Ndom
        wrt_info=nf90_inq_varid(ncid, 'cDOCfrac_c', varid)
        wrt_info=nf90_put_var(ncid, varid, cDOCfrac_c(:,ng),            &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'cDOCfrac_c', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RtUVR_DIC', varid)
        wrt_info=nf90_put_var(ncid, varid, RtUVR_DIC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RtUVR_DIC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RtUVR_DOC', varid)
        wrt_info=nf90_put_var(ncid, varid, RtUVR_DIC(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RtUVR_DOC', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        start(1)=1
        total(1)=Nfec

        wrt_info=nf90_inq_varid(ncid, 'WF', varid)
        wrt_info=nf90_put_var(ncid, varid, WF(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'WF', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegTbase', varid)
        wrt_info=nf90_put_var(ncid, varid, RegTbase(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegTbase', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegTfac', varid)
        wrt_info=nf90_put_var(ncid, varid, RegTfac(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegTfac', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegCR', varid)
        wrt_info=nf90_put_var(ncid, varid, RegCR(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegCR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegNR', varid)
        wrt_info=nf90_put_var(ncid, varid, RegNR(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegNR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegSR', varid)
        wrt_info=nf90_put_var(ncid, varid, RegSR(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegSR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegPR', varid)
        wrt_info=nf90_put_var(ncid, varid, RegPR(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegPR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'RegFR', varid)
        wrt_info=nf90_put_var(ncid, varid, RegFR(:,ng), start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'RegFR', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

# endif
#endif
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch.
!
        wrt_info=nf90_inq_varid(ncid, 'spherical', varid)
        WRITE (char1,'(l1)') spherical
        wrt_info=nf90_put_var(ncid, varid, char1)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'spherical', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
!
!  Domain Length.
!
        wrt_info=nf90_inq_varid(ncid, 'xl', varid)
        wrt_info=nf90_put_var(ncid, varid, xl)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'xl', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'el', varid)
        wrt_info=nf90_put_var(ncid, varid, el)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'el', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
        wrt_info=nf90_inq_varid(ncid, 'theta_s', varid)
        wrt_info=nf90_put_var(ncid, varid, theta_s(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'theta_s', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'theta_b', varid)
        wrt_info=nf90_put_var(ncid, varid, theta_b(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'theta_b', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Tcline', varid)
        wrt_info=nf90_put_var(ncid, varid, Tcline(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Tcline', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'hc', varid)
        wrt_info=nf90_put_var(ncid, varid, hc(ng))
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'hc', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
!
!  S-coordinate non-dimensional independent variables.
!
        wrt_info=nf90_inq_varid(ncid, 's_rho', varid)
        start(1)=1
        total(1)=N(ng)
        wrt_info=nf90_put_var(ncid, varid, SCALARS(ng)%sc_r(:),         &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 's_rho', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 's_w', varid)
        start(1)=1
        total(1)=N(ng)+1
        wrt_info=nf90_put_var(ncid, varid, SCALARS(ng)%sc_w(0:),        &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 's_w', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
!
!  S-coordinate non-dimensional stretching curves.
!
        wrt_info=nf90_inq_varid(ncid, 'Cs_r', varid)
        start(1)=1
        total(1)=N(ng)
        wrt_info=nf90_put_var(ncid, varid, SCALARS(ng)%Cs_r(:),         &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Cs_r', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF

        wrt_info=nf90_inq_varid(ncid, 'Cs_w', varid)
        start(1)=1
        total(1)=N(ng)+1
        wrt_info=nf90_put_var(ncid, varid, SCALARS(ng)%Cs_w(0:),        &
     &                        start, total)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Cs_w', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
#endif
!
!  User generic parameters.
!
        IF (Nuser.gt.0) THEN
          wrt_info=nf90_inq_varid(ncid, 'user', varid)
          start(1)=1
          total(1)=Nuser
          wrt_info=nf90_put_var(ncid, varid, user, start, total)
          IF (wrt_info.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'user', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
#ifdef STATIONS
!
!  Stations positions.
!
        IF (ncid.eq.ncSTAid(ng)) THEN
          wrt_info=nf90_inq_varid(ncid, 'Ipos', varid)
          start(1)=1
          total(1)=Nstation(ng)
          wrt_info=nf90_put_var(ncid, varid, SCALARS(ng)%SposX(:),      &
     &                          start, total)
          IF (wrt_info.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'Ipos', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'Jpos', varid)
          start(1)=1
          total(1)=Nstation(ng)
          wrt_info=nf90_put_var(ncid, varid,SCALARS(ng)%SposY(:),       &
     &                          start, total)
          IF (wrt_info.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'Jpos', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
#endif
      END IF
#ifdef NO_WRITE_GRID
      IF (ncid.eq.ncSTAid(ng)) THEN
#else
      IF (ncid.ne.ncFLTid(ng)) THEN
#endif
# if !(defined SED_MORPH && defined SEDIMENT)
!
!  Bathymetry.
!
        scale=1.0_r8
        wrt_info=nf90_inq_varid(ncid, 'h', varid)
        IF (ncid.eq.ncSTAid(ng)) THEN
#  ifdef STATIONS
          CALL extract_sta2d (ng, model, Cgrid, varid, r2dvar,          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, GRID(ng)%h,                        &
     &                        Nstation(ng), SCALARS(ng)%SposX,          &
     &                        SCALARS(ng)%SposY, wrk)
          IF (wrtThread) THEN
            start(1)=1
            total(1)=Nstation(ng)
            wrt_info=nf90_put_var(ncid, varid, wrk, start, total)
          END IF
#  endif
        ELSE
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
#  endif
     &                         GRID(ng) % h(LBi,LBj))
        END IF
        IF (wrtThread.and.(wrt_info.ne.nf90_noerr)) THEN
          WRITE (stdout,10) 'h', TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
# endif
!
!  Coriolis parameter.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'f', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % f(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'f', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'pm', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % pm(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'pm', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'pn', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % pn(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'pn', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'lon_rho', varid)
          IF (ncid.eq.ncSTAid(ng)) THEN
# ifdef STATIONS
            CALL extract_sta2d (ng, model, Cgrid, varid, r2dvar,        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%lonr,                   &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              start(1)=1
              total(1)=Nstation(ng)
              wrt_info=nf90_put_var(ncid, varid, wrk, start, total)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % lonr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf90_noerr)) THEN
            WRITE (stdout,10) 'lon_rho', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'lat_rho', varid)
          IF (ncid.eq.ncSTAid(ng)) THEN
# ifdef STATIONS
            CALL extract_sta2d (ng, model, Cgrid, varid, r2dvar,        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%latr,                   &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              start(1)=1
              total(1)=Nstation(ng)
              wrt_info=nf90_put_var(ncid, varid, wrk, start, total)
            END IF 
# endif
          ELSE
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % latr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf90_noerr)) THEN
            WRITE (stdout,10) 'lat_rho', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        ELSE
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'x_rho', varid)
          IF (ncid.eq.ncSTAid(ng)) THEN
# ifdef STATIONS
            CALL extract_sta2d (ng, model, Cgrid, varid, r2dvar,        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%xr,                     &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              start(1)=1
              total(1)=Nstation(ng)
              wrt_info=nf90_put_var(ncid, varid, wrk, start, total)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % xr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf90_noerr)) THEN
            WRITE (stdout,10) 'x_rho', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'y_rho', varid)
          IF (ncid.eq.ncSTAid(ng)) THEN
# ifdef STATIONS
            CALL extract_sta2d (ng, model, Cgrid, varid, r2dvar,        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%yr,                     &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              start(1)=1
              total(1)=Nstation(ng)
              wrt_info=nf90_put_var(ncid, varid, wrk, start, total)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % rmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % yr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf90_noerr)) THEN
            WRITE (stdout,10) 'y_rho', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'lon_u', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % umask(LBi,LBj),             &
# endif
     &                           GRID(ng) % lonu(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'lon_u', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF

          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'lat_u', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % umask(LBi,LBj),             &
# endif
     &                           GRID(ng) % latu(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'lat_u', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF
        ELSE
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'x_u', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % umask(LBi,LBj),             &
# endif
     &                           GRID(ng) % xu(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'x_u', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF

          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'y_u', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % umask(LBi,LBj),             &
# endif
     &                           GRID(ng) % yu(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'y_u', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'lon_v', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % vmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % lonv(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'lon_v', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF

          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'lat_v', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % vmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % latv(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'lat_v', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF
        ELSE
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'x_v', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % vmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % xv(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'x_v', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF

          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'y_v', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % vmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % yv(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'y_v', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'lon_psi', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % pmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % lonp(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'lon_p', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF

          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'lat_psi', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % vmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % latp(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'lat_p', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF
        ELSE
          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'x_psi', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % pmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % xp(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'x_psi', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF

          IF (ncid.ne.ncSTAid(ng)) THEN
            scale=1.0_r8
            wrt_info=nf90_inq_varid(ncid, 'y_psi', varid)
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
# ifdef MASKING
     &                           GRID(ng) % pmask(LBi,LBj),             &
# endif
     &                           GRID(ng) % yp(LBi,LBj))
            IF (wrt_info.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'y_psi', TRIM(ncname)
              END IF
              exit_flag=3
              ioerror=wrt_info
              RETURN
            END IF
          END IF
        END IF

# ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (spherical) THEN
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'angle', varid)
          IF (ncid.eq.ncSTAid(ng)) THEN
#  ifdef STATIONS
            CALL extract_sta2d (ng, model, Cgrid, varid, r2dvar,        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, GRID(ng)%angler,                 &
     &                          Nstation(ng), SCALARS(ng)%SposX,        &
     &                          SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              start(1)=1
              total(1)=Nstation(ng)
              wrt_info=nf90_put_var(ncid, varid, wrk, start, total)
            END IF
#  endif
          ELSE
            wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask(LBi,LBj),             &
#  endif
     &                           GRID(ng) % angler(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf90_noerr)) THEN
            WRITE (stdout,10) 'angle', TRIM(ncname)
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
# endif
# ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'mask_rho', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
     &                         GRID(ng) % rmask(LBi,LBj),               &
     &                         GRID(ng) % rmask(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'mask_rho', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'mask_u', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
     &                         GRID(ng) % umask(LBi,LBj),               &
     &                         GRID(ng) % umask(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'mask_u', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'mask_v', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
     &                         GRID(ng) % vmask(LBi,LBj),               &
     &                         GRID(ng) % vmask(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'mask_v', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'mask_psi', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
     &                         GRID(ng) % pmask(LBi,LBj),               &
     &                         GRID(ng) % pmask(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'mask_psi', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
!
!  Adjoint sensitivity spatial scope mask at RHO-, U-, and V-points.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf90_inq_varid(ncid, 'scope_rho', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
#  endif
     &                         GRID(ng) % Rscope(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'scope_rho', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'scope_u', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % umask(LBi,LBj),               &
#  endif
     &                         GRID(ng) % Uscope(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'scope_u', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF

          wrt_info=nf90_inq_varid(ncid, 'scope_v', varid)
          wrt_info=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,       &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % vmask(LBi,LBj),               &
#  endif
     &                         GRID(ng) % Vscope(LBi,LBj))
          IF (wrt_info.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) 'scope_v', TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=wrt_info
            RETURN
          END IF
        END IF
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (wrtThread) THEN
        wrt_info=nf90_sync(ncid)
        IF (wrt_info.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=3
          ioerror=wrt_info
          RETURN
        END IF
      END IF

  10  FORMAT (/,' WRT_INFO - error while writing variable: ',a,/,       &
     &        12x,'into NetCDF file: ',a)
  20  FORMAT (/,' WRT_INFO - unable to synchronize NetCDF to ',         &
     &        'disk:',/,12x,a)

      RETURN
      END FUNCTION wrt_info
