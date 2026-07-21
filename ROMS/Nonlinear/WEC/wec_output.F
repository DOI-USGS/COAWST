#include "cppdefs.h"
      MODULE wec_output_mod

#if defined WEC || defined WEC_VF
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes Waves Effect on Currents variables into  !
!  output NetCDF files.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef AVERAGES
      USE mod_average
# endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
#ifdef INWAVE_MODEL
      USE mod_inwave_vars
#endif
!
      USE def_var_mod,     ONLY : def_var
# ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
#  ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta3d
#  endif
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod,       ONLY : scale_omega
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: wec_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: wec_def_pio
# endif
# ifdef STATIONS
      PUBLIC :: wec_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: wec_def_station_pio
#  endif
# endif
      PUBLIC :: wec_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: wec_wrt_pio
# endif
# ifdef STATIONS
      PUBLIC :: wec_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: wec_wrt_station_pio
#  endif
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_def_nf90 (ng, model, ldef, VarOut, S,              &
     &                         t2dgrd, u2dgrd, v2dgrd,                  &
     &                         t3dgrd, u3dgrd, v3dgrd, w3dgrd,          &
     &                         a3degrd, u3degrd, v3degrd, t3degrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
      integer, intent(in), optional :: t3dgrd(:), u3dgrd(:), v3dgrd(:)
      integer, intent(in), optional :: w3dgrd(:)
      integer, intent(in), optional :: a3degrd(:), u3degrd(:), v3degrd(:)
      integer, intent(in), optional :: t3degrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, nvd3, nvd4, status
!
      real(r8) :: Aval(6)
!
# ifdef ADJOINT
      character (len=21)  :: Prefix
# else
      character (len=13)  :: Prefix
# endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Waves Effect on Currents output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
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
!  Set long name prefix string.
!
# ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
# else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
# endif
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

# ifdef WEC
!
!  Define 2D total Waves Effect on Currents U-stress.
!
        IF (VarOut(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU2rs))
          ELSE
            Vinfo( 2)=Vname(2,idU2rs)
          END IF
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idU2rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2rs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU2rs),      &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D total Waves Effect on Currents V-stress.
!
        IF (VarOut(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV2rs))
          ELSE
            Vinfo( 2)=Vname(2,idV2rs)
          END IF
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idV2rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2rs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV2rs),      &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Stokes U-velocity.
!
        IF (VarOut(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU2Sd))
          ELSE
            Vinfo( 2)=Vname(2,idU2Sd)
          END IF
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idU2Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU2Sd),      &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Stokes V-velocity.
!
        IF (VarOut(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV2Sd))
          ELSE
            Vinfo( 2)=Vname(2,idV2Sd)
          END IF
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idV2Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV2Sd),      &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define 3D total Waves Effect on Currents U-stress.
!
        IF (VarOut(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU3rs))
          ELSE
            Vinfo( 2)=Vname(2,idU3rs)
          END IF
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idU3rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3rs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU3rs),      &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D total Waves Effect on Currents V-stress.
!
        IF (VarOut(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV3rs))
          ELSE
            Vinfo( 2)=Vname(2,idV3rs)
          END IF
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idV3rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3rs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV3rs),      &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes U-velocity.
!
        IF (VarOut(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU3Sd))
          ELSE
            Vinfo( 2)=Vname(2,idU3Sd)
          END IF
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idU3Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU3Sd),      &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes V-velocity.
        IF (VarOut(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV3Sd))
          ELSE
            Vinfo( 2)=Vname(2,idV3Sd)
          END IF
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idV3Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV3Sd),      &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes omega-velocity.
!
        IF (VarOut(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3Sd))
          ELSE
            Vinfo( 2)=Vname(2,idW3Sd)
          END IF
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idW3Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3Sd,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idW3Sd),      &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes 'true' W-velocity.
!
        IF (VarOut(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3St))
          ELSE
            Vinfo( 2)=Vname(2,idW3St)
          END IF
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idW3St)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3St,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idW3St),      &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
# ifdef WEC_VF
!
!  Define Waves Effect on Currents quasi-static sea level adjustment.
!
        IF (VarOut(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWztw))
          ELSE
            Vinfo( 2)=Vname(2,idWztw)
          END IF
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWztw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWztw,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWztw),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents quasi-static pressure.
!
        IF (VarOut(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWqsp))
          ELSE
            Vinfo( 2)=Vname(2,idWqsp)
          END IF
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWqsp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWqsp,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWqsp),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents Bernoulli head.
!
        IF (VarOut(idWbeh,ng)) THEN
          Vinfo( 1)=Vname(1,idWbeh)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWbeh))
          ELSE
            Vinfo( 2)=Vname(2,idWbeh)
          END IF
          Vinfo( 3)=Vname(3,idWbeh)
          Vinfo(14)=Vname(4,idWbeh)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWbeh)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbeh,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWbeh),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if (defined BOTTOM_STREAMING && defined WEC_VF) || \
      defined WAV_COUPLING
!
!  Define wave dissipation due to bottom friction.
!
        IF (VarOut(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdif))
          ELSE
            Vinfo( 2)=Vname(2,idWdif)
          END IF
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdif,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdif),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
     defined WAVES_DISS     || defined WDISS_INWAVE
#  if defined DISSIP_BREAK_DIR
!
!  Define wave dissipation due to depth limited breaking in x- dir.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdix)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdix))
          ELSE
            Vinfo( 2)=Vname(2,idWdix)
          END IF
          Vinfo( 3)=Vname(3,idWdix)
          Vinfo(14)=Vname(4,idWdix)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdix)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdix,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdix),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to depth limited breaking in y- dir.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiy)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdiy))
          ELSE
            Vinfo( 2)=Vname(2,idWdiy)
          END IF
          Vinfo( 3)=Vname(3,idWdiy)
          Vinfo(14)=Vname(4,idWdiy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdiy)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiy,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdiy),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to whitecapping in x- dir. 
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWcpx)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWcpx))
          ELSE
            Vinfo( 2)=Vname(2,idWcpx)
          END IF
          Vinfo( 3)=Vname(3,idWcpx)
          Vinfo(14)=Vname(4,idWcpx)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWcpx)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWcpx,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWcpx),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to whitecapping in y- dir. 
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWcpy)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWcpy))
          ELSE
            Vinfo( 2)=Vname(2,idWcpy)
          END IF
          Vinfo( 3)=Vname(3,idWcpy)
          Vinfo(14)=Vname(4,idWcpy)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWcpy)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWcpy,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWcpy),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define wave dissipation due to breaking.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdib)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdib))
          ELSE
            Vinfo( 2)=Vname(2,idWdib)
          END IF
          Vinfo( 3)=Vname(3,idWdib)
          Vinfo(14)=Vname(4,idWdib)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdib)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdib,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdib),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to whitecapping.
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdiw))
          ELSE
            Vinfo( 2)=Vname(2,idWdiw)
          END IF
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdiw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiw,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdiw),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
     defined INWAVE_MODEL
!
!  Define percent wave breaking.
!
        IF (VarOut(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWbrk))
          ELSE
            Vinfo( 2)=Vname(2,idWbrk)
          END IF
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWbrk)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbrk,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWbrk),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# ifdef WEC_ROLLER
!
!  Define wave roller dissipation.
!
        IF (VarOut(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdis))
          ELSE
            Vinfo( 2)=Vname(2,idWdis)
          END IF
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdis,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdis),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# ifdef ROLLER_RENIERS
!
!  Define roller wave action density.
!
        IF (VarOut(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWrol))
          ELSE
            Vinfo( 2)=Vname(2,idWrol)
          END IF
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWrol)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWrol,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWrol),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# ifdef INWAVE_MODEL
!
!  Define 3D Action density wave energy.
!
        IF (VarOut(idACen,ng)) THEN
          Vinfo( 1)=Vname(1,idACen)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idACen))
          ELSE
            Vinfo( 2)=Vname(2,idACen)
          END IF
          Vinfo( 3)=Vname(3,idACen)
          Vinfo(14)=Vname(4,idACen)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idACen)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dinw,r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACen),      &
     &                   NF_FOUT, nvd4, a3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D wave celerity xi-direction.
!
        IF (VarOut(idACcx,ng)) THEN
          Vinfo( 1)=Vname(1,idACcx)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idACcx))
          ELSE
            Vinfo( 2)=Vname(2,idACcx)
          END IF
          Vinfo( 3)=Vname(3,idACcx)
          Vinfo(14)=Vname(4,idACcx)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idACcx)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(c3uinw,r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACcx),      &
     &                   NF_FOUT, nvd4, u3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D wave celerity eta-direction.
!
        IF (VarOut(idACcy,ng)) THEN
          Vinfo( 1)=Vname(1,idACcy)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idACcy))
          ELSE
            Vinfo( 2)=Vname(2,idACcy)
          END IF
          Vinfo( 3)=Vname(3,idACcy)
          Vinfo(14)=Vname(4,idACcy)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idACcy)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(c3vinw,r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idACcy),       &
     &                   NF_FOUT, nvd4,v3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D wave celerity theta-direction.
!
        IF (VarOut(idACct,ng)) THEN
          Vinfo( 1)=Vname(1,idACct)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idACct))
          ELSE
            Vinfo( 2)=Vname(2,idACct)
          END IF
          Vinfo( 3)=Vname(3,idACct)
          Vinfo(14)=Vname(4,idACct)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idACct)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(c3tinw,r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idACct),       &
     &                   NF_FOUT, nvd4, t3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Peak period.
!
        IF (VarOut(idACtp,ng)) THEN
          Vinfo( 1)=Vname(1,idACtp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idACtp))
          ELSE
            Vinfo( 2)=Vname(2,idACtp)
          END IF
          Vinfo( 3)=Vname(3,idACtp)
          Vinfo(14)=Vname(4,idACtp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idACtp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dinw,r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idACtp),       &
     &                   NF_FOUT, nvd4, a3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave number.
!
        IF (VarOut(idACkc,ng)) THEN
          Vinfo( 1)=Vname(1,idACkc)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idACkc))
          ELSE
            Vinfo( 2)=Vname(2,idACkc)
          END IF
          Vinfo( 3)=Vname(3,idACkc)
          Vinfo(14)=Vname(4,idACkc)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idACkc)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dinw,r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idACkc),       &
     &                   NF_FOUT, nvd4, a3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined WAVE_PARTITION
!
!  Define HsPT01 Wave height of partition 01 (m).
!
!       IF (VarOut(idHs01,ng)) THEN
          Vinfo( 1)=Vname(1,idHs01)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHs01))
          ELSE
            Vinfo( 2)=Vname(2,idHs01)
          END IF
          Vinfo( 3)=Vname(3,idHs01)
          Vinfo(14)=Vname(4,idHs01)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idHs01)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHs01,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idHs01),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define HsPT02 Wave height of partition 02 (m).
!
!       IF (VarOut(idHs02,ng)) THEN
          Vinfo( 1)=Vname(1,idHs02)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHs02))
          ELSE
            Vinfo( 2)=Vname(2,idHs02)
          END IF
          Vinfo( 3)=Vname(3,idHs02)
          Vinfo(14)=Vname(4,idHs02)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idHs02)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHs02,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idHs02),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define TpPT01 Relative peak period of partition 01 (s).
!
!       IF (VarOut(idTp01,ng)) THEN
          Vinfo( 1)=Vname(1,idTp01)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTp01))
          ELSE
            Vinfo( 2)=Vname(2,idTp01)
          END IF
          Vinfo( 3)=Vname(3,idTp01)
          Vinfo(14)=Vname(4,idTp01)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTp01)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTp01,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idTp01),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define TpPT02 Relative peak period of partition 02 (s).
!
!       IF (VarOut(idTp02,ng)) THEN
          Vinfo( 1)=Vname(1,idTp02)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTp02))
          ELSE
            Vinfo( 2)=Vname(2,idTp02)
          END IF
          Vinfo( 3)=Vname(3,idTp02)
          Vinfo(14)=Vname(4,idTp02)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTp02)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTp02,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idTp02),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define WlPT01 Average wave length of partition 01 (m).
!
!       IF (VarOut(idWl01,ng)) THEN
          Vinfo( 1)=Vname(1,idWl01)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWl01))
          ELSE
            Vinfo( 2)=Vname(2,idWl01)
          END IF
          Vinfo( 3)=Vname(3,idWl01)
          Vinfo(14)=Vname(4,idWl01)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWl01)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWl01,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idWl01),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define WlPT02 Average wave length of partition 02 (m).
!
!       IF (VarOut(idWl02,ng)) THEN
          Vinfo( 1)=Vname(1,idWl02)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWl02))
          ELSE
            Vinfo( 2)=Vname(2,idWl02)
          END IF
          Vinfo( 3)=Vname(3,idWl02)
          Vinfo(14)=Vname(4,idWl02)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWl02)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWl02,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idWl02),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define DrPT01 Average wave direction of partition 01 (deg).
!
!       IF (VarOut(idDr01,ng)) THEN
          Vinfo( 1)=Vname(1,idDr01)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDr01))
          ELSE
            Vinfo( 2)=Vname(2,idDr01)
          END IF
          Vinfo( 3)=Vname(3,idDr01)
          Vinfo(14)=Vname(4,idDr01)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idDr01)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDr01,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idDr01),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
!
!  Define DrPT02 Average wave direction of partition 02 (deg).
!
!       IF (VarOut(idDr02,ng)) THEN
          Vinfo( 1)=Vname(1,idDr02)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDr02))
          ELSE
            Vinfo( 2)=Vname(2,idDr02)
          END IF
          Vinfo( 3)=Vname(3,idDr02)
          Vinfo(14)=Vname(4,idDr02)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idDr02)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDr02,ng),r8)
          status=def_var(ng, iNLM, S(ng)%ncid, S(ng)%Vid(idDr02),       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!       END IF
# endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize local logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
# ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            S(ng)%Vid(idU2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            S(ng)%Vid(idV2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            S(ng)%Vid(idU2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            S(ng)%Vid(idV2Sd)=var_id(i)
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            S(ng)%Vid(idU3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            S(ng)%Vid(idV3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            S(ng)%Vid(idU3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            S(ng)%Vid(idV3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            S(ng)%Vid(idW3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            S(ng)%Vid(idW3St)=var_id(i)
#  endif
# endif
# ifdef WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            S(ng)%Vid(idWztw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            S(ng)%Vid(idWqsp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            S(ng)%Vid(idWbeh)=var_id(i)
# endif
# if (defined BOTTOM_STREAMING && defined WEC_VF) || \
      defined WAV_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            S(ng)%Vid(idWdif)=var_id(i)
# endif
# if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
     defined WAVES_DISS     || defined WDISS_INWAVE
#  if defined DISSIP_BREAK_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdix)=.TRUE.
            S(ng)%Vid(idWdix)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdiy)=.TRUE.
            S(ng)%Vid(idWdiy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWcpx)=.TRUE.
            S(ng)%Vid(idWcpx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWcpy)=.TRUE.
            S(ng)%Vid(idWcpy)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            S(ng)%Vid(idWdib)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            S(ng)%Vid(idWdiw)=var_id(i)
# endif
# if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
     defined INWAVE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            S(ng)%Vid(idWbrk)=var_id(i)
# endif
# ifdef WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            S(ng)%Vid(idWdis)=var_id(i)
# endif
# ifdef ROLLER_RENIERS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            S(ng)%Vid(idWrol)=var_id(i)
# endif
# ifdef INWAVE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACen))) THEN
            got_var(idACen)=.TRUE.
            S(ng)%Vid(idACen)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACcx))) THEN
            got_var(idACcx)=.TRUE.
            S(ng)%Vid(idACcx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACcy))) THEN
            got_var(idACcy)=.TRUE.
            S(ng)%Vid(idACcy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACct))) THEN
            got_var(idACct)=.TRUE.
            S(ng)%Vid(idACct)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACtp))) THEN
            got_var(idACtp)=.TRUE.
            S(ng)%Vid(idACtp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACkc))) THEN
            got_var(idACkc)=.TRUE.
            S(ng)%Vid(idACkc)=var_id(i)
# endif
# if defined WAVE_PARTITION
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHs01))) THEN
            got_var(idHs01)=.TRUE.
            S(ng)%Vid(idHs01)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHs02))) THEN
            got_var(idHs02)=.TRUE.
            S(ng)%Vid(idHs02)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTp01))) THEN
            got_var(idTp01)=.TRUE.
            S(ng)%Vid(idTp01)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTp02))) THEN
            got_var(idTp02)=.TRUE.
            S(ng)%Vid(idTp02)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWl01))) THEN
            got_var(idWl01)=.TRUE.
            S(ng)%Vid(idWl01)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWl02))) THEN
            got_var(idWl02)=.TRUE.
            S(ng)%Vid(idWl02)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDr01))) THEN
            got_var(idDr01)=.TRUE.
            S(ng)%Vid(idDr01)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDr02))) THEN
            got_var(idDr02)=.TRUE.
            S(ng)%Vid(idDr02)=var_id(i)
# endif
          END IF
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef WEC
        IF (.not.got_var(idU2rs).and.VarOut(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.VarOut(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2Sd).and.VarOut(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.VarOut(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idU3Sd).and.VarOut(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.VarOut(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.VarOut(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.VarOut(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.VarOut(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.VarOut(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef WEC_VF
        IF (.not.got_var(idWztw).and.VarOut(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.VarOut(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.VarOut(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbeh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if (defined BOTTOM_STREAMING && defined WEC_VF) || \
      defined WAV_COUPLING
        IF (.not.got_var(idWdif).and.VarOut(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
     defined WAVES_DISS     || defined WDISS_INWAVE
#  if defined DISSIP_BREAK_DIR
        IF (.not.got_var(idWdix).and.VarOut(idWdix,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdix)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiy).and.VarOut(idWdiy,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdiy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWcpx).and.VarOut(idWcpx,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWcpx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWcpy).and.VarOut(idWcpy,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWcpy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idWdib).and.VarOut(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.VarOut(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
     defined INWAVE_MODEL
        IF (.not.got_var(idWbrk).and.VarOut(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WEC_ROLLER
        IF (.not.got_var(idWdis).and.VarOut(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef ROLLER_RENIERS
        IF (.not.got_var(idWrol).and.VarOut(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef INWAVE_MODEL
        IF (.not.got_var(idACen).and.VarOut(idACen,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACcx).and.VarOut(idACcx,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACcx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACcy).and.VarOut(idACcy,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACcy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACct).and.VarOut(idACct,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACct)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACtp).and.VarOut(idACtp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACtp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACkc).and.VarOut(idACkc,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACkc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined WAVE_PARTITION
        IF (.not.got_var(idHs01).and.VarOut(idHs01,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idHs01)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHs02).and.VarOut(idHs02,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idHs02)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTp01).and.VarOut(idTp01,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTp01)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTp02).and.VarOut(idTp02,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTp02)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWl01).and.VarOut(idWl01,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWl01)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWl02).and.VarOut(idWl02,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWl02)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDr01).and.VarOut(idDr01,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idDr01)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDr02).and.VarOut(idDr02,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idDr02)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
      END IF QUERY
!
  10  FORMAT (/,' WEC_DEF_NF90 - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE wec_def_nf90

# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE wec_def_station_nf90 (ng, model, ldef, VarOut, S,      &
     &                                 pgrd, rgrd, igrd, cgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: pgrd(:), rgrd(:)
      integer, intent(in), optional :: igrd(:), cgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, status
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Waves Effect on Currents output stations variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
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

#  ifdef WEC
!
!  Define 2D Stokes U-velocity.
!
        IF (VarOut(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          Vinfo( 2)=Vname(2,idU2Sd)
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU2Sd),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Stokes V-velocity.
!
        IF (VarOut(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          Vinfo( 2)=Vname(2,idV2Sd)
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV2Sd),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D total Waves Effect on Currents U-stress.
!
        IF (VarOut(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          Vinfo( 2)=Vname(2,idU2rs)
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU2rs),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D total Waves Effect on Currents V-stress.
!
        IF (VarOut(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          Vinfo( 2)=Vname(2,idV2rs)
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV2rs),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SOLVE3D
!
!  Define 3D Stokes U-velocity.
!
        IF (VarOut(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          Vinfo( 2)=Vname(2,idU3Sd)
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU3Sd),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes V-velocity.
!
        IF (VarOut(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          Vinfo( 2)=Vname(2,idV3Sd)
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV3Sd),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes omega-velocity.
!
        IF (VarOut(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          Vinfo( 2)=Vname(2,idW3Sd)
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idW3Sd),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D "true" Stokes W-velocity.
!
        IF (VarOut(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          Vinfo( 2)=Vname(2,idW3St)
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idW3St),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D total Waves Effect on Currents U-stress.
!
        IF (VarOut(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          Vinfo( 2)=Vname(2,idU3rs)
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idU3rs),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D total Waves Effect on Current V-stress.
!
        IF (VarOut(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          Vinfo( 2)=Vname(2,idV3rs)
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idV3rs),      &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
#  ifdef WEC_VF
!
!  Define Waves Effect on Currents quasi-static sea level adjustment.
!
        IF (VarOut(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          Vinfo( 2)=Vname(2,idWztw)
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWztw),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents quasi-static pressure.
!
        IF (VarOut(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          Vinfo( 2)=Vname(2,idWqsp)
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWqsp),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents Bernoulli head.
!
        IF (VarOut(idWbeh,ng)) THEN
          Vinfo( 1)=Vname(1,idWbeh)
          Vinfo( 2)=Vname(2,idWbeh)
          Vinfo( 3)=Vname(3,idWbeh)
          Vinfo(14)=Vname(4,idWbeh)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWbeh),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
!
!  Define wave dissipation due to bottom friction.
!
        IF (VarOut(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          Vinfo( 2)=Vname(2,idWdif)
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdif),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
#   if defined DISSIP_BREAK_DIR
!
!  Define wave dissipation due to breaking in x-direction.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdix)
          Vinfo( 2)=Vname(2,idWdix)
          Vinfo( 3)=Vname(3,idWdix)
          Vinfo(14)=Vname(4,idWdix)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdix),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to breaking in y-direction.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiy)
          Vinfo( 2)=Vname(2,idWdiy)
          Vinfo( 3)=Vname(3,idWdiy)
          Vinfo(14)=Vname(4,idWdiy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdiy),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to white capping in x-direction.
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWcpx)
          Vinfo( 2)=Vname(2,idWcpx)
          Vinfo( 3)=Vname(3,idWcpx)
          Vinfo(14)=Vname(4,idWcpx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWcpx),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to white capping in y-direction.
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWcpy)
          Vinfo( 2)=Vname(2,idWcpy)
          Vinfo( 3)=Vname(3,idWcpy)
          Vinfo(14)=Vname(4,idWcpy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWcpy),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define wave dissipation due to breaking.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdib)
          Vinfo( 2)=Vname(2,idWdib)
          Vinfo( 3)=Vname(3,idWdib)
          Vinfo(14)=Vname(4,idWdib)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdib),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to whitecapping.
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          Vinfo( 2)=Vname(2,idWdiw)
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdiw),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
      defined INWAVE_MODEL
!
!  Define percent wave breaking.
!
        IF (VarOut(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          Vinfo( 2)=Vname(2,idWbrk)
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWbrk),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef WEC_ROLLER
!
!  Define wave roller dissipation.
!
        IF (VarOut(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          Vinfo( 2)=Vname(2,idWdis)
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdis),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef ROLLER_RENIERS
!
!  Define wave roller action density.
!
        IF (VarOut(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          Vinfo( 2)=Vname(2,idWrol)
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWrol),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef INWAVE_MODEL
!
!  Define 3D Action density wave energy.
!
        IF (VarOut(idACen,ng)) THEN
          Vinfo( 1)=Vname(1,idACen)
          Vinfo( 2)=Vname(2,idACen)
          Vinfo( 3)=Vname(3,idACen)
          Vinfo(14)=Vname(4,idACen)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACen),      &
     &                   NF_FOUT, 3, igrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D wave celerity xi-direction.
!
        IF (VarOut(idACcx,ng)) THEN
          Vinfo( 1)=Vname(1,idACcx)
          Vinfo( 2)=Vname(2,idACcx)
          Vinfo( 3)=Vname(3,idACcx)
          Vinfo(14)=Vname(4,idACcx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACcx),      &
     &                   NF_FOUT, 3, igrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D wave celerity eta-direction.
!
        IF (VarOut(idACcx,ng)) THEN
          Vinfo( 1)=Vname(1,idACcy)
          Vinfo( 2)=Vname(2,idACcy)
          Vinfo( 3)=Vname(3,idACcy)
          Vinfo(14)=Vname(4,idACcy)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACcy),      &
     &                   NF_FOUT, 3, igrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D wave celerity theta-direction.
!
        IF (VarOut(idACct,ng)) THEN
          Vinfo( 1)=Vname(1,idACct)
          Vinfo( 2)=Vname(2,idACct)
          Vinfo( 3)=Vname(3,idACct)
          Vinfo(14)=Vname(4,idACct)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACct),      &
     &                   NF_FOUT, 3, cgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Peak period.
!
        IF (VarOut(idACtp,ng)) THEN
          Vinfo( 1)=Vname(1,idACtp)
          Vinfo( 2)=Vname(2,idACtp)
          Vinfo( 3)=Vname(3,idACtp)
          Vinfo(14)=Vname(4,idACtp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACtp),      &
     &                   NF_FOUT, 3, igrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave number.
!
        IF (VarOut(idACkc,ng)) THEN
          Vinfo( 1)=Vname(1,idACkc)
          Vinfo( 2)=Vname(2,idACkc)
          Vinfo( 3)=Vname(3,idACkc)
          Vinfo(14)=Vname(4,idACkc)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idACkc),      &
     &                   NF_FOUT, 3, igrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
#  ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            S(ng)%Vid(idU2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            S(ng)%Vid(idV2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            S(ng)%Vid(idU2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            S(ng)%Vid(idV2rs)=var_id(i)
#   ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            S(ng)%Vid(idU3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            S(ng)%Vid(idV3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            S(ng)%Vid(idW3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            S(ng)%Vid(idW3St)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            S(ng)%Vid(idU3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            S(ng)%Vid(idV3rs)=var_id(i)
#   endif
#  endif
#  ifdef WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            S(ng)%Vid(idWztw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            S(ng)%Vid(idWqsp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            S(ng)%Vid(idWbeh)=var_id(i)
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            S(ng)%Vid(idWdif)=var_id(i)
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
#   if defined DISSIP_BREAK_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdix)=.TRUE.
            S(ng)%Vid(idWdix)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdiy)=.TRUE.
            S(ng)%Vid(idWdiy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWcpx)=.TRUE.
            S(ng)%Vid(idWcpx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWcpy)=.TRUE.
            S(ng)%Vid(idWcpy)=var_id(i)
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            S(ng)%Vid(idWdib)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            S(ng)%Vid(idWdiw)=var_id(i)
#  endif
#  if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
      defined INWAVE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            S(ng)%Vid(idWbrk)=var_id(i)
#  endif
#  ifdef WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            S(ng)%Vid(idWdis)=var_id(i)
#  endif
#  ifdef ROLLER_RENIERS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            S(ng)%Vid(idWrol)=var_id(i)
#  endif
#  ifdef INWAVE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACen))) THEN
            got_var(idACen)=.TRUE.
            S(ng)%Vid(idACen)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACcx))) THEN
            got_var(idACcx)=.TRUE.
            S(ng)%Vid(idACcx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACcy))) THEN
            got_var(idACcy)=.TRUE.
            S(ng)%Vid(idACcy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACct))) THEN
            got_var(idACct)=.TRUE.
            S(ng)%Vid(idACct)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACtp))) THEN
            got_var(idACtp)=.TRUE.
            S(ng)%Vid(idACtp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idACkc))) THEN
            got_var(idACkc)=.TRUE.
            S(ng)%Vid(idACkc)=var_id(i)
#  endif
          END IF
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef WEC
        IF (.not.got_var(idU2Sd).and.Sout(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Sout(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Sout(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Sout(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SOLVE3D
        IF (.not.got_var(idU3Sd).and.Sout(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Sout(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.Sout(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.Sout(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Sout(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Sout(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef WEC_VF
        IF (.not.got_var(idWztw).and.Sout(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.Sout(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.Sout(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbeh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
        IF (.not.got_var(idWdif).and.Sout(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
#   if defined DISSIP_BREAK_DIR
        IF (.not.got_var(idWdix).and.Sout(idWdix,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdix)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiy).and.Sout(idWdiy,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdiy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWcpx).and.Sout(idWcpx,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWcpx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWcpy).and.Sout(idWcpy,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWcpy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idWdib).and.Sout(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.Sout(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
      defined INWAVE_MODEL
        IF (.not.got_var(idWbrk).and.Sout(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WEC_ROLLER
        IF (.not.got_var(idWdis).and.Sout(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ROLLER_RENIERS
        IF (.not.got_var(idWrol).and.Sout(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef INWAVE_MODEL
        IF (.not.got_var(idACen).and.Sout(idACen,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACcx).and.Sout(idACcx,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACcx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACcy).and.Sout(idACcy,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACcy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACct).and.Sout(idACct,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACct)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACtp).and.Sout(idACtp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACtp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idACkc).and.Sout(idACkc,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idACkc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
      END IF QUERY
!
  10  FORMAT (/,' WEC_DEF_STATION_NF90 - unable to find variable: ',    &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE wec_def_station_nf90
# endif
!
!*********************************************************************** 
      SUBROUTINE wec_wrt_nf90 (ng, model, tile,                         &
    &                          LBi, UBi, LBj, UBj,                      &
    &                          VarOut, S)
!*********************************************************************** 
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Linstataneous
!
      integer :: gfactor, gtype, status
!
      real(dp) :: scale

# ifdef SOLVE3D
!
      real(r8), allocatable :: Wr3d(:,:,:)
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents output variables into specified
!  output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!  Set instantaneous fields.
!
      IF ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                              &
     &    (S(ng)%ncid.eq.QCK(ng)%ncid)) THEN
        Linstataneous=.TRUE.
      ELSE
        Linstataneous=.FALSE.                  ! time-averged fiels
      END IF

# ifdef WEC
!
!  Write out 2D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU2rs,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idU2rs,             &
     &                       S(ng)%Vid(idU2rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       MIXING(ng) % rustr2d)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idU2rs,             &
     &                       S(ng)%Vid(idU2rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       AVERAGE(ng) % avgu2rs)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D total Waves Effect on Currents V-stress.
!
      IF (VarOut(idV2rs,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idV2rs,             &
     &                       S(ng)%Vid(idV2rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       MIXING(ng) % rvstr2d)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idV2rs,             &
     &                       S(ng)%Vid(idV2rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       AVERAGE(ng) % avgv2rs)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Stokes U-velocity.
!
      IF (VarOut(idU2Sd,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idU2sd,             &
     &                       S(ng)%Vid(idU2sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       OCEAN(ng) % ubar_stokes)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idU2sd,             &
     &                       S(ng)%Vid(idU2sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       AVERAGE(ng) % avgu2Sd)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Stokes V-velocity.
!
      IF (VarOut(idV2Sd,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idV2sd,             &
     &                       S(ng)%Vid(idV2sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       OCEAN(ng) % vbar_stokes)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idV2sd,             &
     &                       S(ng)%Vid(idV2sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       AVERAGE(ng) % avgv2Sd)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU3rs,ng)) THEN
        scale=rho0
        gtype=gfactor*u3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idU3rs,             &
     &                       S(ng)%Vid(idU3rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       MIXING(ng) % rustr3d)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idU3rs,             &
     &                       S(ng)%Vid(idU3rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask,                          &
#    endif
     &                       AVERAGE(ng) % avgu3rs)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D total Waves Effect on Currents V-stress.
!
      IF (VarOut(idV3rs,ng)) THEN
        scale=rho0
        gtype=gfactor*v3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idV3rs,             &
     &                       S(ng)%Vid(idV3rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       MIXING(ng) % rvstr3d)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idV3rs,             &
     &                       S(ng)%Vid(idV3rs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#    endif
     &                       AVERAGE(ng) % avgv3rs)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes U-velocity.
!
      IF (VarOut(idU3Sd,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idU3Sd,             &
     &                       S(ng)%Vid(idU3Sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       OCEAN(ng) % u_stokes)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idU3Sd,             &
     &                       S(ng)%Vid(idU3Sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask,                          &
#    endif
     &                       AVERAGE(ng) % avgu3Sd)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes V-velocity.
!
      IF (VarOut(idV3Sd,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idV3Sd,             &
     &                       S(ng)%Vid(idV3Sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       OCEAN(ng) % v_stokes)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idV3Sd,             &
     &                       S(ng)%Vid(idV3Sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#    endif
     &                       AVERAGE(ng) % avgv3Sd)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes omega-velocity.
!
      IF (VarOut(idW3Sd,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        IF (Linstataneous) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
            Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
          END IF
          CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),     &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % W_stokes,                       &
     &                      Wr3d)
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idW3Sd,             &
     &                       S(ng)%Vid(idW3Sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       Wr3d)
          deallocate (Wr3d)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idW3Sd,             &
     &                       S(ng)%Vid(idW3Sd),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgw3d)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes 'true' W-velocity (m/s)
!
      IF (VarOut(idW3St,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idW3St,             &
     &                       S(ng)%Vid(idW3St),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       OCEAN(ng) % wstvel)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idW3St,             &
     &                       S(ng)%Vid(idW3St),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgW3St)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3St)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef WEC_VF
!
!  Write out Waves Effect on Currents quasi-static sea level adjustment.
!
      IF (VarOut(idWztw,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWztw,             &
     &                       S(ng)%Vid(idWztw),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % zetaw)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWztw,             &
     &                       S(ng)%Vid(idWztw),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWztw)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWztw)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Waves Effect on Currents quasi-static pressure.
!
      IF (VarOut(idWqsp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWqsp,             &
     &                       S(ng)%Vid(idWqsp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % qsp)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWqsp,             &
     &                       S(ng)%Vid(idWqsp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWqsp)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWqsp)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Waves Effect on Currents Bernoulli head.
!
      IF (VarOut(idWbeh,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWbeh,             &
     &                       S(ng)%Vid(idWbeh),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % bh)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWbeh,             &
     &                       S(ng)%Vid(idWbeh),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWbeh)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbeh)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if (defined BOTTOM_STREAMING && defined WEC_VF) || \
      defined WAV_COUPLING
!
!  Write out wave dissipation due to bottom friction.
!
      IF (VarOut(idWdif,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdif,             &
     &                       S(ng)%Vid(idWdif),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Dissip_fric)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdif,             &
     &                       S(ng)%Vid(idWdif),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWdif)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdif)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
     defined WAVES_DISS     || defined WDISS_INWAVE
#  if defined DISSIP_BREAK_DIR
!
!  Write out wave dissipation due to breaking in x-direction.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdix,             &
     &                       S(ng)%Vid(idWdix),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_breakx)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdix,             &
     &                       S(ng)%Vid(idWdix),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWdix)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdix)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave dissipation due to breaking in y-direction.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdiy,             &
     &                       S(ng)%Vid(idWdiy),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_breaky)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdiy,             &
     &                       S(ng)%Vid(idWdiy),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdiy)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdiy)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave dissipation due to whitecapping in x-direction.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWcpx,             &
     &                       S(ng)%Vid(idWcpx),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_wcapx)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWcpx,             &
     &                       S(ng)%Vid(idWcpx),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWcpx)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWcpx)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave dissipation due to whitecapping in y-direction.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWcpy,             &
     &                       S(ng)%Vid(idWcpy),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_wcapy)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWcpy,             &
     &                       S(ng)%Vid(idWcpy),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWcpy)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWcpy)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out wave dissipation due to breaking.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdib,             &
     &                       S(ng)%Vid(idWdib),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       FORCES(ng) % Dissip_break)
# ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdib,             &
     &                       S(ng)%Vid(idWdib),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgWdib)
# endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdib)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave dissipation due to whitecapping.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdiw,             &
     &                       S(ng)%Vid(idWdiw),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       FORCES(ng) % Dissip_wcap)
# ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdiw,             &
     &                       S(ng)%Vid(idWdiw),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgWdiw)
# endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdiw)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined ROLLER_SVENDSEN || defined WW3_MODEL || \
     defined INWAVE_MODEL || defined SWAN_MODEL
!
!  Write out percent wave breaking.
!
      IF (VarOut(idWbrk,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWbrk,             &
     &                       S(ng)%Vid(idWbrk),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Wave_break)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWbrk,             &
     &                       S(ng)%Vid(idWbrk),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWbrk)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbrk)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WEC_ROLLER
!
!  Write out wave roller dissipation.
!
      IF (VarOut(idWdis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdis,             &
     &                       S(ng)%Vid(idWdis),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Dissip_roller)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdis,             &
     &                       S(ng)%Vid(idWdis),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWdis)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdis)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ROLLER_RENIERS
!
!  Write out roller wave action density.
!
      IF (VarOut(idWrol,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWrol,             &
     &                       S(ng)%Vid(idWrol),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % rollA)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWrol,             &
     &                       S(ng)%Vid(idWrol),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWrol)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWrol)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef INWAVE_MODEL
!
!  Write out 3D Action density wave energy.
!
      IF (VarOut(idACen,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACen,             &
     &                       S(ng)%Vid(idACen),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       WAVEP(ng) % AC(:,:,:,NOUT))
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACen,             &
     &                       S(ng)%Vid(idACen),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#   ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#   endif
     &                       AVERAGE(ng) % avgACen)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACen)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D wave celerity in xi direction.
!
      IF (VarOut(idACcx,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACcx,             &
     &                       S(ng)%Vid(idACcx),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#  ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#  endif
     &                       WAVEP(ng) % cx)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACcx,             &
     &                       S(ng)%Vid(idACcx),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#   ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#   endif
     &                       AVERAGE(ng) % avgACcx)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACcx)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D wave celerity in eta direction.
!
      IF (VarOut(idACcy,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACcy,             &
     &                       S(ng)%Vid(idACcy),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#  ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#  endif
     &                       WAVEP(ng) % cy)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACcy,             &
     &                       S(ng)%Vid(idACcy),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       AVERAGE(ng) % avgACcy)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACcy)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D wave celerity in theta direction.
!
      IF (VarOut(idACct,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACct,             &
     &                       S(ng)%Vid(idACct),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND+1, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       WAVEP(ng) % ct)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACct,             &
     &                       S(ng)%Vid(idACct),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND+1, scale,        &
#   ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#   endif
     &                       AVERAGE(ng) % avgACct)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACct)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave peak period.
!
      IF (VarOut(idACtp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACtp,             &
     &                       S(ng)%Vid(idACtp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       WAVEP(ng) % Ta)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACtp,             &
     &                       S(ng)%Vid(idACtp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#   ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#   endif
     &                       AVERAGE(ng) % avgACtp)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACtp)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave peak period.
!
      IF (VarOut(idACkc,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACkc,             &
     &                       S(ng)%Vid(idACkc),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       WAVEP(ng) % kwc)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idACkc,             &
     &                       S(ng)%Vid(idACkc),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, ND, scale,          &
#   ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#   endif
     &                       AVERAGE(ng) % avgACkc)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACkc)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
  10  FORMAT (/," WEC_WRT_NF90 - error while writing variable '",a,     &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE wec_wrt_nf90

# ifdef STATIONS
!
!*********************************************************************** 
     SUBROUTINE wec_wrt_station_nf90 (ng, model, tile,                  &
    &                                 LBi, UBi, LBj, UBj,               &
    &                                 VarOut, S)
!*********************************************************************** 
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposR, NposW
      integer :: i, k, np, status
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
#  ifdef SOLVE3D
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
#  endif
#  ifdef INWAVE_MODEL
      real(r8), dimension(Nstation(ng)*(ND)) :: ista
      real(r8), dimension(Nstation(ng)*(ND)) :: dsta
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents output variables into specified
!  stations output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
#  ifdef SOLVE3D
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
#  endif
      END DO

#  ifdef WEC
!
!  Write out 2D Stokes U-velocity.
!
      IF (VarOut(idU2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idU2Sd, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % ubar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idU2Sd)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idU2Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Stokes V-velocity.
!
      IF (VarOut(idV2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idV2Sd, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % vbar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idV2Sd)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idV2Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idU2rs, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rustr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idU2rs)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idU2rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D total Waves Effect on Currents V-stress.
!
      IF (VarOut(idV2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idV2rs, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rvstr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idV2rs)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idV2rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef SOLVE3D
!
!  Write out 3D Stokes U-velocity.
!
      IF (VarOut(idU3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idU3Sd, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % u_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idU3Sd)), rsta,              &
     &                        (/1,1,S(ng)%Rindex/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idU3Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Stokes V-velocity.
!
      IF (VarOut(idV3Sd,ng)) THEN
        scale=1.0_dp
       CALL extract_sta3d (ng, model, Cgrid, idV3Sd, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % v_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idV3Sd)), rsta,              &
     &                        (/1,1,S(ng)%Rindex/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idV3Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Stokes omega-velocity.
!
      IF (VarOut(idW3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idW3Sd, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % W_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idW3Sd)), rsta,              &
     &                        (/1,1,S(ng)%Rindex/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idW3Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D "true" Stokes W-velocity.
!
      IF (VarOut(idW3St,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idW3St, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % wstvel,                  &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idW3St)), rsta,              &
     &                        (/1,1,S(ng)%Rindex/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idW3St))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, model, Cgrid, idU3rs, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rustr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idU3rs)), rsta,              &
     &                        (/1,1,S(ng)%Rindex/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idU3rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D total Waves Effect on Current V-stress.
!
      IF (VarOut(idV3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, model, Cgrid, idV3rs, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rvstr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idV3rs)), rsta,              &
     &                        (/1,1,S(ng)%Rindex/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idV3rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#  endif
#  ifdef WEC_VF
!
!  Write out Waves Effect on Currents quasi-static sea level adjustment.
!
      IF (VarOut(idWztw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWztw, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % zeta(:,:,KOUT),          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWztw)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWztw))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out Waves Effect on Currents quasi-static pressure.
!
      IF (VarOut(idWqsp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWqsp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % qsp,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWqsp)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWqsp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out Waves Effect on Currents Bernoulli head.
!
      IF (VarOut(idWbeh,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWbeh, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % bh,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWbeh)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWbeh))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
!
!  Write out wave dissipation due to bottom friction.
!
      IF (VarOut(idWdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdif, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_fric,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdif)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
#   if defined DISSIP_BREAK_DIR
!
!  Write out wave dissipation due to breaking in x-direction.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdix, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_breakx,          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdix)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdix))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to breaking in y-direction.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdiy, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_breaky,          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdiy)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdiy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to whitecapping in x-direction.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWcpx, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_wcapx,           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWcpx)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWcpx))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to whitecapping in y-direction.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWcpy, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_wcapy,           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWcpx)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWcpy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
!
!  Write out wave dissipation due to breaking.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdib, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_break,           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdib)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdib))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to whitecapping.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdiw, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_wcap,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdiw)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdiw))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef ROLLER_SVENDSEN || defined WW3_MODEL || \
   defined INWAVE_MODEL || defined SWAN_MODEL
!
!  Write out percent wave breaking.
!
      IF (VarOut(idWbrk,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWbrk, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_break,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWbrk)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWbrk))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WEC_ROLLER
!
!  Write out wave roller dissipation.
!
      IF (VarOut(idWdis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdis, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_roller,          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdis)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdis))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef ROLLER_RENIERS
!
!  Write out wave roller action density.
!
      IF (VarOut(idWrol,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWrol, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % rollA,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWrol)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWrol))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef INWAVE_MODEL_notyet
!
!  Write out 3D Action density wave energy.
!
      IF (VarOut(idACen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idACen, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, ND,                  &
     &                      scale, WAVEP(ng)%AC(:,:,:,NOUT),            &
     &                      NposR, XposR, YposR, ZposR, ista)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idACen)), ista,              &
     &                        (/1,S(ng)%Rindex/),(/ND,Nstation(ng),1/), &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idACen))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Action wave celerity xi direction.
!
      IF (VarOut(idACcx,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idACcx, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, ND,                  &
     &                      scale, WAVEP(ng)%cx,                        &
     &                      NposR, XposR, YposR, ZposR, ista)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idACcx)), ista,              &
     &                        (/1,S(ng)%Rindex/),(/ND,Nstation(ng),1/), &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idACcx))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Action wave celerity eta direction.
!
      IF (VarOut(idACcy,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idACcy, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, ND,                  &
     &                      scale, WAVEP(ng)%cy,                        &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idACcy)), ista,              &
     &                        (/1,S(ng)%Rindex/),(/ND,Nstation(ng),1/), &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idACcy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Action wave celerity theta direction.
!
      IF (VarOut(idACct,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idACct, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, ND+1,                &
     &                      scale, WAVEP(ng)%ct,                        &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idACct)), dsta,              &
     &                      (/1,S(ng)%Rindex/),(/ND+1,Nstation(ng),1/), &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idACct))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Action wave peak period.
!
      IF (VarOut(idACtp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idACtp, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, ND,                  &
     &                      scale, WAVEP(ng)%Ta,                        &
     &                      NposR, XposR, YposR, ZposR, ista)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idACct)), ista,              &
     &                        (/1,S(ng)%Rindex/),(/ND,Nstation(ng),1/), &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idACtp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Action wave peak period.
!
      IF (VarOut(idACkc,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idACkc, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, ND,                  &
     &                      scale, WAVEP(ng)%kwc,                       &
     &                      NposR, XposR, YposR, ZposR, ista)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idACkc)), ista,              &
     &                        (/1,S(ng)%Rindex/),(/ND,Nstation(ng),1/), &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idACkc))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
!
      RETURN
      END SUBROUTINE wec_wrt_station_nf90
# endif

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wec_def_pio (ng, model, ldef, VarOut, S,               &
     &                        t2dgrd, u2dgrd, v2dgrd,                   &
     &                        t3dgrd, u3dgrd, v3dgrd, w3dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
      integer, intent(in), optional :: t3dgrd(:), u3dgrd(:), v3dgrd(:)
      integer, intent(in), optional :: w3dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, nvd3, nvd4, status
!
      real(r8) :: Aval(6)
!
#  ifdef ADJOINT
      character (len=21)  :: Prefix
#  else
      character (len=13)  :: Prefix
#  endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Waves Effect on Currents output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
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
!  Set long name prefix string.
!
#  ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
#  else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
#  endif
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

#  ifdef WEC
!
!  Define 2D total Waves Effect on Currents U-stress.
!
        IF (VarOut(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU2rs))
          ELSE
            Vinfo( 2)=Vname(2,idU2rs)
          END IF
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idU2rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2rs,ng),r8)
          S(ng)%pioVar(idU2rs)%dkind=PIO_FOUT
          S(ng)%pioVar(idU2rs)%gtype=u2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idU2rs)%vd,                       &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D total Waves Effect on Currents V-stress.
!
        IF (VarOut(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV2rs))
          ELSE
            Vinfo( 2)=Vname(2,idV2rs)
          END IF
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idV2rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2rs,ng),r8)
          S(ng)%pioVar(idV2rs)%dkind=PIO_FOUT
          S(ng)%pioVar(idV2rs)%gtype=v2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idV2rs)%vd,                       &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Stokes U-velocity.
!
        IF (VarOut(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU2Sd))
          ELSE
            Vinfo( 2)=Vname(2,idU2Sd)
          END IF
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idU2Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
          S(ng)%pioVar(idU2Sd)%dkind=PIO_FOUT
          S(ng)%pioVar(idU2Sd)%gtype=u2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idU2Sd)%vd,                       &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Stokes V-velocity.
!
        IF (VarOut(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV2Sd))
          ELSE
            Vinfo( 2)=Vname(2,idV2Sd)
          END IF
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idV2Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
          S(ng)%pioVar(idV2Sd)%dkind=PIO_FOUT
          S(ng)%pioVar(idV2Sd)%gtype=v2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idV2Sd)%vd,                       &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SOLVE3D
!
!  Define 3D total Waves Effect on Currents U-stress.
!
        IF (VarOut(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU3rs))
          ELSE
            Vinfo( 2)=Vname(2,idU3rs)
          END IF
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idU3rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3rs,ng),r8)
          S(ng)%pioVar(idU3rs)%dkind=PIO_FOUT
          S(ng)%pioVar(idU3rs)%gtype=u3dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idU3rs)%vd,                       &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D total Waves Effect on Currents V-stress.
!
        IF (VarOut(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV3rs))
          ELSE
            Vinfo( 2)=Vname(2,idV3rs)
          END IF
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idV3rs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3rs,ng),r8)
          S(ng)%pioVar(idV3rs)%dkind=PIO_FOUT
          S(ng)%pioVar(idV3rs)%gtype=v3dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idV3rs)%vd,                       &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes U-velocity.
!
        IF (VarOut(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU3Sd))
          ELSE
            Vinfo( 2)=Vname(2,idU3Sd)
          END IF
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idU3Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
          S(ng)%pioVar(idU3Sd)%dkind=PIO_FOUT
          S(ng)%pioVar(idU3Sd)%gtype=u3dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idU3Sd)%vd,                       &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes V-velocity.
!
        IF (VarOut(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV3Sd))
          ELSE
            Vinfo( 2)=Vname(2,idV3Sd)
          END IF
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idV3Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
          S(ng)%pioVar(idV3Sd)%dkind=PIO_FOUT
          S(ng)%pioVar(idV3Sd)%gtype=v3dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idV3Sd)%vd,                       &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes omega-velocity.
!
        IF (VarOut(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3Sd))
          ELSE
            Vinfo( 2)=Vname(2,idW3Sd)
          END IF
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idW3Sd)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3Sd,ng),r8)
          S(ng)%pioVar(idW3Sd)%dkind=PIO_FOUT
          S(ng)%pioVar(idW3Sd)%gtype=w3dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idW3Sd)%vd,                       &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes W-velocity.
!
        IF (VarOut(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3St))
          ELSE
            Vinfo( 2)=Vname(2,idW3St)
          END IF
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idW3St)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3St,ng),r8)
          S(ng)%pioVar(idW3St)%dkind=PIO_FOUT
          S(ng)%pioVar(idW3St)%gtype=w3dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idW3St)%vd,                       &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
#  ifdef WEC_VF
!
!  Define Waves Effect on Currents quasi-static sea level adjustment.
!
        IF (VarOut(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWztw))
          ELSE
            Vinfo( 2)=Vname(2,idWztw)
          END IF
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWztw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWztw,ng),r8)
          S(ng)%pioVar(idWztw)%dkind=PIO_FOUT
          S(ng)%pioVar(idWztw)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWztw)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents quasi-static pressure.
!
        IF (VarOut(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWqsp))
          ELSE
            Vinfo( 2)=Vname(2,idWqsp)
          END IF
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWqsp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWqsp,ng),r8)
          S(ng)%pioVar(idWqsp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWqsp)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWqsp)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents Bernoulli head.
!
        IF (VarOut(idWbeh,ng)) THEN
          Vinfo( 1)=Vname(1,idWbeh)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWbeh))
          ELSE
            Vinfo( 2)=Vname(2,idWbeh)
          END IF
          Vinfo( 3)=Vname(3,idWbeh)
          Vinfo(14)=Vname(4,idWbeh)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWbeh)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbeh,ng),r8)
          S(ng)%pioVar(idWbeh)%dkind=PIO_FOUT
          S(ng)%pioVar(idWbeh)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWbeh)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
!
!  Define wave dissipation due to bottom friction.
!
        IF (VarOut(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdif))
          ELSE
            Vinfo( 2)=Vname(2,idWdif)
          END IF
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdif,ng),r8)
          S(ng)%pioVar(idWdif)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdif)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdif)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
!
!  Define wave dissipation due to breaking.
!
        IF (VarOut(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdib)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdib))
          ELSE
            Vinfo( 2)=Vname(2,idWdib)
          END IF
          Vinfo( 3)=Vname(3,idWdib)
          Vinfo(14)=Vname(4,idWdib)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdib)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdib,ng),r8)
          S(ng)%pioVar(idWdib)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdib)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdib)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to whitecapping.
!
        IF (VarOut(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdiw))
          ELSE
            Vinfo( 2)=Vname(2,idWdiw)
          END IF
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdiw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiw,ng),r8)
          S(ng)%pioVar(idWdiw)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdiw)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdiw)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef ROLLER_SVENDSEN || defined WW3_MODEL
!
!  Define percent wave breaking.
!
        IF (VarOut(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWbrk))
          ELSE
            Vinfo( 2)=Vname(2,idWbrk)
          END IF
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWbrk)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbrk,ng),r8)
          S(ng)%pioVar(idWbrk)%dkind=PIO_FOUT
          S(ng)%pioVar(idWbrk)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWbrk)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef WEC_ROLLER
!
!  Define wave roller dissipation.
!
        IF (VarOut(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdis))
          ELSE
            Vinfo( 2)=Vname(2,idWdis)
          END IF
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdis,ng),r8)
          S(ng)%pioVar(idWdis)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdis)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdis)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef ROLLER_RENIERS
!
!  Define roller wave action density.
!
        IF (VarOut(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWrol))
          ELSE
            Vinfo( 2)=Vname(2,idWrol)
          END IF
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWrol)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWrol,ng),r8)
          S(ng)%pioVar(idWrol)%dkind=PIO_FOUT
          S(ng)%pioVar(idWrol)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWrol)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Waves Effect on Currents variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
#  ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            S(ng)%pioVar(idU2rs)%vd=var_desc(i)
            S(ng)%pioVar(idU2rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idU2rs)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            S(ng)%pioVar(idV2rs)%vd=var_desc(i)
            S(ng)%pioVar(idV2rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idV2rs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            S(ng)%pioVar(idU2Sd)%vd=var_desc(i)
            S(ng)%pioVar(idU2Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idU2Sd)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            S(ng)%pioVar(idV2Sd)%vd=var_desc(i)
            S(ng)%pioVar(idV2Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idV2Sd)%gtype=r2dvar
#   ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            S(ng)%pioVar(idU3rs)%vd=var_desc(i)
            S(ng)%pioVar(idU3rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idU3rs)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            S(ng)%pioVar(idV3rs)%vd=var_desc(i)
            S(ng)%pioVar(idV3rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idV3rs)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            S(ng)%pioVar(idU3Sd)%vd=var_desc(i)
            S(ng)%pioVar(idU3Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idU3Sd)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            S(ng)%pioVar(idV3Sd)%vd=var_desc(i)
            S(ng)%pioVar(idV3Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idV3Sd)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            S(ng)%pioVar(idW3Sd)%vd=var_desc(i)
            S(ng)%pioVar(idW3Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idW3Sd)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            S(ng)%pioVar(idW3St)%vd=var_desc(i)
            S(ng)%pioVar(idW3St)%dkind=PIO_FOUT
            S(ng)%pioVar(idW3St)%gtype=r2dvar
#   endif
#  endif
#  ifdef WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            S(ng)%pioVar(idWztw)%vd=var_desc(i)
            S(ng)%pioVar(idWztw)%dkind=PIO_FOUT
            S(ng)%pioVar(idWztw)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            S(ng)%pioVar(idWqsp)%vd=var_desc(i)
            S(ng)%pioVar(idWqsp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWqsp)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            S(ng)%pioVar(idWbeh)%vd=var_desc(i)
            S(ng)%pioVar(idWbeh)%dkind=PIO_FOUT
            S(ng)%pioVar(idWbeh)%gtype=r2dvar
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            S(ng)%pioVar(idWdif)%vd=var_desc(i)
            S(ng)%pioVar(idWdif)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdif)%gtype=r2dvar
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            S(ng)%pioVar(idWdib)%vd=var_desc(i)
            S(ng)%pioVar(idWdib)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdib)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            S(ng)%pioVar(idWdiw)%vd=var_desc(i)
            S(ng)%pioVar(idWdiw)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdiw)%gtype=r2dvar
#  endif
#  ifdef ROLLER_SVENDSEN || defined WW3_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            S(ng)%pioVar(idWbrk)%vd=var_desc(i)
            S(ng)%pioVar(idWbrk)%dkind=PIO_FOUT
            S(ng)%pioVar(idWbrk)%gtype=r2dvar
#  endif
#  ifdef WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            S(ng)%pioVar(idWdis)%vd=var_desc(i)
            S(ng)%pioVar(idWdis)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdis)%gtype=r2dvar
#  endif
#  ifdef ROLLER_RENIERS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            S(ng)%pioVar(idWrol)%vd=var_desc(i)
            S(ng)%pioVar(idWrol)%dkind=PIO_FOUT
            S(ng)%pioVar(idWrol)%gtype=r2dvar
#  endif
          END IF
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef WEC
        IF (.not.got_var(idU2rs).and.VarOut(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.VarOut(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2Sd).and.VarOut(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.VarOut(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SOLVE3D
        IF (.not.got_var(idU3Sd).and.VarOut(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.VarOut(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.VarOut(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.VarOut(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.VarOut(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.VarOut(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef WEC_VF
        IF (.not.got_var(idWztw).and.VarOut(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.VarOut(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.VarOut(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbeh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
        IF (.not.got_var(idWdif).and.VarOut(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
        IF (.not.got_var(idWdib).and.VarOut(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.VarOut(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ROLLER_SVENDSEN || defined WW3_MODEL
        IF (.not.got_var(idWbrk).and.VarOut(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WEC_ROLLER
        IF (.not.got_var(idWdis).and.VarOut(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ROLLER_RENIERS
        IF (.not.got_var(idWrol).and.VarOut(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
      END IF QUERY
!
  10  FORMAT (/,' WEC_DEF_PIO  - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE wec_def_pio
!
!*********************************************************************** 
      SUBROUTINE wec_wrt_pio (ng, model, tile,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        VarOut, S)
!*********************************************************************** 
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Linstataneous
!
      integer :: status
!
      real(dp) :: scale
!
#  ifdef SOLVE3D
!
      real(r8), allocatable :: Wr3d(:,:,:)
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents output variables into specified
!  output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set instantaneous fields.
!
      IF ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                              &
     &    (S(ng)%ncid.eq.QCK(ng)%ncid)) THEN
        Linstataneous=.TRUE.
      ELSE
        Linstataneous=.FALSE.                  ! time-averged fiels
      END IF

#  ifdef WEC
!
!  Write out 2D total Waves Effect on Currents U-stress.
!
      IF (Aout(idU2rs,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idU2rs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idU2rs,          &
     &                       S(ng)%pioVar(idU2rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       MIXING(ng) % rustr2d)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idU2rs,          &
     &                       S(ng)%pioVar(idU2rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % umask,                          &
#    endif
     &                       AVERAGE(ng) % avgu2rs)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D total Waves Effect on Currents V-stress.
!
      IF (Aout(idV2rs,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idV2rs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idV2rs,          &
     &                       S(ng)%pioVar(idV2rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       MIXING(ng) % rvstr2d)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idV2rs,          &
     &                       S(ng)%pioVar(idV2rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#    endif
     &                       AVERAGE(ng) % avgv2rs)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Stokes U-velocity.
!
      IF (Aout(idU2Sd,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idU2Sd)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idU2Sd,          &
     &                       S(ng)%pioVar(idU2Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       OCEAN(ng) % ubar_stokes)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idU2Sd,          &
     &                       S(ng)%pioVar(idU2Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % umask,                          &
#    endif
     &                       AVERAGE(ng) % avgu2Sd)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Stokes V-velocity.
!
      IF (Aout(idV2Sd,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idV2Sd)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idV2Sd,          &
     &                       S(ng)%pioVar(idV2Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       OCEAN(ng) % vbar_stokes)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idV2Sd,          &
     &                       S(ng)%pioVar(idV2Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#    endif
     &                       AVERAGE(ng) % avgv2Sd)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SOLVE3D
!
!  Write out 3D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU3rs,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idU3rs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idU3rs,          &
     &                       S(ng)%pioVar(idU3rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask,                          &
#    endif
     &                       MIXING(ng) % rustr3d)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idU3rs,          &
     &                       S(ng)%pioVar(idU3rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#     ifdef MASKING
     &                       GRID(ng) % umask,                          &
#     endif
     &                       AVERAGE(ng) % avgu3rs)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D total Waves Effect on Currents V-stress.
!
      IF (VarOut(idV3rs,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idV3rs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idV3rs,          &
     &                       S(ng)%pioVar(idV3rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#    endif
     &                       MIXING(ng) % rvstr3d)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idV3rs,          &
     &                       S(ng)%pioVar(idV3rs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#     ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#     endif
     &                       AVERAGE(ng) % avgv3rs)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3rs)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes U-velocity.
!
      IF (VarOut(idU3Sd,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idU3Sd)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idU3Sd,          &
     &                       S(ng)%pioVar(idU3Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask,                          &
#    endif
     &                       OCEAN(ng) % u_stokes)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idU3Sd,          &
     &                       S(ng)%pioVar(idU3Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#     ifdef MASKING
     &                       GRID(ng) % umask,                          &
#     endif
     &                       AVERAGE(ng) % avgu3Sd)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes V-velocity.
!
      IF (VarOut(idV3Sd,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idV3Sd)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idV3Sd,          &
     &                       S(ng)%pioVar(idV3Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#    endif
     &                       OCEAN(ng) % v_stokes)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idV3Sd,          &
     &                       S(ng)%pioVar(idV3Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#     ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#     endif
     &                       AVERAGE(ng) % avgv3Sd)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes omega-velocity.
!
      IF (VarOut(idW3Sd,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idW3Sd)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        IF (Linstataneous) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
            Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
          END IF
          CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),     &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % W_stokes,                       &
     &                      Wr3d)
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idW3Sd,          &
     &                       S(ng)%pioVar(idW3Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       Wr3d)
          deallocate (Wr3d)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idW3Sd,          &
     &                       S(ng)%pioVar(idW3Sd),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#     ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#     endif
     &                       AVERAGE(ng) % avgw3d)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3Sd)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Stokes 'true' W-velocity (m/s).
!
      IF (VarOut(idW3St,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idW3St)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idW3St,          &
     &                       S(ng)%pioVar(idW3St),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       OCEAN(ng) % wstvel)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idW3St,          &
     &                       S(ng)%pioVar(idW3St),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#     ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#     endif
     &                       AVERAGE(ng) % avgW3St)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3St)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   ifdef WEC_VF
!
!  Write out Waves Effect on Currents quasi-static sea level adjustment.
!
      IF (VarOut(idWztw,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWztw)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWztw,          &
     &                       S(ng)%pioVar(idWztw),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       OCEAN(ng) % zetaw)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWztw,          &
     &                       S(ng)%pioVar(idWztw),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#     ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#     endif
     &                       AVERAGE(ng) % avgWztw)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWztw)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Waves Effect on Currents quasi-static pressure.
!
      IF (VarOut(idWqsp,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWqsp)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWqsp,          &
     &                       S(ng)%pioVar(idWqsp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &

     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       OCEAN(ng) % qsp)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWqsp,          &
     &                       S(ng)%pioVar(idWqsp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &

     &                       LBi, UBi, LBj, UBj, scale,                 &
#     ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#     endif
     &                       AVERAGE(ng) % avgWqsp)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWqsp)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Waves Effect on Currents Bernoulli head.
!
      IF (VarOut(idWbeh,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWbeh)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWbeh,          &
     &                       S(ng)%pioVar(idWbeh),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       OCEAN(ng) % bh)
#    ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWbeh,          &
     &                       S(ng)%pioVar(idWbeh),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#     ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#     endif
     &                       AVERAGE(ng) % avgWbeh)
#    endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbeh)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
#  if (defined BOTTOM_STREAMING && defined WEC_VF) || \
       defined WAV_COUPLING
!
!  Write out wave dissipation due to bottom friction.
!
      IF (VarOut(idWdif,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        IF (S(ng)%pioVar(idWdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdif,          &
     &                       S(ng)%pioVar(idWdif),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_fric)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdif,          &
     &                       S(ng)%pioVar(idWdif),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdif)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdif)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
      defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
      defined WAVES_DISS     || defined WDISS_INWAVE
!
!  Write out wave dissipation due to breaking.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWdib)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdib,          &
     &                       S(ng)%pioVar(idWdib),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_break)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdib,          &
     &                       S(ng)%pioVar(idWdib),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdib)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdib)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wave dissipation due to whitecapping.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWdiw)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdiw,          &
     &                       S(ng)%pioVar(idWdiw),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_wcap)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdiw,          &
     &                       S(ng)%pioVar(idWdiw),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdiw)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdiw)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef ROLLER_SVENDSEN || defined WW3_MODEL || \
   defined INWAVE_MODEL || defined SWAN_MODEL
!
!  Write out percent wave breaking.
!
      IF (VarOut(idWbrk,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWbrk)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWbrk,          &
     &                       S(ng)%pioVar(idWbrk),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Wave_break)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWbrk,          &
     &                       S(ng)%pioVar(idWbrk),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWbrk)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbrk)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WEC_ROLLER
!
!  Write out wave roller dissipation.
!
      IF (VarOut(idWdis,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWdis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdis,          &
     &                       S(ng)%pioVar(idWdis),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dissip_roller)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdis,          &
     &                       S(ng)%pioVar(idWdis),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdis)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdis)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef ROLLER_RENIERS
!
!  Write out roller wave action density.
!
      IF (VarOut(idWrol,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWrol)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWrol,          &
     &                       S(ng)%pioVar(idWrol),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % rollA)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWrol,          &
     &                       S(ng)%pioVar(idWrol),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWrol)
#   endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWrol)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
  10  FORMAT (/," WEC_WRT_PIO  - error while writing variable '",a,     &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE wec_wrt_pio

#  ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE wec_def_station_pio (ng, model, ldef, VarOut, S,      &
     &                                pgrd, rgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: pgrd(:), rgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, status
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Waves Effect on Currents output stations variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
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

#   ifdef WEC
!
!  Define 2D Stokes U-velocity.
!
        IF (Sout(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          Vinfo( 2)=Vname(2,idU2Sd)
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idU2Sd)%dkind=PIO_FOUT
          STA(ng)%pioVar(idU2Sd)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idU2Sd)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Stokes V-velocity.
!
        IF (Sout(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          Vinfo( 2)=Vname(2,idV2Sd)
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idV2Sd)%dkind=PIO_FOUT
          STA(ng)%pioVar(idV2Sd)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idV2Sd)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D total Waves Effect on Currents U-stress.
!
        IF (Sout(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          Vinfo( 2)=Vname(2,idU2rs)
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idU2rs)%dkind=PIO_FOUT
          STA(ng)%pioVar(idU2rs)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idU2rs)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D total Waves Effect on Currents V-stress.
!
        IF (Sout(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          Vinfo( 2)=Vname(2,idV2rs)
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idV2rs)%dkind=PIO_FOUT
          STA(ng)%pioVar(idV2rs)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idV2rs)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#    ifdef SOLVE3D
!
!  Define 3D Stokes U-velocity.
!
        IF (Sout(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          Vinfo( 2)=Vname(2,idU3Sd)
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idU3Sd)%dkind=PIO_FOUT
          STA(ng)%pioVar(idU3Sd)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idU3Sd)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes V-velocity.
!
        IF (Sout(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          Vinfo( 2)=Vname(2,idV3Sd)
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idV3Sd)%dkind=PIO_FOUT
          STA(ng)%pioVar(idV3Sd)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idV3Sd)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Stokes omega-velocity.
!
        IF (Sout(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          Vinfo( 2)=Vname(2,idW3Sd)
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idW3Sd)%dkind=PIO_FOUT
          STA(ng)%pioVar(idW3Sd)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idW3Sd)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D 'true' Stokes W-velocity.
!
        IF (Sout(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          Vinfo( 2)=Vname(2,idW3St)
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idW3St)%dkind=PIO_FOUT
          STA(ng)%pioVar(idW3St)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idW3St)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D total Waves Effect on Currents U-stress.
!
        IF (Sout(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          Vinfo( 2)=Vname(2,idU3rs)
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idU3rs)%dkind=PIO_FOUT
          STA(ng)%pioVar(idU3rs)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idU3rs)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D total Waves Effect on Current V-stress.
!
        IF (Sout(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          Vinfo( 2)=Vname(2,idV3rs)
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idV3rs)%dkind=PIO_FOUT
          STA(ng)%pioVar(idV3rs)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idV3rs)%vd,                     &
     &                   PIO_FOUT, 3, rgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#    endif
#   endif
#   ifdef WEC_VF
!
!  Define Waves Effect on Currents quasi-static sea level adjustment.
!
        IF (Sout(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          Vinfo( 2)=Vname(2,idWztw)
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWztw)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWztw)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWztw)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents quasi-static pressure.
!
        IF (Sout(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          Vinfo( 2)=Vname(2,idWqsp)
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWqsp)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWqsp)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWqsp)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Waves Effect on Currents Bernoulli head.
!
        IF (Sout(idWbeh,ng)) THEN
          Vinfo( 1)=Vname(1,idWbeh)
          Vinfo( 2)=Vname(2,idWbeh)
          Vinfo( 3)=Vname(3,idWbeh)
          Vinfo(14)=Vname(4,idWbeh)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWbeh)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWbeh)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWbeh)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if (defined BOTTOM_STREAMING && defined WEC_VF) || \
        defined WAV_COUPLING
!
!  Define wave dissipation due to bottom friction.
!
        IF (Sout(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          Vinfo( 2)=Vname(2,idWdif)
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWdif)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWdif)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWdif)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
       defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
       defined WAVES_DISS     || defined WDISS_INWAVE
!
!  Define wave dissipation due to breaking.
!
        IF (Sout(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          Vinfo( 2)=Vname(2,idWdis)
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWdis)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWdis)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWdis)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave dissipation due to whitecapping.
!
        IF (Sout(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          Vinfo( 2)=Vname(2,idWdiw)
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWdiw)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWdiw)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWdiw)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef ROLLER_SVENDSEN || defined WW3_MODEL
!
!  Define percent wave breaking.
!
        IF (Sout(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          Vinfo( 2)=Vname(2,idWbrk)
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWbrk)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWbrk)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWbrk)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef WEC_ROLLER
!
!  Define wave roller dissipation.
!
        IF (Sout(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          Vinfo( 2)=Vname(2,idWdis)
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWdis)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWdis)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWdis)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef ROLLER_RENIERS
!
!  Define wave roller action density.
!
        IF (Sout(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          Vinfo( 2)=Vname(2,idWrol)
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
          STA(ng)%pioVar(idWrol)%dkind=PIO_FOUT
          STA(ng)%pioVar(idWrol)%gtype=0
!
          status=def_var(ng, iNLM, STA(ng)%pioFile,                     &
     &                   STA(ng)%pioVar(idWrol)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  stations variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
#   ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            S(ng)%pioVar(idU2Sd)%vd=var_desc(i)
            S(ng)%pioVar(idU2Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idU2Sd)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            S(ng)%pioVar(idV2Sd)%vd=var_desc(i)
            S(ng)%pioVar(idV2Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idV2Sd)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            S(ng)%pioVar(idU2rs)%vd=var_desc(i)
            S(ng)%pioVar(idU2rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idU2rs)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            S(ng)%pioVar(idV2rs)%vd=var_desc(i)
            S(ng)%pioVar(idV2rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idV2rs)%gtype=0
#    ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            S(ng)%pioVar(idU3Sd)%vd=var_desc(i)
            S(ng)%pioVar(idU3Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idU3Sd)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            S(ng)%pioVar(idV3Sd)%vd=var_desc(i)
            S(ng)%pioVar(idV3Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idV3Sd)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            S(ng)%pioVar(idW3Sd)%vd=var_desc(i)
            S(ng)%pioVar(idW3Sd)%dkind=PIO_FOUT
            S(ng)%pioVar(idW3Sd)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            S(ng)%pioVar(idW3St)%vd=var_desc(i)
            S(ng)%pioVar(idW3St)%dkind=PIO_FOUT
            S(ng)%pioVar(idW3St)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            S(ng)%pioVar(idU3rs)%vd=var_desc(i)
            S(ng)%pioVar(idU3rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idU3rs)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            S(ng)%pioVar(idV3rs)%vd=var_desc(i)
            S(ng)%pioVar(idV3rs)%dkind=PIO_FOUT
            S(ng)%pioVar(idV3rs)%gtype=0
#    endif
#   endif
#   ifdef WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            S(ng)%pioVar(idWztw)%vd=var_desc(i)
            S(ng)%pioVar(idWztw)%dkind=PIO_FOUT
            S(ng)%pioVar(idWztw)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            S(ng)%pioVar(idWqsp)%vd=var_desc(i)
            S(ng)%pioVar(idWqsp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWqsp)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            S(ng)%pioVar(idWbeh)%vd=var_desc(i)
            S(ng)%pioVar(idWbeh)%dkind=PIO_FOUT
            S(ng)%pioVar(idWbeh)%gtype=0
#   endif
#   if (defined BOTTOM_STREAMING && defined WEC_VF) || \
        defined WAV_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            S(ng)%pioVar(idWdif)%vd=var_desc(i)
            S(ng)%pioVar(idWdif)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdif)%gtype=0
#   endif
#   if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
       defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
       defined WAVES_DISS     || defined WDISS_INWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            S(ng)%pioVar(idWdib)%vd=var_desc(i)
            S(ng)%pioVar(idWdib)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdib)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            S(ng)%pioVar(idWdiw)%vd=var_desc(i)
            S(ng)%pioVar(idWdiw)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdiw)%gtype=0
#   endif
#   ifdef ROLLER_SVENDSEN || defined WW3_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            S(ng)%pioVar(idWbrk)%vd=var_desc(i)
            S(ng)%pioVar(idWbrk)%dkind=PIO_FOUT
            S(ng)%pioVar(idWbrk)%gtype=0
#   endif
#   ifdef WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            S(ng)%pioVar(idWdis)%vd=var_desc(i)
            S(ng)%pioVar(idWdis)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdis)%gtype=0
#   endif
#   ifdef ROLLER_RENIERS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            S(ng)%pioVar(idWrol)%vd=var_desc(i)
            S(ng)%pioVar(idWrol)%dkind=PIO_FOUT
            S(ng)%pioVar(idWrol)%gtype=0
#   endif
          END IF
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
#   ifdef WEC
        IF (.not.got_var(idU2Sd).and.Sout(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Sout(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Sout(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Sout(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef SOLVE3D
        IF (.not.got_var(idU3Sd).and.Sout(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Sout(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.Sout(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.Sout(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Sout(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Sout(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#   ifdef WEC_VF
        IF (.not.got_var(idWztw).and.Sout(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.Sout(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.Sout(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbeh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if (defined BOTTOM_STREAMING && defined WEC_VF) || \
        defined WAV_COUPLING
        IF (.not.got_var(idWdif).and.Sout(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
       defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
       defined WAVES_DISS     || defined WDISS_INWAVE
        IF (.not.got_var(idWdib).and.Sout(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.Sout(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef ROLLER_SVENDSEN || defined WW3_MODEL
        IF (.not.got_var(idWbrk).and.Sout(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WEC_ROLLER
        IF (.not.got_var(idWdis).and.Sout(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef ROLLER_RENIERS
        IF (.not.got_var(idWrol).and.Sout(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
      END IF QUERY
!
  10  FORMAT (/,' WEC_DEF_STATION_PIO - unable to find variable: ',     &
     &        a,2x,' in stations NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE wec_def_station_pio
!
!*********************************************************************** 
      SUBROUTINE wec_wrt_station_pio (ng, model, tile,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                VarOut, S)
!*********************************************************************** 
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposR, NposW
      integer :: i, k, np, status
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
#   ifdef SOLVE3D
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
#   endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wec_wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents output variables into specified
!  stations output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#   ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#   else
      Cgrid=.FALSE.
#   endif
!
!  Set positions for generic extraction routine.
!
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
#   ifdef SOLVE3D
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
#   endif
      END DO

#   ifdef WEC
!
!  Write out 2D Stokes U-velocity.
!
      IF (VarOut(idU2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idU2Sd, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % ubar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idU2Sd)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idU2Sd)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Stokes V-velocity.
!
      IF (VarOut(idV2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idV2Sd, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % vbar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idV2Sd)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idV2Sd)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idU2rs, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rustr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idU2rs)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idU2rs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D total Waves Effect on Currents V-stress.
!
      IF (VarOut(idV2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idV2rs, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rvstr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idV2rs)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idV2rs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#    ifdef SOLVE3D
!
!  Write out 3D Stokes U-velocity.
!
      IF (VarOut(idU3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idU3Sd, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % u_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idU3Sd)), rsta,          &
     &                            (/1,1,S(ng)%Rindex/),                 &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idU3Sd)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Stokes V-velocity.
!
      IF (VarOut(idV3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idV3Sd, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % v_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idV3Sd)), rsta,          &
     &                            (/1,1,S(ng)%Rindex/),                 &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idV3Sd)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Stokes omega-velocity.
!
      IF (VarOut(idW3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idW3Sd, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % W_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idW3Sd)), rsta,          &
     &                            (/1,1,S(ng)%Rindex/),                 &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idW3Sd)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D "true" Stokes W-velocity.
!
      IF (VarOut(idW3St,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idW3St, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % wstvel,                  &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idW3St)), rsta,          &
     &                            (/1,1,S(ng)%Rindex/),                 &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idW3St)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D total Waves Effect on Currents U-stress.
!
      IF (VarOut(idU3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, model, Cgrid, idU3rs, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rustr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idU3rs)), rsta,          &
     &                            (/1,1,S(ng)%Rindex/),                 &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idU3rs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D total Waves Effect on Current V-stress.
!
      IF (VarOut(idV3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, model, Cgrid, idV3rs, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rvstr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idV3rs)), rsta,          &
     &                            (/1,1,S(ng)%Rindex/),                 &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idV3rs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#    endif
#   endif
#   ifdef WEC_VF
!
!  Write out Waves Effect on Currents quasi-static sea level adjustment.
!
      IF (VarOut(idWztw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWztw, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % zeta(:,:,KOUT),          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWztw)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWztw)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out Waves Effect on Currents quasi-static pressure.
!
      IF (VarOut(idWqsp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWqsp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % qsp,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWqsp)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWqsp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out Waves Effect on Currents Bernoulli head.
!
      IF (VarOut(idWbeh,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWbeh, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % bh,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWbeh)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWbeh)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if (defined BOTTOM_STREAMING && defined WEC_VF) || \
        defined WAV_COUPLING
!
!  Write out wave dissipation due to bottom friction.
!
      IF (VarOut(idWdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdif, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_fric,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWdif)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWdif)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined TKE_WAVEDISS   || defined WAV_COUPLING   || \
       defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \ 
       defined WAVES_DISS     || defined WDISS_INWAVE
!
!  Write out wave dissipation due to breaking.
!
      IF (VarOut(idWdib,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdib, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_break,           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWdib)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWdib)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to whitecapping.
!
      IF (VarOut(idWdiw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdiw, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_wcap,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWdiw)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWdiw)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   ifdef ROLLER_SVENDSEN || defined WW3_MODEL || \
    defined INWAVE_MODEL || defined SWAN_MODEL
!
!  Write out percent wave breaking.
!
      IF (VarOut(idWbrk,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWbrk, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_break,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWbrk)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWbrk)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   ifdef WEC_ROLLER
!
!  Write out wave roller dissipation.
!
      IF (VarOut(idWdis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWdis, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_roller,          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWdis)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWdis)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   ifdef ROLLER_RENIERS
!
!  Write out wave roller action density.
!
      IF (VarOut(idWrol,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWrol, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % rollA,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWrol)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWrol)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
!
      RETURN
      END SUBROUTINE wec_wrt_station_pio
#  endif
# endif
#endif
!
      END MODULE wec_output_mod
