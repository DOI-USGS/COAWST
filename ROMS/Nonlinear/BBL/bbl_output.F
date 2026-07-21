#include "cppdefs.h"
      MODULE bbl_output_mod

#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes Bottom Boundary Layer model variables    !
!  into output NetCDF files.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef AVERAGES
      USE mod_average
# endif
# ifdef BBL_MODEL
      USE mod_bbl
# endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
# ifdef BBL_MODEL
      USE mod_sedbed
# endif
# ifdef SEDIMENT
      USE mod_sediment
# endif
      USE mod_stepping
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
      PUBLIC :: bbl_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: bbl_def_pio
# endif
# ifdef STATIONS
      PUBLIC :: bbl_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: bbl_def_station_pio
#  endif
# endif
      PUBLIC :: bbl_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: bbl_wrt_pio
# endif
# ifdef STATIONS
      PUBLIC :: bbl_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: bbl_wrt_station_pio
#  endif
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bbl_def_nf90 (ng, model, ldef, VarOut, S,              &
     &                         t2dgrd, u2dgrd, v2dgrd, r3degrd)
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
      integer, intent(in), optional :: r3degrd(:)
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
     &  __FILE__//", bbl_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Bottom Boundary Layer (BBL) and Waves output variables.
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

# if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
     defined WAV_COUPLING
!
!  Define wind-induced bottom orbital velocity.
!
        IF (VarOut(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWorb))
          ELSE
            Vinfo( 2)=Vname(2,idWorb)
          END IF
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWorb)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWorb,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWorb),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef BBL_MODEL
!
!  Define bottom U-current stress.
!
        IF (VarOut(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbrs))
          ELSE
            Vinfo( 2)=Vname(2,idUbrs)
          END IF
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUbrs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbrs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbrs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-current stress.
!
        IF (VarOut(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbrs))
          ELSE
            Vinfo( 2)=Vname(2,idVbrs)
          END IF
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVbrs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbrs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbrs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (VarOut(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbws))
          ELSE
            Vinfo( 2)=Vname(2,idUbws)
          END IF
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUbws)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbws,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbws),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (VarOut(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbws))
          ELSE
            Vinfo( 2)=Vname(2,idVbws)
          END IF
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVbws)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbws,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbws),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (VarOut(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbcs))
          ELSE
            Vinfo( 2)=Vname(2,idUbcs)
          END IF
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUbcs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbcs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbcs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (VarOut(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbcs))
          ELSE
            Vinfo( 2)=Vname(2,idVbcs)
          END IF
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVbcs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbcs,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbcs),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wave and current bottom stress magnitude.
!
        IF (VarOut(idUVwc,ng)) THEN
          Vinfo( 1)=Vname(1,idUVwc)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUVwc))
          ELSE
            Vinfo( 2)=Vname(2,idUVwc)
          END IF
          Vinfo( 3)=Vname(3,idUVwc)
          Vinfo(14)=Vname(4,idUVwc)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUVwc)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUVwc,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUVwc),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (VarOut(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbot))
          ELSE
            Vinfo( 2)=Vname(2,idUbot)
          END IF
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUbot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbot,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbot),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (VarOut(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbot))
          ELSE
            Vinfo( 2)=Vname(2,idVbot)
          END IF
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVbot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbot,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbot),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (VarOut(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbur))
          ELSE
            Vinfo( 2)=Vname(2,idUbur)
          END IF
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUbur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbur,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbur),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (VarOut(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbvr))
          ELSE
            Vinfo( 2)=Vname(2,idVbvr)
          END IF
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVbvr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbvr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbvr),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef UV_KIRBY
!
!  Define U-velocity from Kirby and Chen.
!
        IF (VarOut(idUwav,ng)) THEN
          Vinfo( 1)=Vname(1,idUwav)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUwav))
          ELSE
            Vinfo( 2)=Vname(2,idUwav)
          END IF
          Vinfo( 3)=Vname(3,idUwav)
          Vinfo(14)=Vname(4,idUwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUwav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUwav,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUwav),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define V-velocity from Kirby and Chen.
!
        IF (VarOut(idVwav,ng)) THEN
          Vinfo( 1)=Vname(1,idVwav)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVwav))
          ELSE
            Vinfo( 2)=Vname(2,idVwav)
          END IF
          Vinfo( 3)=Vname(3,idVwav)
          Vinfo(14)=Vname(4,idVwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVwav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVwav,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVwav),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# ifdef UV_BANIHASHEMI
!
!  Define U-velocity from Banihashemi.
!
        IF (VarOut(idUwav,ng)) THEN
          Vinfo( 1)=Vname(1,idUwav)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUwav))
          ELSE
            Vinfo( 2)=Vname(2,idUwav)
          END IF
          Vinfo( 3)=Vname(3,idUwav)
          Vinfo(14)=Vname(4,idUwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUwav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUwav,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUwav),      &
     &                   NF_FOUT, nvd4, r3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define V-velocity from Kirby and Chen.
!
        IF (VarOut(idVwav,ng)) THEN
          Vinfo( 1)=Vname(1,idVwav)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVwav))
          ELSE
            Vinfo( 2)=Vname(2,idVwav)
          END IF
          Vinfo( 3)=Vname(3,idVwav)
          Vinfo(14)=Vname(4,idVwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVwav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVwav,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVwav),      &
     &                   NF_FOUT, nvd4, r3degrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_HEIGHT
!
!  Define wind-induced significant wave height.
!
        IF (VarOut(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWamp))
          ELSE
            Vinfo( 2)=Vname(2,idWamp)
          END IF
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWamp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWamp,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWamp),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

!
!  Define wind-induced significant wave height squared.
!
        IF (VarOut(idWam2,ng)) THEN
          Vinfo( 1)=Vname(1,idWam2)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWam2))
          ELSE
            Vinfo( 2)=Vname(2,idWam2)
          END IF
          Vinfo( 3)=Vname(3,idWam2)
          Vinfo(14)=Vname(4,idWam2)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWam2)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWam2,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWam2),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_LENGTH
!
!  Define wind-induced mean wavelength.
!
        IF (VarOut(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWlen))
          ELSE
            Vinfo( 2)=Vname(2,idWlen)
          END IF
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWlen)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlen,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWlen),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_LENGTHP
!
!  Define wind-induced peak wave wavelength.
!
        IF (VarOut(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWlep))
          ELSE
            Vinfo( 2)=Vname(2,idWlep)
          END IF
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWlep)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlep,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWlep),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_DIR
!
!  Define wind-induced mean wave direction.
!
        IF (VarOut(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdir))
          ELSE
            Vinfo( 2)=Vname(2,idWdir)
          END IF
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdir)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdir,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdir),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_DIRP
!
!  Define wind-induced peak wave direction.
!
        IF (VarOut(idWdip,ng)) THEN
          Vinfo( 1)=Vname(1,idWdip)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdip))
          ELSE
            Vinfo( 2)=Vname(2,idWdip)
          END IF
          Vinfo( 3)=Vname(3,idWdip)
          Vinfo(14)=Vname(4,idWdip)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdip)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdip,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdip),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_TOP_PERIOD
!
!  Define wind-induced surface wave period.
!
        IF (VarOut(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWptp))
          ELSE
            Vinfo( 2)=Vname(2,idWptp)
          END IF
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWptp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWptp,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWptp),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef WAVES_BOT_PERIOD
!
!  Define wind-induced bottom wave period.
!
        IF (VarOut(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWpbt))
          ELSE
            Vinfo( 2)=Vname(2,idWpbt)
          END IF
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWpbt)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWpbt,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWpbt),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# endif

# ifdef WAVES_DSPR
!
!  Define waves directional spreading.
!
        IF (VarOut(idWvds,ng)) THEN
          Vinfo( 1)=Vname(1,idWvds)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvds))
          ELSE
            Vinfo( 2)=Vname(2,idWvds)
          END IF
          Vinfo( 3)=Vname(3,idWvds)
          Vinfo(14)=Vname(4,idWvds)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWvds)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvds,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWvds),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define waves spectrum peakeness.
!
        IF (VarOut(idWvqp,ng)) THEN
          Vinfo( 1)=Vname(1,idWvqp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvqp))
          ELSE
            Vinfo( 2)=Vname(2,idWvqp)
          END IF
          Vinfo( 3)=Vname(3,idWvqp)
          Vinfo(14)=Vname(4,idWvqp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWvqp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvqp,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWvqp),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
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
!  Bottom Boundary Layer model variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
# if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
     defined WAV_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            S(ng)%Vid(idWorb)=var_id(i)
# endif
# ifdef BBL_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            S(ng)%Vid(idUbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            S(ng)%Vid(idVbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            S(ng)%Vid(idUbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            S(ng)%Vid(idVbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            S(ng)%Vid(idUbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            S(ng)%Vid(idVbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            S(ng)%Vid(idUVwc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            S(ng)%Vid(idUbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            S(ng)%Vid(idVbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            S(ng)%Vid(idUbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            S(ng)%Vid(idVbvr)=var_id(i)
# endif
# if defined UV_KIRBY || defined UV_BANIHASHEMI
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwav))) THEN
            got_var(idUwav)=.TRUE.
            S(ng)%Vid(idUwav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwav))) THEN
            got_var(idVwav)=.TRUE.
            S(ng)%Vid(idVwav)=var_id(i)
# endif
# ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            S(ng)%Vid(idWamp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWam2))) THEN
            got_var(idWam2)=.TRUE.
            S(ng)%Vid(idWam2)=var_id(i)
# endif
# ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            S(ng)%Vid(idWlen)=var_id(i)
# endif
# ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            S(ng)%Vid(idWlep)=var_id(i)
# endif
# ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            S(ng)%Vid(idWdir)=var_id(i)
# endif
# ifdef WAVES_DIRP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdip))) THEN
            got_var(idWdip)=.TRUE.
            S(ng)%Vid(idWdip)=var_id(i)
# endif
# ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            S(ng)%Vid(idWptp)=var_id(i)
# endif
# ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            S(ng)%Vid(idWpbt)=var_id(i)
# endif
# ifdef WAVES_DSPR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvds))) THEN
            got_var(idWvds)=.TRUE.
            S(ng)%Vid(idWvds)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvqp))) THEN
            got_var(idWvqp)=.TRUE.
            S(ng)%Vid(idWvqp)=var_id(i)
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
# if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
     defined WAV_COUPLING
        IF (.not.got_var(idWorb).and.VarOut(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.VarOut(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.VarOut(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.VarOut(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.VarOut(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.VarOut(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.VarOut(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVwc).and.VarOut(idUVwc,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUVwc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.VarOut(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.VarOut(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.VarOut(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.VarOut(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined UV_KIRBY || defined UV_BANIHASHEMI
        IF (.not.got_var(idUwav).and.VarOut(idUwav,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwav).and.VarOut(idVwav,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.VarOut(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWam2).and.VarOut(idWam2,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWam2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.VarOut(idWlen,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.VarOut(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.VarOut(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_DIRP
        IF (.not.got_var(idWdip).and.VarOut(idWdip,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdip)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.VarOut(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.VarOut(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_DSPR
        IF (.not.got_var(idWvds).and.VarOut(idWvds,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvds)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvqp).and.VarOut(idWvqp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvqp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
      END IF QUERY
!
  10  FORMAT (/,' BBL_DEF_NF90 - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE bbl_def_nf90

# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE bbl_def_station_nf90 (ng, model, ldef, VarOut, S,      &
     &                                 pgrd, rgrd)
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
     &  __FILE__//", bbl_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment output stations variables.
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

#  ifdef WAVES_UB
!
!  Define wind-induced wave bottom orbital velocity.
!
        IF (VarOut(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          Vinfo( 2)=Vname(2,idWorb)
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWorb),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef BBL_MODEL
!
!  Define current-induced, bottom U-current stress.
!
        IF (VarOut(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          Vinfo( 2)=Vname(2,idUbrs)
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbrs),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define current-induced, bottom V-current stress.
!
        IF (VarOut(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          Vinfo( 2)=Vname(2,idVbrs)
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbrs),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (VarOut(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          Vinfo( 2)=Vname(2,idUbws)
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbws),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (VarOut(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          Vinfo( 2)=Vname(2,idVbws)
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbws),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (VarOut(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          Vinfo( 2)=Vname(2,idUbcs)
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbcs),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (VarOut(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          Vinfo( 2)=Vname(2,idVbcs)
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbcs),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (VarOut(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          Vinfo( 2)=Vname(2,idUbot)
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbot),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (VarOut(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          Vinfo( 2)=Vname(2,idVbot)
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbot),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (VarOut(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          Vinfo( 2)=Vname(2,idUbur)
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idUbur),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (VarOut(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          Vinfo( 2)=Vname(2,idVbvr)
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idVbvr),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_HEIGHT
!
!  Define wind-induced significant wave height.
!
        IF (VarOut(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          Vinfo( 2)=Vname(2,idWamp)
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWamp),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_LENGTH
!
!  Define wind-induced mean wavelenght.
!
        IF (VarOut(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          Vinfo( 2)=Vname(2,idWlen)
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWlen),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_LENGTHP
!
!  Define wind-induced peak wave wavelength.
!
        IF (VarOut(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          Vinfo( 2)=Vname(2,idWlep)
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWlep),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_DIR
!
!  Define wind-induced mean wave direction.
!
        IF (VarOut(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          Vinfo( 2)=Vname(2,idWdir)
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdir),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_DIRP
!
!  Define wind-induced peak wave direction.
!
        IF (VarOut(idWdip,ng)) THEN
          Vinfo( 1)=Vname(1,idWdip)
          Vinfo( 2)=Vname(2,idWdip)
          Vinfo( 3)=Vname(3,idWdip)
          Vinfo(14)=Vname(4,idWdip)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdip),      &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_TOP_PERIOD
!
!  Define wind-induced surface wave period.
!
        IF (VarOut(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          Vinfo( 2)=Vname(2,idWptp)
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWptp),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  endif

#  ifdef WAVES_BOT_PERIOD
!
!  Define wind-induced bottom wave period.
!
        IF (VarOut(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          Vinfo( 2)=Vname(2,idWpbt)
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWpbt),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  endif

#  ifdef WAVES_DSPR
!
!  Define waves directional spreading.
!
        IF (VarOut(idWvds,ng)) THEN
          Vinfo( 1)=Vname(1,idWvds)
          Vinfo( 2)=Vname(2,idWvds)
          Vinfo( 3)=Vname(3,idWvds)
          Vinfo(14)=Vname(4,idWvds)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWvds),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave spectrum peakedness.
!
        IF (VarOut(idWvqp,ng)) THEN
          Vinfo( 1)=Vname(1,idWvqp)
          Vinfo( 2)=Vname(2,idWvqp)
          Vinfo( 3)=Vname(3,idWvqp)
          Vinfo(14)=Vname(4,idWvqp)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWvqp),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
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
!  Bottom Boundary Layer model variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
#  ifdef WAVES_UB
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            S(ng)%Vid(idWorb)=var_id(i)
#  endif
#  ifdef BBL_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            S(ng)%Vid(idUbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            S(ng)%Vid(idVbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            S(ng)%Vid(idUbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            S(ng)%Vid(idVbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            S(ng)%Vid(idUbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            S(ng)%Vid(idVbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            S(ng)%Vid(idUbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            S(ng)%Vid(idVbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            S(ng)%Vid(idUbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            S(ng)%Vid(idVbvr)=var_id(i)
#  endif
#  ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            S(ng)%Vid(idWamp)=var_id(i)
#  endif
#  ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            S(ng)%Vid(idWlen)=var_id(i)
#  endif
#  ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            S(ng)%Vid(idWlep)=var_id(i)
#  endif
#  ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            S(ng)%Vid(idWdir)=var_id(i)
#  endif
#  ifdef WAVES_DIRP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdip))) THEN
            got_var(idWdip)=.TRUE.
            S(ng)%Vid(idWdip)=var_id(i)
#  endif
#  ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            S(ng)%Vid(idWptp)=var_id(i)
#  endif
#  ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            S(ng)%Vid(idWpbt)=var_id(i)
#  endif
#  ifdef WAVES_DSPR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvds))) THEN
            got_var(idWvds)=.TRUE.
            S(ng)%Vid(idWvds)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvqp))) THEN
            got_var(idWvqp)=.TRUE.
            S(ng)%Vid(idWvqp)=var_id(i)
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
#  ifdef WAVES_UB
        IF (.not.got_var(idWorb).and.VarOut(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.VarOut(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.VarOut(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.VarOut(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.VarOut(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.VarOut(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.VarOut(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.VarOut(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.VarOut(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.VarOut(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.VarOut(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.VarOut(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.VarOut(idWlen,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.VarOut(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.VarOut(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_DIRP
        IF (.not.got_var(idWdip).and.VarOut(idWdip,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdip)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.VarOut(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.VarOut(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_DSPR
        IF (.not.got_var(idWvds).and.VarOut(idWvds,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvds)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvqp).and.VarOut(idWvqp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvqp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
      END IF QUERY
!
  10  FORMAT (/,' BBL_DEF_STATION_NF90 - unable to find variable:',     &
     &        1x,a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE bbl_def_station_nf90
# endif
!
!***********************************************************************
      SUBROUTINE bbl_wrt_nf90 (ng, model, tile,                         &
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
!
      real(r8), allocatable :: wrk2d(:,:)
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", bbl_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Bottom Boundary Layer model output variables into specified
!  output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
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

# if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
     defined WAV_COUPLING
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (VarOut(idWorb,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWorb,             &
     &                       S(ng)%Vid(idWorb),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Uwave_rms)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWorb,             &
     &                       S(ng)%Vid(idWorb),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWorb)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbrs,             &
     &                       S(ng)%Vid(idUbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % bustrc)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbrs,             &
     &                       S(ng)%Vid(idUbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUbrs)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbrs,             &
     &                       S(ng)%Vid(idVbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % bvstrc)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbrs,             &
     &                       S(ng)%Vid(idVbrs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVbrs)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbws,             &
     &                       S(ng)%Vid(idUbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % bustrw)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbws,             &
     &                       S(ng)%Vid(idUbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUbws)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbws,             &
     &                       S(ng)%Vid(idVbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % bvstrw)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbws,             &
     &                       S(ng)%Vid(idVbws),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVbws)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbcs,             &
     &                       S(ng)%Vid(idUbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % bustrcwmax)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbcs,             &
     &                       S(ng)%Vid(idUbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUbcs)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbcs,             &
     &                       S(ng)%Vid(idVbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % bvstrcwmax)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbcs,             &
     &                       S(ng)%Vid(idVbcs),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVbcs)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wave and current bottom stress magnitude.
!
      IF (VarOut(idUVwc,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          IF (.not.allocated(wrk2d)) THEN
            allocate ( wrk2d(LBi:UBi, LBj:UBj) )
            wrk2d(LBi:UBi,LBj:UBj)=0.0_r8
          END IF
          wrk2d=SQRT(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+             &
     &               BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
!
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUVwc,             &
     &                       S(ng)%Vid(idUVwc),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       wrk2d)
          deallocate (wrk2d)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUVwc,             &
     &                       S(ng)%Vid(idUVwc),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUVwc)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity at RHO-points.
!
      IF (VarOut(idUbot,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbot,             &
     &                       S(ng)%Vid(idUbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % Ubot)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbot,             &
     &                       S(ng)%Vid(idUbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUbot)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity at RHO-points
!
      IF (VarOut(idVbot,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbot,             &
     &                       S(ng)%Vid(idVbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % Vbot)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbot,             &
     &                       S(ng)%Vid(idVbot),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVbot)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-velocity above bed at RHO-points.
!
      IF (VarOut(idUbur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbur,             &
     &                       S(ng)%Vid(idUbur),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % Ur)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUbur,             &
     &                       S(ng)%Vid(idUbur),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUbur)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-velocity above bed at RHO-points.
!
      IF (VarOut(idVbvr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbvr,             &
     &                       S(ng)%Vid(idVbvr),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       BBL(ng) % Vr)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVbvr,             &
     &                       S(ng)%Vid(idVbvr),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVbvr)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# if defined UV_KIRBY
!
!  Write out U-velocity from Kirby and Chen.
!
      IF (VarOut(idUwav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUwav,             &
     &                       S(ng)%Vid(idUwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) %uwave )
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idUwav,             &
     &                       S(ng)%Vid(idUwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUwav)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out V-velocity from Kirby and Chen.
!
      IF (VarOut(idVwav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVwav,             &
     &                       S(ng)%Vid(idVwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) %vwave )
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idVwav,             &
     &                       S(ng)%Vid(idVwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgVwav)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# if defined UV_BANIHASHEMI
!
!  Write out U-velocity from Banihashemi.
!
      IF (VarOut(idUwav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idUwav,             &
     &                       S(ng)%Vid(idUwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, MSCs, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) %uwave )
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idUwav,             &
     &                       S(ng)%Vid(idUwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, MSCs,  scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUwav)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out V-velocity from Banihashemi.
!
      IF (VarOut(idVwav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idVwav,             &
     &                       S(ng)%Vid(idVwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, MSCS, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) %vwave )
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idVwav,             &
     &                       S(ng)%Vid(idVwav),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, MSCS, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgVwav)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_HEIGHT
!
!  Write out wind-induced signiticant wave height.
!
      IF (VarOut(idWamp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWamp,             &
     &                       S(ng)%Vid(idWamp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Hwave)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWamp,             &
     &                       S(ng)%Vid(idWamp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWamp)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

!
!  Write out wind-induced significant wave height squared.
!
      IF (VarOut(idWam2,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          IF (.not.allocated(wrk2d)) THEN
            allocate ( wrk2d(LBi:UBi,LBj:UBj) )
            wrk2d(LBi:UBi,LBj:UBj)=0.0_r8
          END IF
          wrk2d=FORCES(ng)%Hwave*FORCES(ng)%Hwave
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWam2,             &
     &                       S(ng)%Vid(idWam2),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       wrk2d)
          deallocate (wrk2d)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWam2,             &
     &                       S(ng)%Vid(idWam2),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWam2)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelength.
!
      IF (VarOut(idWlen,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWlen,             &
     &                       S(ng)%Vid(idWlen),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Lwave)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWlen,             &
     &                       S(ng)%Vid(idWlen),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWlen)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (VarOut(idWlep,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWlep,             &
     &                       S(ng)%Vid(idWlep),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Lwavep)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWlep,             &
     &                       S(ng)%Vid(idWlep),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWlep)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlep)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (VarOut(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdir,             &
     &                       S(ng)%Vid(idWdir),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Dwave)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdir,             &
     &                       S(ng)%Vid(idWdir),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWdir)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (VarOut(idWdip,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdip,             &
     &                       S(ng)%Vid(idWdip),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Dwavep)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdip,             &
     &                       S(ng)%Vid(idWdip),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWdip)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdip)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (VarOut(idWptp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWptp,             &
     &                       S(ng)%Vid(idWptp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Pwave_top)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWptp,             &
     &                       S(ng)%Vid(idWptp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWptp)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (VarOut(idWpbt,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWpbt,             &
     &                       S(ng)%Vid(idWpbt),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Pwave_bot)
#  ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWpbt,             &
     &                       S(ng)%Vid(idWpbt),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgWpbt)
#  endif
        END IF
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), S(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef WAVES_DSPR
!
!  Write out waves directional spreading.
!
      IF (VarOut(idWvds,ng)) THEN
        IF (Linstataneous) THEN
          scale=1.0_dp
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWvds,             &
     &                       S(ng)%Vid(idWvds),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Wave_ds)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idWvds)), S(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out waves spectrum peakeness.
!
      IF (VarOut(idWvqp,ng)) THEN
        IF (Linstataneous) THEN
          scale=1.0_dp
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idWvqp,             &
     &                       S(ng)%Vid(idWvqp),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % Wave_qp)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idWvqp)), S(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
# endif
!
  10  FORMAT (/," BBL_WRT_NF90 - error while writing variable '",       &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE bbl_wrt_nf90

# ifdef STATIONS
!
!***********************************************************************
     SUBROUTINE bbl_wrt_station_nf90 (ng, model, tile,                  &
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
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", bbl_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment output variables into specified stations NetCDF
!  file.
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
#  endif
      END DO

#  ifdef WAVES_UB
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (VarOut(idWorb,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWorb, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Uwave_rms,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWorb)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWorb))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (VarOut(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idUbrs)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idUbrs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (VarOut(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idVbrs)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idVbrs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (VarOut(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idUbws)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idUbws))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (VarOut(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idVbws)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idVbws))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (VarOut(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idUbcs)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idUbcs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (VarOut(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idVbcs)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idVbcs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (VarOut(idUbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idUbot)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idUbot))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (VarOut(idVbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idVbot)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idVbot))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (VarOut(idUbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idUbur)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idUbur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (VarOut(idVbvr,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbvr, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idVbvr)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idVbvr))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_HEIGHT
!
!  Write out wind-induced significant wave height.
!
      IF (VarOut(idWamp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWamp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Hwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWamp)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWamp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelenght.
!
      IF (VarOut(idWlen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlen, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWlen)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWlen))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (VarOut(idWlep,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlep, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWlep)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWlep))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (VarOut(idWdir,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdir, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdir)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdir))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (VarOut(idWdip,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdip, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWdip)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWdip))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (VarOut(idWptp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWptp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_top,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWptp)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWptp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (VarOut(idWpbt,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWpbt, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_bot,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWpbt)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWpbt))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  if defined WAVES_DSPR
!
!  Write out waves directional spreading.
!
      IF (VarOut(idWvds,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvds, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_ds,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWvds)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWvds))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave spectrum peakedness.
!
      IF (VarOut(idWvqp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvqp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_qp,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idWvqp)), psta,              &
     &                        (/1,S(ng)%Rindex/), (/Nstation(ng),1/),   &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idWvqp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
!
      RETURN
      END SUBROUTINE bbl_wrt_station_nf90
# endif

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE bbl_def_pio (ng, model, ldef, VarOut, S,               &
     &                        t2dgrd, u2dgrd, v2dgrd)
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
     &  __FILE__//", bbl_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define Bottom Boundary Layer model output variables.
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

#  if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
      defined WAV_COUPLING
!
!  Define wind-induced bottom orbital velocity.
!
        IF (VarOut(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWorb))
          ELSE
            Vinfo( 2)=Vname(2,idWorb)
          END IF
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWorb)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWorb,ng),r8)
          S(ng)%pioVar(idWorb)%dkind=PIO_FOUT
          S(ng)%pioVar(idWorb)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWorb)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef BBL_MODEL
!
!  Define bottom U-current stress.
!
        IF (VarOut(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbrs))
          ELSE
            Vinfo( 2)=Vname(2,idUbrs)
          END IF
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUbrs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbrs,ng),r8)
          S(ng)%pioVar(idUbrs)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbrs)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbrs)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-current stress.
!
        IF (VarOut(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbrs))
          ELSE
            Vinfo( 2)=Vname(2,idVbrs)
          END IF
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVbrs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbrs,ng),r8)
          S(ng)%pioVar(idVbrs)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbrs)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbrs)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (VarOut(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbws))
          ELSE
            Vinfo( 2)=Vname(2,idUbws)
          END IF
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUbws)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbws,ng),r8)
          S(ng)%pioVar(idUbws)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbws)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbws)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (VarOut(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbws))
          ELSE
            Vinfo( 2)=Vname(2,idVbws)
          END IF
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVbws)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbws,ng),r8)
          S(ng)%pioVar(idVbws)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbws)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbws)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (VarOut(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbcs))
          ELSE
            Vinfo( 2)=Vname(2,idUbcs)
          END IF
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUbcs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbcs,ng),r8)
          S(ng)%pioVar(idUbcs)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbcs)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbcs)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (VarOut(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbcs))
          ELSE
            Vinfo( 2)=Vname(2,idVbcs)
          END IF
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVbcs)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbcs,ng),r8)
          S(ng)%pioVar(idVbcs)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbcs)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbcs)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wave and current bottom stress magnitude.
!
        IF (VarOut(idUVwc,ng)) THEN
          Vinfo( 1)=Vname(1,idUVwc)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUVwc))
          ELSE
            Vinfo( 2)=Vname(2,idUVwc)
          END IF
          Vinfo( 3)=Vname(3,idUVwc)
          Vinfo(14)=Vname(4,idUVwc)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUVwc)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUVwc,ng),r8)
          S(ng)%pioVar(idUVwc)%dkind=PIO_FOUT
          S(ng)%pioVar(idUVwc)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUVwc)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (VarOut(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbot))
          ELSE
            Vinfo( 2)=Vname(2,idUbot)
          END IF
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUbot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbot,ng),r8)
          S(ng)%pioVar(idUbot)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbot)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbot)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (VarOut(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbot))
          ELSE
            Vinfo( 2)=Vname(2,idVbot)
          END IF
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVbot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbot,ng),r8)
          S(ng)%pioVar(idVbot)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbot)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbot)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (VarOut(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbur))
          ELSE
            Vinfo( 2)=Vname(2,idUbur)
          END IF
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUbur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbur,ng),r8)
          S(ng)%pioVar(idUbur)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbur)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbur)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (VarOut(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbvr))
          ELSE
            Vinfo( 2)=Vname(2,idVbvr)
          END IF
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVbvr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbvr,ng),r8)
          S(ng)%pioVar(idVbvr)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbvr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbvr)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined UV_KIRBY && defined AVERAGES
!
!  Define U-velocity from Kirby and Chen.
!
        IF (VarOut(idUwav,ng)) THEN
          Vinfo( 1)=Vname(1,idUwav)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUwav))
          Vinfo( 3)=Vname(3,idUwav)
          Vinfo(14)=Vname(4,idUwav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUwav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUwav,ng),r8)
          AVG(ng)%pioVar(idUwav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUwav)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUwav)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define V-velocity from Kirby and Chen.
!
        IF (VarOut(idVwav,ng)) THEN
          Vinfo( 1)=Vname(1,idVwav)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVwav))
          Vinfo( 3)=Vname(3,idVwav)
          Vinfo(14)=Vname(4,idVwav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVwav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVwav,ng),r8)
          AVG(ng)%pioVar(idVwav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVwav)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVwav)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_HEIGHT
!
!  Define wind-induced significant wave height.
!
        IF (VarOut(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWamp))
          ELSE
            Vinfo( 2)=Vname(2,idWamp)
          END IF
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWamp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWamp,ng),r8)
          S(ng)%pioVar(idWamp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWamp)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWamp)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef AVERAGES
!
!  Define wind-induced significant wave height squared.
!
        IF (VarOut(idWam2,ng)) THEN
          Vinfo( 1)=Vname(1,idWam2)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWam2))
          Vinfo( 3)=Vname(3,idWam2)
          Vinfo(14)=Vname(4,idWam2)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idWam2)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWam2,ng),r8)
          AVG(ng)%pioVar(idWam2)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idWam2)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWam2)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif

#  ifdef WAVES_LENGTH
!
!  Define wind-induced mean wavelength.
!
        IF (VarOut(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWlen))
          ELSE
            Vinfo( 2)=Vname(2,idWlen)
          END IF
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWlen)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlen,ng),r8)
          S(ng)%pioVar(idWlen)%dkind=PIO_FOUT
          S(ng)%pioVar(idWlen)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWlen)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_LENGTHP
!
!  Define wind-induced peak wave wavelength.
!
        IF (VarOut(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWlep))
          ELSE
            Vinfo( 2)=Vname(2,idWlep)
          END IF
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWlep)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlep,ng),r8)
          S(ng)%pioVar(idWlep)%dkind=PIO_FOUT
          S(ng)%pioVar(idWlep)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWlep)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_DIR
!
!  Define wind-induced mean wave direction.
!
        IF (VarOut(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdir))
          ELSE
            Vinfo( 2)=Vname(2,idWdir)
          END IF
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdir)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdir,ng),r8)
          S(ng)%pioVar(idWdir)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdir)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdir)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_DIRP
!
!  Define wind-induced peak wave direction.
!
        IF (VarOut(idWdip,ng)) THEN
          Vinfo( 1)=Vname(1,idWdip)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdip))
          ELSE
            Vinfo( 2)=Vname(2,idWdip)
          END IF
          Vinfo( 3)=Vname(3,idWdip)
          Vinfo(14)=Vname(4,idWdip)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdip)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdip,ng),r8)
          S(ng)%pioVar(idWdip)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdip)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdip)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_TOP_PERIOD
!
!  Define wind-induced surface wave period.
!
        IF (VarOut(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWptp))
          ELSE
            Vinfo( 2)=Vname(2,idWptp)
          END IF
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWptp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWptp,ng),r8)
          S(ng)%pioVar(idWptp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWptp)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWptp)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef WAVES_BOT_PERIOD
!
!  Define wind-induced bottom wave period.
!
        IF (VarOut(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWpbt))
          ELSE
            Vinfo( 2)=Vname(2,idWpbt)
          END IF
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWpbt)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWpbt,ng),r8)
          S(ng)%pioVar(idWpbt)%dkind=PIO_FOUT
          S(ng)%pioVar(idWpbt)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWpbt)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  endif

#  ifdef WAVES_DSPR
!
!  Define waves directional spreading.
!
        IF (VarOut(idWvds,ng)) THEN
          Vinfo( 1)=Vname(1,idWvds)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvds))
          ELSE
            Vinfo( 2)=Vname(2,idWvds)
          END IF
          Vinfo( 3)=Vname(3,idWvds)
          Vinfo(14)=Vname(4,idWvds)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWvds)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvds,ng),r8)
          S(ng)%pioVar(idWvds)%dkind=PIO_FOUT
          S(ng)%pioVar(idWvds)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWvds)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define waves spectrum peakeness.
!
        IF (VarOut(idWvqp,ng)) THEN
          Vinfo( 1)=Vname(1,idWvqp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvqp))
          ELSE
            Vinfo( 2)=Vname(2,idWvqp)
          END IF
          Vinfo( 3)=Vname(3,idWvqp)
          Vinfo(14)=Vname(4,idWvqp)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWvqp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvqp,ng),r8)
          S(ng)%pioVar(idWvqp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWvqp)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWvqp)%vd,                       &
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
!  Bottom Boundary Layer model variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
#  if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
      defined WAV_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            S(ng)%pioVar(idWorb)%vd=var_desc(i)
            S(ng)%pioVar(idWorb)%dkind=PIO_FOUT
            S(ng)%pioVar(idWorb)%gtype=r2dvar
#  endif
#  ifdef BBL_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            S(ng)%pioVar(idUbrs)%vd=var_desc(i)
            S(ng)%pioVar(idUbrs)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbrs)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            S(ng)%pioVar(idVbrs)%vd=var_desc(i)
            S(ng)%pioVar(idVbrs)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbrs)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            S(ng)%pioVar(idUbws)%vd=var_desc(i)
            S(ng)%pioVar(idUbws)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbws)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            S(ng)%pioVar(idVbws)%vd=var_desc(i)
            S(ng)%pioVar(idVbws)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbws)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            S(ng)%pioVar(idUbcs)%vd=var_desc(i)
            S(ng)%pioVar(idUbcs)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbcs)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            S(ng)%pioVar(idVbcs)%vd=var_desc(i)
            S(ng)%pioVar(idVbcs)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbcs)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            S(ng)%pioVar(idUVwc)%vd=var_desc(i)
            S(ng)%pioVar(idUVwc)%dkind=PIO_FOUT
            S(ng)%pioVar(idUVwc)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            S(ng)%pioVar(idUbot)%vd=var_desc(i)
            S(ng)%pioVar(idUbot)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbot)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            S(ng)%pioVar(idVbot)%vd=var_desc(i)
            S(ng)%pioVar(idVbot)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbot)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            S(ng)%pioVar(idUbur)%vd=var_desc(i)
            S(ng)%pioVar(idUbur)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            S(ng)%pioVar(idVbvr)%vd=var_desc(i)
            S(ng)%pioVar(idVbvr)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbvr)%gtype=r2dvar
#  endif
#  if defined UV_KIRBY && defined AVERAGES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwav))) THEN
            got_var(idUwav)=.TRUE.
            S(ng)%pioVar(idUwav)%vd=var_desc(i)
            S(ng)%pioVar(idUwav)%dkind=PIO_FOUT
            S(ng)%pioVar(idUwav)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwav))) THEN
            got_var(idVwav)=.TRUE.
            S(ng)%pioVar(idVwav)%vd=var_desc(i)
            S(ng)%pioVar(idVwav)%dkind=PIO_FOUT
            S(ng)%pioVar(idVwav)%gtype=r2dvar
#  endif
#  ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            S(ng)%pioVar(idWamp)%vd=var_desc(i)
            S(ng)%pioVar(idWamp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWamp)%gtype=r2dvar
#   ifdef AVERAGES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWam2))) THEN
            got_var(idWam2)=.TRUE.
            AVG(ng)%pioVar(idWam2)%vd=var_desc(i)
            AVG(ng)%pioVar(idWam2)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idWam2)%gtype=r2dvar
#   endif
#  endif
#  ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            S(ng)%pioVar(idWlen)%vd=var_desc(i)
            S(ng)%pioVar(idWlen)%dkind=PIO_FOUT
            S(ng)%pioVar(idWlen)%gtype=r2dvar
#  endif
#  ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            S(ng)%pioVar(idWlep)%vd=var_desc(i)
            S(ng)%pioVar(idWlep)%dkind=PIO_FOUT
            S(ng)%pioVar(idWlep)%gtype=r2dvar
#  endif
#  ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            S(ng)%pioVar(idWdir)%vd=var_desc(i)
            S(ng)%pioVar(idWdir)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdir)%gtype=r2dvar
#  endif
#  ifdef WAVES_DIRP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdip))) THEN
            got_var(idWdip)=.TRUE.
            S(ng)%pioVar(idWdip)%vd=var_desc(i)
            S(ng)%pioVar(idWdip)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdip)%gtype=r2dvar
#  endif
#  ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            S(ng)%pioVar(idWptp)%vd=var_desc(i)
            S(ng)%pioVar(idWptp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWptp)%gtype=r2dvar
#  endif
#  ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            S(ng)%pioVar(idWpbt)%vd=var_desc(i)
            S(ng)%pioVar(idWpbt)%dkind=PIO_FOUT
            S(ng)%pioVar(idWpbt)%gtype=r2dvar
#  endif
#  ifdef WAVES_DSPR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvds))) THEN
            got_var(idWvds)=.TRUE.
            S(ng)%pioVar(idWvds)%vd=var_desc(i)
            S(ng)%pioVar(idWvds)%dkind=PIO_FOUT
            S(ng)%pioVar(idWvds)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvqp))) THEN
            got_var(idWvqp)=.TRUE.
            S(ng)%pioVar(idWvqp)%vd=var_desc(i)
            S(ng)%pioVar(idWvqp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWvqp)%gtype=r2dvar
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
#  if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
      defined WAV_COUPLING
        IF (.not.got_var(idWorb).and.VarOut(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.VarOut(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.VarOut(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.VarOut(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.VarOut(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.VarOut(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.VarOut(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVwc).and.VarOut(idUVwc,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUVwc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.VarOut(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.VarOut(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.VarOut(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.VarOut(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined UV_KIRBY && defined AVERAGES
        IF (.not.got_var(idUwav).and.VarOut(idUwav,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwav).and.VarOut(idVwav,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.VarOut(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef AVERAGES
        IF (.not.got_var(idWam2).and.Aout(idWam2,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWam2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.VarOut(idWlen,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.VarOut(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.VarOut(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_DIRP
        IF (.not.got_var(idWdip).and.VarOut(idWdip,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdip)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.VarOut(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.VarOut(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WAVES_DSPR
        IF (.not.got_var(idWvds).and.VarOut(idWvds,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvds)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvqp).and.VarOut(idWvqp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvqp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
      END IF QUERY
!
  10  FORMAT (/,' BBL_DEF_PIO  - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE bbl_def_pio
!
!***********************************************************************
      SUBROUTINE bbl_wrt_pio (ng, model, tile,                          &
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
      real(r8), allocatable :: wrk2d(:,:)
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", bbl_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment output variables into specified NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set instantaneous fields.
!
      IF ((S(ng)%ncid.eq.S(ng)%ncid).or.                              &
     &    (S(ng)%ncid.eq.QCK(ng)%ncid)) THEN
        Linstataneous=.TRUE.
      ELSE
        Linstataneous=.FALSE.                  ! time-averged fiels
      END IF

#  if defined BBL_MODEL || defined BEDLOAD_VANDERA || \
      defined WAV_COUPLING
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (VarOut(idWorb,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWorb)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWorb,          &
     &                       S(ng)%pioVar(idWorb),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Uwave_rms)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWorb,          &
     &                       S(ng)%pioVar(idWorb),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWorb)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbrs,ng)) THEN
        scale=-rho0
        IF (S(ng)%pioVar(idUbrs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbrs,          &
     &                       S(ng)%pioVar(idUbrs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % bustrc)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbrs,          &
     &                       S(ng)%pioVar(idUbrs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgUbrs)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbrs,ng)) THEN
        scale=-rho0
        IF (S(ng)%pioVar(idVbrs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbrs,          &
     &                       S(ng)%pioVar(idVbrs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % bvstrc)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbrs,          &
     &                       S(ng)%pioVar(idVbrs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgVbrs)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbws,ng)) THEN
        scale=rho0
        IF (S(ng)%pioVar(idUbws)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbws,          &
     &                       S(ng)%pioVar(idUbws),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % bustrw)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbws,          &
     &                       S(ng)%pioVar(idUbws),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgUbws)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbws,ng)) THEN
        scale=rho0
        IF (S(ng)%pioVar(idVbws)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbws,          &
     &                       S(ng)%pioVar(idVbws),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % bvstrw)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbws,          &
     &                       S(ng)%pioVar(idVbws),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgVbws)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-stress at RHO-points.
!
      IF (VarOut(idUbcs,ng)) THEN
        scale=rho0
        IF (S(ng)%pioVar(idUbcs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbcs,          &
     &                       S(ng)%pioVar(idUbcs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % bustrcwmax)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbcs,          &
     &                       S(ng)%pioVar(idUbcs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgUbcs)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-stress at RHO-points.
!
      IF (VarOut(idVbcs,ng)) THEN
        scale=rho0
        IF (S(ng)%pioVar(idVbcs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbcs,          &
     &                       S(ng)%pioVar(idVbcs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % bvstrcwmax)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbcs,          &
     &                       S(ng)%pioVar(idVbcs),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgVbcs)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wave and current bottom stress magnitude.
!
      IF (VarOut(idUVwc,ng)) THEN
        scale=rho0
        IF (S(ng)%pioVar(idUVwc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          IF (.not.allocated(wrk2d)) THEN
            allocate ( wrk2d(LBi:UBi, LBj:UBj) )
            wrk2d(LBi:UBi,LBj:UBj)=0.0_r8
          END IF
          wrk2d=SQRT(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+             &
     &               BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
!
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUVwc,          &
     &                       S(ng)%pioVar(idUVwc),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       wrk2d)
          deallocate (wrk2d)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUVwc,          &
     &                       S(ng)%pioVar(idUVwc),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgUVwc)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity at RHO-points.
!
      IF (VarOut(idUbot,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idUbot)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbot,          &
     &                       S(ng)%pioVar(idUbot),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % Ubot)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbot,          &
     &                       S(ng)%pioVar(idUbot),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgUbot)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity at RHO-points
!
      IF (VarOut(idVbot,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idVbot)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbot,          &
     &                       S(ng)%pioVar(idVbot),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % Vbot)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbot,          &
     &                       S(ng)%pioVar(idVbot),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgVbot)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-velocity above bed at RHO-points.
!
      IF (VarOut(idUbur,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idUbur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbur,          &
     &                       S(ng)%pioVar(idUbur),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % Ur)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUbur,          &
     &                       S(ng)%pioVar(idUbur),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgUbur)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-velocity above bed at RHO-points.
!
      IF (VarOut(idVbvr,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idVbvr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbvr,          &
     &                       S(ng)%pioVar(idVbvr),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       BBL(ng) % Vr)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVbvr,          &
     &                       S(ng)%pioVar(idVbvr),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgVbvr)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  if defined UV_KIRBY && defined AVERAGES
!
!  Write out U-velocity from Kirby and Chen.
!
      IF (VarOut(idUwav,ng)) THEN
        IF (.not.Linstataneous) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idUwav)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUwav,          &
     &                       S(ng)%pioVar(idUwav),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgUwav)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUwav)), S(ng)%Rindex,    &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out V-velocity from Kirby and Chen.
!
      IF (VarOut(idVwav,ng)) THEN
        IF (.not.Linstataneous) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idVwav)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idVwav,          &
     &                       S(ng)%pioVar(idVwav),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVwav)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVwav)), S(ng)%Rindex,    &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  endif

#  ifdef WAVES_HEIGHT
!
!  Write out wind-induced sifnificant wave height.
!
      IF (VarOut(idWamp,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWamp)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWamp,          &
     &                       S(ng)%pioVar(idWamp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Hwave)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWamp,          &
     &                       S(ng)%pioVar(idWamp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWamp)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef AVERAGES
!
!  Write out wind-induced significant wave height squared.
!
      IF (VarOut(idWam2,ng)) THEN
        IF (.not.Linstataneous) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idWam2)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWam2,          &
     &                       S(ng)%pioVar(idWam2),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWam2)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idWam2)), AVG(ng)%Rindex,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#   endif
#  endif

#  ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelength.
!
      IF (VarOut(idWlen,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWlen)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWlen,          &
     &                       S(ng)%pioVar(idWlen),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Lwave)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWlen,          &
     &                       S(ng)%pioVar(idWlen),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWlen)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (VarOut(idWlep,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWlen)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWlep,          &
     &                       S(ng)%pioVar(idWlep),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Lwavep)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWlep,          &
     &                       S(ng)%pioVar(idWlep),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWlep)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlep)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (VarOut(idWdir,ng)) THEN
        scale=rad2deg
        IF (S(ng)%pioVar(idWdir)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdir,          &
     &                       S(ng)%pioVar(idWdir),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dwave)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdir,          &
     &                       S(ng)%pioVar(idWdir),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdir)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (VarOut(idWdip,ng)) THEN
        scale=rad2deg
        IF (S(ng)%pioVar(idWdir)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdip,          &
     &                       S(ng)%pioVar(idWdip),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Dwavep)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdip,          &
     &                       S(ng)%pioVar(idWdip),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWdip)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdip)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (VarOut(idWptp,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWptp)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWptp,          &
     &                       S(ng)%pioVar(idWptp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Pwave_top)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWptp,          &
     &                       S(ng)%pioVar(idWptp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWptp)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (VarOut(idWpbt,ng)) THEN
        scale=1.0_dp
        IF (S(ng)%pioVar(idWpbt)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        IF (Linstataneous) THEN
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWpbt,          &
     &                       S(ng)%pioVar(idWpbt),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Pwave_bot)
#   ifdef AVERAGES
        ELSE
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWpbt,          &
     &                       S(ng)%pioVar(idWpbt),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       AVERAGE(ng) % avgWpbt)
#   endif
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef WAVES_DSPR
!
!  Write out waves directional spreading.
!
      IF (VarOut(idWvds,ng)) THEN
        IF (Linstataneous) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idWvds)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWvds,          &
     &                       S(ng)%pioVar(idWvds),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Wave_ds)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idWvds)), S(ng)%Rindex,    &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out waves spectrum peakeness.
!
      IF (VarOut(idWvqp,ng)) THEN
        IF (Linstataneous) THEN
          scale=1.0_dp
          IF (S(ng)%pioVar(idWvds)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWvqp,          &
     &                       S(ng)%pioVar(idWvqp),                      &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % Wave_qp)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idWvqp)), S(ng)%Rindex,    &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  endif
!
  10  FORMAT (/," BBL_WRT_PIO  - error while writing variable '",       &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE bbl_wrt_pio

#  ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE bbl_def_station_pio (ng, model, ldef, VarOut, S,       &
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
     &  __FILE__//", bbl_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment output stations variables.
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

#   ifdef WAVES_UB
!
!  Define wind-induced wave bottom orbital velocity.
!
        IF (VarOut(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          Vinfo( 2)=Vname(2,idWorb)
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWorb)%dkind=PIO_FOUT
          S(ng)%pioVar(idWorb)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWorb)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef BBL_MODEL
!
!  Define bottom U-current stress.
!
        IF (VarOut(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          Vinfo( 2)=Vname(2,idUbrs)
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idUbrs)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbrs)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbrs)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-current stress.
!
        IF (VarOut(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          Vinfo( 2)=Vname(2,idVbrs)
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idVbrs)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbrs)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbrs)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (VarOut(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          Vinfo( 2)=Vname(2,idUbws)
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idUbws)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbws)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbws)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (VarOut(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          Vinfo( 2)=Vname(2,idVbws)
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idVbws)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbws)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbws)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom U-wave stress.
!
        IF (VarOut(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          Vinfo( 2)=Vname(2,idUbcs)
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idUbcs)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbcs)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbcs)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define maximum wind and current, bottom V-wave stress.
!
        IF (VarOut(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          Vinfo( 2)=Vname(2,idVbcs)
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idVbcs)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbcs)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbcs)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (VarOut(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          Vinfo( 2)=Vname(2,idUbot)
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idUbot)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbot)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbot)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (VarOut(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          Vinfo( 2)=Vname(2,idVbot)
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idVbot)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbot)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbot)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (VarOut(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          Vinfo( 2)=Vname(2,idUbur)
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idUbur)%dkind=PIO_FOUT
          S(ng)%pioVar(idUbur)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUbur)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (VarOut(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          Vinfo( 2)=Vname(2,idVbvr)
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idVbvr)%dkind=PIO_FOUT
          S(ng)%pioVar(idVbvr)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idVbvr)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_HEIGHT
!
!  Define wind-induced significant wave height.
!
        IF (VarOut(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          Vinfo( 2)=Vname(2,idWamp)
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWamp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWamp)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWamp)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_LENGTH
!
!  Define wind-induced mean wavelenght.
!
        IF (VarOut(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          Vinfo( 2)=Vname(2,idWlen)
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWlen)%dkind=PIO_FOUT
          S(ng)%pioVar(idWlen)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWlen)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_LENGTHP
!
!  Define wind-induced peak wave wavelength.
!
        IF (VarOut(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          Vinfo( 2)=Vname(2,idWlep)
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWlep)%dkind=PIO_FOUT
          S(ng)%pioVar(idWlep)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWlep)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_DIR
!
!  Define wind-induced wave direction.
!
        IF (VarOut(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          Vinfo( 2)=Vname(2,idWdir)
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWdir)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdir)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWdir)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_DIRP
!
!  Define wind-induced peak wave direction.
!
        IF (VarOut(idWdip,ng)) THEN
          Vinfo( 1)=Vname(1,idWdip)
          Vinfo( 2)=Vname(2,idWdip)
          Vinfo( 3)=Vname(3,idWdip)
          Vinfo(14)=Vname(4,idWdip)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWdip)%dkind=PIO_FOUT
          S(ng)%pioVar(idWdip)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWdip)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_TOP_PERIOD
!
!  Define wind-induced surface wave period.
!
        IF (VarOut(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          Vinfo( 2)=Vname(2,idWptp)
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWptp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWptp)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWptp)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_BOT_PERIOD
!
!  Define wind-induced bottom wave period.
!
        IF (VarOut(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          Vinfo( 2)=Vname(2,idWpbt)
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWpbt)%dkind=PIO_FOUT
          S(ng)%pioVar(idWpbt)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWpbt)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef WAVES_DSPR
!
!  Define waves directional spreading.
!
        IF (VarOut(idWvds,ng)) THEN
          Vinfo( 1)=Vname(1,idWvds)
          Vinfo( 2)=Vname(2,idWvds)
          Vinfo( 3)=Vname(3,idWvds)
          Vinfo(14)=Vname(4,idWvds)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWvds)%dkind=PIO_FOUT
          S(ng)%pioVar(idWvds)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWvds)%vd,                     &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define wave spectrum peakedness.
!
        IF (VarOut(idWvqp,ng)) THEN
          Vinfo( 1)=Vname(1,idWvqp)
          Vinfo( 2)=Vname(2,idWvqp)
          Vinfo( 3)=Vname(3,idWvqp)
          Vinfo(14)=Vname(4,idWvqp)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idWvqp)%dkind=PIO_FOUT
          S(ng)%pioVar(idWvqp)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idWvqp)%vd,                     &
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
#   ifdef WAVES_UB
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            S(ng)%pioVar(idWorb)%vd=var_desc(i)
            S(ng)%pioVar(idWorb)%dkind=PIO_FOUT
            S(ng)%pioVar(idWorb)%gtype=0
#   endif
#   ifdef BBL_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            S(ng)%pioVar(idUbrs)%vd=var_desc(i)
            S(ng)%pioVar(idUbrs)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbrs)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            S(ng)%pioVar(idVbrs)%vd=var_desc(i)
            S(ng)%pioVar(idVbrs)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbrs)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            S(ng)%pioVar(idUbws)%vd=var_desc(i)
            S(ng)%pioVar(idUbws)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbws)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            S(ng)%pioVar(idVbws)%vd=var_desc(i)
            S(ng)%pioVar(idVbws)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbws)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            S(ng)%pioVar(idUbcs)%vd=var_desc(i)
            S(ng)%pioVar(idUbcs)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbcs)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            S(ng)%pioVar(idVbcs)%vd=var_desc(i)
            S(ng)%pioVar(idVbcs)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbcs)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            S(ng)%pioVar(idUbot)%vd=var_desc(i)
            S(ng)%pioVar(idUbot)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbot)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            S(ng)%pioVar(idVbot)%vd=var_desc(i)
            S(ng)%pioVar(idVbot)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbot)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            S(ng)%pioVar(idUbur)%vd=var_desc(i)
            S(ng)%pioVar(idUbur)%dkind=PIO_FOUT
            S(ng)%pioVar(idUbur)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            S(ng)%pioVar(idVbvr)%vd=var_desc(i)
            S(ng)%pioVar(idVbvr)%dkind=PIO_FOUT
            S(ng)%pioVar(idVbvr)%gtype=0
#   endif
#   ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            S(ng)%pioVar(idWamp)%vd=var_desc(i)
            S(ng)%pioVar(idWamp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWamp)%gtype=0
#   endif
#   ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            S(ng)%pioVar(idWlen)%vd=var_desc(i)
            S(ng)%pioVar(idWlen)%dkind=PIO_FOUT
            S(ng)%pioVar(idWlen)%gtype=0
#   endif
#   ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            S(ng)%pioVar(idWlep)%vd=var_desc(i)
            S(ng)%pioVar(idWlep)%dkind=PIO_FOUT
            S(ng)%pioVar(idWlep)%gtype=0
#   endif
#   ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            S(ng)%pioVar(idWdir)%vd=var_desc(i)
            S(ng)%pioVar(idWdir)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdir)%gtype=0
#   endif
#   ifdef WAVES_DIRP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdip))) THEN
            got_var(idWdip)=.TRUE.
            S(ng)%pioVar(idWdip)%vd=var_desc(i)
            S(ng)%pioVar(idWdip)%dkind=PIO_FOUT
            S(ng)%pioVar(idWdip)%gtype=0
#   endif
#   ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            S(ng)%pioVar(idWptp)%vd=var_desc(i)
            S(ng)%pioVar(idWptp)%dkind=PIO_FOUT
            S(ng)%pioVar(idWptp)%gtype=0
#   endif
#   ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            S(ng)%pioVar(idWpbt)%vd=var_desc(i)
            S(ng)%pioVar(idWpbt)%dkind=PIO_FOUT
            S(ng)%pioVar(idWpbt)%gtype=0
#   endif
          END IF
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
#   ifdef BBL_MODEL
        IF (.not.got_var(idUbrs).and.VarOut(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.VarOut(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.VarOut(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.VarOut(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.VarOut(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.VarOut(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.VarOut(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.VarOut(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.VarOut(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.VarOut(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_UB
        IF (.not.got_var(idWorb).and.VarOut(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.VarOut(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_LENGTH
        IF (.not.got_var(idWlen).and.VarOut(idWlen,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlen)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.VarOut(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.VarOut(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_DIRP
        IF (.not.got_var(idWdip).and.VarOut(idWdip,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdip)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.VarOut(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.VarOut(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef WAVES_DSPR
        IF (.not.got_var(idWvds).and.VarOut(idWvds,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvds)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvqp).and.VarOut(idWvqp,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWvqp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
      END IF QUERY
!
  10  FORMAT (/,' BBL_DEF_STATION_PIO - unable to find variable: ',     &
     &        a,2x,' in stations NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE bbl_def_station_pio
!
!***********************************************************************
      SUBROUTINE bbl_wrt_station_pio (ng, model, tile,                  &
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
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
#   endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", bbl_wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment output variables into specified stations NetCDF
!  file.
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
#   endif
      END DO
!
#   ifdef WAVES_UB
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (VarOut(idWorb,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWorb, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Ub_swan,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWorb)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWorb)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (VarOut(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idUbrs)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idUbrs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (VarOut(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idVbrs)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idVbrs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (VarOut(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idUbws)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idUbws)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (VarOut(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idVbws)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idVbws)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (VarOut(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idUbcs)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idUbcs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (VarOut(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idVbcs)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idVbcs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (VarOut(idUbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idUbot)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idUbot)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (VarOut(idVbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idVbot)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idVbot)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (VarOut(idUbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idUbur)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idUbur)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (VarOut(idVbvr,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbvr, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idVbvr)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idVbvr)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_HEIGHT
!
!  Write out wind-induced significant wave height.
!
      IF (VarOut(idWamp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWamp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Hwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWamp)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWamp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelenght.
!
      IF (VarOut(idWlen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlen, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWlen)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWlen)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (VarOut(idWlep,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlep, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWlep)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWlep)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (VarOut(idWdir,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdir, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWdir)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWdir)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (VarOut(idWdip,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdip, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWdip)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWdip)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (VarOut(idWptp,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWptp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_top,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWptp)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWptp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (VarOut(idWpbt,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWpbt, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_bot,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWpbt)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWpbt)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   if defined WAVES_DSPR
!
!  Write out waves directional spreading.
!
      IF (VarOut(idWvds,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvds, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_ds,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWvds)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWvds)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave spectrum peakedness.
!
      IF (VarOut(idWvqp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvqp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_qp,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idWvqp)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idWvqp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
!
      RETURN
      END SUBROUTINE bbl_wrt_station_pio
#  endif
# endif
#endif
!
      END MODULE bbl_output_mod
