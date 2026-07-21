#include "cppdefs.h"
      MODULE vegetation_output_mod

#if defined VEGETATION && defined SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes submerged aquatic vegetation variables   !
!  into output NetCDF files.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
      USE mod_vegetation
!
      USE def_var_mod,     ONLY : def_var
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: vegetation_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: vegetation_def_pio
# endif
      PUBLIC :: vegetation_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: vegetation_wrt_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_def_nf90 (ng, model, ldef, VarOut, S,       &
     &                                sapgrd, t2dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: sapgrd(:), t2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lrestart, got_var(NV)
!
      integer :: i, itrc, j, nvd3, nvd4, status, vtype
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
     &  __FILE__//", vegetation_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define submerged aquatic vegetation output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        Lrestart=.TRUE.
        Vtype=NF_FRST
      ELSE
        Lrestart=.FALSE.
        Vtype=NF_FOUT
      END IF
!
!  Defined vegetation variables.
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

# if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Define vegetation plant properties (density, diameter, height, and
!  thickness).
!
        DO i=1,NVEGP
          IF (Lrestart.or.VarOut(idvprp(i),ng)) THEN
            Vinfo( 1)=Vname(1,idvprp(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idvprp(i)))
            ELSE
              Vinfo( 2)=Vname(2,idvprp(i))
            END IF
            Vinfo( 3)=Vname(3,idvprp(i))
            Vinfo(14)=Vname(4,idvprp(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idvprp(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idvprp(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idvprp(i)), &
     &                     vtype, nvd4, sapgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif

# if defined VEG_STREAMING
!
!  Define wave dissipation due to vegetation (W/m2).
!
        IF (Lrestart.or.VarOut(idWdvg,ng)) THEN
          Vinfo( 1)=Vname(1,idWdvg)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idWdvg))
          ELSE
            Vinfo( 2)=Vname(2,idWdvg)
          END IF
          Vinfo( 3)=Vname(3,idWdvg)
          Vinfo(14)=Vname(4,idWdvg)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWdvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdvg,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idWdvg),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define spectral drag coefficient due to waves and vegetation.
!
        IF (Lrestart.or.VarOut(idCdvg,ng)) THEN
          Vinfo( 1)=Vname(1,idCdvg)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idCdvg))
          ELSE
            Vinfo( 2)=Vname(2,idCdvg)
          END IF
          Vinfo( 3)=Vname(3,idCdvg)
          Vinfo(14)=Vname(4,idCdvg)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idCdvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idCdvg,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idCdvg),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# ifdef MARSH_DYNAMICS
!
!  Define marsh masking (1=marsh cell or 0=non-march cell), stored at
!  time intervals. It is independent from the wetting and drying masks.
!
        IF (Lrestart.or.VarOut(idTims,ng)) THEN
          Vinfo( 1)=Vname(1,idTims)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTims))
          ELSE
            Vinfo( 2)=Vname(2,idTims)
          END IF
          Vinfo( 3)=Vname(3,idTims)
          Vinfo(14)=Vname(4,idTims)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTims)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTims,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTims),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname ,     &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef MARSH_WAVE_THRUST
!
!  Define lateral total wave thrust from all directions acting on the
!  marsh cells.
!
        IF (Lrestart.or.VarOut(idTtot,ng)) THEN
          Vinfo( 1)=Vname(1,idTtot)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTtot))
          ELSE
            Vinfo( 2)=Vname(2,idTtot)
          END IF
          Vinfo( 3)=Vname(3,idTtot)
          Vinfo(14)=Vname(4,idTtot)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTtot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTtot,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTtot),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef MARSH_SED_EROSION
!
!  Define sediment tracer flux out from marsh cells for each class (kg).
!
        DO i=1,NST
          IF (Lrestart.or.VarOut(idTmfo(i),ng)) THEN
            Vinfo( 1)=Vname(1,idTmfo(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idTmfo(i)))
            ELSE
              Vinfo( 2)=Vname(2,idTmfo(i))
            END IF
            Vinfo( 3)=Vname(3,idTmfo(i))
            Vinfo(14)=Vname(4,idTmfo(i))
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,idTmfo(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTmfo(i),ng))
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmfo(i)), &
     &                     vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#    ifdef MARSH_RETREAT
!
!  Define amount of marshland retreat from every directions (m).
!
        IF (Lrestart.or.VarOut(idTmmr,ng)) THEN
          Vinfo( 1)=Vname(1,idTmmr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmmr))
          ELSE
            Vinfo( 2)=Vname(2,idTmmr)
          END IF
          Vinfo( 3)=Vname(3,idTmmr)
          Vinfo(14)=Vname(4,idTmmr)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idTmmr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmmr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmmr),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#    endif
#   endif
#  endif

#  ifdef MARSH_TIDAL_RANGE
!
!  Define tidal range (m), for a given frequency, leading to marshland
!  growth.
!
        IF (.not.Lrestart.and.VarOut(idTmtr,ng)) THEN
          Vinfo( 1)=Vname(1,idTmtr)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmtr))
          ELSE
            Vinfo( 2)=Vname(2,idTmtr)
          END IF
          Vinfo( 3)=Vname(3,idTmtr)
          Vinfo(14)=Vname(4,idTmtr)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTmtr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmtr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmtr),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined MARSH_VERT_GROWTH
!
!  Define amount of mean high water (m) in marsh cells for a given tidal
!  frequency.
!
        IF (.not.Lrestart.and.VarOut(idTmhw,ng)) THEN
          Vinfo( 1)=Vname(1,idTmhw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmhw))
          ELSE
            Vinfo( 2)=Vname(2,idTmhw)
          END IF
          Vinfo( 3)=Vname(3,idTmhw)
          Vinfo(14)=Vname(4,idTmhw)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTmhw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmhw,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmhw),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define amount of mean low water (m) in marsh cells for a given tidal
!  frequency.
!
        IF (.not.Lrestart.and.VarOut(idTmlw,ng)) THEN
          Vinfo( 1)=Vname(1,idTmlw)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmlw))
          ELSE
            Vinfo( 2)=Vname(2,idTmlw)
          END IF
          Vinfo( 3)=Vname(3,idTmlw)
          Vinfo(14)=Vname(4,idTmlw)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTmlw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmlw,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmlw),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define amount of marsh biomass peak production (kg/m2).
!
        IF (.not.Lrestart.and.VarOut(idTmbp,ng)) THEN
          Vinfo( 1)=Vname(1,idTmbp)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmbp))
          ELSE
            Vinfo( 2)=Vname(2,idTmbp)
          END IF
          Vinfo( 3)=Vname(3,idTmbp)
          Vinfo(14)=Vname(4,idTmbp)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTmbp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmbp,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmbp),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define  marsh vertical growth rate (m/year).
!
        IF (.not.Lrestart.and.VarOut(idTmvg,ng)) THEN
          Vinfo( 1)=Vname(1,idTmvg)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmvg))
          ELSE
            Vinfo( 2)=Vname(2,idTmvg)
          END IF
          Vinfo( 3)=Vname(3,idTmvg)
          Vinfo(14)=Vname(4,idTmvg)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTmvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmvg,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmvg),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define marsh vertical accretion (m) due to biomass production.
!
        IF (.not.Lrestart.and.VarOut(idTmvt,ng)) THEN
          Vinfo( 1)=Vname(1,idTmvt)
          IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmvt))
          ELSE
            Vinfo( 2)=Vname(2,idTmvt)
          END IF
          Vinfo( 3)=Vname(3,idTmvt)
          Vinfo(14)=Vname(4,idTmvt)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTmvt)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmvt,ng),r8)
          status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idTmvt),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
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
!  submerged aquatic vegetation variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
# if defined VEG_STREAMING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdvg))) THEN
            got_var(idWdvg)=.TRUE.
            S(ng)%Vid(idWdvg)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCdvg))) THEN
            got_var(idCdvg)=.TRUE.
            S(ng)%Vid(idCdvg)=var_id(i)
# ifdef MARSH_DYNAMICS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTims))) THEN
            got_var(idTims)=.TRUE.
            S(ng)%Vid(idTims)=var_id(i)
#  ifdef MARSH_WAVE_THRUST
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTtot))) THEN
            got_var(idTtot)=.TRUE.
            S(ng)%Vid(idTtot)=var_id(i)
#   if defined MARSH_SED_EROSION && defined MARSH_RETREAT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmmr))) THEN
            got_var(idTmmr)=.TRUE.
            S(ng)%Vid(idTmmr)=var_id(i)
#   endif
#  endif
#  ifdef MARSH_TIDAL_RANGE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmtr))) THEN
            got_var(idTmtr)=.TRUE.
            S(ng)%Vid(idTmtr)=var_id(i)
#  endif
#  if defined MARSH_VERT_GROWTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmhw))) THEN
            got_var(idTmhw)=.TRUE.
            S(ng)%Vid(idTmhw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmlw))) THEN
            got_var(idTmlw)=.TRUE.
            S(ng)%Vid(idTmlw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmbp))) THEN
            got_var(idTmbp)=.TRUE.
            S(ng)%Vid(idTmbp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmvg))) THEN
            got_var(idTmvg)=.TRUE.
            S(ng)%Vid(idTmvg)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmvt))) THEN
            got_var(idTmvt)=.TRUE.
            S(ng)%Vid(idTmvt)=var_id(i)
#  endif
# endif
          END IF
# if defined VEG_DRAG || defined VEG_BIOMASS
          DO itrc=1,NVEGP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idvprp(itrc)))) THEN
              got_var(idvprp(itrc))=.TRUE.
              S(ng)%Vid(idvprp(itrc))=var_id(i)
            END IF
          END DO
# endif
#  if defined MARSH_DYNAMICS && defined MARSH_WAVE_THRUST
#   ifdef MARSH_SED_EROSION
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmfo(itrc)))) THEN
              got_var(idTmfo(itrc))=.TRUE.
            END IF
          END DO
#   endif
#  endif
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
# if defined VEG_STREAMING
        IF (.not.got_var(idWdvg).and.                                   &
     &      (Lrestart.or.VarOut(idWdvg,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdvg)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idCdvg).and.                                   &
     &      (Lrestart.or.VarOut(idCdvg,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idCdvg)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef MARSH_DYNAMICS
        IF (.not.got_var(idTims).and.                                   &
     &      (Lrestart.or.VarOut(idTims,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTims)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef MARSH_WAVE_THRUST
        IF (.not.got_var(idTtot).and.                                   &
     &      (Lrestart.or.VarOut(idTtot,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTtot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   if defined MARSH_SED_EROSION && defined MARSH_RETREAT
        IF (.not.got_var(idTmmr).and.                                   &
     &      (Lrestart.or.VarOut(idTmmr,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmmr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef MARSH_TIDAL_RANGE
        IF (.not.got_var(idTmtr).or.                                    &
     &      (.not.Lrestart.and.VarOut(idTmtr,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmtr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined MARSH_VERT_GROWTH
        IF (.not.got_var(idTmhw).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmhw,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmhw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmlw).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmlw,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmlw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmbp).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmbp,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmbp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmvg).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmvg,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmvg)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmvt).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmvt,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmvt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# if defined VEG_DRAG || defined VEG_BIOMASS
        DO i=1,NVEGP
          IF (.not.got_var(idvprp(i)).and.                              &
     &        (Lrestart.or.VarOut(idvprp(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idvprp(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# if defined MARSH_DYNAMICS && defined MARSH_WAVE_THRUST
#  ifdef MARSH_SED_EROSION
        DO itrc=1,NST
          IF (.not.got_var(idTmfo(itrc)).and.                           &
     &        (Lrestart.or.VarOut(idTmfo(itrc),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmfo(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif

      END IF QUERY
!
  10  FORMAT (/,' VEGETATION_DEF_NF90 - unable to find variable: ',     &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE vegetation_def_nf90
!
!***********************************************************************
      SUBROUTINE vegetation_wrt_nf90 (ng, model, tile,                  &
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
      logical :: Lrestart
!
      integer :: gfactor, gtype, i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", vegetation_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out submerged aquatic vegetation output variables into
!  specified output NetCDF file.
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
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        Lrestart=.TRUE.
      ELSE
        Lrestart=.FALSE.
      END IF

# if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Write out vegetation plant properties (density, diameter, height, and
!  thickness).
!
      DO i=1,NVEGP
        IF (Lrestart.or.VarOut(idvprp(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*a3dvar
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idvprp(i),          &
     &                       S(ng)%Vid(idvprp(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, NVEG, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       VEG(ng) % plant(:,:,:,i),                  &
     &                       SetFillVal= .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idvprp(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif

# ifdef VEG_STREAMING
!
!  Write out wave dissipation due to vegetation (W/m2)
!
      IF (Lrestart.or.VarOut(idWdvg,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idWdvg,               &
     &                     S(ng)%Vid(idWdvg),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     VEG(ng)%Dissip_veg)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdvg)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out spectral drag coefficient due to waves and vegetation.
!
      IF (Lrestart.or.VarOut(idCdvg,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idCdvg,               &
     &                     S(ng)%Vid(idCdvg),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     VEG(ng)%Cdwave_veg)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idCdvg)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef MARSH_DYNAMICS
!
!  Write out marsh masking (1=marsh cell or 0=non-march cell), stored at
!  time intervals. It is independent from the wetting and drying masks.
!
      IF (Lrestart.or.VarOut(idTims,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTims,               &
     &                     S(ng)%Vid(idTims),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     VEG(ng)%marsh_mask)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTims)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef MARSH_WAVE_THRUST
!
!  Write out lateral total wave thrust from all directions acting on the
!  marsh cells.
!
      IF (Lrestart.or.VarOut(idTtot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTtot,               &
     &                     S(ng)%Vid(idTtot),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%Thrust_total)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTtot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   if defined MARSH_SED_EROSION
!
!  Write out sediment tracer flux out from marsh cells for each class
!  (kg).
!
      DO i=1,NST
        IF (Lrestart.or.VarOut(idTmfo(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmfo(i),          &
     &                       S(ng)%Vid(idTmfo(i)),                      &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       VEG(ng) % marsh_flux_out(:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTmfo(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#    if defined MARSH_RETREAT
!
!  Write out amount of marshland retreat from every directions (m).
!
      IF (Lrestart.or.VarOut(idTmmr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmmr,               &
     &                     S(ng)%Vid(idTmmr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     VEG(ng)%marsh_retreat)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmmr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#  endif

#  if defined MARSH_TIDAL_RANGE
!
!  Write out tidal range (m), for a given frequency, leading to
!  marshland growth.
!
      IF (.not.Lrestart.and.VarOut(idTmtr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmtr,               &
     &                     S(ng)%Vid(idTmtr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_tidal_range)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmtr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  if defined MARSH_VERT_GROWTH
!
!  Write out amount of mean high water (m) in marsh cells for a given
!  tidal frequency.
!
      IF (.not.Lrestart.and.VarOut(idTmhw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmhw,               &
     &                     S(ng)%Vid(idTmhw),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_high_water)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmhw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out amount of mean low water (m) in marsh cells for a given
!  tidal frequency.
!
      IF (.not.Lrestart.and.VarOut(idTmlw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmlw,               &
     &                     S(ng)%Vid(idTmlw),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_low_water)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmlw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out amount of marsh biomass peak production (kg/m2).
!
      IF (.not.Lrestart.and.VarOut(idTmbp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmbp,               &
     &                     S(ng)%Vid(idTmbp),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_biomass_peak)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmbp)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out marsh vertical growth rate (m/year).
!
      IF (.not.Lrestart.and.VarOut(idTmvg,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmvg,               &
     &                     S(ng)%Vid(idTmvg),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_vert_rate)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmvg)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out marsh vertical accretion (m) due to biomass production.
!
      IF (.not.Lrestart.and.VarOut(idTmvt,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idTmvt,               &
     &                     S(ng)%Vid(idTmvt),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_accret)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                   __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmvt)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
  10  FORMAT (/," VEGETATION_WRT_NF90 - error while writing variable '",&
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE vegetation_wrt_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE vegetation_def_pio (ng, model, ldef, VarOut, S,        &
     &                               sapgrd, t2dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: sapgrd(:), t2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lrestart, got_var(NV)
!
      integer :: FileH, i, itrc, j, nvd3, nvd4, status, vtype
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
     &  __FILE__//", vegetation_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define submerged aquatic vegetation output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
      FileH=ABS(S(ng)%pioFile%fh)
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        Lrestart=.TRUE.
        Vtype=PIO_FRST
      ELSE
        Lrestart=.FALSE.
        Vtype=PIO_FOUT
      END IF
!
!  Defined vegetation variables.
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

#  if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Define vegetation module properties (density, diameter, height, and
!  thickness)
!
        DO i=1,NVEGP
          IF (Lrestart.or.VarOut(idvprp(i),ng)) THEN
            Vinfo( 1)=Vname(1,idvprp(i))
            IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idvprp(i)))
            ELSE
              Vinfo( 2)=Vname(2,idvprp(i))
            END IF
            Vinfo( 3)=Vname(3,idvprp(i))
            Vinfo(14)=Vname(4,idvprp(i))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idvprp(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idvprp(i),ng),r8)
            S(ng)%pioVar(idvprp(i))%dkind=vtype
            S(ng)%pioVar(idvprp(i))%gtype=a3dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idvprp(i))%vd,                  &
     &                     vtype, nvd4, sapgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif

#  if defined VEG_STREAMING
!
!  Define wave dissipation due to vegetation (W/m2).
!
        IF (Lrestart.or.VarOut(idWdvg,ng)) THEN
          Vinfo( 1)=Vname(1,idWdvg)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idWdvg))
          ELSE
            Vinfo( 2)=Vname(2,idWdvg)
          END IF
          Vinfo( 3)=Vname(3,idWdvg)
          Vinfo(14)=Vname(4,idWdvg)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWdvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdvg,ng),r8)
          S(ng)%pioVar(idWdvg)%dkind=vtype
          S(ng)%pioVar(idWdvg)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idWdvg)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define spectral drag coefficient due to waves and vegetation.
!
        IF (Lrestart.or.VarOut(idCdvg,ng)) THEN
          Vinfo( 1)=Vname(1,idCdvg)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idCdvg))
          ELSE
            Vinfo( 2)=Vname(2,idCdvg)
          END IF
          Vinfo( 3)=Vname(3,idCdvg)
          Vinfo(14)=Vname(4,idCdvg)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idCdvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idCdvg,ng),r8)
          S(ng)%pioVar(idCdvg)%dkind=vtype
          S(ng)%pioVar(idCdvg)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idCdvg)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  ifdef MARSH_DYNAMICS
!
!  Define marsh masking (1=marsh cell or 0=non-march cell), stored at
!  time intervals. It is independent from the wetting and drying masks.
!
        IF (Lrestart.or.VarOut(idTims,ng)) THEN
          Vinfo( 1)=Vname(1,idTims)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTims))
          ELSE
            Vinfo( 2)=Vname(2,idTims)
          END IF
          Vinfo( 3)=Vname(3,idTims)
          Vinfo(14)=Vname(4,idTims)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTims)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTims,ng),r8)
          S(ng)%pioVar(idTims)%dkind=vtype
          S(ng)%pioVar(idTims)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTims)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname ,     &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef MARSH_WAVE_THRUST
!
!  Define lateral total wave thrust from all directions acting on the
!  marsh cells.
!
        IF (Lrestart.or.VarOut(idTtot,ng)) THEN
          Vinfo( 1)=Vname(1,idTtot)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTtot))
          ELSE
            Vinfo( 2)=Vname(2,idTtot)
          END IF
          Vinfo( 3)=Vname(3,idTtot)
          Vinfo(14)=Vname(4,idTtot)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTtot)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTtot,ng),r8)
          S(ng)%pioVar(idTtot)%dkind=vtype
          S(ng)%pioVar(idTtot)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTtot)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#    ifdef MARSH_SED_EROSION
!
!  Define sediment tracer flux out from marsh cells for each class (kg).
!
        DO i=1,NST
          IF (Lrestart.or.VarOut(idTmfo(i),ng)) THEN
            Vinfo( 1)=Vname(1,idTmfo(i))
            IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idTmfo(i)))
            ELSE
              Vinfo( 2)=Vname(2,idTmfo(i))
            END IF
            Vinfo( 3)=Vname(3,idTmfo(i))
            Vinfo(14)=Vname(4,idTmfo(i))
            Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#     endif
            Vinfo(21)=Vname(6,idTmfo(i))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTmfo(i),ng))
            S(ng)%pioVar(idTmfo(i))%dkind=vtype
            S(ng)%pioVar(idTmfo(i))%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idTmfo(i))%vd,                  &
     &                     vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#     ifdef MARSH_RETREAT
!
!  Define amount of marshland retreat from every directions (m).
!
        IF (Lrestart.or.VarOut(idTmmr,ng)) THEN
          Vinfo( 1)=Vname(1,idTmmr)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmmr))
          ELSE
            Vinfo( 2)=Vname(2,idTmmr)
          END IF
          Vinfo( 3)=Vname(3,idTmmr)
          Vinfo(14)=Vname(4,idTmmr)
          Vinfo(16)=Vname(1,idtime)
#      if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#      endif
          Vinfo(21)=Vname(6,idTmmr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmmr,ng),r8)
          S(ng)%pioVar(idTmmr)%dkind=vtype
          S(ng)%pioVar(idTmmr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmmr)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#     endif
#    endif
#   endif

#   ifdef MARSH_TIDAL_RANGE
!
!  Define tidal range (m), for a given frequency, leading to marshland
!  growth.
!
        IF (.not.Lrestart.and.VarOut(idTmtr,ng)) THEN
          Vinfo( 1)=Vname(1,idTmtr)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmtr))
          ELSE
            Vinfo( 2)=Vname(2,idTmtr)
          END IF
          Vinfo( 3)=Vname(3,idTmtr)
          Vinfo(14)=Vname(4,idTmtr)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmtr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmtr,ng),r8)
          S(ng)%pioVar(idTmtr)%dkind=vtype
          S(ng)%pioVar(idTmtr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmtr)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   if defined MARSH_VERT_GROWTH
!
!  Define amount of mean high water (m) in marsh cells for a given tidal
!  frequency.
!
        IF (.not.Lrestart.and.VarOut(idTmhw,ng)) THEN
          Vinfo( 1)=Vname(1,idTmhw)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmhw))
          ELSE
            Vinfo( 2)=Vname(2,idTmhw)
          END IF
          Vinfo( 3)=Vname(3,idTmhw)
          Vinfo(14)=Vname(4,idTmhw)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmhw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmhw,ng),r8)
          S(ng)%pioVar(idTmhw)%dkind=vtype
          S(ng)%pioVar(idTmhw)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmhw)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define amount of mean low water (m) in marsh cells for a given tidal
!  frequency.
!
        IF (.not.Lrestart.and.VarOut(idTmlw,ng)) THEN
          Vinfo( 1)=Vname(1,idTmlw)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmlw))
          ELSE
            Vinfo( 2)=Vname(2,idTmlw)
          END IF
          Vinfo( 3)=Vname(3,idTmlw)
          Vinfo(14)=Vname(4,idTmlw)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmlw)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmlw,ng),r8)
          S(ng)%pioVar(idTmlw)%dkind=vtype
          S(ng)%pioVar(idTmlw)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmlw)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define amount of marsh biomass peak production (kg/m2).
!
        IF (.not.Lrestart.and.VarOut(idTmbp,ng)) THEN
          Vinfo( 1)=Vname(1,idTmbp)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmbp))
          ELSE
            Vinfo( 2)=Vname(2,idTmbp)
          END IF
          Vinfo( 3)=Vname(3,idTmbp)
          Vinfo(14)=Vname(4,idTmbp)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmbp)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmbp,ng),r8)
          S(ng)%pioVar(idTmbp)%dkind=vtype
          S(ng)%pioVar(idTmbp)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmbp)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define  marsh vertical growth rate (m/year).
!
        IF (.not.Lrestart.and.VarOut(idTmvg,ng)) THEN
          Vinfo( 1)=Vname(1,idTmvg)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmvg))
          ELSE
            Vinfo( 2)=Vname(2,idTmvg)
          END IF
          Vinfo( 3)=Vname(3,idTmvg)
          Vinfo(14)=Vname(4,idTmvg)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmvg)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmvg,ng),r8)
          S(ng)%pioVar(idTmvg)%dkind=vtype
          S(ng)%pioVar(idTmvg)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmvg)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define marsh vertical accretion (m) due to biomass production.
!
        IF (.not.Lrestart.and.VarOut(idTmvt,ng)) THEN
          Vinfo( 1)=Vname(1,idTmvt)
          IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTmvt))
          ELSE
            Vinfo( 2)=Vname(2,idTmvt)
          END IF
          Vinfo( 3)=Vname(3,idTmvt)
          Vinfo(14)=Vname(4,idTmvt)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTmvt)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTmvt,ng),r8)
          S(ng)%pioVar(idTmvt)%dkind=vtype
          S(ng)%pioVar(idTmvt)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idTmvt)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
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
!  Initialize local logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  submerged aquatic vegetation variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
#  if defined VEG_STREAMING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdvg))) THEN
            got_var(idWdvg)=.TRUE.
            S(ng)%pioVar(idWdvg)%vd=var_desc(i)
            S(ng)%pioVar(idWdvg)%dkind=vtype
            S(ng)%pioVar(idWdvg)%gtype=r2dvar
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idCdvg))) THEN
            got_var(idCdvg)=.TRUE.
            S(ng)%pioVar(idCdvg)%vd=var_desc(i)
            S(ng)%pioVar(idCdvg)%dkind=vtype
            S(ng)%pioVar(idCdvg)%gtype=r2dvar
#  ifdef MARSH_DYNAMICS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTims))) THEN
            got_var(idTims)=.TRUE.
            S(ng)%pioVar(idTims)%vd=var_desc(i)
            S(ng)%pioVar(idTims)%dkind=vtype
            S(ng)%pioVar(idTims)%gtype=r2dvar
#   ifdef MARSH_WAVE_THRUST
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTtot))) THEN
            got_var(idTtot)=.TRUE.
            S(ng)%pioVar(idTtot)%vd=var_desc(i)
            S(ng)%pioVar(idTtot)%dkind=vtype
            S(ng)%pioVar(idTtot)%gtype=r2dvar
#    if defined MARSH_SED_EROSION && defined MARSH_RETREAT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmmr))) THEN
            got_var(idTmmr)=.TRUE.
            S(ng)%pioVar(idTmmr)%vd=var_desc(i)
            S(ng)%pioVar(idTmmr)%dkind=vtype
            S(ng)%pioVar(idTmmr)%gtype=r2dvar
#    endif
#   endif
#   ifdef MARSH_TIDAL_RANGE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmtr))) THEN
            got_var(idTmtr)=.TRUE.
            S(ng)%pioVar(idTmtr)%vd=var_desc(i)
            S(ng)%pioVar(idTmtr)%dkind=vtype
            S(ng)%pioVar(idTmtr)%gtype=r2dvar
#   endif
#   if defined MARSH_VERT_GROWTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmhw))) THEN
            got_var(idTmhw)=.TRUE.
            S(ng)%pioVar(idTmhw)%vd=var_desc(i)
            S(ng)%pioVar(idTmhw)%dkind=vtype
            S(ng)%pioVar(idTmhw)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmlw))) THEN
            got_var(idTmlw)=.TRUE.
            S(ng)%pioVar(idTmlw)%vd=var_desc(i)
            S(ng)%pioVar(idTmlw)%dkind=vtype
            S(ng)%pioVar(idTmlw)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmbp))) THEN
            got_var(idTmbp)=.TRUE.
            S(ng)%pioVar(idTmbp)%vd=var_desc(i)
            S(ng)%pioVar(idTmbp)%dkind=vtype
            S(ng)%pioVar(idTmbp)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmvg))) THEN
            got_var(idTmvg)=.TRUE.
            S(ng)%pioVar(idTmvg)%vd=var_desc(i)
            S(ng)%pioVar(idTmvg)%dkind=vtype
            S(ng)%pioVar(idTmvg)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmvt))) THEN
            got_var(idTmvt)=.TRUE.
            S(ng)%pioVar(idTmvt)%vd=var_desc(i)
            S(ng)%pioVar(idTmvt)%dkind=vtype
            S(ng)%pioVar(idTmvt)%gtype=r2dvar
#   endif
#  endif
          END IF
#  if defined VEG_DRAG || defined VEG_BIOMASS
          DO itrc=1,NVEGP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idvprp(itrc)))) THEN
              got_var(idvprp(itrc))=.TRUE.
              S(ng)%pioVar(idvprp(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idvprp(itrc))%dkind=vtype
              S(ng)%pioVar(idvprp(itrc))%gtype=a3dvar
            END IF
          END DO
#  endif
#  if defined MARSH_DYNAMICS && defined MARSH_WAVE_THRUST
#   ifdef MARSH_SED_EROSION
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTmfo(itrc)))) THEN
              got_var(idTmfo(itrc))=.TRUE.
              S(ng)%pioVar(idTmfo(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idTmfo(itrc))%dkind=vtype
              S(ng)%pioVar(idTmfo(itrc))%gtype=r2dvar
            END IF
          END DO
#   endif      
#  endif   
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
#  if defined VEG_STREAMING
        IF (.not.got_var(idWdvg).and.                                   &
     &      (Lrestart.or.VarOut(idWdvg,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idWdvg)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idCdvg).and.                                   &
     &      (Lrestart.or.VarOut(idCdvg,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idCdvg)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef MARSH_DYNAMICS
        IF (.not.got_var(idTims).and.                                   &
     &      (Lrestart.or.VarOut(idTims,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTims)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef MARSH_WAVE_THRUST
        IF (.not.got_var(idTtot).and.                                   &
     &      (Lrestart.or.VarOut(idTtot,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTtot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    if defined MARSH_SED_EROSION && defined MARSH_RETREAT
        IF (.not.got_var(idTmmr).and.                                   &
     &      (Lrestart.or.VarOut(idTmmr,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmmr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#   ifdef MARSH_TIDAL_RANGE
        IF (.not.got_var(idTmtr).or.                                    &
     &      (.not.Lrestart.and.VarOut(idTmtr,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmtr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined MARSH_VERT_GROWTH
        IF (.not.got_var(idTmhw).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmhw,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmhw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmlw).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmlw,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmlw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmbp).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmbp,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmbp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmvg).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmvg,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmvg)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTmvt).and.                                   &
     &      (.not.Lrestart.and.VarOut(idTmvt,ng))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmvt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  if defined VEG_DRAG || defined VEG_BIOMASS
        DO i=1,NVEGP
          IF (.not.got_var(idvprp(i)).and.                              &
     &        (Lrestart.or.VarOut(idvprp(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idvprp(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined MARSH_DYNAMICS && defined MARSH_WAVE_THRUST
#   ifdef MARSH_SED_EROSION
        DO itrc=1,NST
          IF (.not.got_var(idTmfo(itrc)).and.                           &
     &        (Lrestart.or.VarOut(idTmfo(itrc),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idTmfo(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif
#  endif

      END IF QUERY
!
  10  FORMAT (/,' VEGETATION_DEF_PIO - unable to find variable: ',      &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE vegetation_def_pio
!
!***********************************************************************
      SUBROUTINE vegetation_wrt_pio (ng, model, tile,                   &
    &                                LBi, UBi, LBj, UBj,                &
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
      logical :: Lrestart
!
      integer :: FileH, gfactor, gtype, i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", vegetation_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out submerged aquatic vegetation output variables into
!  specified NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
      FileH=ABS(S(ng)%pioFile%fh)
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
      gfactor=-1
#  else
      gfactor=1
#  endif
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        Lrestart=.TRUE.
      ELSE
        Lrestart=.FALSE.
      END IF

#  if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Write out vegetation plant properties (density, diameter, height, and
!  thickness).
!
      DO i=1,NVEGP
        IF (Lrestart.or.VarOut(idvprp(i),ng)) THEN
          scale=1.0_r8
          IF (S(ng)%pioVar(idvprp(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_a3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_a3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idvprp(i),       &
     &                       S(ng)%pioVar(idvprp(i)),                   &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, NVEG, scale,        &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       VEG(ng) % plant(:,:,:,i),                  &
     &                       SetFillVal= .FALSE.)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idvprp(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  ifdef VEG_STREAMING
!
!  Write out wave dissipation due to vegetation (W/m2)
!
      IF (Lrestart.or.VarOut(idWdvg,ng)) THEN
        scale=rho0                          ! W m /kg to W/m2
        IF (S(ng)%pioVar(idWdvg)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idWdvg,            &
     &                     S(ng)%pioVar(idWdvg),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%Dissip_veg)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdvg)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out spectral drag coefficient due to waves and vegetation.
!
      IF (Lrestart.or.VarOut(idCdvg,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idCdvg)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idCdvg,            &
     &                     S(ng)%pioVar(idCdvg),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%Cdwave_veg)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idCdvg)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef MARSH_DYNAMICS
!
!  Write out marsh masking (1=marsh cell or 0=non-march cell), stored at
!  time intervals. It is independent from the wetting and drying masks.
!
      IF (Lrestart.or.VarOut(idTims,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTims)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTims,            &
     &                     S(ng)%pioVar(idTims),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     VEG(ng)%marsh_mask)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTims)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef MARSH_WAVE_THRUST
!
!  Write out lateral total wave thrust from all directions acting on the
!  marsh cells.
!
      IF (Lrestart.or.VarOut(idTtot,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTtot)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTtot,            &
     &                     S(ng)%pioVar(idTtot),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%Thrust_total)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTtot)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#    if defined MARSH_SED_EROSION
!
!  Write out sediment tracer flux out from marsh cells for each class
!  (kg).
!
      DO i=1,NST
        IF (Lrestart.or.VarOut(idTmfo(i),ng)) THEN
          scale=1.0_r8
          IF (S(ng)%pioVar(idTmfo(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmfo(i),       &
     &                       S(ng)%pioVar(idTmfo(i)),                   &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#     ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#     endif
     &                       VEG(ng) % marsh_flux_out(:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTmfo(i))), S(ng)%Rindex, &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#     if defined MARSH_RETREAT
!
!  Write out amount of marshland retreat from every directions (m).
!
      IF (Lrestart.or.VarOut(idTmmr,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmmr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmmr,            &
     &                     S(ng)%pioVar(idTmmr),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#      ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#      endif
     &                     VEG(ng)%marsh_retreat)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmmr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#     endif
#    endif
#   endif

#   if defined MARSH_TIDAL_RANGE
!
!  Write out tidal range (m), for a given frequency, leading to
!  marshland growth.
!
      IF (.not.Lrestart.and.VarOut(idTmtr,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmtr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmtr,            &
     &                     S(ng)%pioVar(idTmtr),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_tidal_range)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmtr)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif

#   if defined MARSH_VERT_GROWTH
!
!  Write out amount of mean high water (m) in marsh cells for a given
!  tidal frequency.
!
      IF (.not.Lrestart.and.VarOut(idTmhw,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmhw)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmhw,            &
     &                     S(ng)%pioVar(idTmhw),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_high_water)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmhw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out amount of mean low water (m) in marsh cells for a given
!  tidal frequency.
!
      IF (.not.Lrestart.and.VarOut(idTmlw,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmlw)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmlw,            &
     &                     S(ng)%pioVar(idTmlw),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_low_water)
        IF (FoundError(status, PIO_noerr, __LINE__, __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmlw)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out amount of marsh biomass peak production (kg/m2).
!
      IF (.not.Lrestart.and.VarOut(idTmbp,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmbp)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmbp,            &
     &                     S(ng)%pioVar(idTmbp),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_biomass_peak)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmbp)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out marsh vertical growth rate (m/year).
!
      IF (.not.Lrestart.and.VarOut(idTmvg,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmvg)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmvg,            &
     &                     S(ng)%pioVar(idTmvg),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_vert_rate)
        IF (FoundError(status, PIO_noerr, __LINE__, __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmvg)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out marsh vertical accretion (m) due to biomass production.
!
      IF (.not.Lrestart.and.VarOut(idTmvt,ng)) THEN
        scale=1.0_r8
        IF (S(ng)%pioVar(idTmvt)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idTmvt,            &
     &                     S(ng)%pioVar(idTmvt),                        &
     &                     S(ng)%Rindex,                                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     VEG(ng)%marsh_accret)
        IF (FoundError(status, PIO_noerr, __LINE__, __FILE__)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmvt)), S(ng)%Rindex,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
!
  10  FORMAT (/," VEGETATION_WRT_PIO - error while writing variable '", &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE vegetation_wrt_pio
# endif
!
#endif
      END MODULE vegetation_output_mod
