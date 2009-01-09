#include "cppdefs.h"
      SUBROUTINE get_state (ng, model, msg, ncname, IniRec, Tindex)
!
!svn $Id: get_state.F 1607 2008-07-27 03:03:34Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested model state from specified NetCDF   !
!  file. It is usually used to read initial conditions.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     msg        Message index for Mstate.                             !
!     ncname     Nonlinear initial conditions NetCDF file name.        !
!     IniRec     Nonlinear initial conditions time record to read.     !
!     Tindex     State variable time index to intialize.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
#endif
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
#if defined GLS_MIXING || defined MY25_MIXING || defined LMD_MIXING || \
    defined FOUR_DVAR
      USE mod_mixing
#endif
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
      USE mod_stepping
      USE mod_strings
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex

      integer, intent(inout) :: IniRec

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Perfect2D, Perfect3D, foundit

      logical :: find_string

      logical, dimension(NV) :: get_var, have_var

      integer :: LBi, UBi, LBj, UBj
      integer :: IDmod, InpRec, gtype, i, itrc, lstr, lend
      integer :: Nrec, mySize, ncINPid, nvatts, nvdim, status, varid
      integer :: Vsize(4), start(4), total(4)

      integer :: nf_fread2d
      integer :: nf_fread3d
#ifdef SOLVE3D
      integer :: nf_fread4d
#endif

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, INPtime, Tmax, scale, time_scale

      real(r8), allocatable :: myTime(:)

      character (len=6 ) :: string
      character (len=14) :: t_code
      character (len=15) :: attnam, tvarnam
      character (len=40) :: tunits
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Set model identification string.
!
      IF (model.eq.iNLM.or.(model.eq.0)) THEN
        string=' NLM: '                       ! nonlinear model, restart
        IDmod=iNLM
      ELSE IF (model.eq.iTLM) THEN
        string=' TLM: '                       ! tangent linear model
        IDmod=iTLM
      ELSE IF (model.eq.iRPM) THEN
        string=' RPM: '                       ! representer model
        IDmod=iRPM
      ELSE IF (model.eq.iADM) THEN
        string=' ADM: '                       ! adjoint model
        IDmod=iADM
      ELSE IF (model.eq.5) THEN
        string=' NRM: '                       ! correlations norm factor
        IDmod=iNLM
      ELSE IF (model.eq.6) THEN
        string=' STD: '                       ! standard deviations
        IDmod=iNLM
      ELSE IF (model.eq.7) THEN
        string=' FRC: '                       ! impulse forcing
        IDmod=iNLM
      END IF
!
!  Set switch to process variables for nonlinear model perfect restart.
!
      Perfect2D=.FALSE.
      Perfect3D=.FALSE.
#ifdef PERFECT_RESTART
      IF ((model.eq.0).and.(nrrec(ng).ne.0)) THEN
        Perfect2D=.TRUE.
        Perfect3D=.TRUE.
      END IF
#endif
      PerfectRST(ng)=Perfect2d.or.Perfect3d
!
!  Determine variables to read.
!
      CALL checkvars (ng, IDmod, ncname, string, Nrec, NV, tvarnam,     &
     &                get_var, have_var)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!-----------------------------------------------------------------------
!  Open input NetCDF file and check time variable.
!-----------------------------------------------------------------------
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Open input NetCDF file.
!
      IF (InpThread) THEN
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncINPid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) string, TRIM(ncname)
          exit_flag=2
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, IDmod, exit_flag, 1)
#endif
      IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the input time variable.
!
      CALL netcdf_inq_var (ng, IDmod, ncname, TRIM(tvarnam), varid,     &
     &                     nvdim, nvatts)
      IF (exit_flag.ne.NoError) RETURN
!
!  Allocate input time variable and read its value(s).  Recall that
!  input time variable is a one-dimensional array with one or several
!  values.
!
      mySize=var_Dsize(1)
      IF (.not.allocated(myTime)) allocate (myTime(mySize))
      CALL netcdf_get_fvar (ng, IDmod, ncname, TRIM(tvarnam), myTime)
      IF (exit_flag.ne.NoError) RETURN
!
!  If using the latest time record from input NetCDF file as the
!  initialization record, assign input time.
!
      IF (LastRec(ng)) THEN
        Tmax=-1.0_r8
        DO i=1,mySize
          IF (myTime(i).gt.Tmax) THEN
            Tmax=myTime(i)
            IniRec=i
          END IF
        END DO
        INPtime=Tmax
        InpRec=IniRec
      ELSE
        IF ((IniRec.ne.0).and.(IniRec.gt.mySize)) THEN
          IF (Master)  WRITE (stdout,40) string, IniRec, TRIM(ncname),   &
     &                                   mySize
          exit_flag=2
          RETURN
        END IF
        IF (IniRec.ne.0) THEN
          InpRec=IniRec
        ELSE
          InpRec=1
        END IF
        INPtime=myTime(InpRec)
      END IF
      IF (allocated(myTime)) deallocate ( myTime )
!
!  Set input time scale by looking at the "units" attribute.
!
      time_scale=0.0_r8
      DO i=1,nvatts
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            time_scale=day2sec
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            time_scale=1.0_r8
          END IF
        END IF
      END DO
      IF (time_scale.gt.0.0_r8) THEN
        INPtime=INPtime*time_scale
      END IF
!
!  Set starting time index and time clock in days.
!
      IF ((model.eq.0).or.(model.eq.iNLM).or.                           &
     &    (model.eq.iTLM).or.(model.eq.iRPM)) THEN
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.                   &
     &      (INPtime.ne.(dstart*day2sec))) THEN
          INPtime=dstart*day2sec
        END IF
        time(ng)=INPtime
        tdays(ng)=time(ng)*sec2day
        ntstart(ng)=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
        IF (ntstart(ng).lt.1) ntstart(ng)=1
        IF (PerfectRST(ng)) THEN
          ntfirst(ng)=1
        ELSE
          ntfirst(ng)=ntstart(ng)
        END IF
      ELSE IF (model.eq.iADM) THEN
        IF (msg.ne.1) THEN
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
        END IF
        ntstart(ng)=ntimes(ng)+1
        ntend(ng)=1
        ntfirst(ng)=ntend(ng)
      END IF
      CALL time_string (time(ng), time_code(ng))
!
!  Over-write "IniRec" to the actual initial record processed.
!
      IF (model.eq.iNLM) THEN
        IniRec=InpRec
      END IF
!
!  Report information.
!
      lstr=SCAN(ncname,'/',BACK=.TRUE.)+1
      lend=LEN_TRIM(ncname)
      IF (Master) THEN
        CALL time_string (INPtime, t_code)
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,80) string, TRIM(Mstate(msg)), t_code, Nrun,    &
     &                      ncname(lstr:lend), InpRec, Tindex
        ELSE
          WRITE (stdout,90) string, TRIM(Mstate(msg)), t_code,          &
     &                      ncname(lstr:lend), InpRec, Tindex
        END IF
      END IF

#ifdef NONLINEAR
!
!-----------------------------------------------------------------------
!  Read in nonlinear state variables. If applicable, read in perfect
!  restart variables.
!-----------------------------------------------------------------------
!
      NLM_STATE: IF ((model.eq.iNLM).or.(model.eq.0)) THEN

# ifdef PERFECT_RESTART
!
!  Read in time-stepping indices.
!
        IF (InpThread.and.((model.eq.0).and.(nrrec(ng).ne.0))) THEN
          start(1)=InpRec
          total(1)=1
#  ifdef SOLVE3D
          status=nf90_inq_varid(ncINPid, 'nstp', varid)
          status=nf90_get_var(ncINPid, varid, nstp(ng:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, 'nstp', InpRec, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
          status=nf90_inq_varid(ncINPid, 'nrhs', varid)
          status=nf90_get_var(ncINPid, varid, nrhs(ng:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, 'nrhs', InpRec, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
          status=nf90_inq_varid(ncINPid, 'nnew', varid)
          status=nf90_get_var(ncINPid, varid, nnew(ng:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, 'nnew', InpRec, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# endif
          status=nf90_inq_varid(ncINPid, 'kstp', varid)
          status=nf90_get_var(ncINPid, varid, kstp(ng:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, 'kstp', InpRec, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
          status=nf90_inq_varid(ncINPid,'krhs',varid)
          status=nf90_get_var(ncINPid, varid, krhs(ng:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, 'krhs', InpRec, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
          status=nf90_inq_varid(ncINPid, 'knew', varid)
          status=nf90_get_var(ncINPid, varid, knew(ng:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, 'knew', InpRec, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
# endif
# if defined SEDIMENT && defined SED_MORPH
!
!  Read in time-evolving bathymetry (m).
!
        IF (get_var(idbath)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idbath)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#  endif
     &                      GRID(ng) % h(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idbath)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idbath)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in nonlinear free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          IF (Perfect2d) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2d) THEN
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, 3,          &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        OCEAN(ng) % zeta(LBi,LBj,1))
          ELSE
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
# endif
     &                        OCEAN(ng) % zeta(LBi,LBj,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of free-surface.
!
        IF (get_var(idRzet).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRzet)),   &
     &                        varid)
          gtype=var_flag(varid)*r3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, 2,            &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % rzeta(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRzet)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRzet)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 2D momentum component (m/s) in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*u3dvar
          ELSE
            gtype=var_flag(varid)*u2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, 3,          &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
# endif
     &                        OCEAN(ng) % ubar(LBi,LBj,1))
          ELSE
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
# endif
     &                        OCEAN(ng) % ubar(LBi,LBj,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 2D momentum component in the XI-direction.
!
        IF (get_var(idRu2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu2d)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, 2,            &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % rubar(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRu2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRu2d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 2D momentum component (m/s) in the ETA-direction.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*v3dvar
          ELSE
            gtype=var_flag(varid)*v2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, 3,          &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
# endif
     &                        OCEAN(ng) % vbar(LBi,LBj,1))
          ELSE
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
# endif
     &                        OCEAN(ng) % vbar(LBi,LBj,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS 2D momentum component in the ETA-direction.
!
        IF (get_var(idRv2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv2d)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, 2,            &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % rvbar(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRv2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRv2d)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifdef ADJUST_WSTRESS
!
!  Read in surface U-momentum stress.
!
        IF (get_var(idUsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % ustr(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in surface V-momentum stress.
!
        IF (get_var(idVsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % vstr(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in nonlinear 3D momentum component (m/s) in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          IF (Perfect3D) THEN
            status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,&
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % u(LBi,LBj,1,1))
          ELSE
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % u(LBi,LBj,1,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 3D momentum component in the XI-direction.
!
        IF (get_var(idRu3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu3d)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,  &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % ru(LBi,LBj,0,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRu3d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRu3d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 3D momentum component (m/s) in the ETA-direction.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          IF (Perfect3D) THEN
            status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,&
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % v(LBi,LBj,1,1))
          ELSE
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % v(LBi,LBj,1,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 3D momentum component in the ETA-direction.
!
        IF (get_var(idRv3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv3d)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,  &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % rv(LBi,LBj,0,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRv3d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRv3d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            IF (Perfect3D) THEN
              status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec,      &
     &                          gtype, Vsize,                           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask(LBi,LBj),              &
#  endif
     &                          OCEAN(ng) % t(LBi,LBj,1,1,itrc))
            ELSE
              status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec,      &
     &                          gtype, Vsize,                           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask(LBi,LBj),              &
#  endif
     &                          OCEAN(ng) % t(LBi,LBj,1,Tindex,itrc))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            IF (itrc.eq.itemp) THEN
              scale=1.0_r8/(rho0*Cp)        ! W/m2 to Celsius m/s
            ELSE
              scale=1.0_r8
            END IF
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),  &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        FORCES(ng) % tflux(LBi,LBj,1,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF             
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTsur(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING || defined LMD_MIXING
!
!  Read in vertical viscosity.
!
        IF (have_var(idVvis)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvis)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, Fmin,Fmax,                            &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % AKv(LBi,LBj,0))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvis)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvis)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in temperature vertical diffusion.
!
        IF (have_var(idTdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idTdif)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, Fmin,Fmax,                            &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % AKt(LBi,LBj,0,itemp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idTdif)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idTdif)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in salinity vertical diffusion.
!
        IF (have_var(idSdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idSdif)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, Fmin,Fmax,                            &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % AKt(LBi,LBj,0,isalt))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idSdif)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idSdif)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Read in turbulent kinetic energy.
!
        IF (get_var(idMtke).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idMtke)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,  &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % tke(LBi,LBj,0,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idMtke)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idMtke)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in turbulent kinetic energy time length scale.
!
        IF (get_var(idMtls).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idMtls)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread4d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,  &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % gls(LBi,LBj,0,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idMtls)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idMtls)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in vertical mixing turbulent length scale.
!
        IF (get_var(idVmLS).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVmLS)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % Lscale(LBi,LBj,0))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVmLS)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVmLS)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in turbulent kinetic energy vertical diffusion coefficient.
!
        IF (get_var(idVmKK).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVmKK)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#   endif
     &                      MIXING(ng) % Akk(LBi,LBj,0))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVmKK)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVmKK)), Fmin, Fmax
            END IF
          END IF
        END IF
#   ifdef GLS_MIXING
!
!  Read in turbulent length scale vertical diffusion coefficient.
!
        IF (get_var(idVmKP).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVmKP)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#    ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#    endif
     &                      MIXING(ng) % Akp(LBi,LBj,0))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVmKP)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVmKP)), Fmin, Fmax
            END IF
          END IF
        END IF
#   endif
#  endif
#  ifdef SEDIMENT
!
!  Read in nonlinear sediment fraction of each size class in each bed
!  layer.
!
        DO i=1,NST
          IF (get_var(idfrac(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idfrac(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % bed_frac(LBi,LBj,1,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idfrac(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
!  Read in nonlinear sediment mass of each size class in each bed layer.
!
          IF (get_var(idBmas(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBmas(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % bed_mass(LBi,LBj,1,Tindex,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idBmas(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
!
!  Read in nonlinear sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (get_var(idSbed(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idSbed(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % bed(LBi,LBj,1,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idSbed(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#   ifdef BEDLOAD
!
!  Read in nonlinear sediment fraction of bed load.
!
        DO i=1,NST
          IF (get_var(idUbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idUbld(i))), varid)
            gtype=var_flag(varid)*u2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#    endif
     &                        OCEAN(ng) % bedldu(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idUbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#    endif
     &                        OCEAN(ng) % bedldv(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idVbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#   endif
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Read in nonlinear sediment properties in exposed bed layer.
!
        DO i=1,MBOTP
          IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBott(i))), varid)
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % bottom(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idBott(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idBott(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
# endif
      END IF NLM_STATE
#endif

#if defined TANGENT || defined TL_IOMS
!
!-----------------------------------------------------------------------
!  Read in tangent linear state variables.
!-----------------------------------------------------------------------
!
      TLM_STATE: IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
!
!  Read in tangent linear free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % tl_zeta(LBi,LBj,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tangent linear 2D momentum component (m/s) in the
!  XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % tl_ubar(LBi,LBj,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifdef ADJUST_WSTRESS
!
!  Read in tangent linear surface U-momentum stress.
!
        IF (get_var(idUsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % tl_ustr(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tangent linear surface V-momentum stress.
!
        IF (get_var(idVsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % tl_vstr(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in tangent linear 2D momentum component (m/s) in the
!  ETA-direction.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % tl_vbar(LBi,LBj,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef SOLVE3D
!
!  Read in tangent linear 3D momentum component (m/s) in the
!  XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % tl_u(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tangent linear 3D momentum component (m/s) in the
!  ETA-direction.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % tl_v(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tangent linear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % tl_t(LBi,LBj,1,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF             
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#  ifdef ADJUST_STFLUX
!
!  Read in tangent linear surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            IF (itrc.eq.itemp) THEN
              scale=1.0_r8/(rho0*Cp)        ! W/m2 to Celsius m/s
            ELSE
              scale=1.0_r8
            END IF
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),  &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        FORCES(ng) % tl_tflux(LBi,LBj,1,          &
     &                                              Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF             
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTsur(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
#  ifdef SEDIMENT
!
!  Read in tangent linear sediment fraction of each size class in each
!  bed layer.
!
        DO i=1,NST
          IF (get_var(idfrac(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idfrac(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % tl_bed_frac(LBi,LBj,1,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idfrac(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
!  Read in tangent linear sediment mass of each size class in each
!  bed layer.
!
          IF (get_var(idBmas(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBmas(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng)%tl_bed_mass(LBi,LBj,1,Tindex,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idBmas(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
!
!  Read in tangent linear sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (get_var(idSbed(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idSbed(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % tl_bed(LBi,LBj,1,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idSbed(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#   ifdef BEDLOAD
!
!  Read in tangent linear sediment fraction of bed load.
!
        DO i=1,NST
          IF (get_var(idUbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idUbld(i))), varid)
            gtype=var_flag(varid)*u2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#    endif
     &                        OCEAN(ng) % tl_bedldu(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idUbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#    endif
     &                        OCEAN(ng) % tl_bedldv(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idVbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#   endif
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Read in tangent linear sediment properties in exposed bed layer.
!
        DO i=1,MBOTP
          IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBott(i))), varid)
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % tl_bottom(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idBott(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idBott(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
# endif
      END IF TLM_STATE
#endif

#ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Read in adjoint state variables.
!-----------------------------------------------------------------------
!
      ADM_STATE: IF (model.eq.iADM) THEN
!
!  Read in adjoint free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % ad_zeta(LBi,LBj,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in adjoint 2D momentum component (m/s) in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % ad_ubar(LBi,LBj,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in adjoint 2D momentum component (m/s) in the ETA-direction.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % ad_vbar(LBi,LBj,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifdef ADJUST_WSTRESS
!
!  Read in adjoint linear surface U-momentum stress.
!
        IF (get_var(idUsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % ad_ustr(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in adjoint linear surface V-momentum stress.
!
        IF (get_var(idVsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % ad_vstr(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in adjoint 3D momentum component (m/s) in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % ad_u(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in adjoint 3D momentum component (m/s) in the ETA-direction.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % ad_v(LBi,LBj,1,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in adjoint tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % ad_t(LBi,LBj,1,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#  ifdef ADJUST_STFLUX
!
!  Read in adjoint surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            IF (itrc.eq.itemp) THEN
              scale=1.0_r8/(rho0*Cp)        ! W/m2 to Celsius m/s
            ELSE
              scale=1.0_r8
            END IF
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),  &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        FORCES(ng) % ad_tflux(LBi,LBj,1,          &
     &                                              Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF             
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTsur(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif

#  ifdef SEDIMENT
!
!  Read in adjoint sediment fraction of each size class in each bed
!  layer.
!
        DO i=1,NST
          IF (get_var(idfrac(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idfrac(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % ad_bed_frac(LBi,LBj,1,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idfrac(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
!  Read in adjoint sediment mass of each size class in each bed layer.
!
          IF (get_var(idBmas(i))) THEN
            foundit=find_string(var_name, n_var,
     &                          TRIM(Vname(1,idBmas(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng)%ad_bed_mass(LBi,LBj,1,Tindex,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idBmas(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
!
!  Read in adjoint sediment properties in each bed layer.
!
        DO i=1,MBEDP
          IF (get_var(idSbed(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idSbed(i))), varid)
            gtype=var_flag(varid)*b3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, Nbed,       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % ad_bed(LBi,LBj,1,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idSbed(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#   ifdef BEDLOAD
!
!  Read in adjoint sediment fraction of bed load.
!
        DO i=1,NST
          IF (get_var(idUbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idUbld(i))), varid)
            gtype=var_flag(varid)*u2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask(LBi,LBj),                &
#    endif
     &                        OCEAN(ng) % ad_bedldu(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idUbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask(LBi,LBj),                &
#    endif
     &                        OCEAN(ng) % ad_bedldv(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idVbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#   endif
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Read in adjoint sediment properties in exposed bed layer.
!
        DO i=1,MBOTP
          IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBott(i))), varid)
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        OCEAN(ng) % ad_bottom(LBi,LBj,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idBott(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idBott(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
# endif
      END IF ADM_STATE
#endif

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in background-error covariance normalization (nondimensional)
!  factors.
!-----------------------------------------------------------------------
!
      NRM_STATE: IF (model.eq.5) THEN
!
!  Read in free-surface norm.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % b_zeta(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D momentum norm in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % b_ubar(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D momentum norm in the ETA-direction.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % b_vbar(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_WSTRESS
!
!  Read in surface U-momentum stress norm.
!
        IF (get_var(idUsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % b_sustr(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in surface V-momentum stress norm.
!
        IF (get_var(idVsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % b_svstr(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in 3D momentum norm in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % b_u(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 3D momentum norm in the ETA-direction.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % b_v(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tracer type variables norm.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % b_t(LBi,LBj,1,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  ifdef ADJUST_STFLUX
!
!  Read in surface tracer flux norm.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        FORCES(ng) % b_stflx(LBi,LBj,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTsur(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
# endif
      END IF NRM_STATE
#endif

#if defined FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in backgound-error standard deviation factors.
!-----------------------------------------------------------------------
!
      STD_STATE: IF (model.eq.6) THEN
!
!  Read in free-surface standard deviation.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % e_zeta(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D momentum standard deviation in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % e_ubar(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D momentum standard deviation in the ETA-direction.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % e_vbar(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_WSTRESS
!
!  Read in surface U-momentum stress standard deviation.
!
        IF (get_var(idUsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % e_sustr(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in surface V-momentum stress standard deviation.
!
        IF (get_var(idVsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      FORCES(ng) % e_svstr(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in 3D momentum standard deviation in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % e_u(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 3D momentum standard deviation in the ETA-direction.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % e_v(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tracer type variables standard deviation.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % e_t(LBi,LBj,1,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  ifdef ADJUST_STFLUX
!
!  Read in surface tracer flux standard deviations.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r2dvar
            IF (itrc.eq.itemp) THEN
              scale=1.0_r8/(rho0*Cp)        ! W/m2 to Celsius m/s
            ELSE
              scale=1.0_r8
            END IF
            status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj,                &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#   endif
     &                        FORCES(ng) % e_stflx(LBi,LBj,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTsur(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
# endif
!
!  Read in convolution horizontal diffusion coefficients.
!
        IF (have_var(idKhor)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idKhor)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, KhMin(ng), KhMax(ng),                 &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      MIXING(ng) % Kh(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idKhor)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idKhor)), KhMin(ng),      &
     &                           KhMax(ng)
            END IF
          END IF
        END IF

# ifdef SOLVE3D
!
!  Read in convolution vertical diffusion coefficient.
!
        IF (have_var(idKver)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idKver)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                      Fscl, KvMin(ng), KvMax(ng),                 &
#  ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
#  endif
     &                      MIXING(ng) % Kv(LBi,LBj,0))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idKver)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idKver)), KvMin(ng),      &
     &                           KvMax(ng)
            END IF
          END IF
        END IF
# endif
      END IF STD_STATE
#endif

#if defined IMPULSE
!
!-----------------------------------------------------------------------
!  Read in adjoint model or tangent linear model impulse forcing terms.
!-----------------------------------------------------------------------
!
      FRC_STATE: IF (model.eq.7) THEN
!
!  Set number of records available.
!
        NrecFrc(ng)=Nrec
!
!  Read in next impulse forcing time to process.
!
        IF (InpThread) THEN
          start(1)=InpRec
          total(1)=1
          status=nf90_inq_varid(ncINPid, TRIM(tvarnam), varid)
          status=nf90_get_var(ncINPid, varid, FrcTime(ng:),             &
     &                        start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) string, TRIM(tvarnam), InpRec,            &
     &                        TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, IDmod, exit_flag, 1)
# endif
        IF (exit_flag.ne.NoError) RETURN
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, IDmod, FrcTime(ng), 1)
# endif
!
!  Read in free-surface forcing.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % f_zeta(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D momentum forcing in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % f_ubar(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D momentum forcing in the ETA-direction.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      OCEAN(ng) % f_vbar(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef SOLVE3D
!
!  Read in 3D momentum forcing in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % f_u(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 3D momentum norm in the ETA-direction.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype,   &
     &                      Vsize, LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
#  endif
     &                      OCEAN(ng) % f_v(LBi,LBj,1))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tracer type variables norm.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncINPid, varid, InpRec, gtype, &
     &                        Vsize, LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask(LBi,LBj),                &
#  endif
     &                        OCEAN(ng) % f_t(LBi,LBj,1,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
# endif
      END IF FRC_STATE
#endif
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      status=nf90_close(ncINPid)
!
  20  FORMAT (/,a,'GET_STATE - unable to open input NetCDF file: ',a)
  30  FORMAT (/,a,'GET_STATE - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,18x,'in input NetCDF file: ',a)
  40  FORMAT (/,a,'GET_STATE - requested input time record = ',i3,/,    &
     &        18x,'not found in input NetCDF: ',a,/,                    &
     &        18x,'number of available records = ',i3)
  50  FORMAT (/,a,'GET_STATE - error while reading attribute: ',a)
  60  FORMAT (/,a,'GET_STATE - error while inquiring attributes',       &
     &          ' for variable: ',a)
  70  FORMAT (/,a,'GET_STATE - cannot inquire about time variable',/,   &
     &        18x,'in input NetCDF file: ',a)
  80  FORMAT (/,a,'GET_STATE - ',a,t62,'t = ',a,                        &
     &        /,19x,'(Iter=',i4.4,', File: ',a,', Rec=',i4.4,           &
     &        ', Index=',i1,')')
  90  FORMAT (/,a,'GET_STATE - ',a,t62,'t = ',a,                        &
     &        /,19x,'(File: ',a,', Rec=',i4.4,', Index=',i1,')')
 100  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
      RETURN
      END SUBROUTINE get_state
