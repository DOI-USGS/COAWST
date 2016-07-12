#include "cppdefs.h"
      SUBROUTINE get_state (ng, model, msg, ncname, IniRec, Tindex)
!
!svn $Id: get_state.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!     msg        Message index for StateMsg.                           !
!     ncname     Nonlinear initial conditions NetCDF file name.        !
!     IniRec     Nonlinear initial conditions time record to read.     !
!     Tindex     State variable time index to intialize.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#if defined ADJUST_BOUNDARY
      USE mod_boundary
#endif
      USE mod_grid
      USE mod_iounits
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
#endif
# if defined PERFECT_RESTART && defined ICE_MODEL
      USE mod_forces
# endif
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
#if defined ESTUARYBGC && defined SAV_BIOMASS
      USE mod_biology
#endif 
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
#if defined VEGETATION 
      USE mod_vegarr
      USE mod_vegetation
#endif 
#if defined ICE_MODEL || defined CICE_MODEL
      USE mod_ice
#endif
#ifdef NCEP_FLUXES
      USE mod_forces
#endif
      USE mod_stepping
      USE mod_strings
#ifdef FILTERED_RST
      USE mod_filter, ONLY: FILN
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
#endif
#ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
# ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
# endif
#endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE nf_fread3d_mod, ONLY : nf_fread3d
#ifdef SOLVE3D
      USE nf_fread4d_mod, ONLY : nf_fread4d
#endif
      USE strings_mod, ONLY : find_string
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
      logical :: Perfect2D, Perfect3D, foundit, read_ice
#if defined ADJUST_BOUNDARY || \
    defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
      logical :: get_adjust
#endif
      logical, dimension(NV) :: get_var, have_var

      integer :: LBi, UBi, LBj, UBj
#ifdef ADJUST_BOUNDARY
      integer :: IorJ, LBij, UBij
#endif
      integer :: IDmod, InpRec, gtype, i, ifield, itrc, lstr, lend
      integer :: Nrec, mySize, ncINPid, nvatts, nvdim, status, varid
      integer :: Vsize(4), start(4), total(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, INPtime, Tmax, scale, time_scale

      real(r8), allocatable :: TimeVar(:)

      character (len=6 ) :: string
      character (len=14) :: t_code
      character (len=15) :: attnam, tvarnam
      character (len=40) :: tunits

#ifdef ADJUST_BOUNDARY
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
      IorJ=IOBOUNDS(ng)%IorJ
#endif
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Set model identification string.
!
      IF (model.eq.iNLM.or.(model.eq.0)) THEN
        string=' NLM: '                    ! nonlinear model, restart
        IDmod=iNLM
      ELSE IF (model.eq.iTLM) THEN
        string=' TLM: '                    ! tangent linear model
        IDmod=iTLM
      ELSE IF (model.eq.iRPM) THEN
        string=' RPM: '                    ! representer model
        IDmod=iRPM
      ELSE IF (model.eq.iADM) THEN
        string=' ADM: '                    ! adjoint model
        IDmod=iADM
      ELSE IF (model.eq.5) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! model or initial conditions
      ELSE IF (model.eq.6) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! model or initial conditions
      ELSE IF (model.eq.7) THEN
        string=' FRC: '                    ! impulse forcing
        IDmod=iNLM
      ELSE IF (model.eq.8) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.9) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.10) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.11) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.12) THEN
        string=' NLM: '                    ! tangent linear forcing and
        IDmod=iNLM                         ! obc increments
      ELSE IF (model.eq.14) THEN
        string=' TLM: '                    ! tangent model error forcing
        IDmod=iTLM
      END IF
!
!  Set switch to process variables for nonlinear model perfect restart.
!
      Perfect2D=.FALSE.
      Perfect3D=.FALSE.
#ifdef PERFECT_RESTART
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        Perfect2D=.TRUE.
        Perfect3D=.TRUE.
      END IF
#endif
      PerfectRST(ng)=Perfect2D.or.Perfect3D
!
!  Determine variables to read.
!
      CALL checkvars (ng, model, ncname, string, Nrec, NV, tvarnam,     &
     &                get_var, have_var)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
      SourceFile='get_state.F'
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
      CALL netcdf_open (ng, IDmod, ncname, 0, ncINPid)
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF

#ifdef NO_LBC_ATT
!
!  Lateral boundary conditions attribute not checked in restart file.
!
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        IF (Master) WRITE (stdout,20) string, 'NLM_LBC', TRIM(ncname)
      END IF
#else
!
!  If restart, read in lateral boundary conditions global attribute
!  from restart file and check keyword strings with structure vlues
!  for consistency.
!
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        CALL lbc_getatt (ng, model, ncINPid, ncname, 'NLM_LBC', LBC)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
!
!  Inquire about the input time variable.
!
      CALL netcdf_inq_var (ng, IDmod, ncname,                           &
     &                     MyVarName = TRIM(tvarnam),                   &
     &                     VarID = varid,                               &
     &                     nVarDim =  nvdim,                            &
     &                     nVarAtt = nvatts)
      IF (exit_flag.ne.NoError) RETURN
!
!  Allocate input time variable and read its value(s).  Recall that
!  input time variable is a one-dimensional array with one or several
!  values.
!
      mySize=var_Dsize(1)
      IF (.not.allocated(TimeVar)) allocate (TimeVar(mySize))
      CALL netcdf_get_fvar (ng, IDmod, ncname, TRIM(tvarnam), TimeVar)
      IF (exit_flag.ne.NoError) RETURN
!
!  If using the latest time record from input NetCDF file as the
!  initialization record, assign input time.
!
      IF (LastRec(ng)) THEN
        Tmax=-1.0_r8
        DO i=1,mySize
          IF (TimeVar(i).gt.Tmax) THEN
            Tmax=TimeVar(i)
            IniRec=i
          END IF
        END DO
        INPtime=Tmax
        InpRec=IniRec
      ELSE
        IF ((IniRec.ne.0).and.(IniRec.gt.mySize)) THEN
          IF (Master)  WRITE (stdout,30) string, IniRec, TRIM(ncname),  &
     &                                   mySize
          exit_flag=2
          RETURN
        END IF
        IF (IniRec.ne.0) THEN
          InpRec=IniRec
        ELSE
          InpRec=1
        END IF
        INPtime=TimeVar(InpRec)
      END IF
      IF (allocated(TimeVar)) deallocate ( TimeVar )
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
!  Set starting time index and time clock in days.  Notice that the
!  global time variables and indices are only over-written when
!  processing initial conditions (msg = 1).
!
      IF ((model.eq.0).or.(model.eq.iNLM).or.                           &
     &    (model.eq.iTLM).or.(model.eq.iRPM)) THEN
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(dstart*day2sec))) THEN
          INPtime=dstart*day2sec
        END IF
        IF (msg.eq.1) THEN            ! processing initial conditions
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
          ntstart(ng)=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
          IF (ntstart(ng).lt.1) ntstart(ng)=1
          IF (PerfectRST(ng)) THEN
            ntfirst(ng)=1
          ELSE
            ntfirst(ng)=ntstart(ng)
          END IF
        END IF
#ifdef WEAK_CONSTRAINT
        IF (msg.eq.4) THEN
          ForceTime(ng)=time(ng)
        END IF
#endif
      ELSE IF (model.eq.iADM) THEN
        IF ((msg.eq.1).and.(INPtime.eq.0.0_r8)) THEN
          INPtime=time(ng)
        ELSE IF (msg.ne.1) THEN
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
#ifdef FILTERED
      filtindx = InpRec
#endif
!
!  Set current input time, io_time .  Notice that the model time,
!  time(ng), is reset above.  This is a THREADPRIVATE variable in
!  shared-memory and this routine is only processed by the MASTER
!  thread since it is an I/O routine. Therefore, we need to update
!  time(ng) somewhere else in a parallel region. This will be done
!  with io_time variable.
!
      io_time=INPtime
!
!  Report information.
!
      lstr=SCAN(ncname,'/',BACK=.TRUE.)+1
      lend=LEN_TRIM(ncname)
      IF (Master) THEN
        CALL time_string (INPtime, t_code)
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,40) string, TRIM(StateMsg(msg)), t_code, ng,    &
     &                      Nrun, ncname(lstr:lend), InpRec, Tindex
        ELSE
          WRITE (stdout,50) string, TRIM(StateMsg(msg)), t_code, ng,    &
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
      NLM_STATE: IF ((model.eq.iNLM).or.(model.eq.0).or.                &
     &               (model.eq.13)) THEN

# ifdef PERFECT_RESTART
!
!  Read in time-stepping indices.
!
!        IF ((model.eq.0).and.(nrrec(ng).ne.0)) THEN
        IF (nrrec(ng).ne.0) THEN
#  ifdef SOLVE3D
          CALL netcdf_get_ivar (ng, IDmod, ncname, 'nstp',              &
     &                          nstp(ng:),                              &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'nrhs',              &
     &                          nrhs(ng:),                              &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'nnew',              &
     &                          nnew(ng:),                              &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN
#  endif
          CALL netcdf_get_ivar (ng, IDmod, ncname, 'kstp',              &
     &                          kstp(ng:),                              &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN


          CALL netcdf_get_ivar (ng, IDmod, ncname, 'krhs',              &
     &                          krhs(ng:),                              &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'knew',              &
     &                          knew(ng:),                              &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN
# ifdef ICE_MODEL
          CALL netcdf_get_ivar (ng, IDmod, ncname, 'linew',             &
     &                          linew(ng:),                             &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'liold',             &
     &                          liold(ng:),                             &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'liunw',             &
     &                          liunw(ng:),                             &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'liuol',             &
     &                          liuol(ng:),                             &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'lienw',             &
     &                          lienw(ng:),                             &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_get_ivar (ng, IDmod, ncname, 'lieol',             &
     &                          lieol(ng:),                             &
     &                          ncid = ncINPid,                         &
     &                          start = (/InpRec/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) RETURN
# endif
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
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idbath), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      GRID(ng) % h)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idbath)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idbath)), Fmin, Fmax
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
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idFsur), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 3,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        OCEAN(ng) % zeta)
          ELSE
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idFsur), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        OCEAN(ng) % zeta(:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
# if defined WET_DRY   && defined PERFECT_RESTART
!
!  Read in wet/dry mask in RHO-points.
!
        IF (get_var(idRwet)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRwet)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                        Vname(1,idRwet), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % rmask_wet)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRwet)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRwet)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in wet/dry mask in U-points.
!
        IF (get_var(idUwet)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUwet)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                        Vname(1,idUwet), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % umask_wet)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUwet)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUwet)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in wet/dry mask in V-points.
!
        IF (get_var(idVwet)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVwet)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                        Vname(1,idVwet), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % vmask_wet)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVwet)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVwet)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in nonlinear RHS of free-surface.
!
        IF (get_var(idRzet).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRzet)),   &
     &                        varid)
          gtype=var_flag(varid)*r3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRzet), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % rzeta)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRzet)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRzet)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 2D U-momentum component (m/s).
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 3,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
     &                        OCEAN(ng) % ubar)
          ELSE
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
     &                        OCEAN(ng) % ubar(:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 2D U-momentum component.
!
        IF (get_var(idRu2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu2d)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRu2d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
# endif
     &                      OCEAN(ng) % rubar)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRu2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRu2d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 2D U-momentum component (m/s).
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 3,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
     &                        OCEAN(ng) % vbar)
          ELSE
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask,                         &
# endif
     &                        OCEAN(ng) % vbar(:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS 2D U-momentum component.
!
        IF (get_var(idRv2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv2d)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRv2d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask,                           &
# endif
     &                      OCEAN(ng) % rvbar)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRv2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRv2d)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifdef SOLVE3D
!
!  Read in nonlinear 3D U-momentum component (m/s).
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          IF (Perfect3D) THEN
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        OCEAN(ng) % u)
          ELSE
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        OCEAN(ng) % u(:,:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 3D U-momentum component.
!
        IF (get_var(idRu3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu3d)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread4d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRu3d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % ru)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRu3d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRu3d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 3D V-momentum component (m/s).
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          IF (Perfect3D) THEN
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        OCEAN(ng) % v)
          ELSE
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        OCEAN(ng) % v(:,:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 3D V-momentum component.
!
        IF (get_var(idRv3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv3d)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread4d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRv3d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % rv)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRv3d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRv3d)), Fmin, Fmax
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
              status=nf_fread4d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          OCEAN(ng) % t(:,:,:,:,itrc))
            ELSE
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          OCEAN(ng) % t(:,:,:,Tindex,itrc))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO

#  ifdef NEMURO_SED1
!
!  Read in bio sediment variables.
!
        IF (get_var(idPONsed)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idPONsed)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idPONsed), varid,                   &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % PONsed)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idPONsed)),        &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idPONsed)),               &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idOPALsed)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idOPALsed)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idOPALsed), varid,                  &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % OPALsed)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idOPALsed)),       &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idOPALsed)),              &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idDENITsed)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idDENITsed)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idDENITsed), varid,                 &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % DENITsed)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idDENITsed)),      &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idDENITsed)),             &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idPONbur)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idPONbur)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idPONbur), varid,                   &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % PON_burial)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idPONbur)),        &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idPONbur)),               &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idOPALbur)) THEN
          foundit=find_string(var_name, n_var,                          &
     &                          TRIM(Vname(1,idOPALbur)), varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idOPALbur), varid,                  &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      OCEAN(ng) % OPAL_burial)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idOPALbur)),       &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idOPALbur)),              &
     &                           Fmin, Fmax
            END IF
          END IF
        END IF
#  endif

#  if defined GLS_MIXING || defined MY25_MIXING || defined LMD_MIXING
!
!  Read in vertical viscosity.
!
        IF (have_var(idVvis)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvis)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvis), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % AKv)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvis)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvis)), Fmin, Fmax
            END IF
          END IF
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % AKv)
#   endif
        END IF
!
!  Read in temperature vertical diffusion.
!
        IF (have_var(idTdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idTdif)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idTdif), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % AKt(:,:,:,itemp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idTdif)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idTdif)), Fmin, Fmax
            END IF
          END IF
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % AKt(:,:,:,itemp))
#   endif
        END IF
#   ifdef SALINITY
!
!  Read in salinity vertical diffusion.
!
        IF (have_var(idSdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idSdif)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idSdif), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
#    ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#    endif
     &                      MIXING(ng) % AKt(:,:,:,isalt))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idSdif)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idSdif)), Fmin, Fmax
            END IF
          END IF
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % AKt(:,:,:,isalt))
#   endif
        END IF
#   endif
#  endif
#  if defined LMD_SKPP
!
!  Read in Hsbl
!
        IF (have_var(idHsbl).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idHsbl)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idHsbl), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      MIXING(ng) % Hsbl)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idHsbl)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idHsbl)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
#  if defined LMD_BKPP
!
!  Read in Hbbl
!
        IF (have_var(idHbbl).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idHbbl)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idHbbl), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      MIXING(ng) % Hbbl)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idHbbl)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idHbbl)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
#  if defined LMD_NONLOCAL && defined PERFECT_RESTART
!
!  Read in Ghats
!
      DO itrc=1,NAT
        IF (have_var(idGhat(itrc))) THEN
          foundit=find_string(var_name, n_var,                          &
     &                        TRIM(Vname(1,idGhat(itrc))), varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idGhat(itrc)), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % Ghats(:,:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idGhat(itrc))),    &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idGhat(itrc))), Fmin, Fmax
            END IF
          END IF
#   ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % Ghats(:,:,:,itrc))
#   endif
        END IF
      END DO
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Read in turbulent kinetic energy.
!
        IF (get_var(idMtke).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idMtke)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread4d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idMtke), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % tke)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idMtke)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idMtke)), Fmin, Fmax
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
          status=nf_fread4d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idMtls), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % gls)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idMtls)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idMtls)), Fmin, Fmax
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
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVmLS), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % Lscale)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVmLS)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVmLS)), Fmin, Fmax
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
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVmKK), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#   ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#   endif
     &                      MIXING(ng) % Akk)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVmKK)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVmKK)), Fmin, Fmax
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
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVmKP), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#    ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#    endif
     &                      MIXING(ng) % Akp)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVmKP)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVmKP)), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idfrac(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % bed_frac(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idfrac(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idBmas(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % bed_mass(:,:,:,Tindex,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBmas(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSbed(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % bed(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idSbed(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
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
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbld(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask,                         &
#    endif
     &                        SEDBED(ng) % bedldu(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbld(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#    endif
     &                        SEDBED(ng) % bedldv(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#   endif
#  endif
!
#if defined ESTUARYBGC && defined SAV_BIOMASS
!
!  Read in AGB and BGB values
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idsagb)),     & 
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idsagb), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                 OCEAN(ng)%AGB)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,idsagb )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,idsagb)), Fmin, Fmax 
          END IF
        END IF
!
! BGB 
!  
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idsbgb)),     & 
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idsbgb), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &               OCEAN(ng)%BGB)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,idsbgb )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,idsbgb)), Fmin, Fmax 
          END IF
        END IF
!
! DINsed
!  
        foundit=find_string(var_name, n_var, TRIM(Vname(1,iddins)),     & 
     &                        varid)
        gtype=var_flag(varid)*r3dvar
        status=nf_fread3d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,iddins), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,  1, N(ng),             &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &               OCEAN(ng)%DINsed)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,iddins )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,iddins)), Fmin, Fmax 
          END IF
        END IF
!
#  endif 
!
#  ifdef VEGETATION 
#   if defined VEG_DRAG || defined VEG_BIOMASS
!
!  Read in plant properties for each vegetation type 
!
        DO i=1,NVEGP
          IF (get_var(idvprp(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idvprp(i))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idvprp(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, NVEG,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        VEG(ng) % plant(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idvprp(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idvprp(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO 
#   endif 
#   ifdef MARSH_WAVE_THRUST
!
!  Read in masking for determining marsh interface 
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTims)),     & 
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTims), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                 VEG(ng) % marsh_mask)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,60) string, TRIM(Vname(1,idTims )),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN 
            WRITE (stdout,70) TRIM(Vname(2,idTims)), Fmin, Fmax 
          END IF
        END IF
!
#   endif 
#  endif 
!
#  if defined SEDIMENT || defined BBL_MODEL
!
!  Read in nonlinear sediment properties in exposed bed layer.
!
        DO i=1,MBOTP
          IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBott(i))), varid)
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idBott(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % bottom(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBott(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBott(i))), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
#if defined BEST_NPZ && defined BERING_10K
    
!  Read in benthic initial conditions.
        DO itrc=1,1 !,NBeT(ng)
          IF (get_var(idBvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idBvar(itrc))), varid)
            gtype=var_flag(varid)*r2dvar
            IF (Perfect2D) THEN
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idBvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, 3,               &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          OCEAN(ng) % bt(:,:,1,:,itrc))
            ELSE
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idBvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          OCEAN(ng) % bt(:,:,1,Tindex,itrc))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
#  ifdef ICE_MODEL
        read_ice = .true.
#   ifdef ANA_ICE
        IF (nrrec(ng).eq.0) THEN
          read_ice = .false.
        END IF
#   endif
!
!  Read in 2D ice momentum component (m/s) in the XI-direction.
!
      IF (read_ice) THEN
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idUice)),     &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*u3dvar
          ELSE
            gtype=var_flag(varid)*u2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUice), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
     &                        ICE(ng) % ui)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idUice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                ICE(ng) % ui(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idUice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idUice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D ice momentum component (m/s) in the ETA-direction.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idVice)),     &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*v3dvar
          ELSE
            gtype=var_flag(varid)*v2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVice), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % vi)
          ELSE
        gtype=var_flag(varid)*v2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idVice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                ICE(ng) % vi(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idVice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idVice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice concentration.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idAice)),     &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idAice), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % ai)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idAice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % ai(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idAice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idAice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice average thickness.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idHice)),     &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idHice), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % hi)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idHice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % hi(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idHice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idHice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in snow average thickness.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idHsno)),     &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idHsno), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % hsn)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idHsno), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % hsn(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idHsno)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idHsno)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice/snow surface temperature.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTice)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % tis)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idTice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idTice)), Fmin, Fmax
          END IF
        END IF
#ifdef MELT_PONDS
!
!  Read in surface water fraction (on ice).
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idApond)),    &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idApond), varid,                  &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % apond)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idApond), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % apond(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idApond)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idApond)), Fmin, Fmax
          END IF
        END IF
!
!  Read in surface water thickness (on ice).
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idHpond)),    &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idHpond), varid,                  &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % hpond)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idHpond), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % hpond(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idHpond)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idHpond)), Fmin, Fmax
          END IF
        END IF
#endif
!
!  Read in ice age.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idAgeice)),   &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idAgeice), varid,                 &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % ageice)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idAgeice), varid, InpRec, gtype,          &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % ageice(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idAgeice)),          &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idAgeice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice internal temperature.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTimid)),    &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTimid), varid,                  &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % ti)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTimid), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % ti(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idTimid)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idTimid)), Fmin, Fmax
          END IF
        END IF
!
!  Read in internal ice stress component 11
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSig11)),    &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSig11), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % sig11)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSig11), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % sig11(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idSig11)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idSig11)), Fmin, Fmax
          END IF
        END IF
!
!  Read in internal ice stress component 22
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSig22)),    &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSig22), varid,                  &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % sig22)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSig22), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % sig22(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idSig22)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idSig22)), Fmin, Fmax
          END IF
        END IF
!
!  Read in internal ice stress component 12
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSig12)),    &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSig12), varid,                  &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 2,                 &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ICE(ng) % sig12)
          ELSE
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSig12), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % sig12(:,:,Tindex))
          END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idSig12)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idSig12)), Fmin, Fmax
          END IF
        END IF
#ifndef ICE_MOM_BULK
!
!  Read in ice-water friction velocity.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTauiw)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTauiw), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % utau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idTauiw)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idTauiw)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice-water momentum transfer coefficient.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idChuiw)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idChuiw), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % chu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idChuiw)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idChuiw)), Fmin, Fmax
          END IF
        END IF
#endif
!
!  Read in salinity of molecular sub-layer under ice
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idS0mk)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idS0mk), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % s0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idS0mk)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idS0mk)), Fmin, Fmax
          END IF
        END IF
!
!  Read in temperature of molecular sub-layer under ice
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idT0mk)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idT0mk), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                ICE(ng) % t0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,60) string, TRIM(Vname(1,idT0mk)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(2,idT0mk)), Fmin, Fmax
          END IF
        END IF
#   ifdef PERFECT_RESTART
!
!  Read in surface tracers flux.
!
        IF (Perfect2D) THEN
          DO itrc=1,NAT
            IF (get_var(idTsur(itrc))) THEN
              foundit=find_string(var_name, n_var,                      &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
              gtype=var_flag(varid)*r2dvar
              scale=1.0_r8
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTsur(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, Fmin, Fmax,                      &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
     &                          FORCES(ng)% stflx_save(:,:,itrc))
              IF (status.ne.nf90_noerr) THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),&
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
                  WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//       &
     &                               ', saved stflx', Fmin, Fmax
                END IF
              END IF
            END IF
          END DO
!
!  Read in surface U-momentum stress.
!
          IF (get_var(idUsms)) THEN
            foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)), &
     &                          varid)
            gtype=var_flag(varid)*u2dvar
            scale=1.0_r8
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUsms), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % umask,                         &
#  endif
     &                        FORCES(ng) % sustr_save)
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUsms)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idUsms))//               &
     &                             ', saved sustr', Fmin, Fmax
              END IF
            END IF
          END IF
!
!  Read in surface V-momentum stress.
!
          IF (get_var(idVsms)) THEN
            foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)), &
     &                          varid)
            gtype=var_flag(varid)*v2dvar
            scale=1.0_r8
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVsms), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#  endif
     &                        FORCES(ng) % svstr_save)
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVsms)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idVsms))//               &
     &                             ', saved svstr', Fmin, Fmax
              END IF
            END IF
          END IF
        END IF
#   endif
      END IF
!
#  endif
#  ifdef NCEP_FLUXES
!
!  Read in NCEP gustiness squared
!
        IF (get_var(idWg2d)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idWg2d)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idWg2d), varid, InpRec, gtype,            &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % wg2_d)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idWg2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idWg2d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in NCEP momentum transfer coefficient
!
        IF (get_var(idCdd)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idCdd)),    &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idCdd), varid, InpRec, gtype,             &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % cd_d)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idCdd)), InpRec,   &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idCdd)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in NCEP sensible heat transfer coefficient
!
        IF (get_var(idChd)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idChd)),    &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idChd), varid, InpRec, gtype,             &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % ch_d)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idChd)), InpRec,   &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idChd)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in NCEP latent heat transfer coefficient
!
        IF (get_var(idCed)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idCed)),    &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idCed), varid, InpRec, gtype,             &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % ce_d)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idCed)), InpRec,   &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idCed)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in model gustiness squared
!
        IF (get_var(idWg2m)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idWg2m)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idWg2m), varid, InpRec, gtype,            &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % wg2_m)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idWg2m)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idWg2m)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in model momentum transfer coefficient
!
        IF (get_var(idCdm)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idCdm)),    &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idCdm), varid, InpRec, gtype,             &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % cd_m)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idCdm)), InpRec,   &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idCdm)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in model sensible heat transfer coefficient
!
        IF (get_var(idChm)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idChm)),    &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idChm), varid, InpRec, gtype,             &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % ch_m)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idChm)), InpRec,   &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idChm)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in model latent heat transfer coefficient
!
        IF (get_var(idCem)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idCem)),    &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idCem), varid, InpRec, gtype,             &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % ce_m)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idCem)), InpRec,   &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idCem)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in NCEP air density
!
        IF (get_var(idRhoa)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRhoa)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                Vname(1,idRhoa), varid, InpRec, gtype,            &
     &                      Vsize, LBi, UBi, LBj, UBj,                  &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      FORCES(ng) % rhoa_n)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idRhoa)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idRhoa)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 4,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % rhoa_n,                      &
     &                        FORCES(ng) % cd_d,                        &
     &                        FORCES(ng) % ch_d,                        &
     &                        FORCES(ng) % ce_d)
          CALL mp_exchange2d (ng, MyRank, IDmod, 3,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % cd_m,                        &
     &                        FORCES(ng) % ch_m,                        &
     &                        FORCES(ng) % ce_m)
          CALL mp_exchange2d (ng, MyRank, IDmod, 2,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % wg2_d,                       &
     &                        FORCES(ng) % wg2_m)
# endif

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

# if defined ADJUST_BOUNDARY || \
    defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
        IF (inner.eq.0.and.model.eq.iRPM) THEN
          get_adjust=.FALSE.
        ELSE
          get_adjust=.TRUE.
        END IF
# endif
!
!  Read in tangent linear free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idFsur), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % tl_zeta(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
        IF (get_var(idSbry(isFsur)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, r2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_zeta_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in tangent linear 2D U-momentum component (m/s).
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
# endif
     &                      OCEAN(ng) % tl_ubar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
        IF (get_var(idSbry(isUbar)).and.get_adjust.and.                 &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, u2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_ubar_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in tangent linear 2D V-momentum component.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask,                           &
# endif
     &                      OCEAN(ng) % tl_vbar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D V-momentum component open boundaries.
!
        IF (get_var(idSbry(isVbar)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, v2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_vbar_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef ADJUST_WSTRESS
!
!  Read in tangent linear surface U-momentum stress.
!
        IF (get_var(idUsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          scale=1.0_r8
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      FORCES(ng) % tl_ustr(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUsms)),          &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUsms))//                 &
     &                          ', adjusted tl_ustr', Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tangent linear surface V-momentum stress.
!
        IF (get_var(idVsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          scale=1.0_r8
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      FORCES(ng) % tl_vstr(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVsms)),          &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVsms))//                 &
     &                          ', adjusted tl_vstr', Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in tangent linear 3D U-momentum component.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % tl_u(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
        IF (get_var(idSbry(isUvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, u3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_u_obc(:,:,:,:,       &
     &                                                   Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
!
!  Read in tangent linear 3D V-momentum component.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % tl_v(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
        IF (get_var(idSbry(isVvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, v3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_v_obc(:,:,:,:,       &
     &                                                   Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
!
!  Read in tangent linear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTvar(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % tl_t(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D tracers open boundaries.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.get_adjust.and.         &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)), &
     &                          varid)
            status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,          &
     &                             Vname(1,ifield), varid,              &
     &                             InpRec, r3dvar,                      &
     &                             LBij, UBij, 1, N(ng), Nbrec(ng),     &
     &                             Fscl, Fmin, Fmax,                    &
     &                             BOUNDARY(ng) % tl_t_obc(:,:,:,:,     &
     &                                                     Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Read in tangent linear surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.get_adjust.and.                 &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            scale=1.0_r8
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTsur(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nfrec(ng),         &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        FORCES(ng)% tl_tflux(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//         &
     &                            ', adjusted tl_tflux', Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idfrac(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % tl_bed_frac(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idfrac(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idBmas(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % tl_bed_mass(:,:,:,Tindex,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBmas(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSbed(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % tl_bed(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idSbed(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
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
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbld(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask,                         &
#    endif
     &                        SEDBED(ng) % tl_bedldu(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbld(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#    endif
     &                        SEDBED(ng) % tl_bedldv(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
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
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idBott(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % tl_bottom(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBott(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBott(i))), Fmin, Fmax
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
!  Read in adjoint free-surface.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idFsur), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % ad_zeta(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_BOUNDARY
!
!  Read in adjoint free-surface open boundaries.
!
        IF (get_var(idSbry(isFsur)).and.                                &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, r2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % ad_zeta_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in adjoint 2D U-momentum component.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
# endif
     &                      OCEAN(ng) % ad_ubar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D adjoint U-momentum component open boundaries.
!
        IF (get_var(idSbry(isUbar)).and.                                &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, u2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % ad_ubar_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
!
!  Read in adjoint 2D V-momentum component.
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask,                           &
# endif
     &                      OCEAN(ng) % ad_vbar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D V-momentum component open boundaries.
!
        IF (get_var(idSbry(isVbar)).and.                                &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, v2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % ad_vbar_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef ADJUST_WSTRESS
!
!  Read in adjoint linear surface U-momentum stress.
!
        IF (get_var(idUsms)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          scale=1.0_r8
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      FORCES(ng) % ad_ustr(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUsms))//                 &
     &                          ', adjusted ad_ustr', Fmin, Fmax
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
          scale=1.0_r8
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      FORCES(ng) % ad_vstr(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVsms))//                 &
     &                          ', adjusted ad_vstr', Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in adjoint 3D U-momentum component.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % ad_u(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
#  ifdef ADJUST_BOUNDARY
!
!  Read in adjoint 3D U-momentum component open boundaries.
!
        IF (get_var(idSbry(isUvel)).and.                                &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, u3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % ad_u_obc(:,:,:,:,       &
     &                                                   Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
!
!  Read in adjoint 3D V-momentum component.
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % ad_v(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
        IF (get_var(idSbry(isVvel)).and.                                &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, v3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % ad_v_obc(:,:,:,:,       &
     &                                                   Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
#  endif
!
!  Read in adjoint tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTvar(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % ad_t(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  ifdef ADJUST_BOUNDARY
!
!  Read in adjoint 3D tracers open boundaries.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.                        &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)), &
     &                          varid)
            status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,          &
     &                             Vname(1,ifield), varid,              &
     &                             InpRec, r3dvar,                      &
     &                             LBij, UBij, 1, N(ng), Nbrec(ng),     &
     &                             Fscl, Fmin, Fmax,                    &
     &                             BOUNDARY(ng) % ad_t_obc(:,:,:,:,     &
     &                                                     Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Read in adjoint surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            scale=1.0_r8
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTsur(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nfrec(ng),         &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        FORCES(ng) % ad_tflux(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//         &
     &                            ', adjusted ad_tflux', Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idfrac(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % ad_bed_frac(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idfrac(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idBmas(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % ad_bed_mass(:,:,:,Tindex,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBmas(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idSbed(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % ad_bed(:,:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idSbed(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
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
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbld(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % umask,                         &
#    endif
     &                        SEDBED(ng) % ad_bedldu(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idUbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
              END IF
            END IF
          END IF
!
          IF (get_var(idVbld(i))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idVbld(i))), varid)
            gtype=var_flag(varid)*v2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbld(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % vmask,                         &
#    endif
     &                        SEDBED(ng) % ad_bedldv(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idVbld(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
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
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idBott(i)), varid,                &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
     &                        SEDBED(ng) % ad_bottom(:,:,i))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idBott(i))),     &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idBott(i))), Fmin, Fmax
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
!  Read in error covariance normalization (nondimensional) factor.
!-----------------------------------------------------------------------
!
      NRM_STATE: IF ((model.eq. 5).or.                                  &
     &               (model.eq.10).or.                                  &
     &               (model.eq.11)) THEN
!
!  Read in free-surface normalization factor.
!
        IF (get_var(idFsur).and.(model.eq.5)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idFsur), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % b_zeta(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % b_zeta(:,:,Tindex))
# endif
        END IF
!
!  Read in 2D U-momentum component normalization factor.
!
        IF (get_var(idUbar).and.(model.eq.5)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
# endif
     &                      OCEAN(ng) % b_ubar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % b_ubar(:,:,Tindex))
# endif
        END IF
!
!  Read in 2D V-momentum component normalization factor.
!
        IF (get_var(idVbar).and.(model.eq.5)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                       LBi, UBi, LBj, UBj,                        &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask,                           &
# endif
     &                      OCEAN(ng) % b_vbar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % b_vbar(:,:,Tindex))
# endif
        END IF
# ifdef SOLVE3D
!
!  Read in 3D U-momentum component normalization factor.
!
        IF (get_var(idUvel).and.(model.eq.5)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % b_u(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % b_u(:,:,:,Tindex))
#  endif
        END IF
!
!  Read in 3D V-momentum component normalization factor.
!
        IF (get_var(idVvel).and.(model.eq.5)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % b_v(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % b_v(:,:,:,Tindex))
#  endif
        END IF
!
!  Read in tracer type variables normalization factor.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc)).and.(model.eq.5)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTvar(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % b_t(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, MyRank, IDmod, 1,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          OCEAN(ng) % b_t(:,:,:,Tindex,itrc))
#  endif
          END IF
        END DO
# endif
# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries normalization factor.
!
        IF (get_var(idSbry(isFsur)).and.(model.eq.10).and.              &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          BOUNDARY(ng) % b_zeta_obc(LBij:,:),     &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1/),                        &
     &                          total = (/IorJ,4/),                     &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D U-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isUbar)).and.(model.eq.10).and.              &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          BOUNDARY(ng) % b_ubar_obc(LBij:,:),     &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1/),                        &
     &                          total = (/IorJ,4/),                     &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D V-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isVbar)).and.(model.eq.10).and.              &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          BOUNDARY(ng) % b_vbar_obc(LBij:,:),     &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1/),                        &
     &                          total = (/IorJ,4/),                     &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
#  ifdef SOLVE3D
!
!  Read in 3D U-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isUvel)).and.(model.eq.10).and.              &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          BOUNDARY(ng) % b_u_obc(LBij:,:,:),      &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1,1/),                      &
     &                          total = (/IorJ,N(ng),4/),               &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 3D V-momentum component open boundaries normalization factor.
!
        IF (get_var(idSbry(isVvel)).and.(model.eq.10).and.              &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          BOUNDARY(ng) % b_v_obc(LBij:,:,:),      &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1,1/),                      &
     &                          total = (/IorJ,N(ng),4/),               &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 3D tracers open boundaries normalization factor.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.(model.eq.10).and.      &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            CALL netcdf_get_fvar (ng, IDmod, ncname,                    &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            BOUNDARY(ng) % b_t_obc(LBij:,:,:,     &
     &                                                   itrc),         &
     &                            ncid = ncINPid,                       &
     &                            start =(/1,1,1/),                     &
     &                            total =(/IorJ,N(ng),4/),              &
     &                            min_val = Fmin, max_val = Fmax)
            IF (exit_flag.ne.NoError) RETURN
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          Fmin, Fmax
            END IF
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Read in surface U-momentum stress normalization factors.
!
        IF (get_var(idUsms).and.(model.eq.11)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      FORCES(ng) % b_sustr)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % b_sustr)
#  endif
        END IF
!
!  Read in surface V-momentum stress normalization factors.
!
        IF (get_var(idVsms).and.(model.eq.11)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      FORCES(ng) % b_svstr)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % b_svstr)
#  endif
        END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Read in surface tracer flux normalization factors.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.(model.eq.11).and.              &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTsur(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        FORCES(ng) % b_stflx(:,:,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % b_stflx(:,:,itrc))
#  endif
          END IF
        END DO
# endif
      END IF NRM_STATE
#endif
#ifdef FILTERED_RST
      IF (nrrec(ng).ne.0) THEN
        DO i=1,FILN
          CALL get_filter (ng, i)
        END DO
      END IF
#endif

#if defined FOUR_DVAR || (defined HESSIAN_SV && defined BNORM)
!
!-----------------------------------------------------------------------
!  Read in error covariance standard deviation factors.
!-----------------------------------------------------------------------
!
      STD_STATE: IF ((model.eq.6).or.                                   &
     &               (model.eq.8).or.                                   &
     &               (model.eq.9)) THEN
!
!  Read in free-surface standard deviation.
!
        IF (get_var(idFsur).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idFsur), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % e_zeta(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % e_zeta(:,:,Tindex))
# endif
        END IF
!
!  Read in 2D U-momentum component standard deviation.
!
        IF (get_var(idUbar).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
# endif
     &                      OCEAN(ng) % e_ubar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % e_ubar(:,:,Tindex))
# endif
        END IF
!
!  Read in 2D V-momentum component standard deviation.
!
        IF (get_var(idVbar).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask,                           &
# endif
     &                      OCEAN(ng) % e_vbar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % e_vbar(:,:,Tindex))
# endif
        END IF
# ifdef SOLVE3D
!
!  Read in 3D U-momentum component standard deviation.
!
        IF (get_var(idUvel).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % e_u(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % e_u(:,:,:,Tindex))
#  endif
        END IF
!
!  Read in 3D V-momentum standard deviation.
!
        IF (get_var(idVvel).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % e_v(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % e_v(:,:,:,Tindex))
#  endif
        END IF
!
!  Read in tracer type variables standard deviation.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc)).and.(model.eq.6)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTvar(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % e_t(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        OCEAN(ng) % e_t(:,:,:,Tindex,itrc))
#  endif
          END IF
        END DO
# endif
!
!  Read in convolution horizontal diffusion coefficients.
!
        IF (have_var(idKhor).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idKhor)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idKhor), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, KhMin(ng), KhMax(ng),                 &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      MIXING(ng) % Kh)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idKhor)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idKhor)), KhMin(ng),       &
     &                          KhMax(ng)
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % Kh)
# endif
        END IF
# ifdef SOLVE3D
!
!  Read in convolution vertical diffusion coefficient.
!
        IF (have_var(idKver).and.(model.eq.6)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idKver)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idKver), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, KvMin(ng), KvMax(ng),                 &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      MIXING(ng) % Kv)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idKver)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idKver)), KvMin(ng),       &
     &                          KvMax(ng)
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % Kv)
#  endif
        END IF
# endif
# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries standard deviation.
!
        IF (get_var(idSbry(isFsur)).and.(model.eq.8).and.               &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          BOUNDARY(ng) % e_zeta_obc(LBij:,:),     &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1/),                        &
     &                          total = (/IorJ,4/),                     &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D U-momentum component open boundaries standard deviation.
!
        IF (get_var(idSbry(isUbar)).and.(model.eq.8).and.               &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          BOUNDARY(ng) % e_ubar_obc(LBij:,:),     &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1/),                        &
     &                          total = (/IorJ,4/),                     &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D V-momentum component open boundaries standard deviation.
!
        IF (get_var(idSbry(isVbar)).and.(model.eq.8).and.               &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          BOUNDARY(ng) % e_vbar_obc(LBij:,:),     &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1/),                        &
     &                          total = (/IorJ,4/),                     &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
#  ifdef SOLVE3D
!
!  Read in 3D U-momentum component open boundaries standard deviation.
!
        IF (get_var(idSbry(isUvel)).and.(model.eq.8).and.               &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          BOUNDARY(ng) % e_u_obc(LBij:,:,:),      &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1,1/),                      &
     &                          total = (/IorJ,N(ng),4/),               &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 3D V-momentum component open boundaries standard deviation.
!
        IF (get_var(idSbry(isVvel)).and.(model.eq.8).and.               &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          CALL netcdf_get_fvar (ng, IDmod, ncname,                      &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          BOUNDARY(ng) % e_v_obc(LBij:,:,:),      &
     &                          ncid = ncINPid,                         &
     &                          start = (/1,1,1/),                      &
     &                          total = (/IorJ,N(ng),4/),               &
     &                          min_val = Fmin, max_val = Fmax)
          IF (exit_flag.ne.NoError) RETURN
          IF (Master) THEN
            WRITE (stdout,70) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
!
!  Read in 3D tracers open boundaries standard deviation.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.(model.eq.8).and.       &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            CALL netcdf_get_fvar (ng, IDmod, ncname,                    &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            BOUNDARY(ng) % e_t_obc(LBij:,:,:,     &
     &                                                   itrc),         &
     &                            ncid = ncINPid,                       &
     &                            start =(/1,1,1/),                     &
     &                            total =(/IorJ,N(ng),4/),              &
     &                            min_val = Fmin, max_val = Fmax)
            IF (exit_flag.ne.NoError) RETURN
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          Fmin, Fmax
            END IF
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Read in surface U-momentum stress standard deviation.
!
        IF (get_var(idUsms).and.(model.eq.9)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      FORCES(ng) % e_sustr)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUsms)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % e_sustr)
#  endif
        END IF
!
!  Read in surface V-momentum stress standard deviation.
!
        IF (get_var(idVsms).and.(model.eq.9)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v2dvar
          scale=1.0_r8/rho0                 ! N/m2 (Pa) to m2/s2
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      FORCES(ng) % e_svstr)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVsms)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVsms)), Fmin, Fmax
            END IF
          END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % e_svstr)
#  endif
        END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Read in surface tracer flux standard deviations.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.(model.eq.9).and.               &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r2dvar
            IF (itrc.eq.itemp) THEN
              scale=1.0_r8/(rho0*Cp)        ! W/m2 to Celsius m/s
            ELSE
              scale=1.0_r8
            END IF
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
                              Vname(1,idTsur(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        FORCES(ng) % e_stflx(:,:,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        FORCES(ng) % e_stflx(:,:,itrc))
#  endif
          END IF
        END DO
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
        CALL netcdf_get_fvar (ng, IDmod, ncname, TRIM(tvarnam),         &
     &                        FrcTime(ng:),                             &
     &                        ncid = ncINPid,                           &
     &                        start = (/InpRec/),                       &
     &                        total = (/1/))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in free-surface forcing.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idFsur), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % f_zeta)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifndef SOLVE3D
!
!  Read in 2D momentum forcing in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % f_ubar)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
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
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % f_vbar)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in 3D momentum forcing in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % f_u)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
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
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % f_v)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTvar(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % f_t(:,:,:,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
# endif
      END IF FRC_STATE
#endif

#if (defined W4DPSAS || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY)          && \
    (defined ADJUST_BOUNDARY || defined ADJUST_WSTRESS ||\
     defined ADJUST_STFLUX)
!
!-----------------------------------------------------------------------
!  Read in tangent linear forcing corrections.
!-----------------------------------------------------------------------
!
      TLM_FORCING: IF (model.eq.12) THEN
!
!  Set switch to process surface forcing and/or open boundaries during
!  4DVar minimization.
!
        get_adjust=.TRUE.

# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
        IF (get_var(idSbry(isFsur)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, r2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_zeta_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D U-momentum component open boundaries.
!
        IF (get_var(idSbry(isUbar)).and.get_adjust.and.                 &
     &        ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, u2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_ubar_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 2D V-momentum component open boundaries.
!
        IF (get_var(idSbry(isVbar)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread2d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, v2dvar,                        &
     &                           LBij, UBij, Nbrec(ng),                 &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_vbar_obc(:,:,:,      &
     &                                                      Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF

#  ifdef SOLVE3D
!
!  Read in 3D U-momentum component open boundaries.
!
        IF (get_var(idSbry(isUvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, u3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_u_obc(:,:,:,:,       &
     &                                                   Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 3D V-momentum component open boundaries.
!
        IF (get_var(idSbry(isVvel)).and.get_adjust.and.                 &
     &      ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        varid)
          status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,            &
     &                           Vname(1,ifield), varid,                &
     &                           InpRec, v3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           Fscl, Fmin, Fmax,                      &
     &                           BOUNDARY(ng) % tl_v_obc(:,:,:,:,       &
     &                                                   Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,ifield)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in 3D tracers open boundaries.
!
        DO itrc=1,NT(ng)
          IF (get_var(idSbry(isTvar(itrc))).and.get_adjust.and.         &
     &        ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)), &
     &                          varid)
            status=nf_fread3d_bry (ng, IDmod, ncname, ncINPid,          &
     &                             Vname(1,ifield), varid,              &
     &                             InpRec, r3dvar,                      &
     &                             LBij, UBij, 1, N(ng), Nbrec(ng),     &
     &                             Fscl, Fmin, Fmax,                    &
     &                             BOUNDARY(ng) % tl_t_obc(:,:,:,:,     &
     &                                                     Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Read in tangent linear surface U-momentum stress.
!
        IF (get_var(idUsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUsms)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          scale=1.0_r8
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      FORCES(ng) % tl_ustr(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUsms)),          &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUsms))//                 &
     &                          ', adjusted tl_ustr', Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in tangent linear surface V-momentum stress.
!
        IF (get_var(idVsms).and.get_adjust) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVsms)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          scale=1.0_r8
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVsms), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      FORCES(ng) % tl_vstr(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVsms)),          &
     &                          InpRec, TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVsms))//                 &
     &                          ', adjusted tl_vstr', Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Read in tangent linear surface tracers flux.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTsur(itrc)).and.get_adjust.and.                 &
     &        Lstflux(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTsur(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            scale=1.0_r8
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTsur(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, Nfrec(ng),         &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        FORCES(ng)% tl_tflux(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTsur(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTsur(itrc)))//         &
     &                            ', adjusted tl_tflux', Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
# endif
      END IF TLM_FORCING
#endif
!
#if defined TIME_CONV
!
!-----------------------------------------------------------------------
!  Read in tangent linear model error forcing terms.
!-----------------------------------------------------------------------
!
      TCS_STATE: IF (model.eq.14) THEN
!
!  Set number of records available.
!
        NrecFrc(ng)=Nrec
!
!  Read in next impulse forcing time to process.
!
        CALL netcdf_get_fvar (ng, IDmod, ncname, TRIM(tvarnam),         &
     &                        ForceTime(ng:),                             &
     &                        ncid = ncINPid,                           &
     &                        start = (/InpRec/),                       &
     &                        total = (/1/))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in free-surface forcing.
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          gtype=var_flag(varid)*r2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idFsur), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % tl_zeta(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifndef SOLVE3D
!
!  Read in 2D momentum forcing in the XI-direction.
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          gtype=var_flag(varid)*u2dvar
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % tl_ubar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUbar)), Fmin, Fmax
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
          status=nf_fread2d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVbar), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % tl_vbar(:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
# endif
# ifdef SOLVE3D
!
!  Read in 3D momentum forcing in the XI-direction.
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idUvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      OCEAN(ng) % tl_u(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idUvel)), Fmin, Fmax
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
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvel), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      OCEAN(ng) % tl_v(:,:,:,Tindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,60) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idVvel)), Fmin, Fmax
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
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idTvar(itrc)), varid,             &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % tl_t(:,:,:,Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idTvar(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
# endif
      END IF TCS_STATE
#endif
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      CALL netcdf_close (ng, IDmod, ncINPid, ncname, .FALSE.)
!
  10  FORMAT (/,a,'GET_STATE - unable to open input NetCDF file: ',a)
  20  FORMAT (/,a,'GET_STATE - Warning - NetCDF global attribute: ',a,  &
     &        /,18x,'for lateral boundary conditions not checked',      &
     &        /,18x,'in restart file: ',a)
  30  FORMAT (/,a,'GET_STATE - requested input time record = ',i3,/,    &
     &        18x,'not found in input NetCDF: ',a,/,                    &
     &        18x,'number of available records = ',i3)
  40  FORMAT (/,a,'GET_STATE - ',a,t62,'t = ',a,                        &
     &        /,19x,'(Grid ',i2.2,', Iter=',i4.4,', File: ',a,          &
     &        ', Rec=',i4.4,', Index=',i1,')')
  50  FORMAT (/,a,'GET_STATE - ',a,t62,'t = ',a,                        &
     &        /,19x,'(Grid ',i2.2,', File: ',a,', Rec=',i4.4,           &
     &        ', Index=',i1,')')
  60  FORMAT (/,a,'GET_STATE - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,18x,'in input NetCDF file: ',a)
  70  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')

      RETURN
      END SUBROUTINE get_state
