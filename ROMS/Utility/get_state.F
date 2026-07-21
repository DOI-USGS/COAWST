#include "cppdefs.h"
      MODULE get_state_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads in requested model state from specified NetCDF   !
!  file. It is usually used to read initial conditions.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     msg        Message index for StateMsg (string)                   !
!     S          File structure, TYPE(T_IO).                           !
!     IniRec     Time record to read (integer)                         !
!     Tindex     State variable time index to load (integer)           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#if defined ADJUST_BOUNDARY
      USE mod_boundary
#endif
#if defined BIOLOGY
      USE mod_biology
#endif
      USE mod_grid
#if defined ICE_MODEL && defined SOLVE3D
      USE mod_ice
#endif
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
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
      USE mod_stepping
      USE mod_strings
#if defined VEGETATION 
      USE mod_vegetation
#endif 
!
      USE dateclock_mod,      ONLY : time_string
      USE checkvars_mod,      ONLY : checkvars
      USE lbc_mod,            ONLY : lbc_getatt
#ifdef DISTRIBUTE
      USE mp_exchange_mod,    ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod,    ONLY : mp_exchange3d
# endif
#endif
#ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
# ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
# endif
#endif
      USE nf_fread2d_mod,     ONLY : nf_fread2d
      USE nf_fread3d_mod,     ONLY : nf_fread3d
#ifdef SOLVE3D
      USE nf_fread4d_mod,     ONLY : nf_fread4d
#endif
      USE strings_mod,        ONLY : find_string
      USE strings_mod,        ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: get_state
      PRIVATE :: get_state_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: get_state_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_state (ng, model, msg, S, IniRec, Tindex)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex

      integer, intent(inout) :: IniRec
!
      TYPE(T_IO), intent(inout) :: S
!
!  Local variable declarations.
!
      integer :: tile
#ifdef ADJUST_BOUNDARY
      integer :: IorJ, LBij, UBij
#endif
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
      tile=MyRank
#else
      tile=-1
#endif
!
#ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
      IorJ=IOBOUNDS(ng)%IorJ
#endif
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (S%IOtype)
        CASE (io_nf90)
          CALL get_state_nf90 (ng, model, msg, S, IniRec, Tindex,       &
#ifdef ADJUST_BOUNDARY
     &                         IorJ, LBij, UBij,                        &
#endif
     &                         LBi, UBi, LBj, UBj)

#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_state_pio (ng, model, msg, S, IniRec, Tindex,        &
# ifdef ADJUST_BOUNDARY
     &                        IorJ, LBij, UBij,                         &
# endif
     &                        LBi, UBi, LBj, UBj)
#endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) S%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_STATE - Illegal input file type, io_type = ',i0,    &
     &        /,13x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_state
!
!***********************************************************************
      SUBROUTINE get_state_nf90 (ng, model, msg, S, IniRec, Tindex,     &
#ifdef ADJUST_BOUNDARY
     &                           IorJ, LBij, UBij,                      &
#endif
     &                           LBi, UBi, LBj, UBj)

!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex
#ifdef ADJUST_BOUNDARY
      integer, intent(in) :: IorJ, LBij, UBij
#endif
      integer, intent(in) :: LBi, UBi, LBj, UBj

      integer, intent(inout) :: IniRec
!
      TYPE(T_IO), intent(inout) :: S
!
!  Local variable declarations.
!
      logical :: Perfect2D, Perfect3D, foundit
#if defined ADJUST_BOUNDARY || \
    defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
      logical :: get_adjust
#endif
      logical, dimension(NV) :: get_var, have_var
!
      integer :: IDmod, InpRec, gtype, i, ifield, itrc, lstr, lend
      integer :: Nrec, mySize, ncINPid, nvatts, nvdim, status, varid
      integer :: Vsize(4), start(4), total(4)
      integer(i8b) :: Fhash
!
      real(dp), parameter :: Fscl = 1.0_r8

      real(dp) :: INPtime, Tmax, my_dstart, scale, time_scale
      real(r8) :: Fmax, Fmin

      real(dp), allocatable :: TimeVar(:)
!
      character (len=  5) :: string
      character (len= 15) :: Tstring, attnam, tvarnam
      character (len= 22) :: t_code
      character (len= 40) :: tunits
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_state_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
      ncname=TRIM(S%name)
!
!  Set model identification string.
!
      IF (model.eq.iNLM.or.(model.eq.0)) THEN
        string='NLM: '                     ! nonlinear model, restart
        IDmod=iNLM
      ELSE IF (model.eq.iTLM) THEN
        string='TLM: '                     ! tangent linear model
        IDmod=iTLM
      ELSE IF (model.eq.iRPM) THEN
        string='RPM: '                     ! representer model
        IDmod=iRPM
      ELSE IF (model.eq.iADM) THEN
        string='ADM: '                     ! adjoint model
        IDmod=iADM
      ELSE IF (model.eq.5) THEN
        string='NLM: '                     ! surface forcing and
        IDmod=iNLM                         ! OBC increments
      ELSE IF (model.eq.6) THEN
        string='TLM: '                     ! tangent linear error
        IDmod=iTLM                         ! forcing (time covariance)
      ELSE IF (model.eq.7) THEN
        string='FRC: '                     ! impulse forcing
        IDmod=iNLM
      ELSE IF (model.eq.8) THEN
        string='TLM: '                     ! v-space increments
        IDmod=iTLM                         ! I4D-Var
      ELSE IF (model.eq.9) THEN
        string='NLM: '                     ! nonlinear model
        IDmod=iNLM                         ! background state
      ELSE IF (model.eq.10) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! initial conditions
      ELSE IF (model.eq.11) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! model error
      ELSE IF (model.eq.12) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.13) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.14) THEN
        string='NRM: '                     ! normalization factors
        IDmod=iNLM                         ! initial conditions
      ELSE IF (model.eq.15) THEN
        string='NRM: '                     ! normalization factors
        IDmod=iNLM                         ! model error
      ELSE IF (model.eq.16) THEN
        string='NRM: '                     ! normalization factor
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.17) THEN
        string='NRM: '                     ! normalization factor
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.19) THEN
        string='ITL: '                     ! cummulative increments for
        IDmod=iTLM                         ! Incremental Analysis Update
      END IF

#ifdef PROFILE
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, IDmod, 80, __LINE__, MyFile)
#endif
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
!  Open input NetCDF file.
!
      CALL netcdf_open (ng, IDmod, ncname, 0, ncINPid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(ncname)
        RETURN
      END IF
!
!  Determine variables to read.
!
      CALL checkvars (ng, model, ncname, ncINPid, string,               &
     &                Nrec, NV, tvarnam, get_var, have_var)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      SourceFile=MyFile

#if defined DEBUGGING || defined NO_LBC_ATT
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
!
!  Inquire about the input time variable.
!
      CALL netcdf_inq_var (ng, IDmod, ncname,                           &
     &                     ncid = ncINPid,                              &
     &                     MyVarName = TRIM(tvarnam),                   &
     &                     VarID = varid,                               &
     &                     nVarDim =  nvdim,                            &
     &                     nVarAtt = nvatts)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Allocate input time variable and read its value(s).  Recall that
!  input time variable is a one-dimensional array with one or several
!  values.
!
      mySize=var_Dsize(1)
      IF (.not.allocated(TimeVar)) allocate (TimeVar(mySize))
      CALL netcdf_get_time (ng, IDmod, ncname, TRIM(tvarnam),           &
     &                      Rclock%DateNumber, TimeVar,                 &
     &                      ncid = ncINPid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      time_scale=0.0_dp
      DO i=1,nvatts
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            time_scale=day2sec
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            time_scale=1.0_dp
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
#ifdef GENERIC_DSTART
        IF (INItime(ng).lt.0.0_dp) THEN
          my_dstart=dstart                  ! ROMS input script
        ELSE
          my_dstart=INItime(ng)/86400.0_dp  ! NLM IC time is known
        END IF
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(my_dstart*day2sec))) THEN
          INPtime=my_dstart*day2sec
        END IF
#else
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(dstart*day2sec))) THEN
          INPtime=dstart*day2sec
        END IF
#endif
        IF (msg.eq.1) THEN            ! processing initial conditions
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
          ntstart(ng)=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
          IF (ntstart(ng).lt.1) ntstart(ng)=1
          ntend(ng)=ntstart(ng)+ntimes(ng)-1
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
        IF ((10.le.model).and.(model.le.17)) THEN
          t_code=' '              ! time is meaningless for these fields
        ELSE
          CALL time_string (INPtime, t_code)
        END IF
        WRITE (Tstring,'(f15.4)') tdays(ng)
#if defined I4DVAR          || defined I4DVAR_ANA_SENSITIVITY || \
    defined WEAK_CONSTRAINT
        WRITE (stdout,40) string, TRIM(StateMsg(msg)),                  &
     &                    t_code, ng, ', Outer=', outer,                &
     &                    TRIM(ADJUSTL(Tstring)), ncname(lstr:lend),    &
     &                    InpRec, Tindex
#else
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,40) string, TRIM(StateMsg(msg)),                &
     &                      t_code, ng, ', Iter=', Nrun,                &
     &                      TRIM(ADJUSTL(Tstring)),  ncname(lstr:lend), &
     &                      InpRec, Tindex
        ELSE
          WRITE (stdout,50) string, TRIM(StateMsg(msg)),                &
     &                      t_code, ng, TRIM(ADJUSTL(Tstring)),         &
     &                      ncname(lstr:lend), InpRec, Tindex
        END IF
#endif
      END IF

#ifdef NONLINEAR
!
!-----------------------------------------------------------------------
!  Read in nonlinear state variables. If applicable, read in perfect
!  restart variables.
!-----------------------------------------------------------------------
!
      NLM_STATE: IF ((model.eq.iNLM).or.(model.eq.0)) THEN
# include "get_state_nlm_nf90.h"
      END IF NLM_STATE
#endif

#if defined TANGENT || defined TL_IOMS
!
!-----------------------------------------------------------------------
!  Read in tangent linear state variables.
!-----------------------------------------------------------------------
!
      TLM_STATE: IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
# include "get_state_tlm_nf90.h"
      END IF TLM_STATE
#endif

#ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Read in adjoint state variables.
!-----------------------------------------------------------------------
!
      ADM_STATE: IF (model.eq.iADM) THEN
# include "get_state_adm_nf90.h"
      END IF ADM_STATE
#endif

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in error covariance normalization (nondimensional) factors.
!-----------------------------------------------------------------------
!
      NRM_STATE: IF ((model.eq.14).or.                                  &
     &               (model.eq.15).or.                                  &
     &               (model.eq.16).or.                                  &
     &               (model.eq.17)) THEN
# include "get_state_nrm_nf90.h"
      END IF NRM_STATE
#endif

#if defined FOUR_DVAR || (defined HESSIAN_SV && defined BNORM)
!
!-----------------------------------------------------------------------
!  Read in error covariance standard deviation factors.
!-----------------------------------------------------------------------
!
      STD_STATE: IF ((model.eq.10).or.                                  &
     &               (model.eq.11).or.                                  &
     &               (model.eq.12).or.                                  &
     &               (model.eq.13)) THEN
# include "get_state_std_nf90.h"
      END IF STD_STATE
#endif

#if defined FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in generic state variables into work "d_" arrays.
!-----------------------------------------------------------------------
!
      GENERIC_STATE: IF (model.eq.19) THEN
# include "get_state_generic_nf90.h"
      END IF GENERIC_STATE
#endif

#if defined IMPULSE
!
!-----------------------------------------------------------------------
!  Read in adjoint model or tangent linear model impulse forcing terms.
!-----------------------------------------------------------------------
!
      FRC_STATE: IF (model.eq.7) THEN
# include "get_state_frc_nf90.h"
      END IF FRC_STATE
#endif

#if (defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined TL_RBL4DVAR)             && \
    (defined ADJUST_BOUNDARY          || \
     defined ADJUST_STFLUX            || \
     defined ADJUST_WSTRESS)
!
!-----------------------------------------------------------------------
!  Read in tangent linear forcing corrections.
!-----------------------------------------------------------------------
!
      TLM_FORCING: IF (model.eq.5) THEN
# include "get_state_tlm_forcing_nf90.h"
      END IF TLM_FORCING
#endif
!
#if defined TIME_CONV
!
!-----------------------------------------------------------------------
!  Read in tangent linear model error forcing terms used in the time
!  convolutions.
!-----------------------------------------------------------------------
!
      TCS_STATE: IF (model.eq.6) THEN
# include "get_state_tcs_nf90.h"
      END IF TCS_STATE
#endif
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      CALL netcdf_close (ng, IDmod, ncINPid, ncname, .FALSE.)

#ifdef PROFILE
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, IDmod, 80, __LINE__, MyFile)
#endif
!
  10  FORMAT (/,2x,'GET_STATE_NF90 - ',a,'unable to open input NetCDF', &
     &        ' file: ',a)
  20  FORMAT (/,2x,'GET_STATE_NF90 - ',a,'Warning - NetCDF global',     &
     &        ' attribute:',a,                                          &
     &        /,19x,'for lateral boundary conditions not checked',      &
     &        /,19x,'in file: ',a)
  30  FORMAT (/,2x,'GET_STATE_NF90 - ',a,'requested input time',        &
     &        ' record = ',i0,/,19x,'not found in input NetCDF: ',a,/,  &
     &        19x,'number of available records = ',i0)
  40  FORMAT (/,2x,'GET_STATE_NF90   - ',a,a,t75,a,                     &
     &        /,22x,'(Grid ',i2.2,a,i4.4, ', t = ',a,                   &
     &        ', File: ',a, ', Rec=',i4.4,', Index=',i1,')')
  50  FORMAT (/,2x,'GET_STATE_NF90   - ',a,a,t75,a,                     &
     &        /,22x,'(Grid ',i2.2, ', t = ',a,                          &
     &        ', File: ',a,', Rec=',i4.4, ', Index=',i1,')')
  60  FORMAT (/,2x,'GET_STATE_NF90 - ',a,'error while reading',         &
     &        ' variable: ',a,2x,'at time record = ',i0,                &
     &        /,19x,'in input NetCDF file: ',a)
#ifdef CHECKSUM
  70  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,' CheckSum = ',i0,')')
#else
  70  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
#endif
  75  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
  80  FORMAT (/,2x,'GET_STATE_NF90 - ',a,'cannot find variable: ',a,    &
     &        /,19x,'in input NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE get_state_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_state_pio (ng, model, msg, S, IniRec, Tindex,      &
# ifdef ADJUST_BOUNDARY
     &                          IorJ, LBij, UBij,                       &
# endif
     &                          LBi, UBi, LBj, UBj)

!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: IorJ, LBij, UBij
# endif
      integer, intent(in) :: LBi, UBi, LBj, UBj

      integer, intent(inout) :: IniRec
!
      TYPE(T_IO), intent(inout) :: S
!
!  Local variable declarations.
!
      logical :: Perfect2D, Perfect3D, foundit
# if defined ADJUST_BOUNDARY || \
    defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
      logical :: get_adjust
# endif
      logical, dimension(NV) :: get_var, have_var
!
      integer :: IDmod, InpRec, i, ifield, itrc, lstr, lend
      integer :: Nrec, mySize, nvatts, nvdim, status, vindex
      integer :: Vsize(4), start(4), total(4)
      integer(i8b) :: Fhash
!
      real(dp), parameter :: Fscl = 1.0_r8

      real(dp) :: INPtime, Tmax, my_dstart, scale, time_scale
      real(r8) :: Fmax, Fmin

      real(dp), allocatable :: TimeVar(:)
!
      character (len=  5) :: string
      character (len= 15) :: Tstring, attnam, tvarnam
      character (len= 22) :: t_code
      character (len= 40) :: tunits
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_state_pio"
!
      TYPE (IO_Desc_t), pointer :: ioDesc

      TYPE (file_desc_t) :: pioFile
      TYPE (Var_desc_t)  :: pioVar
      TYPE (My_VarDesc)  :: my_pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
      ncname=TRIM(S%name)
!
!  Set model identification string.
!
      IF (model.eq.iNLM.or.(model.eq.0)) THEN
        string='NLM: '                     ! nonlinear model, restart
        IDmod=iNLM
      ELSE IF (model.eq.iTLM) THEN
        string='TLM: '                     ! tangent linear model
        IDmod=iTLM
      ELSE IF (model.eq.iRPM) THEN
        string='RPM: '                     ! representer model
        IDmod=iRPM
      ELSE IF (model.eq.iADM) THEN
        string='ADM: '                     ! adjoint model
        IDmod=iADM
      ELSE IF (model.eq.5) THEN
        string='NLM: '                     ! surface forcing and
        IDmod=iNLM                         ! OBC increments
      ELSE IF (model.eq.6) THEN
        string='TLM: '                     ! tangent linear error
        IDmod=iTLM                         ! forcing (time covariance)
      ELSE IF (model.eq.7) THEN
        string='FRC: '                     ! impulse forcing
        IDmod=iNLM
      ELSE IF (model.eq.8) THEN
        string='TLM: '                     ! v-space increments
        IDmod=iTLM                         ! I4D-Var
      ELSE IF (model.eq.9) THEN
        string='NLM: '                     ! nonlinear model
        IDmod=iNLM                         ! background state
      ELSE IF (model.eq.10) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! initial conditions
      ELSE IF (model.eq.11) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! model error
      ELSE IF (model.eq.12) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.13) THEN
        string='STD: '                     ! standard deviation
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.14) THEN
        string='NRM: '                     ! normalization factors
        IDmod=iNLM                         ! initial conditions
      ELSE IF (model.eq.15) THEN
        string='NRM: '                     ! normalization factors
        IDmod=iNLM                         ! model error
      ELSE IF (model.eq.16) THEN
        string='NRM: '                     ! normalization factor
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.17) THEN
        string='NRM: '                     ! normalization factor
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.19) THEN
        string='ITL: '                     ! cummulative increments for
        IDmod=iTLM                         ! Incremental Analysis Update
      END IF

# ifdef PROFILE
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, IDmod, 80, __LINE__, MyFile)
# endif
!
!  Set switch to process variables for nonlinear model perfect restart.
!
      Perfect2D=.FALSE.
      Perfect3D=.FALSE.
# ifdef PERFECT_RESTART
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        Perfect2D=.TRUE.
        Perfect3D=.TRUE.
      END IF
# endif
      PerfectRST(ng)=Perfect2D.or.Perfect3D
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
!  Open input NetCDF file.
!
      CALL pio_netcdf_open (ng, IDmod, ncname, 0, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(ncname)
        RETURN
      END IF
!
!  Determine variables to read.
!
      CALL checkvars (ng, model, ncname, pioFile, string, Nrec, NV,     &
     &                tvarnam, get_var, have_var)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      SourceFile=MyFile

# if defined DEBUGGING || defined NO_LBC_ATT
!
!  Lateral boundary conditions attribute not checked in restart file.
!
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        IF (Master) WRITE (stdout,20) string, 'NLM_LBC', TRIM(ncname)
      END IF
# else
!
!  If restart, read in lateral boundary conditions global attribute
!  from restart file and check keyword strings with structure vlues
!  for consistency.
!
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        CALL lbc_getatt (ng, model, pioFile, ncname, 'NLM_LBC', LBC)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
!
!  Inquire about the input time variable.
!
      CALL pio_netcdf_inq_var (ng, IDmod, ncname,                       &
     &                         pioFile = pioFile,                       &
     &                         MyVarName = TRIM(tvarnam),               &
     &                         pioVar = pioVar,                         &
     &                         nVarDim =  nvdim,                        &
     &                         nVarAtt = nvatts)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Allocate input time variable and read its value(s).  Recall that
!  input time variable is a one-dimensional array with one or several
!  values.
!
      mySize=var_Dsize(1)
      IF (.not.allocated(TimeVar)) allocate (TimeVar(mySize))
      CALL pio_netcdf_get_time (ng, IDmod, ncname, TRIM(tvarnam),       &
     &                          Rclock%DateNumber, TimeVar,             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      time_scale=0.0_dp
      DO i=1,nvatts
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            time_scale=day2sec
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            time_scale=1.0_dp
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
# ifdef GENERIC_DSTART
        IF (INItime(ng).lt.0.0_dp) THEN
          my_dstart=dstart                  ! ROMS input script
        ELSE
          my_dstart=INItime(ng)/86400.0_dp  ! NLM IC time is known
        END IF
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(my_dstart*day2sec))) THEN
          INPtime=my_dstart*day2sec
        END IF
# else
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(dstart*day2sec))) THEN
          INPtime=dstart*day2sec
        END IF
# endif
        IF (msg.eq.1) THEN            ! processing initial conditions
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
          ntstart(ng)=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
          IF (ntstart(ng).lt.1) ntstart(ng)=1
          ntend(ng)=ntstart(ng)+ntimes(ng)-1
          IF (PerfectRST(ng)) THEN
            ntfirst(ng)=1
          ELSE
            ntfirst(ng)=ntstart(ng)
          END IF
        END IF
# ifdef WEAK_CONSTRAINT
        IF (msg.eq.4) THEN
          ForceTime(ng)=time(ng)
        END IF
# endif
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
        IF ((10.le.model).and.(model.le.17)) THEN
          t_code=' '              ! time is meaningless for these fields
        ELSE
          CALL time_string (INPtime, t_code)
        END IF
        WRITE (Tstring,'(f15.4)') tdays(ng)
# if defined I4DVAR          || defined I4DVAR_ANA_SENSITIVITY || \
    defined WEAK_CONSTRAINT
        WRITE (stdout,40) string, TRIM(StateMsg(msg)),                  &
     &                    t_code, ng, ', Outer=', outer,                &
     &                    TRIM(ADJUSTL(Tstring)), ncname(lstr:lend),    &
     &                    InpRec, Tindex
# else
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,40) string, TRIM(StateMsg(msg)),                &
     &                      t_code, ng, ', Iter=', Nrun,                &
     &                      TRIM(ADJUSTL(Tstring)),  ncname(lstr:lend), &
     &                      InpRec, Tindex
        ELSE
          WRITE (stdout,50) string, TRIM(StateMsg(msg)),                &
     &                      t_code, ng, TRIM(ADJUSTL(Tstring)),         &
     &                      ncname(lstr:lend), InpRec, Tindex
        END IF
# endif
      END IF

# ifdef NONLINEAR
!
!-----------------------------------------------------------------------
!  Read in nonlinear state variables. If applicable, read in perfect
!  restart variables.
!-----------------------------------------------------------------------
!
      NLM_STATE: IF ((model.eq.iNLM).or.(model.eq.0)) THEN
#  include "get_state_nlm_pio.h"
      END IF NLM_STATE
# endif

# if defined TANGENT || defined TL_IOMS
!
!-----------------------------------------------------------------------
!  Read in tangent linear state variables.
!-----------------------------------------------------------------------
!
      TLM_STATE: IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
#  include "get_state_tlm_pio.h"
      END IF TLM_STATE
# endif

# ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Read in adjoint state variables.
!-----------------------------------------------------------------------
!
      ADM_STATE: IF (model.eq.iADM) THEN
#  include "get_state_adm_pio.h"
      END IF ADM_STATE
# endif

# ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in error covariance normalization (nondimensional) factors.
!-----------------------------------------------------------------------
!
      NRM_STATE: IF ((model.eq.14).or.                                  &
     &               (model.eq.15).or.                                  &
     &               (model.eq.16).or.                                  &
     &               (model.eq.17)) THEN
#  include "get_state_nrm_pio.h"
      END IF NRM_STATE
# endif

# if defined FOUR_DVAR || (defined HESSIAN_SV && defined BNORM)
!
!-----------------------------------------------------------------------
!  Read in error covariance standard deviation factors.
!-----------------------------------------------------------------------
!
      STD_STATE: IF ((model.eq.10).or.                                  &
     &               (model.eq.11).or.                                  &
     &               (model.eq.12).or.                                  &
     &               (model.eq.13)) THEN
#  include "get_state_std_pio.h"
      END IF STD_STATE
# endif

# if defined FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in generic state variables into work "d_" arrays.
!-----------------------------------------------------------------------
!
      GENERIC_STATE: IF (model.eq.19) THEN
#  include "get_state_generic_pio.h"
      END IF GENERIC_STATE
# endif

# if defined IMPULSE
!
!-----------------------------------------------------------------------
!  Read in adjoint model or tangent linear model impulse forcing terms.
!-----------------------------------------------------------------------
!
      FRC_STATE: IF (model.eq.7) THEN
#  include "get_state_frc_pio.h"
      END IF FRC_STATE
# endif

# if (defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined TL_RBL4DVAR)             && \
    (defined ADJUST_BOUNDARY          || \
     defined ADJUST_STFLUX            || \
     defined ADJUST_WSTRESS)
!
!-----------------------------------------------------------------------
!  Read in tangent linear forcing corrections.
!-----------------------------------------------------------------------
!
      TLM_FORCING: IF (model.eq.5) THEN
#  include "get_state_tlm_forcing_pio.h"
      END IF TLM_FORCING
# endif
!
# if defined TIME_CONV
!
!-----------------------------------------------------------------------
!  Read in tangent linear model error forcing terms used in the time
!  convolutions.
!-----------------------------------------------------------------------
!
      TCS_STATE: IF (model.eq.6) THEN
#  include "get_state_tcs_pio.h"
      END IF TCS_STATE
# endif
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_close (ng, IDmod, pioFile, ncname, .FALSE.)

# ifdef PROFILE
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, IDmod, 80, __LINE__, MyFile)
# endif
!
  10  FORMAT (/,2x,'GET_STATE_PIO  - ',a,'unable to open input NetCDF', &
     &        ' file: ',a)
  20  FORMAT (/,2x,'GET_STATE_PIO  - ',a,'Warning - NetCDF global',     &
     &        ' attribute:',a,                                          &
     &        /,19x,'for lateral boundary conditions not checked',      &
     &        /,19x,'in file: ',a)
  30  FORMAT (/,2x,'GET_STATE_PIO  - ',a,'requested input time',        &
     &        ' record = ',i0,/,19x,'not found in input NetCDF: ',a,/,  &
     &        19x,'number of available records = ',i0)
  40  FORMAT (/,2x,'GET_STATE_PIO    - ',a,a,t75,a,                     &
     &        /,22x,'(Grid ',i2.2,a,i4.4, ', t = ',a,                   &
     &        ', File: ',a, ', Rec=',i4.4,', Index=',i1,')')
  50  FORMAT (/,2x,'GET_STATE_PIO    - ',a,a,t75,a,                     &
     &        /,22x,'(Grid ',i2.2, ', t = ',a,                          &
     &        ', File: ',a,', Rec=',i4.4, ', Index=',i1,')')
  60  FORMAT (/,2x,'GET_STATE_PIO  - ',a,'error while reading',         &
     &        ' variable: ',a,2x,'at time record = ',i0,                &
     &        /,19x,'in input NetCDF file: ',a)
#ifdef CHECKSUM
  70  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,' CheckSum = ',i0,')')
#else
  70  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
#endif
  75  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
  80  FORMAT (/,2x,'GET_STATE_PIO  - ',a,'cannot find variable: ',a,    &
     &        /,19x,'in input NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE get_state_pio
#endif
      END MODULE get_state_mod
