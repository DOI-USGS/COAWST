#include "cppdefs.h"
      SUBROUTINE check_multifile (ng, model)
!
!svn $Id: check_multifile.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  If applicable, this routine checks input NetCDF multi-files and     !
!  sets several parameters in the file information structure so the    !
!  appropriate file is selected during initialization or restart.      !
!                                                                      !
!  Multi-files are allowed for several input fields. That is, the      !
!  time records for a particular input field can be split into         !
!  several NetCDF files.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod, ONLY : time_string
      USE strings_mod,   ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: i, max_files
!
      real(dp) :: Tfinal
!
      character (len= 1), parameter :: blank = ' '
      character (len=40) :: file_type
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!=======================================================================
!  If applicable, initialize parameters for input multi-files.
!=======================================================================
!
!  Initialize I/O information variables to facilitate to reset the
!  values in iterative algorithms that call the NLM, TLM, RPM, and
!  ADM kernels repetitevely.  Notice that Iinfo(1,:,:) is not reset
!  because it is part of the metadata.
!
      DO i=1,NV
        Cinfo(i,ng)=blank
        Linfo(1,i,ng)=.FALSE.
        Linfo(2,i,ng)=.FALSE.
        Linfo(3,i,ng)=.FALSE.
        Linfo(4,i,ng)=.FALSE.
        Linfo(5,i,ng)=.FALSE.
        Linfo(6,i,ng)=.FALSE.
        Iinfo(2,i,ng)=-1
        Iinfo(3,i,ng)=-1
        Iinfo(4,i,ng)=0
        Iinfo(5,i,ng)=0
        Iinfo(6,i,ng)=0
        Iinfo(7,i,ng)=0
        Iinfo(8,i,ng)=2
        Iinfo(9,i,ng)=0
        Iinfo(10,i,ng)=0
        Finfo(1,i,ng)=0.0_r8
        Finfo(2,i,ng)=0.0_r8
        Finfo(3,i,ng)=0.0_r8
        Finfo(5,i,ng)=0.0_r8
        Finfo(6,i,ng)=0.0_r8
        Finfo(7,i,ng)=0.0_r8
        Finfo(10,i,ng)=1.0_r8
        Fpoint(1,i,ng)=0.0_r8
        Fpoint(2,i,ng)=0.0_r8
        Tintrp(1,i,ng)=0.0_r8
        Tintrp(2,i,ng)=0.0_r8
        Vtime(1,i,ng)=0.0_r8
        Vtime(2,i,ng)=0.0_r8
        ncFRCid(i,ng)=-1
      END DO
!
!  Get initialization time string.
!
      CALL time_string (time(ng), I_code)
!
!  Get final time string for simulation.
!
      IF (model.eq.iADM) THEN
        Tfinal=time(ng)-ntimes(ng)*dt(ng)
      ELSE
!       Tfinal=time(ng)+ntimes(ng)*dt(ng)
        Tfinal=dstart*day2sec+ntimes(ng)*dt(ng)
      END IF
      CALL time_string (tfinal, F_code)
!
!-----------------------------------------------------------------------
!  Input lateral boundary conditions data.
!-----------------------------------------------------------------------
!
      IF (ObcData(ng)) THEN
        file_type='Lateral Boundary'
        max_files=nBCfiles(ng)
        CALL multifile_info_s2d (ng, model, file_type, nBCfiles,        &
     &                           max_files, BRY)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Input climatology data.
!-----------------------------------------------------------------------
!
      IF (CLM_FILE(ng)) THEN
        file_type='Climatology'
        max_files=nCLMfiles(ng)
        CALL multifile_info_s2d (ng, model, file_type, nCLMfiles,       &
     &                           max_files, CLM)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Nonlinear model forward trajectory data.
!-----------------------------------------------------------------------
!
      IF (LreadFWD(ng)) THEN
        CALL multifile_info_s1d (ng, model, 'Forward', FWD)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#if defined TANGENT || defined TL_IOMS
!
!-----------------------------------------------------------------------
!  Adjoint model trajectory data.
!-----------------------------------------------------------------------
!
      IF (LreadTLM(ng)) THEN
        CALL multifile_info_s1d (ng, model, 'Tangent Linear', TLM)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#ifdef FORWARD_FLUXES
!
!-----------------------------------------------------------------------
!  Nonlinear model surface fluxes data.
!-----------------------------------------------------------------------
!
      IF (LreadBLK(ng)) THEN
        CALL multifile_info_s1d (ng, model, 'Nonlinear Fluxes', BLK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#ifndef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Input Source/Sink data (river runoff).
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        file_type='Sources/Sinks Data'
        CALL multifile_info_s1d (ng, model, file_type, SSF)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#ifdef FRC_FILE
!
!-----------------------------------------------------------------------
!  Input forcing data.
!-----------------------------------------------------------------------
!
      IF (LreadFRC(ng)) THEN
        file_type='Forcing'
        max_files=MAXVAL(nFfiles)
        CALL multifile_info_s2d (ng, model, file_type, nFfiles,         &
     &                           max_files, FRC)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
!
      RETURN
      END SUBROUTINE check_multifile
!
      SUBROUTINE multifile_info_s1d (ng, model, file_type, S)
!
!=======================================================================
!                                                                      !
!  Updates multi-file information for the given 1D file structure, S.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     file_type  File structure indentifier (string)                   !
!     S          2D Derived type structure, TYPE(T_IO)                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod, ONLY : time_string
      USE strings_mod,   ONLY : FoundError, lowercase
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: file_type
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lcheck, foundit
      logical :: check_file
!
      integer :: Fcount, Nfiles, i, ifile
!
      real(dp) :: Tfinal, Tmax, Tmin, Tscale
!
      character (len= 22) :: Tmin_code, Tmax_code
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multifile_info_s1d"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Update file structure parameters.
!-----------------------------------------------------------------------
!
!  Set final time of simulation.
!
      IF (model.eq.iADM) THEN
        Tfinal=time(ng)-ntimes(ng)*dt(ng)
      ELSE
!       Tfinal=time(ng)+ntimes(ng)*dt(ng)
        Tfinal=dstart*day2sec+ntimes(ng)*dt(ng)
      END IF
!
!  Set minimum and maximum time values.
!
      Nfiles=S(ng)%Nfiles
      DO ifile=1,Nfiles
        ncname=TRIM(S(ng)%files(ifile))
        foundit=check_file(ng, model, Nfiles, ncname,                   &
     &                     Tmin, Tmax, Tscale,  Lcheck)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        S(ng)%time_min(ifile)=Tmin
        S(ng)%time_max(ifile)=Tmax
      END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
      Fcount=0
      IF (Lcheck) THEN
        IF (model.eq.iADM) THEN
          DO ifile=Nfiles,1,-1
            Tmax=Tscale*S(ng)%time_max(ifile)
            IF (time(ng).le.Tmax) THEN
              Fcount=ifile
            END IF
          END DO
        ELSE
          DO ifile=1,Nfiles
            Tmin=Tscale*S(ng)%time_min(ifile)
            IF (time(ng).ge.Tmin) THEN
              Fcount=ifile
            END IF
          END DO
        END IF
      ELSE
        Fcount=1
      END IF
!
!  Initialize other structure parameters or issue an error if data does
!  not include initalization time.
!
      IF (Fcount.gt.0) THEN
        S(ng)%Fcount=Fcount
        ncname=S(ng)%files(Fcount)
        S(ng)%name=TRIM(ncname)
      ELSE
        IF (Master.and.Lcheck) THEN
          WRITE (stdout,10) TRIM(lowercase(file_type)), I_code
          DO ifile=1,Nfiles
            Tmin=Tscale*S(ng)%time_min(ifile)
            Tmax=Tscale*S(ng)%time_max(ifile)
            CALL time_string (Tmin, Tmin_code)
            CALL time_string (Tmax, Tmax_code)
            WRITE (stdout,20) Tmin_code, Tmax_code,                     &
     &                        TRIM(S(ng)%files(ifile))
          END DO
        END IF
        exit_flag=4
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Check if there is data up to the end of the simulation.
!
      IF (Lcheck) THEN
        IF (model.eq.iADM) THEN
          Tmin=Tscale*S(ng)%time_min(1)
          IF (Tfinal.lt.Tmin) THEN
            CALL time_string (Tmin, Tmin_code)
            IF (Master) THEN
              WRITE (stdout,30) TRIM(file_type)//' (adjoint)',          &
     &                          TRIM(S(ng)%files(1)),                   &
     &                          'first ', Tmin_code, F_code
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        ELSE
          Tmax=Tscale*S(ng)%time_max(Nfiles)
          IF (Tfinal.gt.Tmax) THEN
            CALL time_string (Tmax, Tmax_code)
            IF (Master) THEN
              WRITE (stdout,30) TRIM(file_type),                        &
     &                          TRIM(S(ng)%files(Nfiles)),              &
     &                          'last  ', Tmax_code, F_code
            END IF
            exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
          END IF
        END IF
      END IF
!
 10   FORMAT (/,' MULTIFILE_INFO_S1D - Error while processing ', a,     &
     &        ' multi-files: ',/,22x,'data does not include',           &
     &        ' initialization time = ', a,/)
 20   FORMAT (3x,a,2x,a,5x,a)
 30   FORMAT (/,' MULTIFILE_INFO_S1D - Error while checking input ', a, &
     &        ' file:',/,22x,a,/,22x,                                   &
     &        a,'data time record available is for day: ',a,/,22x,      &
     &        'but data is needed to finish run until day: ',a)
!
      RETURN
      END SUBROUTINE multifile_info_s1d
!
      SUBROUTINE multifile_info_s2d (ng, model, file_type, InpFiles,    &
     &                               idim, S)
!
!=======================================================================
!                                                                      !
!  Updates multi-file information for the given 2D file structure, S.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     file_type  File structure indentifier (string)                   !
!     InpFiles   Number of source input files per grid (integer vector)!
!     idim       Size of structure inner dimension (integer):          !
!                  idim=MAXVAL(InpFiles) in calling routine            !
!     S          2D Derived type structure, TYPE(T_IO)                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod, ONLY : time_string
      USE strings_mod,   ONLY : FoundError, lowercase
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, idim
      integer, intent(in) :: InpFiles(Ngrids)
!
      character (len=*), intent(in) :: file_type
!
      TYPE(T_IO), intent(inout) :: S(idim,Ngrids)
!
!  Local variable declarations.
!
      logical :: Lcheck, foundit
      logical :: check_file
!
      integer :: Fcount, Nfiles, i, ifile
!
      real(dp) :: Tfinal, Tmax, Tmin, Tscale
!
      character (len= 22) :: Tmin_code, Tmax_code
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multifile_info_s2d"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Update file structure parameters.
!-----------------------------------------------------------------------
!
!  Set final time of simulation.
!
      IF (model.eq.iADM) THEN
        Tfinal=time(ng)-ntimes(ng)*dt(ng)
      ELSE
!       Tfinal=time(ng)+ntimes(ng)*dt(ng)
        Tfinal=dstart*day2sec+ntimes(ng)*dt(ng)
      END IF
!
!  Set minimum and maximum time values.
!
      DO i=1,InpFiles(ng)
        Nfiles=S(i,ng)%Nfiles
        DO ifile=1,Nfiles
          ncname=TRIM(S(i,ng)%files(ifile))
          foundit=check_file(ng, model, Nfiles, ncname,                 &
     &                       Tmin, Tmax, Tscale, Lcheck)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          S(i,ng)%time_min(ifile)=Tmin
          S(i,ng)%time_max(ifile)=Tmax
        END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
        Fcount=0
        IF (Lcheck) THEN
          IF (model.eq.iADM) THEN
            DO ifile=Nfiles,1,-1
              Tmax=Tscale*S(i,ng)%time_max(ifile)
              IF (time(ng).le.Tmax) THEN
                Fcount=ifile
              END IF
            END DO
          ELSE
            DO ifile=1,Nfiles
              Tmin=Tscale*S(i,ng)%time_min(ifile)
              IF (time(ng).ge.Tmin) THEN
                Fcount=ifile
              END IF
            END DO
          END IF
        ELSE
          Fcount=1
        END IF
!
!  Initialize other structure parameters or issue an error if data does
!  not include initalization time.
!
        IF (Fcount.gt.0) THEN
          S(i,ng)%Fcount=Fcount
          ncname=S(i,ng)%files(Fcount)
          S(i,ng)%name=TRIM(ncname)
        ELSE
          IF (Master.and.Lcheck) THEN
            WRITE (stdout,10) TRIM(lowercase(file_type)), I_code
            DO ifile=1,Nfiles
              Tmin=Tscale*S(i,ng)%time_min(ifile)
              Tmax=Tscale*S(i,ng)%time_max(ifile)
              CALL time_string (Tmin, Tmin_code)
              CALL time_string (Tmax, Tmax_code)
              WRITE (stdout,20) Tmin_code, Tmax_code,                   &
     &                          TRIM(S(i,ng)%files(ifile))
            END DO
          END IF
          exit_flag=4
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Check if there is data up to the end of the simulation.
!
        IF (Lcheck) THEN
          IF (model.eq.iADM) THEN
            Tmin=Tscale*S(i,ng)%time_min(1)
            IF (Tfinal.lt.Tmin) THEN
              CALL time_string (Tmin, Tmin_code)
              IF (Master) THEN
                WRITE (stdout,30) TRIM(file_type)//' (adjoint)',        &
     &                            TRIM(S(i,ng)%files(1)),               &
     &                            'first ', Tmin_code, F_code
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          ELSE
            Tmax=Tscale*S(i,ng)%time_max(Nfiles)
            IF (Tfinal.gt.Tmax) THEN
              CALL time_string (Tmax, Tmax_code)
              IF (Master) THEN
                WRITE (stdout,30) TRIM(file_type),                      &
     &                            TRIM(S(i,ng)%files(Nfiles)),          &
     &                            'last  ', Tmax_code, F_code
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END IF
      END DO
!
 10   FORMAT (/,' MULTIFILE_INFO_S2D - Error while processing ', a,     &
     &        ' multi-files: ',/,22x,'data does not include',           &
     &        ' initialization time = ', a,/)
 20   FORMAT (3x,a,2x,a,5x,a)
 30   FORMAT (/,' MULTIFILE_INFO_S2D - Error while checking input ', a, &
     &        ' file:',/,22x,a,/,22x,                                   &
     &        a,'data time record available is for day: ',a,/,22x,      &
     &        'but data is needed to finish run until day: ',a)
!
      RETURN
      END SUBROUTINE multifile_info_s2d
!
      FUNCTION check_file (ng, model, Nfiles, ncname, Tmin, Tmax,       &
     &                     Tscale, Lcheck) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This logical function scans the variables of the provided input     !
!  NetCDF for the time record variable and gets its range of values.   !
!  It used elsewhere to determine which input NetCDF multi-file is     !
!  needed for initialization or restart.                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     Nfiles       Number of multi-files.                              !
!     ncname       NetCDF file name to process (string).               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Tmin         Available minimum time variable value.              !
!     Tmax         Available maximum time variable value.              !
!     Tscale       Scale to convert time variable units to seconds     !
!     Lcheck       Switch to indicate that the time range needs to be  !
!                    checked by the calling routine.                   !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    time variable is found or not.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError, lowercase
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: Lcheck
!
      integer, intent(in) :: ng, model, Nfiles
!
      character (*), intent(in) :: ncname
!
      real(dp), intent(out) :: Tmin, Tmax, Tscale
!
!  Local variable declarations.
!
      logical :: Lcycle, Lperpetual, Lspectral, foundit
!
      integer :: Nrec, TvarID, i, ncid, nvdim, nvatt
!
      character (len=40) :: Tunits, TvarName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", check_file"
!
      SourceFile=MyFile
!
!------------------------------------------------------------------------
!  Check if requested time is within the NetCDF file dataset.
!------------------------------------------------------------------------
!
!  Initialize.
!
      foundit=.FALSE.
      Lcheck=.TRUE.
      Lcycle=.FALSE.
      Lperpetual=.FALSE.
      Lspectral =.FALSE.
      Tscale=1.0_r8                        ! seconds
      Tmin=0.0_r8
      Tmax=0.0_r8
!
!  Open NetCDF file for reading.
!
      CALL netcdf_open (ng, model, ncname, 0, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF
!
!  Inquire about all the variables
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Search for the time variable: any 1D array variable with the string
!  'time' in the variable name.
!
      DO i=1,n_var
        IF ((INDEX(TRIM(lowercase(var_name(i))),'time').ne.0).and.      &
     &            (var_ndim(i).eq.1)) THEN
          TvarName=TRIM(var_name(i))
          foundit=.TRUE.
          EXIT
        ELSE IF ((INDEX(TRIM(var_name(i)),'tide_period').ne.0).and.     &
     &            (var_ndim(i).eq.1)) THEN
          TvarName=TRIM(var_name(i))
          foundit=.TRUE.
          Lspectral=.TRUE.          ! we do not need to check tidal data
          EXIT
        END IF
      END DO
      IF (.not.foundit) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(ncname)
        END IF
        exit_flag=4
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Inquire about requested variable.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ncid,                                 &
     &                     MyVarName = TRIM(TvarName),                  &
     &                     VarID = TvarID,                              &
     &                     nVarDim = nvdim,                             &
     &                     nVarAtt = nvatt)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set number of records available and check the 'units' attribute.
!  Also, set output logical switch 'Lcheck' for the calling to check
!  the available data time range. For example, we need to check it
!  there is enough data to finish the simulation.  Notice that for
!  data with 'cycle_length', Lcheck = FALSE.  Also,  Lcheck = FALSE
!  for perpetual time axis: the 'calendar' attribute is 'none' or
!  the number of records in the time dimension is one (Nrec=1).
!
      Nrec=var_Dsize(1)              ! time is a 1D array
      DO i=1,nvatt
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          Tunits=TRIM(var_Achar(i))
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            Tscale=86400.0_r8
          ELSE IF (INDEX(TRIM(var_Achar(i)),'hour').ne.0) THEN
            Tscale=3600.0_r8
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            Tscale=1.0_r8
          END IF
        ELSE IF (TRIM(var_Aname(i)).eq.'calendar') THEN
          IF (((Nrec.eq.1).and.(Nfiles.eq.1)).or.                       &
     &        (INDEX(TRIM(var_Achar(i)),'none').ne.0)) THEN
            Lperpetual=.TRUE.
          END IF
        ELSE IF (TRIM(var_Aname(i)).eq.'cycle_length') THEN
          Lcycle=.TRUE.
        END IF
      END DO
!
!  Turn off the checking of time range if cycling, perpectual, or
!  spectral time axis.
!
      IF (Lcycle.or.Lperpetual.or.Lspectral.or.                         &
     &    ((Nrec.eq.1).and.(Nfiles.eq.1))) THEN
        Lcheck=.FALSE.
      END IF
!
!  Read in time variable minimum and maximun values (input time units).
!
      CALL netcdf_get_time (ng, model, ncname, TvarName,                &
     &                      Rclock%DateNumber, Tmin,                    &
     &                      ncid = ncid,                                &
     &                      start = (/1/),                              &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL netcdf_get_time (ng, model, ncname, TvarName,                &
     &                      Rclock%DateNumber, Tmax,                    &
     &                      ncid = ncid,                                &
     &                      start = (/Nrec/),                           &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Close NetCDF file.
!
      CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
 10   FORMAT (/, ' CHECK_FILE - unable to open grid NetCDF file: ',a)
 20   FORMAT (/, ' CHECK_FILE - unable to find time variable in input', &
     &        ' NetCDF file:', /, 14x, a, /, 14x,                       &
     &        'variable name does not contains the "time" string.')
!
      RETURN
      END FUNCTION check_file
