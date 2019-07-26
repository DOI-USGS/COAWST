#include "cppdefs.h"
      SUBROUTINE check_multifile (ng, model)
!
!svn $Id: check_multifile.F 927 2018-10-16 03:51:56Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
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
      USE mod_parallel
      USE mod_ncparam
      USE mod_iounits
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
      logical :: Lcheck, foundit
      logical :: check_file

      integer :: Fcount, Nfiles, i, ifile, lstr

      real(dp) :: Tfinal, Tmax, Tmin, Tscale

      character (len=1), parameter :: blank = ' '
      character (len= 22) :: F_code, I_code, Tmin_code, Tmax_code
      character (len=256) :: ncname
!
      SourceFile=__FILE__
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
        DO i=1,nBCfiles(ng)
          Nfiles=BRY(i,ng)%Nfiles
          DO ifile=1,Nfiles
            ncname=BRY(i,ng)%files(ifile)
            foundit=check_file(ng, model, ncname, Tmin, Tmax, Tscale,   &
     &                         Lcheck)
            BRY(i,ng)%time_min(ifile)=Tmin
            BRY(i,ng)%time_max(ifile)=Tmax
          END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
          Fcount=0
          IF (Lcheck) THEN
            IF (model.eq.iADM) THEN
              DO ifile=Nfiles,1,-1
                Tmax=Tscale*BRY(i,ng)%time_max(ifile)
                IF (time(ng).le.Tmax) THEN
                  Fcount=ifile
                END IF
              END DO
            ELSE
              DO ifile=1,Nfiles
                Tmin=Tscale*BRY(i,ng)%time_min(ifile)
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
            BRY(i,ng)%Fcount=Fcount
            ncname=BRY(i,ng)%files(Fcount)
            lstr=LEN_TRIM(ncname)
            BRY(i,ng)%name=TRIM(ncname)
            BRY(i,ng)%base=ncname(1:lstr-3)
          ELSE
            IF (Master.and.Lcheck) THEN
              WRITE (stdout,10) 'Lateral Boundary', I_code
              DO ifile=1,Nfiles
                Tmin=Tscale*BRY(i,ng)%time_min(ifile)
                Tmax=Tscale*BRY(i,ng)%time_max(ifile)
                CALL time_string(Tmin, Tmin_code)
                CALL time_string(Tmax, Tmax_code)
                WRITE (stdout,20) Tmin_code, Tmax_code,                 &
     &                            TRIM(BRY(i,ng)%files(ifile))
              END DO
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END IF
!
!  Check if there is boundary data to the end of the simulation.
!
          IF (Lcheck) THEN
            IF (model.eq.iADM) THEN
              Tmin=Tscale*BRY(i,ng)%time_min(1)
              IF (Tfinal.lt.Tmin) THEN
                CALL time_string (Tmin, Tmin_code)
                IF (Master) THEN
                  WRITE (stdout,30) 'Lateral Boundary (adjoint)',       &
     &                              TRIM(BRY(i,ng)%files(1)),           &
     &                              'first ', Tmin_code, F_code
                END IF
                exit_flag=4
              END IF
            ELSE
              Tmax=Tscale*BRY(i,ng)%time_max(Nfiles)
              IF (Tfinal.gt.Tmax) THEN
                CALL time_string(Tmax, Tmax_code)
                IF (Master) THEN
                  WRITE (stdout,30) 'Lateral Boundary',                 &
     &                              TRIM(BRY(i,ng)%files(Nfiles)),      &
     &                              'last  ', Tmax_code, F_code
                END IF
                exit_flag=4
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         __FILE__)) RETURN
              END IF
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Input Source/Sink data (river runoff).
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        Nfiles=SSF(ng)%Nfiles
        DO ifile=1,Nfiles
          ncname=SSF(ng)%files(ifile)
          foundit=check_file(ng, model, ncname, Tmin, Tmax, Tscale,     &
     &                       Lcheck)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
          SSF(ng)%time_min(ifile)=Tmin
          SSF(ng)%time_max(ifile)=Tmax
        END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
        Fcount=0
        IF (Lcheck) THEN
          IF (model.eq.iADM) THEN
            DO ifile=Nfiles,1,-1
              Tmax=Tscale*SSF(ng)%time_max(ifile)
              IF (time(ng).le.Tmax) THEN
                Fcount=ifile
              END IF
            END DO
          ELSE
            DO ifile=1,Nfiles
              Tmin=Tscale*SSF(ng)%time_min(ifile)
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
          SSF(ng)%Fcount=Fcount
          ncname=SSF(ng)%files(Fcount)
          lstr=LEN_TRIM(ncname)
          SSF(ng)%name=TRIM(ncname)
          SSF(ng)%base=ncname(1:lstr-3)
        ELSE
          IF (Master.and.Lcheck) THEN
            WRITE (stdout,10) 'Sources/Sinks Data', I_code
            DO ifile=1,Nfiles
              Tmin=Tscale*SSF(ng)%time_min(ifile)
              Tmax=Tscale*SSF(ng)%time_max(ifile)
              CALL time_string (Tmin, Tmin_code)
              CALL time_string (Tmax, Tmax_code)
              WRITE (stdout,20) Tmin_code, Tmax_code,                   &
     &                          TRIM(SSF(ng)%files(ifile))
            END DO
          END IF
          exit_flag=4
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Check if there is Sources/Sinks data up to the end of the simulation.
!
        IF (Lcheck) THEN
          IF (model.eq.iADM) THEN
            Tmin=Tscale*SSF(ng)%time_min(1)
            IF (Tfinal.lt.Tmin) THEN
              CALL time_string (Tmin, Tmin_code)
              IF (Master) THEN
                WRITE (stdout,30) 'Sources/Sinks Data (adjoint)',       &
     &                            TRIM(SSF(ng)%files(1)),               &
     &                            'first ', Tmin_code, F_code
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__)) RETURN
            END IF
          ELSE
            Tmax=Tscale*SSF(ng)%time_max(Nfiles)
            IF (Tfinal.gt.Tmax) THEN
              CALL time_string (Tmax, Tmax_code)
              IF (Master) THEN
                WRITE (stdout,30) 'Sources/Sinks Data',                 &
     &                            TRIM(SSF(ng)%files(Nfiles)),          &
     &                            'last  ', Tmax_code, F_code
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__)) RETURN
            END IF
          END IF
        END IF
      END IF

#ifdef FRC_FILE
!
!-----------------------------------------------------------------------
!  Input forcing data.
!-----------------------------------------------------------------------
!
      DO i=1,nFfiles(ng)
        Nfiles=FRC(i,ng)%Nfiles
        DO ifile=1,Nfiles
          ncname=FRC(i,ng)%files(ifile)
          foundit=check_file(ng, model, ncname, Tmin, Tmax, Tscale,     &
     &                       Lcheck)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
          FRC(i,ng)%time_min(ifile)=Tmin
          FRC(i,ng)%time_max(ifile)=Tmax
        END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
        Fcount=0
        IF (Lcheck) THEN
          IF (model.eq.iADM) THEN
            DO ifile=Nfiles,1,-1
              Tmax=Tscale*FRC(i,ng)%time_max(ifile)
              IF (time(ng).le.Tmax) THEN
                Fcount=ifile
              END IF
            END DO
          ELSE
            DO ifile=1,Nfiles
              Tmin=Tscale*FRC(i,ng)%time_min(ifile)
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
          FRC(i,ng)%Fcount=Fcount
          ncname=FRC(i,ng)%files(Fcount)
          lstr=LEN_TRIM(ncname)
          FRC(i,ng)%name=TRIM(ncname)
          FRC(i,ng)%base=ncname(1:lstr-3)
        ELSE
          IF (Master.and.Lcheck) THEN
            WRITE (stdout,10) 'Forcing', I_code
            DO ifile=1,Nfiles
              Tmin=Tscale*FRC(i,ng)%time_min(ifile)
              Tmax=Tscale*FRC(i,ng)%time_max(ifile)
              CALL time_string (Tmin, Tmin_code)
              CALL time_string (Tmax, Tmax_code)
              WRITE (stdout,20) Tmin_code, Tmax_code,                   &
     &                          TRIM(FRC(i,ng)%files(ifile))
            END DO
          END IF
          exit_flag=4
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Check if there is forcing data up to the end of the simulation.
!
        IF (Lcheck) THEN
          IF (model.eq.iADM) THEN
            Tmin=Tscale*FRC(i,ng)%time_min(1)
            IF (Tfinal.lt.Tmin) THEN
              CALL time_string (Tmin, Tmin_code)
              IF (Master) THEN
                WRITE (stdout,30) 'Forcing (adjoint)',                  &
     &                            TRIM(FRC(i,ng)%files(1)),             &
     &                            'first ', Tmin_code, F_code
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__)) RETURN
            END IF
          ELSE
            Tmax=Tscale*FRC(i,ng)%time_max(Nfiles)
            IF (Tfinal.gt.Tmax) THEN
              CALL time_string (Tmax, Tmax_code)
              IF (Master) THEN
                WRITE (stdout,30) 'Forcing',                            &
     &                            TRIM(FRC(i,ng)%files(Nfiles)),        &
     &                            'last  ', Tmax_code, F_code
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       __FILE__)) RETURN
            END IF
          END IF
        END IF
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Input climatology data.
!-----------------------------------------------------------------------
!
      IF (CLM_FILE(ng)) THEN
        DO i=1,nCLMfiles(ng)
          Nfiles=CLM(i,ng)%Nfiles
          DO ifile=1,Nfiles
            ncname=CLM(i,ng)%files(ifile)
            foundit=check_file(ng, model, ncname, Tmin, Tmax, Tscale,   &
     &                       Lcheck)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                   __FILE__)) RETURN
            CLM(i,ng)%time_min(ifile)=Tmin
            CLM(i,ng)%time_max(ifile)=Tmax
          END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
          Fcount=0
          IF (Lcheck) THEN
            IF (model.eq.iADM) THEN
              DO ifile=Nfiles,1,-1
                Tmax=Tscale*CLM(i,ng)%time_max(ifile)
                IF (time(ng).le.Tmax) THEN
                  Fcount=ifile
                END IF
              END DO
            ELSE
              DO ifile=1,Nfiles
                Tmin=Tscale*CLM(i,ng)%time_min(ifile)
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
            CLM(i,ng)%Fcount=Fcount
            ncname=CLM(i,ng)%files(Fcount)
            lstr=LEN_TRIM(ncname)
            CLM(i,ng)%name=TRIM(ncname)
            CLM(i,ng)%base=ncname(1:lstr-3)
          ELSE
            IF (Master.and.Lcheck) THEN
              WRITE (stdout,10) 'Climatology', I_code
              DO ifile=1,Nfiles
                Tmin=Tscale*CLM(i,ng)%time_min(ifile)
                Tmax=Tscale*CLM(i,ng)%time_max(ifile)
                CALL time_string(Tmin, Tmin_code)
                CALL time_string(Tmax, Tmax_code)
                WRITE (stdout,20) Tmin_code, Tmax_code,                 &
     &                            TRIM(CLM(i,ng)%files(ifile))
              END DO
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN

          END IF
!
!  Check if there is climatology data to the end of the simulation.
!
          IF (Lcheck) THEN
            IF (model.eq.iADM) THEN
              Tmin=Tscale*CLM(i,ng)%time_min(1)
              IF (Tfinal.lt.Tmin) THEN
                CALL time_string (Tmin, Tmin_code)
                IF (Master) THEN
                  WRITE (stdout,30) 'Climatology (adjoint)',            &
     &                              TRIM(CLM(i,ng)%files(1)),           &
     &                              'first ', Tmin_code, F_code
                END IF
                exit_flag=4
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         __FILE__)) RETURN
              END IF
            ELSE
              Tmax=Tscale*CLM(i,ng)%time_max(Nfiles)
              IF (Tfinal.gt.Tmax) THEN
                CALL time_string(Tmax, Tmax_code)
                IF (Master) THEN
                  WRITE (stdout,30) 'Climatology',                      &
     &                              TRIM(CLM(i,ng)%files(Nfiles)),      &
     &                              Tmax_code, F_code
                END IF
                exit_flag=4
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         __FILE__)) RETURN
              END IF
            END IF
          END IF
        END DO
      END IF

 10   FORMAT (/,' CHECK_MULTIFILE - Error while processing ', a,        &
     &        ' multi-files: ',/,19x,'data does not include',           &
     &        ' initialization time = ', a,/)
 20   FORMAT (3x,a,2x,a,5x,a)
 30   FORMAT (/,' CHECK_MULTIFILE - Error while checking input ', a,    &
     &        ' file:',/,19x,a,/,19x,                                   &
     &        a,'data time record available is for day: ',a,/,19x,      &
     &        'but  data is needed to finish run until day: ',a)

      RETURN
      END SUBROUTINE check_multifile
!
      FUNCTION check_file (ng, model, ncname, Tmin, Tmax, Tscale,       &
     &                     Lcheck) RESULT (foundit)
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
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: Lcheck

      integer, intent(in) :: ng, model

      character (*), intent(in) :: ncname

      real(dp), intent(out) :: Tmin, Tmax, Tscale
!
!  Local variable declarations.
!
      logical :: Lcycle, Lperpetual, Lspectral, foundit

      integer :: Nrec, TvarID, i, ncid, nvdim, nvatt

      character (len=40) :: Tunits, TvarName
!
      SourceFile=__FILE__ // ", check_file"
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
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF
!
!  Inquire about all the variables
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Search for the time variable: any 1D array variable with the string
!  'time' in the variable name.
!
      DO i=1,n_var
        IF ((INDEX(TRIM(var_name(i)),'time').ne.0).and.                 &
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
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
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
          IF ((Nrec.eq.1).or.                                           &
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
      IF (Lcycle.or.Lperpetual.or.Lspectral.or.(Nrec.eq.1)) THEN
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
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
      CALL netcdf_get_time (ng, model, ncname, TvarName,                &
     &                      Rclock%DateNumber, Tmax,                    &
     &                      ncid = ncid,                                &
     &                      start = (/Nrec/),                           &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Close NetCDF file.
!
      CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
 10   FORMAT (/, ' CHECK_FILE - unable to open grid NetCDF file: ',a)
 20   FORMAT (/, ' CHECK_FILE - unable to find time variable in input', &
     &        ' NetCDF file:', /, 14x, a, /, 14x,                       &
     &        'variable name does not contains the "time" string.')

      RETURN
      END FUNCTION check_file
