#include "cppdefs.h"
      SUBROUTINE check_multifile (ng, model)
!
!svn $Id: check_multifile.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      USE mod_iounits
      USE mod_scalars
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

      real(r8) :: Tfinal, Tmax, Tmin, Tscale

      character(len=14 ) :: F_code, I_code, Tmin_code, Tmax_code
      character(len=256) :: ncname
!
!=======================================================================
!  If applicable, initialize parameters for input multi-files.
!=======================================================================
!
!  Get initialization time string.
!
      CALL time_string(time(ng), I_code)
!
!  Get final time string for simulation.
!
      Tfinal=dstart*day2sec+ntimes(ng)*dt(ng)
      CALL time_string(tfinal, F_code)
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
            IF (exit_flag.ne.NoError) RETURN
            BRY(i,ng)%time_min(ifile)=Tmin
            BRY(i,ng)%time_max(ifile)=Tmax
          END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart. The EXIT below is removed because when restarting there is
!  a possibility that the restart time is not included in the first
!  boundary file in the list to avoid getting an IO error.
!
          Fcount=0
          IF (Lcheck) THEN
            DO ifile=1,Nfiles
              Tmin=Tscale*BRY(i,ng)%time_min(ifile)
              IF (time(ng).ge.Tmin) THEN
                Fcount=ifile
!!              EXIT
              END IF
            END DO
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
          END IF
!
!  Check if there is boundary data to the end of the simulation.
!
          IF (Lcheck) THEN
            Tmax=Tscale*BRY(i,ng)%time_max(Nfiles)
            IF (Tfinal.gt.Tmax) THEN
              CALL time_string(Tmax, Tmax_code)
              IF (Master) THEN
                WRITE (stdout,30) 'Lateral Boundary',                   &
     &                            TRIM(BRY(i,ng)%files(Nfiles)),        &
     &                            Tmax_code, F_code
              END IF
            END IF
          END IF
        END DO
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
          IF (exit_flag.ne.NoError) RETURN
          FRC(i,ng)%time_min(ifile)=Tmin
          FRC(i,ng)%time_max(ifile)=Tmax
        END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart. The EXIT below is removed because when restarting there is
!  a possibility that the restart time is not included in the first
!  forcing file in the list to avoid getting an IO error.
!
        Fcount=0
        IF (Lcheck) THEN
          DO ifile=1,Nfiles
            Tmin=Tscale*FRC(i,ng)%time_min(ifile)
            IF (time(ng).ge.Tmin) THEN
              Fcount=ifile
!!            EXIT
            END IF
          END DO
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
              CALL time_string(Tmin, Tmin_code)
              CALL time_string(Tmax, Tmax_code)
              WRITE (stdout,20) Tmin_code, Tmax_code,                   &
     &                          TRIM(FRC(i,ng)%files(ifile))
            END DO
          END IF
          exit_flag=4
        END IF
!
!  Check if there is focing data to the end of the simulation.
!
        IF (Lcheck) THEN
          Tmax=Tscale*FRC(i,ng)%time_max(Nfiles)
          IF (Tfinal.gt.Tmax) THEN
            CALL time_string(Tmax, Tmax_code)
            IF (Master) THEN
              WRITE (stdout,30) 'Forcing',                              &
     &                          TRIM(FRC(i,ng)%files(Nfiles)),          &
     &                          Tmax_code, F_code
            END IF
            exit_flag=4
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
            IF (exit_flag.ne.NoError) RETURN
            CLM(i,ng)%time_min(ifile)=Tmin
            CLM(i,ng)%time_max(ifile)=Tmax
          END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
          Fcount=0
          IF (Lcheck) THEN
            DO ifile=1,Nfiles
              Tmin=Tscale*CLM(i,ng)%time_min(ifile)
              IF (time(ng).ge.Tmin) THEN
                Fcount=ifile
              END IF
            END DO
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
          END IF
!
!  Check if there is climatology data to the end of the simulation.
!
          IF (Lcheck) THEN
            Tmax=Tscale*CLM(i,ng)%time_max(Nfiles)
            IF (Tfinal.gt.Tmax) THEN
              CALL time_string(Tmax, Tmax_code)
              IF (Master) THEN
                WRITE (stdout,30) 'Climatology',                        &
     &                            TRIM(CLM(i,ng)%files(Nfiles)),        &
     &                            Tmax_code, F_code
              END IF
              exit_flag=4
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
     &        'last data time record available is for day: ',a,/,19x,   &
     &        'but data is needed to finish run until day: ',a)

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
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: Lcheck

      integer, intent(in) :: ng, model

      character (*), intent(in) :: ncname

      real(r8), intent(out) :: Tmin, Tmax, Tscale
!
!  Local variable declarations.
!
      logical :: Lcycle, Lperpetual, Lspectral, foundit

      integer :: Nrec, TvarID, i, nvdim, nvatt

      character (len=40) :: Tunits, TvarName
!
      SourceFile='check_multifile.F, check_file'
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
!  Inquire about all the variables
!
      CALL netcdf_inq_var (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
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
          WRITE (stdout,10) TRIM(ncname)
        END IF
        exit_flag=4
      END IF
!
!  Inquire about requested variable.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     MyVarName = TRIM(TvarName),                  &
     &                     VarID = TvarID,                              &
     &                     nVarDim = nvdim,                             &
     &                     nVarAtt = nvatt)
      IF (exit_flag.ne.NoError) RETURN
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
      IF (Lcycle.or.Lperpetual.or.Lspectral) THEN
        Lcheck=.FALSE.
      END IF
!
!  Read in time variable minimum and maximun values (input time units).
!
      CALL netcdf_get_fvar (ng, model, ncname, TvarName,                &
     &                      Tmin,                                       &
     &                      start = (/1/),                              &
     &                      total = (/1/))
      IF (exit_flag.ne.NoError) RETURN
!
      CALL netcdf_get_fvar (ng, model, ncname, TvarName,                &
     &                      Tmax,                                       &
     &                      start = (/Nrec/),                           &
     &                      total = (/1/))
      IF (exit_flag.ne.NoError) RETURN

 10   FORMAT (/, ' CHECK_FILE - unable to find time variable in input', &
     &        ' NetCDF file:', /, 14x, a, /, 14x,                       &
     &        'variable name does not contains the "time" string.')

      RETURN
      END FUNCTION check_file
