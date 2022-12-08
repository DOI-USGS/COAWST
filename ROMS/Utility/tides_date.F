#include "cppdefs.h"
#if !defined ANA_TIDES              && \
    (defined TIDE_GENERATING_FORCES || \
     defined SSH_TIDES              || defined UV_TIDES)
      SUBROUTINE tides_date (ng)
!
!git $Id$
!svn $Id: tides_date.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine checks the input tide reference parameter TIDE_START,  !
!  defined as the time of phase zero when preparing the input forcing  !
!  tidal boundary data. The tide reference time is important and often !
!  ignored ROMS parameter by users.                                    !
!                                                                      !
!  Currently, there are two ways to specify the "tide_start" in an     !
!  application:                                                        !
!                                                                      !
!  (1) Set the "zero_phase_date" variable in the input tidal forcing   !
!      NetCDF file (recommended). It a floating-point variable of the  !
!      form YYYYMMDD.dddd with the following metadata:                 !
!                                                                      !
!  double zero_phase_date                                              !
!    zero_phase_date:long_name = "tidal reference date for zero phase" !
!    zero_phase_date:units = "days as %Y%m%d.%f"                       !
!    zero_phase_date:C_format = "%13.4f"                               !
!    zero_phase_date:FORTRAN_format = "(f13.4)"                        !
!                                                                      !
!  For example, if the tide reference date is 2005-01-01 12:00:00, the !
!  zero_phase_date = 20050101.5000. Again, it corresponds to the zero  !
!  tidal phase when preparing the forcing NetCDF file from the tides   !
!  dataset (say, OTPS).                                                !
!                                                                      !
!  Use "forcing/add_tide_date.m" from the ROMS Matlab repository to    !
!  add the "zero_phase_date" variable to your existing tidal forcing   !
!  NetCDF file. It is highly recommended to use this approach. If such !
!  a variable is found, the TIDE_START value will overwritten below.   !
!                                                                      !
!  Notice that it is possible to have different reference values for   !
!  "zero_phase_date" and ROMS clock defined as seconds from reference  !
!  date (standard input parameter TIME_REF).  If "time_ref" is earlier !
!  than "zero_phase_date", the frequencies (omega) to harmonic terms   !
!  will be negative since they are computed as follows:                !
!                                                                      !
!     tide_start = Rclock%tide_DateNumber(2) -                         !
!                  Rclock%DateNumber(2))                               !
!     omega = 2 * pi * (time - tide_start) / Tperiod                   !
!                                                                      !
!  Notice that "tide_start" (in seconds) is recomputed and the value   !
!  specified in the input standard file is ignored.                    !
!                                                                      !
!  (2) The specify TIDE_START in the ROMS standard input file as days  !
!      from the application reference time (TIME_REF) is used it the   !
!      variable "zero_phase_date" is not found in the NetCDF file. In  !
!      the input standard input file we have:                          !
!                                                                      !
!  TIDE_START =  0.0d0                      ! days                     !
!    TIME_REF =  20050101.5d0               ! yyyymmdd.dd              !
!                                                                      !
!  Usually, tide_start = 0 implies that the zero-phase's tidal forcing !
!  date is the same as the application reference date "time_ref". ROMS,!
!  does not know how to check if it is the case.  Thus, it is assumed  !
!  that  the user was carefull  when configuring his/her application.  !
!                                                                      !
!======================================================================!
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
      USE mod_scalars
!
      USE dateclock_mod, ONLY : datenum, datestr
!
!  Imported variables declarations.
!
      integer, intent(in) :: ng
!
!  Local variables declarations.
!
      logical :: foundit
!
      integer :: iday, ihour, isec, iyear, minute, month
!
      real(dp) :: day, sec, zero_phase_date
!
      character (len=19) :: string
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Check if "zero_phase_date" variable is available in input tidal
!  forcing file.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) THEN
        IF (idTref.eq.0) THEN
          IF (Master) WRITE (stdout,10) 'idTref', TRIM(varname)
          exit_flag=5
          RETURN
          END IF
        END IF
!
        SELECT CASE (TIDE(1)%IOtype)
          CASE (io_nf90)
            CALL netcdf_inq_var (ng, iNLM, TIDE(ng)%name,               &
     &                           MyVarName = Vname(1,idTref),           &
     &                           SearchVar = foundit)
            IF (foundit) THEN
              CALL netcdf_get_fvar (ng, iNLM, TIDE(ng)%name,            &
     &                              Vname(1,idTref),                    &
     &                              zero_phase_date)
            END IF

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_inq_var (ng, iNLM, TIDE(ng)%name,           &
     &                               MyVarName = Vname(1,idTref),       &
     &                               SearchVar = foundit)
            IF (foundit) THEN
              CALL pio_netcdf_get_fvar (ng, iNLM, TIDE(ng)%name,        &
     &                                  Vname(1,idTref),                &
     &                                  zero_phase_date)
            END IF
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  If found "zero_phase_date", decode from YYYYMMDD.dddd and compute
!  datenum. Overwrite given "TIDE_START".
!
        IF (foundit) THEN
          iyear=MAX(1,INT(zero_phase_date*0.0001_dp))
          month=MIN(12,MAX(1,INT((zero_phase_date-                      &
     &                            REAL(iyear*10000,dp))*0.01_dp)))
          day=zero_phase_date-AINT(zero_phase_date*0.01_dp)*100.0_dp
          iday=MAX(1,INT(day))
          sec=(day-AINT(day))*86400.0_dp
          ihour=INT(sec/3600.0_dp)
          minute=INT(MOD(sec,3600.0_dp)/60.0_dp)
          isec=INT(MOD(sec,60.0_dp))
          CALL datenum (Rclock%tide_DateNumber, iyear, month, iday,     &
     &                  ihour, minute, REAL(isec,dp))
          CALL datestr (Rclock%tide_DateNumber(1), .TRUE., string)
!
          IF (Master) THEN
            WRITE (stdout,20) 'zero_phase_date = ', zero_phase_date,    &
     &                        TRIM(TIDE(ng)%name),                      &
     &                        'tide_DateNumber = ',                     &
     &                        Rclock%tide_DateNumber(1), TRIM(string),  &
     &                        'old tide_start  = ', tide_start,         &
     &                        ' (days)',                                &
     &                        'new tide_start  = ',                     &
     &                        Rclock%tide_DateNumber(1)-                &
     &                        Rclock%DateNumber(1), ' (days)'
          END IF
          tide_start=Rclock%tide_DateNumber(1)-Rclock%DateNumber(1)
!
!  Otherwise, compute datenum from "tide_start" and "time_ref".
!
        ELSE
          Rclock%tide_DateNumber(1)=Rclock%DateNumber(1)+               &
     &                              tide_start
          Rclock%tide_DateNumber(2)=Rclock%DateNumber(2)+               &
     &                              tide_start*day2sec
          CALL datestr (Rclock%tide_DateNumber(1), .TRUE., string)
!
          IF (Master) THEN
            WRITE (stdout,30) 'zero_phase_date', TRIM(TIDE(ng)%name),   &
     &                        'given tide_start = ', tide_start,        &
     &                        ' (days)',                                &
     &                        'tide_DateNumber  = ',                    &
     &                        Rclock%tide_DateNumber(1), TRIM(string)
        END IF
      END IF
!
  10  FORMAT (/,' TIDES_DATE - Variable index not yet loaded, ', a,     &
     &        /,14x,'Update your metadata file: ',a)
  20  FORMAT (/,2x,'TIDES_DATE - Checking tidal reference date for ',   &
     &        'zero phase: ',/,17x,a,f13.4,' (read from ',a,')',        &
     &        /,17x,a,f13.4,' (',a,')',/,17x,a,f13.4,a,/,17x,a,f13.4,a)
  30  FORMAT (/,2x,'TIDE_DATE - Checking tidal reference date for ',    &
     &        'zero phase: ',/,17x,'''',a,''' variable not found in: ', &
     &        a,/,17x,a,f13.4,a,/,17x,a,f13.4,' (',a,')')
!
      RETURN
      END SUBROUTINE tides_date
#else
      SUBROUTINE tides_date
      END SUBROUTINE tides_date
#endif
