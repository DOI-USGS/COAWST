#include "cppdefs.h"
      MODULE dateclock_mod
!
!git $Id$
!svn $Id: dateclock.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to manage ROMS date, clocks,  !
!  and calendars:                                                      !
!                                                                      !
!    caldate      Converts current model time (days) to calendar date. !
!                   All the returned variables require keyword syntax  !
!                   since they are optional.                           !
!                                                                      !
!    datenum      Converts requested date (year, month, day, ...) into !
!                   a serial number according to the supported calendar!
!                   options.                                           !
!                                                                      !
!    datestr      Converts date number to date string of the form:     !
!                   YYYY-MM-DD hh:mm:ss.ss                             !
!                                                                      !
!    datevec      Converts a given date number to a date vector. It is !
!                   inverse routine to "datenum".                      !
!                                                                      !
!    day_code     Given (month, day, year) it returns a numerical code !
!                   (0 to 6) for the day of the week.                  !
!                                                                      !
!    get_date     Retuns today date string of the form:                !
!                   DayOfWeak - Month day, year - hh:mm:ss ?M          !
!                                                                      !
!    ref_clock    Sets application time clock/reference and loads its  !
!                   to structure Rclock of TYPE T_CLOCK.               !
!                                                                      !
!    ROMS_clock   Given (year, month, day, hour, minutes, seconds),    !
!                   this routine returns ROMS clock time since         !
!                   initialization from the reference date. It is      !
!                   used when importing fields from coupled models.    !
!                                                                      !
!    time_iso8601 Encodes current model time to am ISO 8601 string     !
!                                                                      !
!    time_string  Encodes current model time to a string.              !
!                                                                      !
!    time_units   Decodes time attributes units.                       !
!                                                                      !
!    yearday      Given (year,month,day) this integer function returns !
!                   the day of the year.                               !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC :: caldate
      PUBLIC :: datenum
      PUBLIC :: datestr
      PUBLIC :: datevec
      PUBLIC :: day_code
      PUBLIC :: get_date
      PUBLIC :: ref_clock
      PUBLIC :: ROMS_clock
      PUBLIC :: time_iso8601
      PUBLIC :: time_string
      PUBLIC :: time_units
      PUBLIC :: yearday
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE caldate (CurrentTime,                                  &
     &                    yy_i, yd_i, mm_i, dd_i, h_i, m_i, s_i,        &
     &                    yd_dp, dd_dp, h_dp, m_dp, s_dp)
!***********************************************************************
!                                                                      !
!  This routine converts current model time (in days) to calendar      !
!  date. All the output arguments require keyword syntax since they    !
!  are all optional.  For Example, to get just the fractional (real)   !
!  day-of-year:                                                        !
!                                                                      !
!     CALL caldate (tdays(ng), yd_dp=yday)                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     CurrentTime   Model current time (real; days)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     yy_i          Year including century (integer; OPTIONAL)         !
!     yd_i          Day of the year (integer; OPTIONAL)                !
!     mm_i          Month of the year, 1=Jan, ... (integer; OPTIONAL)  !
!     dd_i          Day of the month (integer; OPTIONAL)               !
!     h_i           Hour of the day (integer; OPTIONAL)                !
!     m_i           Minutes of the hour, 0 - 59 (integer; OPTIONAL)    !
!     s_i           Seconds of the minute (integer; OPTIONAL)          !
!                                                                      !
!     yd_dp         Day of the year (real, fraction; OPTIONAL)         !
!     dd_dp         Day of the month (real, fraction; OPTIONAL)        !
!     h_dp          Hour of the day (real, fraction; OPTION)           !
!     m_dp          Minutes of the hour (real, fraction; OPTION)       !
!     s_dp          Seconds of the minute (real, fraction; OPTIONAL)   !
!                                                                      !
!  Notice that a calendar obtained by extending backward in time from  !
!  its invention or implementation is called the Proleptic version of  !
!  the calendar.                                                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE round_mod, ONLY : ROUND
!
!  Imported variable declarations.
!
      real(dp), intent(in)  :: CurrentTime
!
      integer,  intent(out), optional :: yy_i
      integer,  intent(out), optional :: yd_i
      integer,  intent(out), optional :: mm_i
      integer,  intent(out), optional :: dd_i
      integer,  intent(out), optional :: h_i
      integer,  intent(out), optional :: m_i
      integer,  intent(out), optional :: s_i
!
      real(dp), intent(out), optional :: yd_dp
      real(dp), intent(out), optional :: dd_dp
      real(dp), intent(out), optional :: h_dp
      real(dp), intent(out), optional :: m_dp
      real(dp), intent(out), optional :: s_dp
!
!  Local variable declarations.
!
      logical :: IsDayUnits

      integer :: MyDay, MyHour, MyMinutes, MySeconds
      integer :: MyMonth, MyYday, MyYear

      real(dp) :: DateNumber, DayFraction, RefDateNumber
      real(dp) :: Hour, Minutes, Seconds
!
!-----------------------------------------------------------------------
!  Get calendar date from model current time (days).
!-----------------------------------------------------------------------
!
      RefDateNumber=Rclock%DateNumber(1)              ! fractional days
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since YYYY-MM-DD hh:mm:ss'.  It is called the Gregorian
!  Calendar or Gregorian Proleptic Calendar.
!
      CALENDAR : IF (INT(time_ref).gt.0) THEN
        DateNumber=RefDateNumber+CurrentTime          ! fractional days
        DayFraction=ABS(DateNumber-AINT(DateNumber))
!
        IsDayUnits=.TRUE.
        CALL datevec (DateNumber, IsDayUnits, MyYear, MyMonth, MyDay,   &
     &                MyHour, MyMinutes, Seconds, Minutes, Hour)
        MyYday=yearday(MyYear, MyMonth, MyDay)
        MySeconds=INT(Seconds)
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 0001-01-01 00:00:00'.  It is used in analytical
!  test cases. It has a year length of 365.2425 days (adapted on
!  15 October 1582 by Gregorian Calendar). It is called the Proleptic
!  Gregorian Calendar.
!
      ELSE IF (INT(time_ref).eq.0) THEN
        DateNumber=RefDateNumber+CurrentTime          ! fractional days
        DayFraction=ABS(DateNumber-AINT(DateNumber))
!
        IsDayUnits=.TRUE.
        CALL datevec (DateNumber, IsDayUnits, MyYear, MyMonth, MyDay,   &
     &                MyHour, MyMinutes, Seconds, Minutes, Hour)
        MyYday=yearday(MyYear, MyMonth, MyDay)
        MySeconds=INT(Seconds)
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 0001-01-01 00:00:00'.  It can be used for
!  climatological solutions. It has a year length of 360 days and
!  every month has 30 days.  It is called the 360_day calendar by
!  numerical modelers.
!
      ELSE IF (INT(time_ref).eq.-1) THEN
        DateNumber=RefDateNumber+CurrentTime          ! fractional days
        DayFraction=ABS(DateNumber-AINT(DateNumber))
!
        IsDayUnits=.TRUE.
        CALL datevec (DateNumber, IsDayUnits, MyYear, MyMonth, MyDay,   &
     &                MyHour, MyMinutes, Seconds, Minutes, Hour)
        MyYday=INT(DateNumber-REAL(MyYear*360,dp)+1)
        MySeconds=INT(Seconds)
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 1968-05-23 00:00:00 GMT'. It is a Truncated Julian
!  day introduced by NASA and primarily used by Astronomers. It has
!  a year length of 365.25 days. It is less used nowadays since the length
!  of the year is 648 seconds less (365.2425) resulting in too many leap
!  years.  So it was corrected after 15 October 1582 and it is now called
!  the Gregorian Calendar.
!
      ELSE IF (INT(time_ref).eq.-2) THEN
        IF (CurrentTime.ge.RefDateNumber) THEN        ! fractional day
          DateNumber=CurrentTime                      ! from origin
        ELSE
          DateNumber=RefDateNumber+CurrentTime        ! fractional days
        END IF                                        ! plus truncation
        DayFraction=ABS(DateNumber-AINT(DateNumber))
!
        IsDayUnits=.TRUE.
        CALL datevec (DateNumber, IsDayUnits, MyYear, MyMonth, MyDay,   &
     &                MyHour, MyMinutes, Seconds, Minutes, Hour)
        MyYday=yearday(MyYear, MyMonth, MyDay)
        MySeconds=INT(Seconds)
      END IF CALENDAR
!
!-----------------------------------------------------------------------
!  Load requested time clock values.
!-----------------------------------------------------------------------
!
      IF (PRESENT(yd_i))  yd_i=MyYday
      IF (PRESENT(yy_i))  yy_i=MyYear
      IF (PRESENT(mm_i))  mm_i=MyMonth
      IF (PRESENT(dd_i))  dd_i=MyDay
      IF (PRESENT(h_i ))  h_i =MyHour
      IF (PRESENT(m_i ))  m_i =MyMinutes
      IF (PRESENT(s_i ))  s_i =MySeconds
!
      IF (PRESENT(yd_dp)) yd_dp=REAL(MyYday,dp)+DayFraction
      IF (PRESENT(dd_dp)) dd_dp=REAL(MyDay,dp)+DayFraction
      IF (PRESENT(h_dp )) h_dp =Hour
      IF (PRESENT(m_dp )) m_dp =Minutes
      IF (PRESENT(s_dp )) s_dp =Seconds
!
      RETURN
      END SUBROUTINE caldate
!
!***********************************************************************
      SUBROUTINE datenum (DateNumber,                                   &
     &                    year, month, day, hour, minutes, seconds)
!***********************************************************************
!                                                                      !
!  Converts requested date (year, month, day, ...) into a serial date  !
!  number according to the supported calendars options:                !
!                                                                      !
!  time_ref = -2            Truncated Julian number (Julian/Gregorian) !
!                             'time-units since 1968-05-23 00:00:00'   !
!  time_ref = -1            360_day calendar (Proleptic Gregorian)     !
!                             'time-units since 0000-12-30 00:00:00'   !
!  time_ref = 0             Proleptic Gregorian calendar               !
!                             'time-units since 0001-01-01 00:00:00'   !
!  time_ref = YYYYMMDD.dd   Gregorian or Proleptic Gregorian calendar  !
!                             'time-units since YYYY-MM-DD hh:mm:ss'   !
!                                                                      !
!  For the Proletic Gregogian calendar, the equations are similar to   !
!  the Matlab function "datenum":                                      !
!                                                                      !
!     Matlab:  datenum(0000,00,00)=0       reference date              !
!              datenum(0000,01,01)=1                                   !
!                                                                      !
!  but for simplicity, the equations coded here have have a different  !
!  origin date (Mar 1, 0000) to facilitate the manipulation of leap    !
!  years (adapted from Gary Katch code, Concordia University, Canada)  !
!  yielding:                                                           !
!                                                                      !
!              datenum(0000,03,01)=0       refecence date: Mar 1, 0000 !
!              datenum(0000,01,01)=-59                                 !
!                                                                      !
!  However, to avoid confusion, an offset of 61 days is added to match !
!  Matlab "datenum" function.  The difference between 0000-00-00 and   !
!  0000-03-01 is 61 days.                                              !
!                                                                      !
!  On 15 October 1582, the Gregorian calendar was introduced with a    !
!  year length of 365.2425 days. This is coded as:                     !
!                                                                      !
!     365 + 1/4 - 1/100 + 1/400   or   365 + 0.25 - 0.01 + 0.0025      !
!                                                                      !
!  which is used to account for leap years.  The base of Mar 1, 0000   !
!  is taken for simplicity since the length of february is not fixed.  !
!                                                                      !
!  Notice that a calendar obtained by extending backward in time from  !
!  its invention or implementation is called the Proleptic version of  !
!  the calendar. For example, the Proleptic Gregorian Calendar extends !
!  backwards the date preceding 15 October 1582 with a year length of  !
!  365.2425 days.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     year         Year including the century (integer)                !
!     month        Month of the year: 1=January, ... (integer)         !
!     day          Day of the month (integer)                          !
!     hour         Hour of the day, 0, ... 23 (integer, OPTIONAL)      !
!     minutes      Minutes of the hour (integer, OPTIONAL)             !
!     seconds      Seconds of the minute (real, OPTIONAL)              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DateNumber   Date number (real 1D array),                        !
!                    DateValue(1) => fractional days                   !
!                    DateValue(2) => fractional seconds                !
!                                                                      !
!=======================================================================
!
      USE mod_scalars, ONLY : time_ref
!
!  Imported variable declarations.
!
      integer, intent(in) :: year, month, day

      integer,  intent(in), optional :: hour
      integer,  intent(in), optional :: minutes

      real(dp), intent(in), optional :: seconds

      real(dp), intent(out), dimension(2) :: DateNumber
!
!  Local variable declarations.
!
      integer, parameter :: offset = 61

      integer :: MyDay, MyHour, MyMinutes, MyMonth, MyYear, y01

      real(dp) :: MySeconds
!
!-----------------------------------------------------------------------
!  Initialize optional arguments.
!-----------------------------------------------------------------------
!
      IF (PRESENT(hour)) THEN
        MyHour=hour
      ELSE
        MyHour=0
      END IF
!
      IF (PRESENT(minutes)) THEN
        MyMinutes=minutes
      ELSE
        MyMinutes=0
      END IF
!
      IF (PRESENT(seconds)) THEN
        MySeconds=seconds
      ELSE
        MySeconds=0.0_dp
      END IF
!
!-----------------------------------------------------------------------
!  Date number for the Julian plus Gregorian correction calendar.
!-----------------------------------------------------------------------
!
!  The origin of the Proleptic Julian Calendar is January 1, 4713 BC
!  (November 24, 4713 BC, in the Proleptic Gregorian Calendar).
!  Although the formal definition of Julian day numbers starts and
!  ends at noon, here Julian day starts and ends at midnight. So it
!  is 12 hours faster (substract 12 hours to agree with formal
!  definition).
!
!     datenum(-4713,11,24) = 0              ! Origin: Nov 24, 4713 BC
!     datenum( 1968,05,23) = 2440000        ! Truncated reference (NASA)
!     datenum( 0000,01,01) = 1721060
!
      CALENDAR : IF (INT(time_ref).eq.-2) THEN
        IF (month.gt.2) THEN
          MyYear=year
          MyMonth=month-3
        ELSE
          MyYear=year-1
          MyMonth=month+9
        END IF
        y01=MyYear/100
        MyYear=MyYear-y01*100
        MyDay=(146097*y01/4) + (1461*MyYear/4) + ((153*MyMonth+2)/5) +  &
     &        day + 1721119
!
!-----------------------------------------------------------------------
!  Date mumber for the 360_day Calendar: the year has a length of 360
!  days and every month has 30 days.
!-----------------------------------------------------------------------
!
!     datenum(0000,01,01) = 0
!     datenum(0001,01,01) = 360
!
      ELSE IF (INT(time_ref).eq.-1) THEN
        MyDay=year*360+(month-1)*30+(day-1)
!
!-----------------------------------------------------------------------
!  Date number for the Gregorian and Gregorian Proleptic Calendar. It
!  has a year length of 365.2425 days (correspoding to the Gregorian
!  Calendar introduced in 15 October 1582).
!-----------------------------------------------------------------------
!
!     datenum(0000,01,01) = 1
!     datenum(0001,01,01) = 367
!
      ELSE
        MyMonth=MOD(month+9, 12)                  ! Mar=0, ..., Feb=11
        MyYear=year-INT(0.1_dp*REAL(MyMonth,dp))  ! if Jan or Feb,
!                                                   substract 1
        MyDay=INT(365.0_dp*REAL(MyYear,dp))+                            &
     &        INT(0.25_dp*REAL(MyYear,dp))-                             &
     &        INT(0.01_dp*REAL(MyYear,dp))+                             &
     &        INT(0.0025_dp*REAL(MyYear,dp))+                           &
     &        INT(0.1_dp*(REAL(MyMonth,dp)*306.0_dp + 5.0_dp))+         &
     &        (day - 1)
!
!  Adjust for Matlab origin 0000-00-00 00:00:00, so we get the same
!  value as their function "datenum". The offset is 61 days.
!
!     datenum(0000,00,00) = 0
!
        IF ((year.eq.0).and.(month.eq.0).and.(day.eq.0)) THEN
          MyDay=0;
        ELSE
          IF (MyDay.lt.0) THEN
            MyDay=MyDay+offset-1
          ELSE
            MyDay=MyDay+offset
          END IF
        END IF
      END IF CALENDAR
!
!-----------------------------------------------------------------------
!  Add fractional day to serial day number (day and seconds).
!-----------------------------------------------------------------------
!
!  Fractional date number (units=day).
!
      DateNumber(1)=REAL(MyDay,dp)+                                     &
     &              REAL(MyHour,dp)/24.0_dp+                            &
     &              REAL(MyMinutes,dp)/1440.0_dp+                       &
     &              MySeconds/86400.0_dp
!
!  Fractional date number (units=second).
!
      DateNumber(2)=REAL(MyDay,dp)*86400.0_dp+                          &
     &              REAL(MyHour,dp)*3600.0_dp+                          &
     &              REAL(MyMinutes,dp)*60.0_dp+                         &
     &              MySeconds

      RETURN
      END SUBROUTINE datenum
!
!***********************************************************************
      SUBROUTINE datestr (DateNumber, IsDayUnits, DateString)
!***********************************************************************
!                                                                      !
!  Converts a given date number as computed by "datenum" to a date     !
!  string. Matlab has similar function.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     DateNumber   Date number (real; scalar) as computed by           !
!                    by "datenum":                                     !
!     IsDayUnits   Date number units (logical):                        !
!                    IsDayUnits = .TRUE.   fractional days             !
!                    IsDayUnits = .FALSE.  frational seconds           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DateSring    Date string (YYYY-MM-DD hh:mm:ss.ss)                !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical,  intent(in) :: IsDayUnits

      real(dp), intent(in) :: DateNumber

      character (len=*), intent(out) :: DateString
!
!  Local variable declarations.
!
      integer :: i, year, month, day, hour, minutes

      real(dp):: F_hour, F_minutes, seconds

      character (len= 5) :: sec_string
      character (len=22) :: string
!
!-----------------------------------------------------------------------
!  Compute date vector from serial date number.
!-----------------------------------------------------------------------
!
      CALL datevec (DateNumber, IsDayUnits, year, month, day, hour,     &
     &              minutes, seconds, F_minutes, F_hour)
!
!-----------------------------------------------------------------------
!  Set date string.
!-----------------------------------------------------------------------
!
!  Encode fractional seconds to a string. Round to one digit.
!
      WRITE (sec_string, '(f5.2)') seconds
      DO i=1,LEN(sec_string)                        ! replace leading
        IF (sec_string(i:i).eq.CHAR(32)) THEN       ! space(s) with
          sec_string(i:i)='0'                       ! zeros(s)
        END IF
      END DO
!
!  Encode date string.
!
      WRITE (string,10) year, month, day, hour, minutes, sec_string
 10   FORMAT (i4.4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',a)
!
      DateString=TRIM(string)
!
      RETURN
      END SUBROUTINE datestr
!
!***********************************************************************
      SUBROUTINE datevec (DateNumber, IsDayUnits,                       &
     &                    year, month, day, hour, minutes, seconds,     &
     &                    F_minutes, F_hour)
!***********************************************************************
!                                                                      !
!  Converts a given date number as computed by "datenum" to a date     !
!  vector (year, month, day, hour, minutes, seconds).  It is the       !
!  inverse routine for "datenum" above.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     DateNumber   Date number (real; scalar) as computed by           !
!                    by "datenum":                                     !
!     IsDayUnits   Date number units (logical):                        !
!                    IsDayUnits = .TRUE.   fractional days             !
!                    IsDayUnits = .FALSE.  frational seconds           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     year         Year including the century (integer; YYYY)          !
!     month        Month of the year: 1=January, ... (integer)         !
!     day          Day of the month (integer)                          !
!     hour         Hour of the day, 0, ... 23 (integer)                !
!     minutes      Minutes of the hour (integer)                       !
!     seconds      Seconds of the minute (real)                        !
!                                                                      !
!     F_minutes    Fractional minutes (real)                           !
!     F_hour       Fractional hours (real)                             !
!                                                                      !
!***********************************************************************
!
      USE mod_scalars, ONLY : Rclock, time_ref
      USE round_mod,   ONLY : ROUND
!
!  Imported variable declarations.
!
      logical,  intent(in) :: IsDayUnits

      real(dp), intent(in) :: DateNumber

      integer,  intent(out) :: year, month, day, hour, minutes

      real(dp), intent(out) :: F_hour, F_minutes, seconds
!
!  Local variable declarations.
!
      logical :: ProlepticJulian = .FALSE.

      integer :: MyDay, MyMonth, MyYear, yday
      integer :: ja, jalpha, jb, jc, jd, jday, je

      integer, parameter :: gregorian = 2299161  ! 15 Oct, 1582 A.D.

      real(dp), parameter :: offset = 61.0_dp

      real(dp) :: CT, DayFraction, MyDateNumber
      real(dp) :: dd, jr, js, mo, yy
!
!-----------------------------------------------------------------------
!  Compute date vector from date number for the Julian with Gregorian
!  Calendar correction.
!-----------------------------------------------------------------------
!
!  If truncated, add reference date number (2440000 days). The origin
!  of the Proleptic Julian Calendar is Jan 1, 4713 BC (that is,
!  Nov 24, 4713 BC in the Proleptic Gregorian Calendar).
!
!  Although the formal definition holds that Julian day starts and ends
!  at noon, here Julian day starts and ends at midnight.
!
!  It is assumed that if input DateNumber is greater or equal to the
!  Gregorian Calendar start, its value is full and not Reduced,
!  Modified, or Truncated.
!
      CALENDAR : IF (INT(time_ref).eq.-2) THEN
        IF (IsDayUnits) THEN
          IF (DateNumber.ge.REAL(gregorian,dp)) THEN
            MyDateNumber=DateNumber
          ELSE
            MyDateNumber=DateNumber+Rclock%DateNumber(1)
          END IF
        ELSE
          IF (DateNumber.ge.(REAL(gregorian,dp)*86400.0_dp)) THEN
            MyDateNumber=DateNumber/86400.0_dp
          ELSE
            MyDateNumber=(DateNumber+Rclock%DateNumber(2))/86400.0_dp
          END IF
        END IF
        DayFraction=ABS(MyDateNumber-AINT(MyDateNumber))
!
        IF (ProlepticJulian) THEN            ! Proleptic Julian Calendar
          jday=INT(MyDateNumber)             ! origin: Jan 1, 4713 BC
          IF (jday.ge.gregorian) THEN
            jalpha=INT(((jday-1867216)-0.25_dp)/36524.25_dp)! Gregorian
            ja=jday+1+jalpha-INT(0.25_dp*REAL(jalpha,dp))   ! correction
          ELSE
            ja=jday
          END IF
          jb=ja+1524
          jc=INT(6680.0_dp+(REAL(jb-2439870,dp)-122.1_dp)/365.25_dp)
          jd=365*jc+INT(0.25_dp*REAL(jc,dp))
          je=INT(REAL(jb-jd,dp)/30.6001_dp)
          day=jb-jd-INT(30.6001_dp*REAL(je,dp))
          month=je-1
          IF (month.gt.12) month=month-12
          year=jc-4715
          IF (month.gt.2) year=year-1
          IF (year .le.0) year=year-1
        ELSE                                    ! Proleptic Gregorian
          jr=FLOOR(MyDateNumber)-1721119.0_dp   ! Calendar, origin:
          js=4.0_dp*jr-1.0_dp                   ! Nov 24, 4713 BC
          yy=FLOOR(js/146097.0_dp)
          jr=js-146097.0_dp*yy
          js=FLOOR(jr*0.25_dp)
          js=4.0_dp*js+3.0_dp
          jr=FLOOR(js/1461.0_dp)
          dd=FLOOR(((js-1461.0_dp*jr)+4.0_dp)*0.25_dp)
          js=5.0_dp*dd-3.0_dp
          mo=FLOOR(js/153.0_dp)
          yy=yy*100.0_dp+jr
!
          IF (mo.lt.10.0_dp) THEN
            year =INT(yy)
            month=INT(mo+3)
          ELSE
            year =INT(yy+1)
            month=INT(mo-9)
          END IF
          day=INT(((js-153.0_dp*mo)+5.0_dp)*0.2_dp)
        END IF
!
        seconds=DayFraction*86400.0_dp
        CT=3.0_dp*EPSILON(seconds)           ! comparison tolerance
        seconds=ROUND(seconds, CT)           ! tolerant round function
        F_hour=seconds/3600.0_dp
        hour=INT(F_hour)
        seconds=ABS(seconds-REAL(hour*3600,dp))
        F_minutes=seconds/60.0_dp
        minutes=INT(F_minutes)
        seconds=ABS(seconds-REAL(minutes*60,dp))
!
!-----------------------------------------------------------------------
!  Compute date vector from date mumber for the 360_day calendar: the
!  year has a length of 360 days and every month has 30 day.
!-----------------------------------------------------------------------
!
      ELSE IF (INT(time_ref).eq.-1) THEN
        DayFraction=ABS(DateNumber-AINT(DateNumber))
!
        IF (IsDayUnits) THEN
          year=INT(DateNumber/360.0_dp)
          yday=INT(DateNumber-REAL(year*360,dp)+1)
        ELSE
          year=INT(DateNumber/31104000.0_dp)              ! 360*86400
          yday=INT((DateNumber-REAL(year*31104000,dp)+1)/86400.0_dp)
        END IF
        month=((yday-1)/30)+1
        day=MOD(yday-1,30)+1
!
        seconds=DayFraction*86400.0_dp
        CT=3.0_dp*EPSILON(seconds)           ! comparison tolerance
        seconds=ROUND(seconds, CT)           ! tolerant round function
        F_hour=seconds/3600.0_dp
        hour=INT(F_hour)
        seconds=ABS(seconds-REAL(hour*3600,dp))
        F_minutes=seconds/60.0_dp
        minutes=INT(F_minutes)
        seconds=ABS(seconds-REAL(F_minutes*60,dp))
!
!-----------------------------------------------------------------------
!  Compute date vector from date number for the Gregorian and Gregorian
!  Proleptic Calendar.
!-----------------------------------------------------------------------
!
      ELSE
        IF (IsDayUnits) THEN                        ! fractional days
          MyDateNumber=DateNumber
        ELSE                                        ! fractional seconds
          MyDateNumber=DateNumber/86400.0_dp
        END IF
        DayFraction=ABS(MyDateNumber-AINT(MyDateNumber))
!
        IF (MyDateNumber.lt.offset) THEN            ! adjust for Matlab
          MyDateNumber=MyDateNumber-offset+1.0_dp   ! zero origin,
        ELSE                                        ! datenum(0,0,0)=0,
          MyDateNumber=MyDateNumber-offset          ! 61 days offset
        ENDIF
!
        MyYear=INT((10000.0_dp*AINT(MyDateNumber)+14780.0_dp)/          &
     &             3652425.0_dp)
        MyDay=INT(MyDateNumber)-                                        &
     &        (INT(365.0_dp*REAL(MyYear,dp))+                           &
     &         INT(0.25_dp*REAL(MyYear,dp))-                            &
     &         INT(0.01_dp*REAL(MyYear,dp))+                            &
     &         INT(0.0025_dp*REAL(MyYear,dp)))
        IF (MyDay.lt.0) THEN                        ! if less than Mar 1
          MyYear=MyYear-1                           ! easy on leap-years
          MyDay=INT(MyDateNumber)-                                      &
     &          (INT(365.0_dp*REAL(MyYear,dp))+                         &
     &           INT(0.25_dp*REAL(MyYear,dp))-                          &
     &           INT(0.01_dp*REAL(MyYear,dp))+                          &
     &           INT(0.0025_dp*REAL(MyYear,dp)))
        END IF
        MyMonth=INT((100.0_dp*REAL(MyDay,dp)+ 52.0_dp)/3060.0_dp)
        month=MOD(MyMonth+2, 12) + 1
        year=MyYear+                                                    &
     &       INT((REAL(MyMonth,dp)+2.0_dp)/12.0_dp)
        day=MyDay-                                                      &
     &      INT(0.1_dp*(REAL(MyMonth,dp)*306.0_dp + 5.0_dp)) + 1
!
!  Fix to match Matlab "datestr" function values with the origin at
!  0000-00-00 00:00:00
!
        IF (DateNumber.eq.0.0_dp) THEN
          year=0
          month=1
          day=0
        END IF
!
!  Convert fraction of a day.
!
        seconds=DayFraction*86400.0_dp
        CT=3.0_dp*EPSILON(seconds)           ! comparison tolerance
        seconds=ROUND(seconds, CT)           ! tolerant round function
!
        F_hour=seconds/3600.0_dp
        hour=INT(F_hour)
        seconds=ABS(seconds-REAL(hour*3600,dp))
        F_minutes=seconds/60.0_dp
        minutes=INT(F_minutes)
        seconds=ABS(seconds-REAL(minutes*60,dp))
      END IF CALENDAR
!
      RETURN
      END SUBROUTINE datevec
!
!***********************************************************************
      SUBROUTINE day_code (month, day, year, code)
!***********************************************************************
!                                                                      !
!  This subroutine computes a code for the day of the week, given      !
!  the date. This code is good for dates after:                        !
!                                                                      !
!                              January 1, 1752 AD                      !
!                                                                      !
!  the year the Gregorian calander was adopted in Britian and the      !
!  American colonies.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     month     The month, 1=January, 2=February, ... (integer).       !
!     day       The day of the month (integer).                        !
!     year      The year, including the century (integer).             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     code      A code for the corresponding day of the week           !
!                 (integer):                                           !
!                 code = 0  =>  Sunday                                 !
!                 code = 1  =>  Monday                                 !
!                 code = 2  =>  Tuesday                                !
!                 code = 3  =>  Wednesday                              !
!                 code = 4  =>  Thursday                               !
!                 code = 5  =>  Friday                                 !
!                 code = 6  =>  Saturday                               !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: month, day, year

      integer, intent(out) :: code
!
!  Local variable declarations.
!
      logical :: leap_flag

      integer, parameter :: base_cen = 1700
      integer, parameter :: base_qcen = 1600
      integer, parameter :: base_qyear = 1748
      integer, parameter :: base_year = 1752
      integer, parameter :: bym1_dec31 = 5
      integer, parameter :: feb_end = 59

      integer :: i, leap, no_day, no_yr, nqy, nyc, nyqc

      integer, dimension(12) :: month_day =                             &
     &         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
!
!-----------------------------------------------------------------------
!  Compute the number of years since the base year, the number of
!  years since the beginning of the base century and the number of
!  years since the beginning of the base 400 year.
!-----------------------------------------------------------------------
!
      no_yr=year-base_year
      nqy=year-base_qyear
      nyc=year-base_cen
      nyqc=year-base_qcen
!
!-----------------------------------------------------------------------
!  Compute the number of leapdays in that time.  Determine if this
!  is a leap year.
!-----------------------------------------------------------------------
!
      leap=nqy/4-nyc/100+nyqc/400
      leap_flag=((MOD(nqy,4).eq.0).and.(MOD(nyc,100).ne.0)).or.         &
     &           (MOD(nyqc,400).eq.0)
!
!-----------------------------------------------------------------------
!  Compute the number of days this year.  The leap year corrections
!  are:
!        Jan. 1 - Feb. 28   Have not had the leap day counted above.
!        Feb.29             Counting leap day twice.
!-----------------------------------------------------------------------
!
      no_day=day
      DO i=1,month-1
        no_day=no_day+month_day(i)
      END DO
      IF (leap_flag.and.(no_day.le.feb_end))  no_day=no_day-1
      IF (leap_flag.and.(month.eq.2).and.(day.eq.29)) no_day=no_day-1
!
!-----------------------------------------------------------------------
!  Compute the total number of days since Jan. 1 of the base year,
!  exclusive of the 364 day per year which represent an even 52
!  weeks.  Actually, only need to do the addition mod 7.
!-----------------------------------------------------------------------
!
      no_day=MOD(no_day,7)+MOD(leap,7)+MOD(no_yr,7)+bym1_dec31
!
!-----------------------------------------------------------------------
!  Get the day of the week code.
!-----------------------------------------------------------------------
!
      code=MOD(no_day,7)
      RETURN
      END SUBROUTINE day_code
!
!***********************************************************************
      SUBROUTINE get_date (date_str)
!***********************************************************************
!                                                                      !
!   This routine gets today's date string. It uses intrinsic fortran   !
!   function "date_and_time" and a 12 hour clock.   The string is of   !
!   the form:                                                          !
!                                                                      !
!                DayOfWeak - Month day, year - hh:mm:ss ?M             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     date_str   Today date string, for example:                       !
!                                                                      !
!                  Friday - February 3, 2017 -  3:40:25 PM             !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      character (len=*), intent(out) :: date_str
!
!  Local variable declarations.
!
      integer :: iyear, imonth, iday, ihour, iminute, isecond
      integer :: Dindex, i, half, len1, len2, len3

      integer, dimension(8) :: values

      integer, dimension(31) :: lday =                                  &
     &          (/ (1,i=1,9), (2,i=1,22) /)

      integer, dimension(12) :: lmonth =                                &
     &          (/ 7, 8, 5, 5, 3, 4, 4, 6, 9, 7, 8, 8 /)

      character (len= 5) :: czone
      character (len= 8) :: cdate
      character (len=10) :: ctime
      character (len=11) :: tstring
      character (len=18) :: today
      character (len=20) :: fmt
      character (len=44) :: dstring

      character (len=3), dimension(0:1) :: ampm =                       &
     &                   (/' AM',' PM'/)

      character (len=9), dimension(0:6) :: day =                        &
     &                   (/ 'Sunday   ','Monday   ','Tuesday  ',        &
     &                      'Wednesday','Thursday ','Friday   ',        &
     &                      'Saturday ' /)

      character (len=9), dimension(12) :: month =                       &
     &                   (/ 'January  ','February ','March    ',        &
     &                      'April    ','May      ','June     ',        &
     &                      'July     ','August   ','September',        &
     &                      'October  ','November ','December ' /)
!
!-----------------------------------------------------------------------
!  Get weekday, date and time in short format, then extract its
!  information.
!-----------------------------------------------------------------------
!
      CALL date_and_time (cdate, ctime, czone, values)
!
      iyear=values(1)            ! 4-digit year
      imonth=values(2)           ! month of the year
      iday=values(3)             ! day of the month
      ihour=values(5)            ! hour of the day, local time
      iminute=values(6)          ! minutes of the hour, local time
      isecond=values(7)          ! seconds of the minute, local time
!
!-----------------------------------------------------------------------
!  Convert from 24 hour clock to 12 hour AM/PM clock.
!-----------------------------------------------------------------------
!
      half=ihour/12
      ihour=ihour-half*12
      IF (ihour.eq.0) ihour=12
      IF (half.eq.2) half=0
!
!-----------------------------------------------------------------------
!  Get index for the day of the week.
!-----------------------------------------------------------------------
!
      CALL day_code (imonth, iday, iyear, Dindex)
!
!-----------------------------------------------------------------------
!  Construct date, time and day of the week output string.
!-----------------------------------------------------------------------
!
      WRITE (fmt,10) lmonth(imonth), lday(iday)
 10   FORMAT ('(a',i1,',1x,i',i1,',1h,,1x,i4)')
      WRITE (today,fmt) month(imonth), iday, iyear
      dstring=day(Dindex)
      WRITE (tstring,20) ihour, iminute, isecond, ampm(half)
 20   FORMAT (i2,':',i2.2,':',i2.2,a3)
!
!  Concatenate date string.
!
      len1=LEN_TRIM(dstring)
      len2=LEN_TRIM(today)
      len3=LEN_TRIM(tstring)
      date_str=TRIM(ADJUSTL(dstring(1:len1)))
      IF (len2.gt.0) THEN
        len1=LEN_TRIM(date_str)
        WRITE (date_str,'(a," - ",a)') TRIM(date_str(1:len1)),          &
     &                                 TRIM(today(1:len2))
      END IF
      IF (len3.gt.0) THEN
        len1=LEN_TRIM(date_str)
        WRITE (date_str,'(a," - ",a)') TRIM(date_str(1:len1)),          &
     &                                 TRIM(tstring(1:len3))
      END IF
      RETURN
      END SUBROUTINE get_date
!
!***********************************************************************
      SUBROUTINE ref_clock (r_time)
!***********************************************************************
!                                                                      !
!  This routine encodes the relative time attribute that gives the     !
!  elapsed interval since a specified reference time.  The "units"     !
!  attribute takes the form "time-unit since reference-time".          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     r_time     Time-reference (real; YYYYMMDD.dd; for example,       !
!                  20020115.5 for 15 Jan 2002, 12:0:0).                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Rclock     The time clock base/reference is loaded into module   !
!                  (mod_scalars.F)  structure:                         !
!                                                                      !
!                  Rclock%yday       => day of the year                !
!                  Rclock%year       => year including century (YYYY)  !
!                  Rclock%month      => month of the year              !
!                  Rclock%day        => day of the month               !
!                  Rclock%hour       => hour of the day (0,...,23)     !
!                  Rclock%minutes    => minutes of the hour            !
!                  Rclock%seconds    => seconds of the minute          !
!                  Rclock%base       => reference date (YYYYMMDD.dd)   !
!                  Rclock%DateNumber => date number, 1: days 2: seconds!
!                  Rclock%string     => attribute (YYYY-MM-DD hh:ss:mm)!
!                  Rclock%calendar   => date calendar                  !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: r_time
!
!  Local variable declarations.
!
      integer :: ifac
      integer :: iday, ihour, isec, iyear, leap, minute, month, yday

      real(dp) :: day, sec

      real(dp), dimension(2) :: DateNumber

      character (len=19) :: string
      character (len=20) :: calendar
!
!-----------------------------------------------------------------------
!  Decode reference time.
!-----------------------------------------------------------------------
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since YYYY-MM-DD hh:mm:ss'.
!
      IF (INT(r_time).gt.0) THEN                 ! day 0: Mar 1, 0000
        calendar='proleptic_gregorian'
        iyear=MAX(1,INT(r_time*0.0001_dp))
        month=MIN(12,MAX(1,INT((r_time-REAL(iyear*10000,dp))*0.01_dp)))
        day=r_time-AINT(r_time*0.01_dp)*100.0_dp
        iday=MAX(1,INT(day))
        sec=(day-AINT(day))*86400.0_dp
        ihour=INT(sec/3600.0_dp)
        minute=INT(MOD(sec,3600.0_dp)/60.0_dp)
        isec=INT(MOD(sec,60.0_dp))
        yday=yearday(iyear, month, iday)
        CALL datenum (DateNumber, iyear, month, iday, ihour, minute,    &
     &                REAL(isec,dp))
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 0001-01-01 00:00:00'. It has a year length of
!  365.2425 days
!
      ELSE IF (INT(r_time).eq.0) THEN            ! day 0: Mar 1, 0000
        calendar='proleptic_gregorian'
        iyear=1
        month=1
        iday=1
        ihour=0
        minute=0
        isec=0
        yday=1
        yday=yearday(iyear, month, iday)
        CALL datenum (DateNumber, iyear, month, iday, ihour, minute,    &
     &                REAL(isec,dp))
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 0001-01-01 00:00:00'.  It has a year length of
!  360 days.
!
!  In this calendar, the time in days is simply:
!
!    Time = year * 360 + (month - 1) * 30 + (day - 1)
!
!  And its inverse
!
!    year  = INT(Time / 360)
!    yday  = INT((Time - year * 360) + 1)
!    month = INT(((yday - 1) / 30) + 1)
!    day   = MOD(yday - 1, 30) + 1
!
!  It assumes that the origin (DayNumber=0) corresponds to 01-Jan-0000.
!  However, historically ROMS assumed that DayNumber=1 corresponded to
!  01-Jan-0000 instead. So, there is one day shift. The equations
!  can be manipulated to give either origin, but it is confusing. The
!  above equations are cleaner and now effective (HGA: 30-Jan-2018). The
!  origin (DayNumber=0) occurs on 01-Jan-0000.
!
!  To guarantee compatibility with previous ROMS solutions with this
!  climatological calendar, the reference date is changed to
!
!  'time-units since 0000-12-30 00:00:00'
!
!  to fix the one date shift because DayNumber=0 on 01-Jan-0000. Anyway,
!  it is a highly idealized calendar used in analytical test cases or
!  climatological solutions.
!
      ELSE IF (INT(r_time).eq.-1) THEN           ! day 0: Jan 1, 0000
        calendar='360_day'
        iyear=0
        month=12
        iday=30
        ihour=0
        minute=0
        isec=0
        yday=360
        DateNumber(1)=359.0_dp
        DateNumber(2)=DateNumber(1)*86400.0_dp
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 1968-05-23 00:00:00 GMT'. It is a Truncated Julian
!  day. It has a year length of 365.25 days.
!
!  The one here is known as the Gregorian Calendar.  Although, it is a
!  minor correction to the Julian Calendar after 15 Oct 1582 with a
!  year length of 365.2425.
!
      ELSE IF (INT(r_time).eq.-2) THEN          ! day 0: Nov 24, 4713 BC
        calendar='proleptic_julian'
        iyear=1968
        month=5
        iday=23
        ihour=0
        minute=0
        isec=0
        yday=yearday(iyear, month, iday)
        DateNumber(1)=2440000.0_dp              ! Truncated offset
        DateNumber(2)=DateNumber(1)*86400.0_dp
      END IF
!
!-----------------------------------------------------------------------
!  Set reference-time string, YYYY-MM-DD hh:mm:ss
!-----------------------------------------------------------------------
!
      WRITE (string,10) iyear, month, iday, ihour, minute, isec
 10   FORMAT (i4.4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)
!
!-----------------------------------------------------------------------
!  Load time reference clock information into structure.
!-----------------------------------------------------------------------
!
      Rclock%yday         =yday
      Rclock%year         =iyear
      Rclock%month        =month
      Rclock%day          =iday
      Rclock%hour         =ihour
      Rclock%minutes      =minute
      Rclock%seconds      =isec
      Rclock%base         =r_time
      Rclock%DateNumber(1)=DateNumber(1)
      Rclock%DateNumber(2)=DateNumber(2)
      Rclock%string       =string
      Rclock%calendar     =TRIM(calendar)
!
      RETURN
      END SUBROUTINE ref_clock
!
!**********************************************************************
      SUBROUTINE ROMS_clock (year, month, day, hour, minutes, seconds,  &
     &                       ClockTime)
!***********************************************************************
!                                                                      !
!  Given any date (year, month, day, hour, minute, second), this       !
!  this routine returns ROMS clock time since initialization in        !
!  seconds from reference date.                                        !
!                                                                      !
!  This clock time is used when importing fields from coupled models.  !
!  It is assumed that coupling applications use Gregorian calendar,    !
!  INT(time_ref) .ge. 0.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     year       The year, including the century (integer)             !
!     month      The month, 1=January, 2=February, ... (integer)       !
!     day        The day of the month (integer)                        !
!     hour       The hour of the day (integer)                         !
!     minute     The minute of the hour (integer)                      !
!     seconds    The seconds of the minute (real)                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ClockTime  ROMS clock time since initialization in seconds       !
!                  from reference time (real)                          !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: year, month, day, hour, minutes

      real(dp), intent(in)  :: seconds
      real(dp), intent(out) :: ClockTime
!
!  Local variable declarations.
!
      real(dp), dimension(2) :: DateNumber
!
!-----------------------------------------------------------------------
!  Compute ROMS clock elapsed time since intialization in seconds from
!  reference time.
!-----------------------------------------------------------------------
!
!  Convert requested date into date number.
!
      CALL datenum (DateNumber, year, month, day,                       &
     &              hour, minutes, seconds)
!
!  Compute ROMS clock elapsed time in seconds.
!
      ClockTime=DateNumber(2)-Rclock%DateNumber(2)
!
      RETURN
      END SUBROUTINE ROMS_clock
!
!***********************************************************************
      SUBROUTINE time_iso8601 (MyTime, date_string)
!***********************************************************************
!                                                                      !
!  This routine encodes current model time in seconds to a date        !
!  string of the form:                                                 !
!                                                                      !
!       YYYY-MM-DDThh:mm:ssZ                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     MyTime        Current model time (seconds)                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     date_string  Current model time date string (20 charactes).      !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: MyTime

      character (len=20), intent(out) :: date_string
!
!  Local variable declarations.
!
      integer :: day, hour, minutes, month, year, seconds
      integer :: i
!
      real(dp) :: Currenttime
!
      character (len=20) :: string
!
!-----------------------------------------------------------------------
!  Encode current model time.
!-----------------------------------------------------------------------
!
!  Convert current model time to calendar date.
!
      CurrentTime=MyTime/86400.0_dp                  ! seconds to days
!
      CALL caldate (CurrentTime,                                        &
     &              yy_i=year,                                          &
     &              mm_i=month,                                         &
     &              dd_i=day,                                           &
     &              h_i=hour,                                           &
     &              m_i=minutes,                                        &
     &              s_i=seconds)
!
!  Encode calendar date into a string.
!
      WRITE (string,10) year, month, day, hour, minutes, seconds
 10   FORMAT (i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,'Z')
!
      date_string=TRIM(string)
!
      RETURN
      END SUBROUTINE time_iso8601
!
!***********************************************************************
      SUBROUTINE time_string (MyTime, date_string)
!***********************************************************************
!                                                                      !
!  This routine encodes current model time in seconds to a date        !
!  string of the form:                                                 !
!                                                                      !
!       YYYY-MM-DD hh:mm:ss.ss                                         !
!                                                                      !
!  The decimal seconds (ss.s) are rounded to the next digit. This      !
!  encoding allows an easy to read reporting time.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     MyTime        Current model time (seconds)                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     date_string  Current model time date string (22 charactes).      !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: MyTime

      character (len=22), intent(out) :: date_string
!
!  Local variable declarations.
!
      integer :: day, hour, minutes, month, year
      integer :: i

      real(dp) :: Currenttime, seconds

      character (len= 5) :: sec_string
      character (len=22) :: string
!
!-----------------------------------------------------------------------
!  Encode current model time.
!-----------------------------------------------------------------------
!
!  Convert current model time to calendar date.
!
      CurrentTime=MyTime/86400.0_dp                  ! seconds to days
!
      CALL caldate (CurrentTime,                                        &
     &              yy_i=year,                                          &
     &              mm_i=month,                                         &
     &              dd_i=day,                                           &
     &              h_i=hour,                                           &
     &              m_i=minutes,                                        &
     &              s_dp=seconds)

!
!  Encode fractional seconds to a string. Round to one digit.
!
      WRITE (sec_string, '(f5.2)') seconds
      DO i=1,LEN(sec_string)                        ! replace leading
        IF (sec_string(i:i).eq.CHAR(32)) THEN       ! space(s) with
          sec_string(i:i)='0'                       ! zeros(s)
        END IF
      END DO
!
!  Encode calendar date into a string.
!
      WRITE (string,10) year, month, day, hour, minutes, sec_string
 10   FORMAT (i4.4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',a)
!
      date_string=TRIM(string)
!
      RETURN
      END SUBROUTINE time_string
!
!***********************************************************************
      SUBROUTINE time_units (Ustring, year, month, day,                 &
     &                       hour, minutes, seconds)
!***********************************************************************
!                                                                      !
!  This routine decodes the time units attribute of the form:          !
!                                                                      !
!    'time-units since YYYY-MM-DD hh:mm:ss'                            !
!    'time-units since YYYY-MM-DD hh:mm:ss.ss'                         !
!                                                                      !
!  and various CF compliant variants.                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     U             Time attribute (string)                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     year          Year including century (integer)                   !
!     month         Month of the year, 1=Jan, ..., 12=Dec (integer)    !
!     day           Day of the month (integer)                         !
!     hour          Hour of the day (integer)                          !
!     minutes       Minutes of the hour, 0 - 59 (integer)              !
!     seconds       Seconds of the minute (real)                       !
!                                                                      !
!  Examples of valid unit attributes:                                  !
!                                                                      !
!     'days since 1900-01-01 00:00:00'                                 !
!     'seconds since 1992-10-8 15:15:42.5 -6'                          !
!     'hours since 1990-1-1 0:0:0'                                     !
!     'days since 1582-10-15 1:30:15'                                  !
!     'days since 1-1-1 0:0:0'                                         !
!     'hour since 1997-4-30 1:5:30.5'                                  !
!     'second since 1961-1-1'                                          !
!     'years since -2000-02-29 00:00:0.000Z'                           !
!     'days since 1-07-15 0:0:0'                                       !
!     'days since 0000-01-01 0:0:0'                                    !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(out) :: year, month, day, hour, minutes
!
      real(dp), intent(out) :: seconds
!
      character (len=*), intent(in) :: Ustring
!
!  Exported variable declarations.
!
      logical :: decode
      integer :: i, iblank, ie, is, iscale, lstr, lvar, nval
      integer :: Schar

      real(dp) :: Rval(10)

      character (len=20)       :: Vstring
      character (LEN(Ustring)) :: Tstring
!
!-----------------------------------------------------------------------
!  Decode time string attribute.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      year=0
      month=0
      day=0
      hour=0
      minutes=0
      seconds=0.0_dp
!
      DO i=1,LEN(Tstring)
        Tstring(i:i)=CHAR(32)                            ! blank space
      END DO
!
!  Replace non-numeric characters with blanks.
!
      Tstring=ADJUSTL(TRIM(Ustring))
      lstr=LEN_TRIM(Tstring)
!
!  Only the following ASCII charactes are unchanged:
!
!    Char  Dec  Control Action
!    ------------------------------
!    SP    32   Space
!    +     43   Plus
!    -     45   Hyphen, dash, minus
!    .     46   Period
!    0     48   Zero
!    1     49   One
!    2     50   Two
!    3     51   Three
!    4     52   Four
!    5     53   Five
!    6     54   Six
!    7     55   Seven
!    8     56   Eight
!    9     57   Nine
!
      DO i=1,lstr
        Schar=ICHAR(Tstring(i:i))
        IF (.not.(((48.le.Schar).and.(Schar.le.57)).or.                 &
     &            (Schar.eq.32).or.(Schar.eq.46))) THEN
           Tstring(i:i)=CHAR(32)                          ! blank space
        END IF
      END DO
      Tstring=ADJUSTL(TRIM(Tstring))
      lstr=LEN_TRIM(Tstring)
!
!  Check for negative year indicating CE, BC or BCE (Common Era, Before
!  Christ or Before Common Era).
!
      IF (INDEX(Ustring, 'since -').gt.0) THEN
        iscale=-1
      ELSE
        iscale=1
      END IF
!
!  Process numrical values.  Since CHAR(45) is retained, take the
!  absolute value except for the first number representing the year.
!  The year is the only numerical value that can be negative (BC or
!  BCE.
!
      is=1
      ie=lstr
      iblank=0
      nval=0
      decode=.FALSE.
      DO i=1,lstr
        IF (Tstring(i:i).eq.CHAR(32)) THEN
          IF (Tstring(i+1:i+1).ne.CHAR(32)) decode=.TRUE.
          iblank=i
        ELSE
          ie=i
        END IF
        IF (decode.or.(i.eq.lstr)) THEN
          nval=nval+1
          Vstring=Tstring(is:ie)
          lvar=LEN_TRIM(Vstring)
          READ (Vstring(1:lvar),*) Rval(nval)
          is=iblank+1
          ie=lvar
          decode=.FALSE.
        END IF
      END DO
!
!  Load values.
!
      DO i=1,nval
        SELECT CASE (i)
          CASE (1)
            year=INT(Rval(i))*iscale
          CASE (2)
            month=INT(Rval(i))
          CASE (3)
            day=INT(Rval(i))
          CASE (4)
            hour=INT(Rval(i))
          CASE (5)
            minutes=INT(Rval(i))
          CASE (6)
            seconds=Rval(i)
        END SELECT
      END DO
      RETURN
      END SUBROUTINE time_units
!
!***********************************************************************
      INTEGER FUNCTION yearday (year, month, day) RESULT (yday)
!***********************************************************************
!                                                                      !
!  Given any date year, month, and day, this function returns the      !
!  day of the year.                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     year       Year including the century (integer; YYYY)            !
!     month      Month of the year: 1=January, ... (integer)           !
!     day        Day of the month (integer)                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     yday       Day of the year (integer)                             !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: year, month, day
!
!  Local variable declarations.
!
      integer :: fac
!
!-----------------------------------------------------------------------
!  Compute day of the year.
!-----------------------------------------------------------------------
!
      IF (((MOD(year,4).eq.0).and.(MOD(year,100).ne.0)).or.             &
     &    (MOD(year,400).eq.0)) THEN
        fac=1                                               ! leap year
      ELSE
        fac=2
      END IF
      yday=INT((275.0*month)/9) - fac*INT((month+9)/12) + day - 30

      RETURN
      END FUNCTION yearday
!
      END MODULE dateclock_mod
