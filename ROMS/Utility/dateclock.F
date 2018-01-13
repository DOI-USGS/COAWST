#include "cppdefs.h"
      MODULE dateclock_mod
!
!svn $Id: dateclock.F 853 2017-07-01 02:24:45Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
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
!                   a serial number. Similar to Matlab "datenum" but   !
!                   the value datenum=0 corresponds to Mar 1, 0000.    !
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
!    time_string  Encodes current model time to a string.              !
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
      PUBLIC :: datevec
      PUBLIC :: day_code
      PUBLIC :: get_date
      PUBLIC :: ref_clock
      PUBLIC :: ROMS_clock
      PUBLIC :: time_string
      PUBLIC :: yearday
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE caldate (CurrentTime,                                  &
     &                    yy_i, yd_i, mm_i, dd_i, h_i, m_i, s_i,        &
     &                    yd_r8, dd_r8, h_r8, m_r8, s_r8)
!***********************************************************************
!                                                                      !
!  This routine converts current model time (in days) to calendar      !
!  date. All the output arguments require keyword syntax since they    !
!  are all optional.  For Example, to get just the fractional (real)   !
!  day-of-year:                                                        !
!                                                                      !
!     CALL caldate (tdays(ng), yd_r8=yday)                             !
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
!     m_i           Minutes of the hour, 1 - 23 (integer; OPTIONAL)    !
!     s_i           Seconds of the minute (integer; OPTIONAL)          !
!                                                                      !
!     yd_r8         Day of the year (real, fraction; OPTIONAL)         !
!     dd_r8         Day of the month (real, fraction; OPTIONAL)        !
!     h_r8          Hour of the day (real, fraction; OPTION)           !
!     m_r8          Minutes of the hour (real, fraction; OPTION)       !
!     s_r8          Seconds of the minute (real, fraction; OPTIONAL)   !
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
      real(r8), intent(in)  :: CurrentTime
!
      integer,  intent(out), optional :: yy_i
      integer,  intent(out), optional :: yd_i
      integer,  intent(out), optional :: mm_i
      integer,  intent(out), optional :: dd_i
      integer,  intent(out), optional :: h_i
      integer,  intent(out), optional :: m_i
      integer,  intent(out), optional :: s_i
!
      real(r8), intent(out), optional :: yd_r8
      real(r8), intent(out), optional :: dd_r8
      real(r8), intent(out), optional :: h_r8
      real(r8), intent(out), optional :: m_r8
      real(r8), intent(out), optional :: s_r8
!
!  Local variable declarations.
!
      logical :: IsDayUnits

      integer, parameter :: gregorian = 2299161  ! 15 Oct, 1582 A.D.

      integer :: MyDay, MyHour, MyMinutes, MySeconds
      integer :: MyMonth, MyYday, MyYear
      integer :: ja, jalpha, jb, jc, jd, jday, je

      real(r8) :: CT, DateNumber
      real(r8) :: DayFraction, Hour, Minutes, Seconds
!
!-----------------------------------------------------------------------
!  Get calendar date from model current time (days).
!-----------------------------------------------------------------------
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since YYYY-MM-DD hh:mm:ss'.  It is called the Gregorian
!  Proleptic Calendar.
!
      IF (INT(time_ref).gt.0) THEN
        DateNumber=Rclock%Dnumber+CurrentTime         ! fractional days
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
!  test cases. It has a year length of 365.2425 days (that is, Gregorian
!  Calendar adapted in 15 October 1582). It is called the Gregorian
!  Proleptic Calendar.
!
      ELSE IF (INT(time_ref).eq.0) THEN
        DateNumber=Rclock%Dnumber+CurrentTime         ! fractional days
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
        DateNumber=Rclock%Dnumber+CurrentTime
        DayFraction=ABS(DateNumber-AINT(DateNumber))
!
        MyYear=INT(DateNumber/360.0_r8)
        MyYday=INT(MOD(DateNumber,360.0_r8))
        IF (MyYday.eq.0) THEN
          MyYday=360            ! last day of current year, MOD(x,360)=0
          MyYear=MyYear-1       ! susbract one, new year tomorrow
        END IF
        IF (MOD(MyYday,30).eq.0) THEN
          MyMonth=MyYday/30
          MyDay=30              ! last day of the month, MOD(x,30)=0
        ELSE
          MyMonth=(MyYday/30)+1
          MyDay=MOD(MyYday,30)
        END IF
!
        Seconds=DayFraction*86400.0_r8
        CT=3.0_r8*EPSILON(Seconds)           ! comparison tolerance
        Seconds=ROUND(Seconds, CT)           ! tolerant round function
        Hour=Seconds/3600.0_r8
        MyHour=INT(Hour)
        Seconds=ABS(Seconds-REAL(MyHour*3600,r8))
        Minutes=Seconds/60.0_r8
        MyMinutes=INT(Minutes)
        Seconds=ABS(Seconds-REAL(MyMinutes*60,r8))
        MySeconds=INT(Seconds)
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 1968-05-23 00:00:00 GMT'. It is a Truncated Julian
!  day introduced by NASA and primarilily used by Astronomers. It has
!  a year length of 365.25 days. It is less used nowadays since the length
!  of the year is 648 seconds less (365.2425) resulting in too many leap
!  years.  So it is correct after 15 October 1582 and it is now called
!  the Gregorian Calendar.
!
      ELSE IF (INT(time_ref).eq.-2) THEN
        DayFraction=ABS(CurrentTime-AINT(CurrentTime))
!
        IF (CurrentTime.ge.Rclock%Dnumber) THEN
          jday=INT(CurrentTime)                 ! Origin: Jan 1, 4713 BC
        ELSE
          jday=INT(Rclock%Dnumber+CurrentTime)  ! Truncated Julian Day
        END IF                                  ! add 2440000 offset
        IF (jday.ge.gregorian) THEN
          jalpha=INT(((jday-1867216)-0.25_r8)/36524.25_r8)  ! Gregorian
          ja=jday+1+jalpha-INT(0.25_r8*REAL(jalpha,r8))     ! correction
        ELSE
          ja=jday
        END IF
        jb=ja+1524
        jc=INT(6680.0_r8+(REAL(jb-2439870,r8)-122.1_r8)/365.25_r8)
        jd=365*jc+INT(0.25_r8*REAL(jc,r8))
        je=INT(REAL(jb-jd,r8)/30.6001_r8)
        MyDay=jb-jd-INT(30.6001_r8*REAL(je,r8))
        MyMonth=je-1
        IF (MyMonth.gt.12) MyMonth=MyMonth-12
        MyYear=jc-4715
        IF (MyMonth.gt.2) MyYear=MyYear-1
        IF (MyYear .le.0) MyYear=MyYear-1
        MyYday=yearday(MyYear, MyMonth, MyDay)
!
        Seconds=DayFraction*86400.0_r8
        CT=3.0_r8*EPSILON(Seconds)           ! comparison tolerance
        Seconds=ROUND(Seconds, CT)           ! tolerant round function
        Hour=Seconds/3600.0_r8
        MyHour=INT(Hour)
        Seconds=ABS(Seconds-REAL(MyHour*3600,r8))
        Minutes=Seconds/60.0_r8
        MyMinutes=INT(Minutes)
        Seconds=ABS(Seconds-REAL(MyMinutes*60,r8))
        MySeconds=INT(Seconds)
      END IF
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
      IF (PRESENT(yd_r8)) yd_r8=REAL(MyYday,r8)+DayFraction
      IF (PRESENT(dd_r8)) dd_r8=REAL(MyDay,r8)+DayFraction
      IF (PRESENT(h_r8 )) h_r8 =Hour
      IF (PRESENT(m_r8 )) m_r8 =Minutes
      IF (PRESENT(s_r8 )) s_r8 =Seconds
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
!  number.  It is similar to Matlab function "datenum", except for the !
!  base date for day 0:                                                !
!                                                                      !
!     Matlab:  datenum(0000,00,00)=0       reference date              !
!              datenum(0000,01,01)=1                                   !
!                                                                      !
!       Here:  datenum(0000,03,01)=0       refecence date: Mar 1, 0000 !
!              datenum(0000,01,01)=-59                                 !
!                                                                      !
!  They just have different reference date.  On 15 October 1582, the   !
!  Gregorian was adapted with a year length of 365.2425 days.  This    !
!  is coded as:                                                        !
!                                                                      !
!     365 + 1/4 - 1/100 + 1/400   or   365 + 0.25 - 0.01 + 0.0025      !
!                                                                      !
!  which is used to account for leap years.  The base of Mar 1, 0000   !
!  is taken for simplicity since the length of february is not fixed.  !
!                                                                      !
!  Although this routine and the Matlab function are not equivalent,   !
!  they yield the same elapsed time between two dates:                 !
!                                                                      !
!     datenum(2017,1,1)-datenum(0000, 3,1) = 736635                    !
!     datenum(2014,3,1)-datenum(1582,10,1) = 158302                    !
!                                                                      !
!  The calendar here is the Gregorian Proleptic Calendar because it    !
!  extends backwards the date preceeding 15 October 1582 with a        !
!  year length of 365.2425 days.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     year         Year including the century (integer)                !
!     month        Month of the year: 1=January, ... (integer)         !
!     day          Day of the month (integer)                          !
!     hour         Hour of the day, 1, ... 23 (integer, OPTIONAL)      !
!     minutes      Minutes of the hour (integer, OPTIONAL)             !
!     seconds      Seconds of the minute (real, OPTIONAL)              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DateNumber   Date number (real 1D array),                        !
!                    DateValue(1) => fractional days                   !
!                    DateValue(2) => fractional seconds                !
!                                                                      !
!  Adapted from Gary Katch, Concordia University, Canada.              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: year, month, day

      integer,  intent(in), optional :: hour
      integer,  intent(in), optional :: minutes

      real(r8), intent(in), optional :: seconds

      real(r8), intent(out), dimension(2) :: DateNumber
!
!  Local variable declarations.
!
      integer :: MyDay, MyHour, MyMinutes, MyMonth, MyYear

      real(r8) :: MySeconds
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
        MySeconds=0.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  Compute fractional serial date number.
!-----------------------------------------------------------------------
!
      MyMonth=MOD(month+9, 12)                    ! Mar=0, ..., Feb=11
      MyYear=year-INT(0.1_r8*REAL(MyMonth,r8))    ! if Jan or Feb,
!                                                   substract 1
      MyDay=INT(365.0_r8*REAL(MyYear,r8))+                              &
     &      INT(0.25_r8*REAL(MyYear,r8))-                               &
     &      INT(0.01_r8*REAL(MyYear,r8))+                               &
     &      INT(0.0025_r8*REAL(MyYear,r8))+                             &
     &      INT(0.1_r8*(REAL(MyMonth,r8)*306.0_r8 + 5.0_r8))+           &
     &      (day - 1)
!
!  Fractional date number (units=day).
!
      DateNumber(1)=REAL(MyDay,r8)+                                     &
     &              REAL(MyHour,r8)/24.0_r8+                            &
     &              REAL(MyMinutes,r8)/1440.0_r8+                       &
     &              MySeconds/86400.0_r8
!
!  Fractional date number (units=second).
!
      DateNumber(2)=REAL(MyDay,r8)*86400.0_r8+                          &
     &              REAL(MyHour,r8)*360.0_r8+                           &
     &              REAL(MyMinutes,r8)*60.0_r8+                         &
     &              MySeconds

      RETURN
      END SUBROUTINE datenum
!
!***********************************************************************
      SUBROUTINE datevec (DateNumber, IsDayUnits,                       &
     &                    year, month, day, hour, minutes, seconds,     &
     &                    F_minutes, F_hour)
!***********************************************************************
!                                                                      !
!  Converts a given date number as computed by "datenum" to a date     !
!  vector (year, month, day, hour, minutes, seconds).  It is the       !
!  inverse routine for "datenum" abobe.  Matlab has similar function.  !
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
!     hour         Hour of the day, 1, ... 23 (integer)                !
!     minutes      Minutes of the hour (integer)                       !
!     seconds      Seconds of the minute (real)                        !
!                                                                      !
!     F_minutes    Fractional minutes (real)                           !
!     F_hour       Fractional hours (real)                             !
!                                                                      !
!  Adapted from Gary Katch, Concordia University, Canada.              !
!                                                                      !
!***********************************************************************
!
      USE round_mod, ONLY : ROUND
!
!  Imported variable declarations.
!
      logical,  intent(in) :: IsDayUnits

      real(r8), intent(in) :: DateNumber

      integer,  intent(out) :: year, month, day, hour, minutes

      real(r8), intent(out) :: F_hour, F_minutes, seconds
!
!  Local variable declarations.
!
      integer :: MyDay, MyMonth, MyYear

      real(r8) :: CT, DayFraction, MyDateNumber
!
!-----------------------------------------------------------------------
!  Initialize according to input units (days or seconds).
!-----------------------------------------------------------------------
!
      IF (IsDayUnits) THEN                          ! fractional days
        MyDateNumber=DateNumber
      ELSE                                          ! fractional seconds
        MyDateNumber=DateNumber/86400.0_r8
      END IF
      DayFraction=ABS(MyDateNumber-AINT(MyDateNumber))
!
!-----------------------------------------------------------------------
!  Compute date vector (year, month, day, hour, minutes, seconds).
!-----------------------------------------------------------------------
!
      MyYear=INT((10000.0_r8*AINT(MyDateNumber)+14780.0_r8)/            &
     &           3652425.0_r8)
      MyDay=INT(MyDateNumber)-                                          &
     &      (INT(365.0_r8*REAL(MyYear,r8))+                             &
     &       INT(0.25_r8*REAL(MyYear,r8))-                              &
     &       INT(0.01_r8*REAL(MyYear,r8))+                              &
     &       INT(0.0025_r8*REAL(MyYear,r8)))
      IF (MyDay.lt.0) THEN                          ! if less than Mar 1
        MyYear=MyYear-1
        MyDay=INT(MyDateNumber)-                                        &
     &        (INT(365.0_r8*REAL(MyYear,r8))+                           &
     &         INT(0.25_r8*REAL(MyYear,r8))-                            &
     &         INT(0.01_r8*REAL(MyYear,r8))+                            &
     &         INT(0.0025_r8*REAL(MyYear,r8)))
      END IF
      MyMonth=INT((100.0_r8*REAL(MyDay,r8)+ 52.0_r8)/3060.0_r8)
      month=MOD(MyMonth+2, 12) + 1
      year=MyYear+                                                      &
     &     INT((REAL(MyMonth,r8)+2.0_r8)/12.0_r8)
      day=MyDay-                                                        &
     &    INT(0.1_r8*(REAL(MyMonth,r8)*306.0_r8 + 5.0_r8)) + 1
!
!  Convert fraction of a day.
!
      seconds=DayFraction*86400.0_r8
      CT=3.0_r8*EPSILON(seconds)             ! comparison tolerance
      seconds=ROUND(seconds, CT)             ! tolerant round function
!
      F_hour=seconds/3600.0_r8
      hour=INT(F_hour)
      seconds=ABS(seconds-REAL(hour*3600,r8))
      F_minutes=seconds/60.0_r8
      minutes=INT(F_minutes)
      seconds=ABS(seconds-REAL(minutes*60,r8))
!
      RETURN
      END SUBROUTINE datevec
!
!***********************************************************************
      SUBROUTINE day_code (month, day, year, code)
!***********************************************************************
!                                                                      !
!  This subroutine computes a code for the day of the week, given      !
!  the date. This code is good for date after:                         !
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
!   This routine gets today date string.  It uses intrinsic fortran    !
!   function "date_and_time" and a 12 hour clock.  The string is of    !
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
!                  Rclock%yday     => day of the year                  !
!                  Rclock%year     => year including century (YYYY)    !
!                  Rclock%month    => month of the year                !
!                  Rclock%day      => day of the month                 !
!                  Rclock%hour     => hour of the day (1,...,23)       !
!                  Rclock%minutes  => minutes of the hour              !
!                  Rclock%seconds  => fractionak seconds of the minute !
!                  Rclock%base     => reference date (YYYYMMDD.dd)     !
!                  Rclock%Dnumber  => date number, fractional days     !
!                  Rclock%Snumber  => date number, fractional seconds  !
!                  Rclock%string   => attribute (YYYY-MM-DD hh:ss:mm)  !
!                  Rclock%calendar => date calendar                    !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: r_time
!
!  Local variable declarations.
!
      integer :: ifac
      integer :: iday, ihour, isec, iyear, leap, minute, month, yday

      real(r8) :: day, sec

      real(r8), dimension(2) :: DateNumber

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
        calendar='gregorian_proleptic'
        iyear=MAX(1,INT(r_time*0.0001_r8))
        month=MIN(12,MAX(1,INT((r_time-REAL(iyear*10000,r8))*0.01_r8)))
        day=r_time-AINT(r_time*0.01_r8)*100.0_r8
        iday=MAX(1,INT(day))
        sec=(day-AINT(day))*86400.0_r8
        ihour=INT(sec/3600.0_r8)
        minute=INT(MOD(sec,3600.0_r8)/60.0_r8)
        isec=INT(MOD(sec,60.0_r8))
        yday=yearday(iyear, month, iday)
        CALL datenum (DateNumber, iyear, month, iday, ihour, minute,    &
     &                REAL(isec,r8))
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 0001-01-01 00:00:00'. It has a year length of
!  365.2425 days
!
      ELSE IF (INT(r_time).eq.0) THEN            ! day 0: Mar 1, 0000
        calendar='gregorian_proleptic'
        iyear=1
        month=1
        iday=1
        ihour=0
        minute=0
        isec=0
        yday=1
        yday=yearday(iyear, month, iday)
        CALL datenum (DateNumber, iyear, month, iday, ihour, minute,    &
     &                REAL(isec,r8))
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 0001-01-01 00:00:00'. It has a year length of
!  360 days.
!
      ELSE IF (INT(r_time).eq.-1) THEN           ! day 0: Jan 1, 0000
        calendar='360_day'
        iyear=1
        month=1
        iday=1
        ihour=0
        minute=0
        isec=0
        yday=MAX(1,MOD(iday,360))
        DateNumber(1)=360.0_r8
        DateNumber(2)=DateNumber(1)*86400.0_r8
!
!  The model clock is the elapsed time since reference time of the form
!  'time-units since 1968-05-23 00:00:00 GMT'. It is a Truncated Julian
!  day. It has a year length of 365.25 days.
!
!  The one here is known as the Gregorian Calendar.  Altough, it is a
!  minor correction to the Julian Calendar after 15 Oct 1582 with a
!  year length of 365.2425.
!
      ELSE IF (INT(r_time).eq.-2) THEN           ! day 0: Jan 1, 4713 BC
        calendar='gregorian'
        iyear=1968
        month=5
        iday=23
        ihour=0
        minute=0
        isec=0
        yday=yearday(iyear, month, iday)
        DateNumber(1)=2440000.0_r8               ! Truncated offset
        DateNumber(2)=DateNumber(1)*86400_r8
      END IF
!
!-----------------------------------------------------------------------
!  Set reference-time string, YYYY-MM-DD hh:mm:ss
!-----------------------------------------------------------------------
!
      WRITE (string,10) iyear, month, iday, ihour, minute, isec
 10   FORMAT (i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)
!
!-----------------------------------------------------------------------
!  Load time reference clock information into structure.
!-----------------------------------------------------------------------
!
      Rclock%yday    =yday
      Rclock%year    =iyear
      Rclock%month   =month
      Rclock%day     =iday
      Rclock%hour    =ihour
      Rclock%minutes =minute
      Rclock%seconds =isec
      Rclock%base    =r_time
      Rclock%Dnumber =DateNumber(1)
      Rclock%Snumber =DateNumber(2)
      Rclock%string  =string
      Rclock%calendar=TRIM(calendar)
!
      RETURN
      END SUBROUTINE ref_clock
!
!**********************************************************************
      SUBROUTINE ROMS_clock (year, month, day, hour, minutes, seconds, &
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

      real(r8), intent(in)  :: seconds
      real(r8), intent(out) :: ClockTime
!
!  Local variable declarations.
!
      real(r8), dimension(2) :: DateNumber
!
!-----------------------------------------------------------------------
!  Compute ROMS clock ellapsed time since intialization in seconds from
!  reference time.
!-----------------------------------------------------------------------
!
!  Convert requested date into date number.
!
      CALL datenum (DateNumber, year, month, day,                      &
     &              hour, minutes, seconds)
!
!  Compute ROMS clock ellapsed time in seconds.
!
      ClockTime=DateNumber(2)-Rclock%Snumber
!
      RETURN
      END SUBROUTINE ROMS_clock
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
!  encoding allows an easy to read when reporting time.                !
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
      real(r8), intent(in) :: MyTime

      character (len=22), intent(out) :: date_string
!
!  Local variable declarations.
!
      integer :: day, hour, minutes, month, year
      integer :: i

      real(r8) :: Currenttime, seconds

      character (len= 5) :: sec_string
      character (len=22) :: string
!
!-----------------------------------------------------------------------
!  Encode current model time.
!-----------------------------------------------------------------------
!
!  Convert current model time to calendar date.
!
      CurrentTime=MyTime/86400.0_r8                  ! seconds to days
!
      CALL caldate (CurrentTime,                                        &
     &              yy_i=year,                                          &
     &              mm_i=month,                                         &
     &              dd_i=day,                                           &
     &              h_i=hour,                                           &
     &              m_i=minutes,                                        &
     &              s_r8=seconds)
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
