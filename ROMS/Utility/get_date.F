#include "cppdefs.h"
      SUBROUTINE get_date (date_str)
!
!svn $Id: get_date.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   This routine gets todays date, day of the week and time called     !
!  (default month & weekday are December & Saturday respectively).     !
!  It uses SUN intrinsic date routine by default.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     date_str   Concatenated string for the day of the week, date     !
!                (month,day,year), and time (12hr clock) of day        !
!                (hour:min:sec).                                       !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
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
      SUBROUTINE day_code (month, day, year, code)
!
!=======================================================================
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
!=======================================================================
!
      USE mod_kinds
!
      implicit none
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

      SUBROUTINE caldate (r_date, time, year, yday, month, day, hour)
!
!=======================================================================
!                                                                      !
!  This routine converts Julian day to calendar date.                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     r_date     Model initialization reference date (vector):         !
!                  r_date(1) => reference date (yyyymmdd.f).           !
!                  r_date(2) => year.                                  !
!                  r_date(3) => year day.                              !
!                  r_date(4) => month.                                 !
!                  r_date(5) => day.                                   !
!                  r_date(6) => hour.                                  !
!                  r_date(7) => minute.                                !
!                  r_date(8) => second.                                !
!     time       Model day (real; days).                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     year      The year, including the century (integer).             !
!     yday      The year day, 1 - 366 (real).                          !
!     month     The month, 1=January, 2=February, ... (integer).       !
!     day       The day of the month (integer).                        !
!     hour      The hour, 1 - 24 (real).                               !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(out) :: year, month, day

      real(r8), dimension(8), intent(in) :: r_date
      real(r8), intent(in) :: time

      real(r8), intent(out) :: yday, hour
!
!  Local variable declarations.
!
      integer, parameter :: gregorian = 2299161

      integer :: leap, ja, jalpha, jb, jc, jd, jday, je

      integer, dimension(13) :: iyd =                                   &
     &         (/ 1,32,60,91,121,152,182,213,244,274,305,335,366 /)

      integer, dimension(13) :: iydl =                                  &
     &         (/ 1,32,61,92,122,153,183,214,245,275,306,336,367 /)

      real(r8) :: rday
!
!-----------------------------------------------------------------------
!  Get calendar day from model time.
!-----------------------------------------------------------------------
!
!  The reference time is a positive date specified at initialization.
!
      IF (INT(r_date(1)).gt.0) THEN
        rday=time+r_date(3)+                                            &
     &       (r_date(6)+(r_date(7)+r_date(8)/60.0_r8)/60.0_r8)/24.0_r8
        year=INT(r_date(2))+INT(rday/365.25_r8)
        yday=MAX(1.0_r8,MOD(rday,365.25_r8))
        leap=MOD(year,4)
        IF (leap.eq.0) THEN
          month=MIN(12,1+INT((yday+1.0_r8)/30.6001_r8))
          day=INT(yday)-iyd(month)+1
        ELSE
          month=MIN(12,1+INT(yday/30.6001_r8))
          IF ((yday.gt.59).and.(yday.le.90)) month=3
          day=INT(yday)-iydl(month)+1
        END IF
        hour=(rday-AINT(rday))*24.0_r8
!
!  The reference time is for a climatological simulation with 365.25
!  days in every year.
!
      ELSE IF (INT(r_date(1)).eq.0) THEN
        rday=time+r_date(3)+                                            &
     &       (r_date(6)+(r_date(7)+r_date(8)/60.0_r8)/60.0_r8)/24.0_r8
        year=INT(r_date(2))+INT(rday/365.25_r8)
        yday=MAX(1.0_r8,MOD(rday,365.25_r8))
        leap=MOD(year,4)
        IF (leap.eq.0) THEN
          month=MIN(12,1+INT((yday+1.0_r8)/30.6001_r8))
          day=INT(yday)-iyd(month)+1
        ELSE
          month=MIN(12,1+INT(yday/30.6001_r8))
          IF ((yday.gt.59).and.(yday.le.90)) month=3
          day=INT(yday)-iydl(month)+1
        END IF
        hour=(rday-AINT(rday))*24.0_r8
!
!  The reference time is for a climatological simulation with 360.0
!  days in every year.
!
      ELSE IF (INT(r_date(1)).eq.-1) THEN
        rday=time+r_date(3)+                                            &
     &       (r_date(6)+(r_date(7)+r_date(8)/60.0_r8)/60.0_r8)/24.0_r8
        year=INT(r_date(2))+INT(rday/360.0_r8)
        yday=MAX(1.0_r8,MOD(rday,360.0_r8))
        month=INT(yday/30.0_r8)
        day=INT(MOD(yday,30.0_r8))
        hour=(rday-AINT(rday))*24.0_r8
!
!  The reference time is in terms of modified Julian days.
!
      ELSE
        IF (time.ge.2440000.0_r8) THEN
          jday=INT(time)
        ELSE
          jday=INT(time)+2440000
        END IF
        hour=(time-AINT(time))*24.0_r8
        IF (jday.ge.gregorian) THEN
          jalpha=INT(((jday-1867216)-0.25_r8)/36524.25_r8)
          ja=jday+1+jalpha-INT(0.25_r8*REAL(jalpha,r8))
        ELSE
          ja=jday
        END IF
        jb=ja+1524
        jc=INT(6680.0_r8+(REAL(jb-2439870,r8)-122.1_r8)/365.25_r8)
        jd=365*jc+INT(0.25_r8*REAL(jc,r8))
        je=INT(REAL(jb-jd,r8)/30.6001_r8)
        day=jb-jd-INT(30.6001_r8*REAL(je,r8))
        month=je-1
        IF (month.gt.12) month=month-12
        year=jc-4715
        IF (month.gt.2) year=year-1
        IF (year.le.0) year=year-1
        leap=MOD(year,4)
        IF (leap.eq.0) THEN
          yday=REAL(iydl(month),r8)+REAL(day,r8)-1.0_r8
        ELSE
          yday=REAL(iyd(month),r8)+REAL(day,r8)-1.0_r8
        END IF
        yday=yday+(time-AINT(time))
      END IF
      RETURN
      END SUBROUTINE caldate

      FUNCTION ref_att (r_time, r_date)
!
!=======================================================================
!                                                                      !
!  This function encodes the relative time attribute that gives the    !
!  elapsed interval since a specified reference time.  The "units"     !
!  attribute takes the form "time-unit since reference-time".          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     r_time     Time-reference (real; %Y%m%d.%f, for example,         !
!                  20020115.5 for 15 Jan 2002, 12:0:0).                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ref_att    Time-reference for "units" attribute (string).        !
!     r_date     Calendar date vector (real):                          !
!                  r_date(1) => reference date (yyyymmdd.f).           !
!                  r_date(2) => year.                                  !
!                  r_date(3) => year day.                              !
!                  r_date(4) => month.                                 !
!                  r_date(5) => day.                                   !
!                  r_date(6) => hour.                                  !
!                  r_date(7) => minute.                                !
!                  r_date(8) => second.                                !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: r_time

      real(r8), dimension(8), intent(out) :: r_date

      character (len=19) :: ref_att
!
!  Local variable declarations.
!
      integer :: iday, ihour, isec, iyear, leap, minute, month

      integer, dimension(13) :: iyd =                                   &
     &         (/ 1,32,60,91,121,152,182,213,244,274,305,335,366 /)

      integer, dimension(13) :: iydl =                                  &
     &         (/ 1,32,61,92,122,153,183,214,245,275,306,336,367 /)

      real(r8) :: day, sec, yday

      character (len=19) :: text
!
!-----------------------------------------------------------------------
!  Decode reference time.
!-----------------------------------------------------------------------
!
      iyear=MAX(1,INT(r_time*0.0001_r8))
      month=MIN(12,MAX(1,INT((r_time-REAL(iyear*10000,r8))*0.01_r8)))
      day=r_time-AINT(r_time*0.01_r8)*100.0_r8
      iday=INT(day)
      sec=(day-AINT(day))*86400.0_r8
      ihour=INT(sec/3600.0_r8)
      minute=INT(MOD(sec,3600.0_r8)/60.0_r8)
      isec=INT(MOD(sec,60.0_r8))
!
!-----------------------------------------------------------------------
!  Get year day.
!-----------------------------------------------------------------------
!
      leap=MOD(iyear,4)
      IF (leap.eq.0) THEN
        yday=REAL(iydl(month),r8)+REAL(iday,r8)-1.0_r8
      ELSE
        yday=REAL(iyd(month),r8)+REAL(iday,r8)-1.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  Build output date vector.
!-----------------------------------------------------------------------
!
      r_date(1)=r_time
      r_date(2)=REAL(iyear,r8)
      r_date(3)=MAX(1.0_r8,yday)
      r_date(4)=REAL(month,r8)
      r_date(5)=MAX(1.0_r8,REAL(iday,r8))
      r_date(6)=REAL(ihour,r8)
      r_date(7)=REAL(minute,r8)
      r_date(8)=REAL(isec,r8)
!
!-----------------------------------------------------------------------
!  Build reference-time string.
!-----------------------------------------------------------------------
!
      WRITE (text,10) iyear, month, iday, ihour, minute, isec
 10   FORMAT (i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)
      ref_att=text
      RETURN
      END FUNCTION ref_att

!----------------------------------------------------------------------
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!---DATE (YEAR,MONTH,DAY).
!---------------------------------------------------------------------
      INTEGER FUNCTION jd (year,month,day)
!
      INTEGER year,month,day,i,j,k
!
      i= year
      j= month
      k= day
!
      jd= k-32075+1461*(i+4800+(j-14)/12)/4+367*(j-2-(j-14)/12*12)      &
     &    /12-3*((I+4900+(j-14)/12)/100)/4
!
      RETURN
      END

      SUBROUTINE time_string (time, time_code)
!
!=======================================================================
!                                                                      !
!  This routine encodes current model time in seconds to a time        !
!  string of the form:                                                 !
!                                                                      !
!       DDDDD HH:MM:SS                                                 !
!                                                                      !
!  This allow a more accurate label when reporting time.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     time       Current model time  (seconds)                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     time_code  Current model time string.                            !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: time

      character(len=14), intent(out) :: time_code
!
!  Local variable declarations.
!
      integer :: iday, ihour, isec, minute

      character (len=14) :: text
!
!-----------------------------------------------------------------------
!  Encode current mode time.
!-----------------------------------------------------------------------
!
      iday=INT(time/86400)
      isec=INT(time-iday*86400.0_r8)
      ihour=isec/3600
      minute=MOD(isec,3600)/60
      isec=MOD(isec,60)
!
      WRITE (text,10) iday, ihour, minute, isec
 10   FORMAT (i5,1x,i2.2,':',i2.2,':',i2.2)
      time_code=text

      RETURN
      END SUBROUTINE time_string
