function [code]=day_code(month, day, year);

%
% DAY_CODE:  Compute day of the week code
%
% [code]=day_code(month, day, year)
%
% This routine computes a code for the day of the week, given the
% date. This code is good for date after:
%
%                              January 1, 1752 AD
%
% the year the Gregorian calander was adopted in Britian and the
% American colonies.
%
% On Input:
%
%    day     The day of the month (integer).
%    month   The month, 1=January, 2=February, ... (integer).
%    year    The year, including the century (integer).
%
% On Output:
%
%    code    A code for the corresponding day of the week (numeric):
%                code = 1  =>  Sunday
%                code = 2  =>  Monday
%                code = 3  =>  Tuesday
%                code = 4  =>  Wednesday
%                code = 5  =>  Thursday
%                code = 6  =>  Friday
%                code = 7  =>  Saturday
%
% Calls:  none
%

% svn $Id: day_code.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
%  Initialize basis data.
%---------------------------------------------------------------------------

month_day = [31,28,31,30,31,30,31,31,30,31,30,31];

base_year = 1752;
feb_end = 59;
bym1_dec31 = 5;
base_cen = 1700;
base_qcen = 1600;
base_qyear = 1748;

%---------------------------------------------------------------------------
%  Make sure dealing with integers.
%---------------------------------------------------------------------------

imonth = round(month);
iday = round(day);
iyear = round(year);

%---------------------------------------------------------------------------
%  Compute the number of years since the base year, the number of years
%  since the beginning of the base century and the number of years since
%  the beginning of the base 400 year.
%---------------------------------------------------------------------------

no_yr = iyear - base_year;
nqy = iyear - base_qyear;
nyc = iyear - base_cen;
nyqc = iyear - base_qcen;

%---------------------------------------------------------------------------
% Compute the number of leapdays in that time.  Determine if this is a
% leap year.
%---------------------------------------------------------------------------

leap = fix(nqy/4) - fix(nyc/100) + fix(nyqc/400);

leap_flag = ((rem(nqy,4) == 0) & (rem(nyc,100) ~= 0)) | ...
            (rem (nyqc, 400) == 0);

%---------------------------------------------------------------------------
% Compute the number of days this year.  The leap year corrections
% are:
%      Jan. 1 - Feb. 28   Haven't had the leap day counted above.
%      Feb.29             Counting leap day twice.
%---------------------------------------------------------------------------

no_day = iday;

for n = 1:(imonth-1),
   no_day = no_day + month_day(n);
end

if ( leap_flag & (no_day <= feb_end) ),  no_day = no_day - 1; end;
if (leap_flag & (imonth==2) & (iday==29)), no_day = no_day - 1; end;

%---------------------------------------------------------------------------
%  Compute the total number of days since Jan. 1 of the base year, exclusive
%  of the 364 day per year which represent an even 52 weeks.  Actually,
%  only need to do the addition mod 7.
%---------------------------------------------------------------------------

no_day = rem(no_day,7) + rem(leap,7) + rem(no_yr,7) + bym1_dec31;

%---------------------------------------------------------------------------
% Get the day of the week code.
%---------------------------------------------------------------------------

code = rem (no_day, 7) + 1;

return
