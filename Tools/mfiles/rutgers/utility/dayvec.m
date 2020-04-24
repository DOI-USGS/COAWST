function [V]=dayvec(dnum)

% DAYVEC: Calculates date from a date number (ROMS Version).
%
% [V]=dayvec(dnum)
%
% Converts a given date number as computed by "daynum" to a date vector
% (year, month, day, hour, minutes, seconds).  It is the inverse function
% to "daynum". It is equivalent the native "datevec" function but computed
% from different equations.
%                                     
% On Input:
%
%    dnum       Day numever computed from function "daynum"
%     
% On Output:
%
%    V          Date structure:
%
%                 V.daynum          input day number
%                 V.yday            day of the year
%                 V.year            year including century
%                 V.month           month of the year
%                 V.day             day of the month
%                 V.hour            hour of the day
%                 V.minutes         minutes of the hour
%                 V.seconds         second of the minute
%
% Adapted from Gary Katch, Concordia University, Canada.
%
%    https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html

% svn $Id: dayvec.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set offset to get the same value as Matlat "datenum" function. Matlab
% origin is 0000-00-00 00:00:00 while for the equations below the origin
% is 0000-03-01 00:00:00.  The difference is 61 days.
  
%offset=0;
 offset=61;

% Use variable precision arithmetic (vpa).

dfrac=vpa(dnum-fix(dnum));

% Substract offset to match Matlab "datestr" values with origin at
% 0000-00-00 00:00:00

if (dnum < offset)
  d=fix(dnum-offset+1);
else
  d=fix(dnum-offset);
end

year = fix((10000*d + 14780)/3652425);
ddd = d - (365*year + fix(year/4) - fix(year/100) + fix(year/400));
if (ddd < 0)
 year = year - 1;
 ddd = d - (365*year + fix(year/4) - fix(year/100) + fix(year/400));
end
mi = fix((100*ddd + 52)/3060);
month = fix(mod(mi + 2,12)) + 1;
year = year + fix((mi + 2)/12);
day = ddd - fix((mi*306 + 5)/10) + 1;

s=(dfrac*86400d0);
s=tround(s, 3*eps(double(s)));
Hour=fix((s/3600d0));
s=abs(s-Hour*3600d0);
Minutes=fix(s/60.0d0);
Seconds=abs(s-Minutes*60d0);

if (mod(year,4) == 0 && mod(year,100) ~= 0 || mod(year,400) == 0)
  fac=1d0;                                           % leap year
  leap='true';
else
  fac=2d0;
  leap='false';
end
yday = fix((275.0*month)/9) - fac*fix((month+9)/12) + day - 30;

if (dnum == 0)
  year=0;               % Fix to match Matlab "datestr" values
  month=1;              % with origin at 0000-00-00 00:00:00
  day=0;
end

mstr = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                       ...
        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

string = [num2str(day, '%2.2i'), '-',                                   ...
          char(mstr(month)), '-',                                       ...
          num2str(year, '%4.4i'), ' ',                                  ...
          num2str(double(Hour), '%2.2i'), ':',                          ...
          num2str(double(Minutes), '%2.2i'), ':',                       ...
          num2str(double(Seconds), '%5.2f')];
	  
V = struct('daynum',[], 'leap_year', [],'yday',[], 'year',[],           ...
	   'month',[], 'day', [], 'hour', [],                           ...
           'minutes',0, 'seconds',0, 'string', []);

V.daynum    = dnum;
V.leap_year = leap;
V.yday      = yday;
V.year      = year;
V.month     = month;
V.day       = day;
V.hour      = double(Hour);
V.minutes   = double(Minutes);
V.seconds   = double(Seconds);
V.string    = string;

return
