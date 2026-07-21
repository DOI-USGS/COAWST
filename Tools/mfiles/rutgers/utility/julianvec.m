function [V]=julianvec(Jday, varargin)

% JULIANVEC: Calculates date from a serial Julian date number
%
% [V] = julianvec(Jday)
% [V] = julianvec(Jday, noon_start)
%
% Converts a given Julian date number as computed by "julian" to a date
% vector (year, month, day, hour, minutes, seconds).  It is the inverse
% function to "julian" as implemented in ROMS.
%
% Although the formal definition of Julian day numbers starts and
% ends at noon, here Julian day starts and ends at midnight. So it
% is 12 hour faster (substract 12 hour to agree with fomal definition).
%
% On Input:
%
%    Jday         Julian day number plus fraction (scalar or vector)
%                   Jday = 0         01 Jan, 4713 BC  (Proleptic Julian) 
%                                    24 Nov, 4713 BC  (Proleptic Gregorian)
%
%    noon_start   Switch to start and end at noon (logical, OPTIONAL)
%                   noon_start = false   (default)
%
% On Output:
%
%    V            Date structure:
%
%                   V.Jday           input day number
%                   V.yday           day of the year
%                   V.year           year including century
%                   V.month          month of the year
%                   V.day            day of the month
%                   V.hour           hour of the day
%                   V.minutes        minutes of the hour
%                   V.seconds        second of the minute
%
% Notice that a calendar obtained by extending backward in time from
% its invention or implementation is called the Proleptic version of
% the calendar. For example, the Proleptic Gregorian Calendar extends
% backwards the date preceeding 15 October 1582 with a year length of
% 365.2425 days instead of 365.25 in the Julian Calendar.
%
% Adapted from Numerical Recipes.
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

GregorianStart=2299161;                            % 15 Oct, 1582 A.D.

% If requested set Julian Days to start and end at noon.

noon_start = false;

switch numel(varargin)
  case 1
    noon_start = varargin{1};
end

if (noon_start)
  Jnum = Jday + 0.5;
else
  Jnum = Jday;
end

% Use variable precision arithmetic (vpa).

dfrac=vpa(Jnum-fix(Jnum));

Jnum=fix(Jnum);

if (Jnum >= GregorianStart)
  jalpha=fix(((Jnum-1867216)-0.25)/36524.25);      % Gregorian
  ja=Jnum+1+jalpha-fix(0.25*jalpha);               % correction
else
  ja=Jnum;
end

jb=ja+1524;
jc=fix(6680+((jb-2439870)-122.1)/365.25);
jd=365*jc+fix(0.25*jc);
je=fix((jb-jd)/30.6001);

day=jb-jd-fix(30.6001*je);

month=je-1;
ind=month > 2;
if (~isempty(ind))
  month(ind)=month(ind)-12;
end

year=jc-4715;
ind=month > 2;
if (~isempty(ind))
  year(ind)=year(ind)-1;
end
if (year <= 0)
  year=year-1;
end

s=(dfrac*86400d0);
s=tround(s, 3*eps(double(s)));
Hour=fix((s/3600d0));
s=abs(s-Hour*3600d0);
Minutes=fix(s/60.0d0);
Seconds=abs(s-Minutes*60d0);

if (mod(year,4) == 0 & mod(year,100) ~= 0 | mod(year,400) == 0)
  fac=1d0;                                           % leap year
  leap='true';
else
  fac=2d0;
  leap='false';
end
yday = fix((275.0*month)/9) - fac*fix((month+9)/12) + day - 30;

mstr = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                       ...
        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

month

if (year < 0)
  ystr=strcat(num2str(abs(year(1)), '%4.4i'), ' BC');
else
  ystr=num2str(abs(year(1)), '%4.4i');
end

string = [num2str(day(1), '%2.2i'), '-',                                ...
          char(mstr(month(1))), '-',                                    ...
          ystr, ' ',                                                    ...
          num2str(double(Hour(1)), '%2.2i'), ':',                       ...
          num2str(double(Minutes(1)), '%2.2i'), ':',                    ...
          num2str(double(Seconds(1)), '%5.2f')];

V = struct('Jday',[], 'leap_year', [],'yday',[], 'year',[],             ...
           'month',[], 'day', [], 'hour', [],                           ...
           'minutes',0, 'seconds',0, 'string', []);

V.Jday      = Jday;
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
