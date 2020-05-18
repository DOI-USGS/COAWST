function [V] = gregorian(Jday, varargin)

%
% GREGORIAN: Converts Julian day number to Proleptic Gregorian date vector
%
% [V] = gregorian(Jday)
% [V] = gregorian(Jday, noon_start)
%
% This function convert serial Julian day number to Proleptic Gregorian
% calendar date.
%  
% Although the formal definition holds that Julian days start and end at
% noon, here by default Julian days start and end at midnight.  The Julian
% is primarily used in astronomy.
%
% On Input:
%
%    Jday         Julian day number plus fraction (scalar or vector)
%                   Jday = 0         24 Nov, 4713 BC  (Proleptic Gregorian)
%                                    01 Jan, 4713 BC  (Proleptic Julian) 
%
%    noon_start   Switch to start and end at noon (logical, OPTIONAL)
%                   noon_start = false   (default)
%
% On Ouput:
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
 
% svn $Id: gregorian.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license               Hernan G. Arango  %
%    See License_ROMS.txt                               Rich Signell      %
%=========================================================================%

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

jd = floor(Jnum) - 1721119;
in = 4*jd - 1;
y  = floor(in/146097);
jd = in - 146097*y;
in = floor(jd/4);
in = 4*in + 3;
jd = floor(in/1461);
d  = floor(((in - 1461*jd) + 4)/4);
in = 5*d - 3;
m  = floor(in/153);
y  = y*100 + jd;

if (m < 10)
 year  = y;
 month = m+3;
else
 year  = y+1;
 month = m-9;
end 
day = fix(((in - 153*m) + 5)/5);

s = dfrac*86400d0;
s = tround(s, 3*eps(double(s)));
Hour = fix((s/3600d0));
s = abs(s-Hour*3600d0);
Minutes = fix(s/60.0d0);
Seconds = abs(s-Minutes*60d0);

if (mod(year,4) == 0 & mod(year,100) ~= 0 | mod(year,400) == 0)
  fac  = 1d0;                                           % leap year
  leap = 'true';
else
  fac  = 2d0;
  leap = 'false';
end
yday = fix((275.0*month)/9) - fac*fix((month+9)/12) + day - 30;

mstr = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                       ...
        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

if (year < 0)
  ystr = strcat(num2str(abs(year(1)), '%4.4i'), ' BC');
else
  ystr = num2str(abs(year(1)), '%4.4i');
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
