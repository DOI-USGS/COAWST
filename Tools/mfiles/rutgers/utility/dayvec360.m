function [V]=dayvec360(d360)

% DAYVEC360: Calculates date from a date number of the 360-day calendar.
%
% [V]=dayvec360(dnum)
%
% Convert 360_day calendar date number to date vector.  The
% perpetual 360_day calendar has a year length of 360 days and
% every month have 30 days.
%
% On Input:
%
%    d360       Day number as computed from function "daynum360"
%     
% On Output:
%
%    V          Date structure:
%
%                 V.day360          input day number
%                 V.yday            day of the year
%                 V.year            year including century
%                 V.month           month of the year
%                 V.day             day of the month
%                 V.hour            hour of the day
%                 V.minutes         minutes of the hour
%                 V.seconds         second of the minute
%
% It assumes that
%
%    d360=0      corresponds to 0000-01-01 00:00:00 (Jan  1, 0000)
%    d360=359    corresponds to 0000-12-30 00:00:00 (Dec 30, 0000)
%    d360=360    corresponds to 0001-01-01 00:00:00 (Jan  1, 0001)
%  

% svn $Id: dayvec360.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% Use variable precision arithmetic (vpa).

dfrac=vpa(d360-fix(d360));
d=fix(d360);

year=fix(d/360);
doy=fix((d-year*360) + 1);

month=fix(((doy-1)/30) + 1);
day=fix(mod(doy-1,30) + 1);

s=(dfrac*86400d0);
s=tround(s, 3*eps(double(s)));
Hour=fix((s/3600d0));
s=s-Hour*3600d0;
Minutes=fix(s/60d0);
Seconds=s-Minutes*60d0;

mstr = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                       ...
        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

string = [num2str(day, '%2.2i'), '-',                                   ...
          char(mstr(max(1,min(month,12)))), '-',                        ...
          num2str(year, '%4.4i'), ' ',                                  ...
          num2str(double(Hour), '%2.2i'), ':',                          ...
          num2str(double(Minutes), '%2.2i'), ':',                       ...
          num2str(double(Seconds), '%5.2f')];


V = struct('day360',[], 'yday',[], 'year',[], 'month',[],               ...
	   'day', [], 'hour', [], 'minutes',0, 'seconds',0,             ...
	   'string',[]);

V.day360  = d360;
V.yday    = doy;
V.year    = year;
V.month   = month;
V.day     = day;
V.hour    = double(Hour);
V.minutes = double(Minutes);
V.seconds = double(Seconds);
V.string  = string;

return
