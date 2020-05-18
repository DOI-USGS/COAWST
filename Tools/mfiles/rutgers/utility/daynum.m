function [d]=daynum(year, month, day, varargin)

% DAYNUM:  Calculates date number from date (ROMS version)
%
% [d]=daynum(year, month, day, hour, minute, second)
%
% Converts requested date (year, month, day, ...) into a serial date
% number. It uses the Proleptic Gregorian Calendar, which extends
% backward the date preceding 15 October 1582 with a year length of
% 365.2425 days. It is similar to native function "datenum" but differs
% on the origin, day 0:
%
%       datenum(0000,00,00) = 0
%        daynum(0000,03,01) = 0             March 1, 0000
%
% This function is used in ROMS to compute:
%
%     time-units since 0001-01-01 00:00:00   (Proleptic Gregorian Calendar)
%
% On Input:
%
%    year       Year including century (YYYY)
%    month      Month, a value ranging between 1-12
%    day        Day of the month
%    hour       Hour of the day: 1, ..., 23 (OPTIONAL)
%    minute     Minutes of the hour (OPTIONAL)
%    second     Seconds of the minute (OPTIONAL)
%
% Adapted from Gary Katch, Concordia University, Canada.
%
%    https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html

% svn $Id: daynum.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Set offset to get the same value as Matlat "datenum" function. Matlab
% origin is 0000-00-00 00:00:00 while for the equations below the origin
% is 0000-03-01 00:00:00.  The difference is 61 days.
  
%offset=0;
 offset=61;

if ((year == 0) && (month == 0) && (day == 0))
  isorigin=true;
else
  isorigin=false;
end

switch numel(varargin)
  case 0
    hour=0;
    minute=0;
    second=0;
  case 1
    hour=varargin{1};
    minute=0;
    second=0;
  case 2
    hour=varargin{1};
    minute=varargin{2};
    second=0;
  case 3    
    hour=varargin{1};
    minute=varargin{2};
    second=varargin{3};
end    

m = (mod(month+9, 12));         % Mar=0, ..., Feb=11

year = year - fix(m/10d0);      % if Jan or Feb, substract 1

g = 365*year + fix(year/4) - fix(year/100) + fix(year/400) +  ...
    fix((m*306 + 5)/10) + (day - 1);

if (isorigin)
  g = 0;
else
  if (g < 0)
    g = g + offset - 1;
  else
    g = g + offset;
  end
end

d = g + (hour/24d0) + (minute/1440d0) + (second/86400d0);

s = g*86400 + hour*360 + minute*60 + second;

return


