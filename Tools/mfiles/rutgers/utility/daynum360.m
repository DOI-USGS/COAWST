function [d]=daynum360(year, month, day, varargin)

% DAYNUM360:  Calculates date number from a 360-day calendar date.
%
% [d]=daynum360(year, month, day, hour, minute, second)
%
% Converts requested date (year, month, day, ...) into a serial date
% number using the 360-day calendar.  It has a year length of 360 days 
% and every month have 30 days.
%
%       datenum360(0000,01,01) = 0             origin: Jan 1, 0000
%        daynum360(0000,12,30) = 359
%
% This function is used in ROMS to compute:
%
%     time-units since 0000-12-30 00:00:00   (360-day Calendar)
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

% svn $Id: daynum360.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
  
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

d = fix(year*360) +                                                     ...
    fix((month-1)*30) +                                                 ...
    fix(day-1) +                                                        ...
    (hour/24d0) +                                                       ...
    (minute/1440d0) +                                                   ...
    (second/86400d0);

return
