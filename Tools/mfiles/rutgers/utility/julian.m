function [Jday]=julian(varargin)

%
% JULIAN:  Converts Gregorian calendar date vector to Julian day numbers
%
% [Jday] = julian([YYY MM DD])
% [Jday] = julian([YYY MM DD hh.hh])
% [Jday] = julian([YYY MM DD hh mm ss])
% [Jday] = julian([YYY MM DD], noon_start)
% [Jday] = julian([YYY MM DD hh.hh], noon_start)
% [Jday] = julian([YYY MM DD hh mm ss], noon_start)
%
% [Jday] = julian(YYYY, MM, DD)
% [Jday] = julian(YYYY, MM, DD, hh.hh)
% [Jday] = julian(YYYY, MM, DD, hh.hh, noon_start)
%
% This function convert Proleptic Gregorian calendar dates vector to
% corresponding serial Julian day numbers.
%  
% Although the formal definition holds that  Julian days start and end at
% noon, here by default Julian days start and end at midnight.
%
% On Input:
%
%    V           Date vector 
%                  [YYYY MM DD]
%                  [YYYY MM DD hh.hh]
%                  [YYYY MM DD hh mm ss]
%    YYYY        Year of the century
%    MM          Month of the year (1-12)
%    DD          Day of the month (1-31)
%    hh          Hour of the day (0-23)
%    hh.hh       Hour of the day plus hour fraction
%    mm          Minutes of the hour (0-59)
%    ss          Seconds of the minute (0-59)
%
%    noon_start  Switch to start and end at noon (logical, OPTIONAL)
%                  noon_start = false   (default)
%
% On Output:
%
%    Jday        Serial Julian day number
%
% Example:
%
%    Jday = julian([-4713 11 24 12 0 0; 1968 5 23 12 0 0], true)
%
%    Jday = [0, 2440000]    Jday = 0         origin:     Nov 24, 4713 BC
%                           Jday = 2440000   truncated:  May 23, 1968
%
% Calls:         hms2h
%

% svn $Id: julian.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license               Hernan G. arango  %
%    See License_ROMS.txt                               Rich Signell      %
%=========================================================================%

noon_start = false;

switch numel(varargin)
  case 1
    V = varargin{1};
    switch (length(V))
      case 3
        YY = V(:,1);
        MM = V(:,2);
        DD = V(:,3);
        hh = 0;
      case 4
        YY = V(:,1);
        MM = V(:,2);
        DD = V(:,3);
        hh = V(:,4);
      case 6
        YY = V(:,1);
        MM = V(:,2);
        DD = V(:,3);
        hh = hms2h(V(:,4), V(:,5), V(:,6));
    end
  case 2
    V = varargin{1};
    noon_start = varargin{2}; 
    switch (length(V))
      case 3
        YY = V(:,1);
        MM = V(:,2);
        DD = V(:,3);
        hh = 0;
      case 4
        YY = V(:,1);
        MM = V(:,2);
        DD = V(:,3);
        hh = V(:,4);
      case 6
        YY = V(:,1);
        MM = V(:,2);
        DD = V(:,3);
        hh = hms2h(V(:,4), V(:,5), V(:,6));
    end
  case 3
    YY = varargin{1}; 
    MM = varargin{2};
    DD = varargin{3};
    hh = 0;
  case 4
    YY = varargin{1}; 
    MM = varargin{2};
    DD = varargin{3};
    hh = varargin{4};
  case 5
    YY = varargin{1}; 
    MM = varargin{2};
    DD = varargin{3};
    hh = varargin{4};
    noon_start = varargin{5}; 
end

%  Compute serial Julian day number.

if (MM > 2)
  yr = YY; 
  mo = MM-3;
else
  yr = YY-1;
  mo = MM+9;
end
  
c    = floor(yr/100);
yr   = yr - c*100;
Jday = floor((146097*c)/4) + floor((1461*yr)/4) +                       ...
       floor((153*mo +2)/5) + DD + 1721119;

% If requested set Julian date number to start and end at noon.

if (noon_start)
  Jday = Jday + (hh-12)/24;
else
  Jday = Jday + hh/24;
end

return
