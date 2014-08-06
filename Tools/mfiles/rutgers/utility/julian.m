function [j]=julian(y,m,d,h)

%
% JULIAN:  Converts Gregorian calendar date to Julian day numbers
%
% [j]=julian(y,m,d,h);   or
% [j]=julian([y m d hour min sec])
%
% This function convert Gregorian calendar dates to corresponding
% Julian day numbers.  Although the formal definition holds that
% Julian days start and end at noon, here Julian days start and
% end at midnight.
%
% In this convention, Julian day 2440000 began at 0000 hours,
% May 23, 1968.
%
% On Input:
%
%    y           Year (scalar, e.g. 1979) or six component vector
%                   computed by the Gregorian function
%    m           Month of the year (1-12)
%    d           Day of the month (1-31)
%    h           Decimal hours (assumed 0 if absent)
%
% On Output:
%
%    j           Julian day number
%
% Calls          hms2h
%

% svn $Id: julian.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                               Rich Signell        %
%===========================================================================%

if nargin==3,
  h=0.;
elseif nargin==1,
  h=hms2h(y(:,4),y(:,5),y(:,6));
  d=y(:,3);
  m=y(:,2);
  y=y(:,1);
end,
mo=m+9;
yr=y-1;
i=(m>2);
mo(i)=m(i)-3;
yr(i)=y(i); 
c = floor(yr/100);
yr = yr - c*100;
j = floor((146097*c)/4) + floor((1461*yr)/4) + ...
floor((153*mo +2)/5) +d +1721119;

% If you want julian days to start and end at noon, replace the following
% line with:
%
% j=j+(h-12)/24;
 
j=j+h/24;

return
