function [gtime]=gregorian(julian);

%
% GREGORIAN:  Converts Julian day number to Gregorian calendar date vector
%
% [gtime]=gregorian(julian);
%
% This function convert Julian day number to Gregorian calendar
% date. Although the formal definition holds that Julian days
% start and end at noon, here Julian days start and end at midnight.
%
% In this convention, Julian day 2440000 began at 0000 hours,
% May 23, 1968.
%
% On Input:
%
%    j           Julian day number
%
% On Ouput:
%
%    gtime       Gregorian date. A six component vector:
%                  gtime=[yyyy month day hours minutes seconds]
%
%                  gtime(1) => year
%                  gtime(2) => month of the year
%                  gtime(1) => day of the month
%                  gtime(1) => hours
%                  gtime(1) => minutes
%                  gtime(1) => seconds
%
% Calls          s2hms
 
% svn $Id: gregorian.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                               Rich Signell        %
%===========================================================================%

julian=julian+5.e-9;    % kludge to prevent roundoff error on seconds

% If you want Julian Days to start at noon, use...
%
% h=rem(julian,1)*24+12;
% i=(h >= 24);
% julian(i)=julian(i)+1;
% h(i)=h(i)-24;

secs=rem(julian,1)*24*3600;

j = floor(julian) - 1721119;
in = 4*j -1;
y = floor(in/146097);
j = in - 146097*y;
in = floor(j/4);
in = 4*in +3;
j = floor(in/1461);
d = floor(((in - 1461*j) +4)/4);
in = 5*d -3;
m = floor(in/153);
d = floor(((in - 153*m) +5)/5);
y = y*100 +j;
mo=m-9;
yr=y+1;
i=(m<10);
mo(i)=m(i)+3;
yr(i)=y(i);
[hour,min,sec]=s2hms(secs);
gtime=[yr(:) mo(:) d(:) hour(:) min(:) sec(:)];

return
