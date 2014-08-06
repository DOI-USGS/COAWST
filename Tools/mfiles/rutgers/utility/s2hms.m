function [hour,min,sec]=s2hms(secs);

%
% S2HMS:  Converts seconds to integer hour, minute, seconds
%
% [hour,min,sec]=s2hms(secs);
%
% This function convert decimal seconds to interer hour, minute, and
% seconds (rounded to the nearest integer).
%
% On Input:
%
%    secs        Decimal seconds
%
% On Output:
%
%    h           Hours (whole number)
%    m           Minutes (whole number)
%    s           Seconds (rounded to nearest integer)
%

% svn $Id: s2hms.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                               Rich Signell        %
%===========================================================================%

sec=round(secs);
hour=floor(sec/3600);
min=floor(rem(sec,3600)/60);
sec=round(rem(sec,60));

return
