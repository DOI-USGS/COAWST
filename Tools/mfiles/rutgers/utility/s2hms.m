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

% svn $Id: s2hms.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                               Rich Signell        %
%===========================================================================%

%sec=round(secs);
%hour=floor(sec/3600);
%min=floor(rem(sec,3600)/60);
%sec=round(rem(sec,60));

sec = secs;
hour = fix(sec / 3600);
sec = sec - 3600*hour;
min = fix(sec / 60);
sec = sec - 60*min;

return
