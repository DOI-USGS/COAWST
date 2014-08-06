function [hours]=hms2h(h,m,s);

%
% HMS2H: converts hours, minutes, and seconds to decimal hours
%
% [hours]=hms2h(h,m,s);   or
% [hours]=hms2h(HHMMSS);
%
% This function converts (hours,minutes,seconds) or HHMMSS to
% decimal hours.
%
% On Input:
%
%    h           Hours or hour-minute-sec (HHMMSS)
%    m           Minutes (OPTIONAL)
%    s           Seconds (OPTIONAL)
%
% On Output:
%
%    hours       Decimal hours
%

% svn $Id: hms2h.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                               Rich Signell        %
%===========================================================================%

if nargin== 1,
  hms=h;
  h=floor(hms/10000);
  ms=hms-h*10000;
  m=floor(ms/100);
  s=ms-m*100;
  hours=h+m/60+s/3600;
else,
  hours=h+(m+s/60)/60;
end,

return
