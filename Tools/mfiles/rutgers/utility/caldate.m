function [Date]=caldate(Julian);

%
% CALDATE:  Converts Julian day numbers to calendar date
%
% [Date]=caldate(Julian);
%
% This function converts Julian day number to calendar date structure.
%
% On Input:
%
%    Julian      Julian day number
%
% On Output:
%
%    Date        Calendar (Gregorian) date (structure array):
%                  Date.year  =>  year
%                  Date.yday  =>  year day
%                  Date.month =>  month
%                  Date.day   =>  day of the month
%                  Date.hour  =>  decimal hours
%                  Date.min   =>  minutes
%                  Date.sec   =>  seconds
%

% svn $Id: caldate.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Check fom modified Julian day.  Modified Julian day beggins at 
% Julian Day 2440000 (0000 hours, May 23, 1968).

JulDay=Julian;
if (length(JulDay) == 1),
 if (JulDay < 2440000 ),
   JulDay=JulDay+2440000;
 end,
else,
  for i=1:length(JulDay),
    if (JulDay(i) < 2440000 ),
      JulDay(i)=JulDay(i)+2440000;
    end,
  end,
end,

% Kludge to prevent roundoff error on seconds.

JulDay=JulDay+5.0e-9;

% Convert to gregorain date.

j=floor(JulDay)-1721119;
in=4*j-1;
y=floor(in/146097);
j=in-146097*y;
in=floor(j/4);
in=4*in+3;
j=floor(in/1461);
d=floor(((in-1461*j)+4)/4);
in=5*d-3;
m=floor(in/153);
d=floor(((in-153*m)+5)/5);
y=y*100+j;
mo=m-9;
yr=y+1;
i=(m<10);
mo(i)=m(i)+3;
yr(i)=y(i);

% Get year day.

iyd =[1 32 60 91 121 152 182 213 244 274 305 335 366];
iydl=[1 32 61 92 122 153 183 214 245 275 306 336 367];
nrec=length(JulDay);
for n=1:nrec,
  if ( rem(yr(n),4) == 0 ),
    yday(n)=iydl(mo(n))+d(n)-1;
  else
    yday(n)=iyd(mo(n))+d(n)-1;
  end,
end,

% Extract hour, minutes and seconds.

secs=rem(JulDay,1)*24*3600;
sec=round(secs);
hour=floor(sec/3600);
min=floor(rem(sec,3600)/60);
sec=round(rem(sec,60));

% Build structure array.

if (length(yr) == 1),
  Date.year =yr;
  Date.yday =yday;
  Date.month=mo;
  Date.day  =d;
  Date.hour =hour; 
  Date.min  =min;
  Date.sec  =sec;
else,
  Date.year =yr(:);
  Date.yday =yday(:);
  Date.month=mo(:);
  Date.day  =d(:);
  Date.hour =hour(:); 
  Date.min  =min(:);
  Date.sec  =sec(:);
end,

return
