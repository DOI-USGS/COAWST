function [doy]=yearday(year, month, day)

% YEARDAY:  Computes the day of the year
%
% [doy]=yearday(year, month, day)
%
% Given the year, month, and day, this function computes the
% day of the year.
%
% On Input:
%
%    year       Year including century (YYYY)
%    month      Month, a value ranging between 1-12.
%    day        Day of the month
%
% Check value:
%
%       Mar 1, 2016:     yeaday(2016,3,1) = 61     (leap year)
%       Mar 1, 2017:     yeaday(2017,3,1) = 60
%
%       Dec 31, 2016:    yeaday(2016,12,31) = 366  (leap year)
%       Dec 31, 2017:    yeaday(2017,12,31) = 365  
%
% Adapted from Gary Katch, Concordia University, Canada.
%
%    https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html

% svn $Id: yearday.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Make sure the input data are whole number.

year  = fix(year);
month = fix(month);
day   = fix(day);

% Determine factor for leap year.

if ((mod(year,4) == 0 && mod(year,100) ~= 0) || mod(year,400) == 0)
  fac = 1;                                                      % leap year
else
  fac = 2;
end

% Compute the day od the year.

doy = fix((275*month)/9) - fac*fix((month+9)/12) + day - 30;

return
