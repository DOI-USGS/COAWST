function [gtime]=gregorian(jd)
% GREGORIAN:  Converts Julian day numbers to Gregorian calendar.
%
% USAGE:      [gtime]=gregorian(jd) 
%
% DESCRIPTION:  Converts decimal Julian days to Gregorian dates using the
%               astronomical convension, but with time zero starting
%               at midnight instead of noon.  In this convention, 
%               Julian day 2440000 begins at 0000 hours, May 23, 1968.
%               The Julian day does not have to be an integer, and with
%               Matlab's double precision, the accuracy of decimal days
%               is about 0.1 milliseconds.  
%
%
% INPUT:   jd  = decimal Julian days 
%     
% OUTPUT:  gtime = six column Gregorian time matrix, where each row is
%                  [yyyy mo da hr mi sec]. 
%                   yyyy = year (e.g., 1979)
%                        mo = month (1-12)
%                           da = day (1-31)
%                              hr = hour (0-23)
%                                 mi = minute (0-59) 
%                                     sec = decimal seconds 
%                example: [1990 12 12 0 0 0] is midnight on Dec 12, 1990.
%
% AUTHOR: Rich Signell  (rsignell@usgs.gov)
%

% Add 0.2 milliseconds before Gregorian calculation to prevent
% roundoff error resulting from math operations on time 
% from occasionally representing midnight as 
% (for example) [1990 11 30 23 59 59.99...] instead of [1990 12 1 0 0 0]);
% If adding a 0.2 ms to time (each time you go back and forth between 
% Julian and Gregorian) bothers you more than the inconvenient representation
% of Gregorian time at midnight you can comment this line out...

      jd=jd+2.e-9;    
                            

%      if you want Julian Days to start at noon...
%      h=rem(jd,1)*24+12;
%      i=(h >= 24);
%      jd(i)=jd(i)+1;
%      h(i)=h(i)-24;

      secs=rem(jd,1)*24*3600;

      j = floor(jd) - 1721119;
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
