function stamp = date_stamp;

%
% DATE_STAMP:  Set current date string
%
% stamp = date_stamp
%
% This function returns a current date stamp.
%
% On Output:
%
%    stamp      date/time stamp of the form:
%                 day of week - month day, year - time
%

% svn $Id: date_stamp.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
%  Define string arrays for day of the week and month.
%---------------------------------------------------------------------------

wkday = ['Sunday   ';'Monday   ';'Tuesday  ';'Wednesday';'Thursday ';...
         'Friday   ';'Saturday '];
lwkday = [6,6,7,9,8,6,8];

month = ['January  ';'February ';'March    ';'April    ';'May      ';...
         'June     ';'July     ';'August   ';'September';'October  ';...
         'November ';'December '];
lmonth = [7,8,5,5,3,4,4,6,9,7,8,8];

ampm = [' AM';' PM'];

%---------------------------------------------------------------------------
%  Get date and time information.
%---------------------------------------------------------------------------

date_info = clock;

year = date_info(1);
imonth = date_info(2);
day = date_info(3);
hour = date_info(4);
minute = date_info(5);
second = date_info(6);

%---------------------------------------------------------------------------
%  Reset hour to 12-hour clock.
%---------------------------------------------------------------------------

half = rem( fix(hour/12),2) + 1;

hour = rem( hour, 12);
if (hour==0), hour = 12; end;

%---------------------------------------------------------------------------
%  Get day of the week code.
%---------------------------------------------------------------------------

iwday = day_code(imonth,day,year);

%---------------------------------------------------------------------------
%  Set date stamp.
%---------------------------------------------------------------------------

stamp = [wkday(iwday,1:lwkday(iwday)),' - ',month(imonth,1:lmonth(imonth)),...
         ' ',num2str(day),', ',num2str(year),' - ',num2str(hour),':',...
         num2str(minute),':',num2str(second),ampm(half,:)];

return
