function outstr=uplow(inpstr)

% UPLOW: Capitalize first letter of a string.
%
% outstr = uplow(inpstr);
%
% Example:
%
% outstr = uplow('my CAPitalized strING')
%
% outstr = 'My Capitalized String'
%

% svn $Id: uplow.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

outstr = lower(inpstr);
ind=regexp([' ' outstr],'(?<=\s+)\S','start')-1;
outstr(ind)=upper(outstr(ind));

return
