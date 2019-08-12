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

% svn $Id: uplow.m 938 2019-01-28 06:35:10Z arango $
%===========================================================================%
%  Copyright (c) 2002-2019 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

outstr = lower(inpstr);
ind=regexp([' ' outstr],'(?<=\s+)\S','start')-1;
outstr(ind)=upper(outstr(ind));

return
