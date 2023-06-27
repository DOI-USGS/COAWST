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

% svn $Id: uplow.m 1156 2023-02-18 01:44:37Z arango $
%===========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

outstr = lower(inpstr);
ind=regexp([' ' outstr],'(?<=\s+)\S','start')-1;
outstr(ind)=upper(outstr(ind));

return
