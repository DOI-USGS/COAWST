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

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                   %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.md                            Hernan G. Arango        %
%===========================================================================%

outstr = lower(inpstr);
ind=regexp([' ' outstr],'(?<=\s+)\S','start')-1;
outstr(ind)=upper(outstr(ind));

return
