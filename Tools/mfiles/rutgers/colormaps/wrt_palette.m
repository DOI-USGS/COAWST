function wrt_palette(cmap, pname, varargin)

% WRT_PALETTE:  Writes colormap palette to ASCII file.
%
% wrt_palette(cmap, pname, ncolors)
%
% Writes color palette for an ASCII file so it can be used in other
% software.
%
% On Input:
%
%    cmap     Mx3 colormap matrix
%    pname    ASCII color palette filename (string)
%    M        Number of colors (integer, OPTIONAL)
%
% Example:
%
%    wrt_palette(viridis, 'viridis_120.dat', 120)
%    wrt_palette(flipud(viridis), 'viridis_240_flip.dat', 240)
%

% svn $Id: wrt_palette.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    M = size(cmap,1);
  case 1
    M = varargin{1};
end

% Interpolate to requested number of colors.

P = size(cmap,1);

if (P ~= M)
  cmap = interp1(1:size(cmap,1), cmap, linspace(1,P,M), 'linear');
end

% Write out (R,G,B) color palette into a file.

fout = fopen(pname,'w');

if (fout < 0)
  error(['Cannot create ' pname '.'])
end

for n = 1:size(cmap,1)
  fprintf (fout, ' %22.16e  %22.16e  %22.16e\n',                        ...
           cmap(n,1), cmap(n,2), cmap(n,3));
end

return
