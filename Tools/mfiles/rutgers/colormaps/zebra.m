function map = zebra(a,n,m);

% ZEBRA:  Banded color palette
%
% map = zebra(nbands, nentries, saturation)
%
% Zebra palette colormap with NBANDS broad bands and NENTRIES rows in
% the color map.
%  
% On Input:
%
%    a         Number of large bands in the palette (OPTIONAL)
%              Default is 4.
%
%    n         Number of entries in the colormap (OPTIONAL)
%              Default is size(get(gcf,'colormap'),1)
%
%    m         Saturation value: 0 > m =< 1.
%              Default is 0.5.
%
% On Output:
%
%    map       A Nx3 colormap matrix
%
% Example:
% 
%   colormap(zebra) 
%
% Reference:
%  
%   Hooker, S. B. et. al, 1995: Detecting Dipole Ring Separatrices with
%     Zebra Palettes, IEEE Transactions on Geosciences and Remote Sensing,
%     vol. 33, 1306-1312.
%

% svn $Id: zebra.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

if nargin < 3
  m = 0.5; % saturation and value go from m to 1
  % don't use m = 0
end
if nargin < 2
   n  = size(get(gcf,'colormap'),1);  % number of entries in the colormap
end
if nargin < 1
   a = 4; % there are this many large bands in the palette
end

x = 0:(n-1);

%hue = exp(-3*x/n);
%sat = m+(1-m)*(0.5*(1+sawtooth(2*pi*x/(n/a))));
%val = m+(1-m)*0.5*(1+cos(2*pi*x/(n/a/2)));

hue = exp(-2*x/n);
sat = m+(1-m)*(0.5*(1+sawtooth(2*pi*x/(n/a))));
val = m+(1-m)*0.5*(1+cos(2*pi*x/(n/a/2)));

map = [hue(:) sat(:) val(:)];
map = hsv2rgb(map);

return
