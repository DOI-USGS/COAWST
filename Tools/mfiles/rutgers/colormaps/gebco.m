function cmap = gebco(varargin)

% GEBCO: 24 color palette for bathymetry/land charts
%
% cmap = gebco(M)
%
% Bathymetry amd elevation colormap.
%
% On Input:
%
%    M        Number of colors (integer, OPTIONAL)
%
% On Ouput:
%
%    cmap     Mx3 colormap matrix
%
% Usage:
%
%    colormap(gebco)
%    colormap(flipud(gebco))
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                   %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.md                            Hernan G. Arango        %
%===========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    M = 24;
  case 1
    M = varargin{1};
end

% Set 24 colormap: 1;13 water, 14:24 land.

cmap = [[0.07059, 0.03922, 0.23137],
        [0.09020, 0.19216, 0.43529],
        [0.07843, 0.35294, 0.54902],
        [0.10588, 0.40784, 0.64314],
        [0.11765, 0.44706, 0.70196],
        [0.11373, 0.54510, 0.76863],
        [0.10588, 0.64706, 0.82353],
        [0.10980, 0.72157, 0.87843],
        [0.10588, 0.80000, 0.92549],
        [0.10588, 0.84706, 0.94510],
        [0.14902, 0.87451, 0.94510],
        [0.19216, 0.90196, 0.92549],
        [0.41176, 0.94902, 0.91373],
        [0.63137, 1.00000, 0.90196],
        [0.76471, 0.81961, 0.31373],
        [0.88627, 0.88235, 0.40000],
        [0.87451, 0.76863, 0.36078],
        [0.82745, 0.69804, 0.31765],
        [0.74118, 0.58824, 0.18431],
        [0.63922, 0.49804, 0.18431],
        [0.60000, 0.46275, 0.16863],
        [0.56078, 0.43137, 0.15294],
        [0.52941, 0.40784, 0.14118],
        [0.45882, 0.34510, 0.11373]];

% Interpolate to requested number of colors.

P = size(cmap,1);

if (P ~= M)
  cmap = interp1(1:size(cmap,1), cmap, linspace(1,P,M), 'linear');
end

return

