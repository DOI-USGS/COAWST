function S = sponge (Gname, varargin)

%
% SPONGE: Set diffusion and viscosity sponge coefficients
%
% S = sponge (Gname, factor, Nfilter, Lplot, Lwrite)
%
% Given a grid application NetCDF file, this function computes enhanced
% viscosity and diffusion scaling variables (visc_factor and diff_factor)
% that can be added input ROMS Grid NetCDF file. These scales are used
% in an application to set sponge areas with larger horizontal mixing
% coefficients for the damping of high frequency noise coming from open
% boundary conditions or nesting. In ROMS, these scales are used as
% follows:
%
%   visc2_r(i,j) = visc2_r(i,j) * visc_factor(i,j)
%   visc4_r(i,j) = visc4_r(i,j) * visc_factor(i,j)
%
%   diff2(i,j,itrc) = diff2(i,j,itrc) * diff_factor(i,j)
%   diff4(i,j,itrc) = diff4(i,j,itrc) * diff_factor(i,j)
%
% where visc_factor and diff_factor are defined at RHO-points. Usually,
% sponges are linearly tapered over several grid points adjacent to the
% open boundaries. Its positive values linearly increases from the inner
% to outer edges of the sponge. At the interior of the grid we can values
% of zero (no mixing) or one (regular mixing). 
%
% On Input:
%
%    Gname        ROMS Grid NetCDF file name (character string)
%
%    factor       Sponge enhancement factor (default: 3.0)
%
%    Nfilter      Number of passes of 3-point 2-D boxcar filter. It
%                   determines sponge width (default: 25)
%
%    Lplot        Switch to plot sponge variables (default: false)
%
%    Lwrite       Switch to add/write sponge variables to GRID NetCDF
%                   file (default: false)
%
% On Output:
%
%    S           Nested grids Contact Points structure (struct array)
%
%                  S.lon_rho               longitude of RHO-points
%                  S.lat_rho               latitude  of RHO-points
%                  S.mask_rho              land/sea mask at RHO-points
%                  S.diff_factor           diffusivity factor
%                  S.visc_factor           viscosity factor
%
% Methodology:
%  
% The generic method (John Wilkin; 2014) is to start out by setting
% perimeter values of a 2-D RHO-points array equal to 1 along open
% boundary "wet" points of the RHO-points mask and have zeroes elsewhere.
% Then, apply several passes of a 3x3 convolution operator to diffuse
% the non-zero perimeter values into the ocean interior. This means that
% where the very edge of the grid is masked land those points do not
% become part of the sponge. Simply setting all points within some range
% next to edge (as was done by the default analytical sponge functions
% in ROMS: ana_hmixcoef.h and now in ana_sponge.h) can introduce elevated
% viscosity/diffusivity in bays and rivers that fall close to the edges
% of the computational domain.
%
  
% svn $Id: sponge.m 742 2014-10-22 18:03:58Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license           Hernan G. Arango      %
%    See License_ROMS.txt                                John Wilkin      %
%=========================================================================%

% Initialize.

factor = 3;
Nfilter = 25;
Lplot = false;
Lwrite = false;

switch numel(varargin)
  case 1
    factor = varargin{1};
  case 2
    factor = varargin{1};
    Nfilter = varargin{2};
  case 3
    factor = varargin{1};
    Nfilter = varargin{2};
    Lplot = varargin{3};
  case 4
    factor = varargin{1};
    Nfilter = varargin{2};
    Lplot = varargin{3};
    Lwrite = varargin{4};
end

%--------------------------------------------------------------------------
% Get ROMS grid structure.
%--------------------------------------------------------------------------

G = get_roms_grid(Gname);

S.lon_rho = G.lon_rho;
S.lat_rho = G.lat_rho;
S.mask_rho = G.mask_rho;

%--------------------------------------------------------------------------
% Compute sponge coefficients
%--------------------------------------------------------------------------

% Put perimeter mask values on edge of zero scratch arrays

P = zeros(size(G.h));
P(:,[1 end]) = G.mask_rho(:,[1 end]);
P([1 end],:) = G.mask_rho([1 end],:);

% Pad ends.

P = P([1 1:end end],[1 1:end end]);

% Trick to force a second row of S=1 inside the perimeter. This will mean
% S remains = 1 on all perimeter points during subsequent smoothing to 
% diffuse the sponge

tmp = conv2(P, ones([5 5]), 'same');
tmp(tmp~=0) = 1;
P = tmp;

% Nfilter passes of a 3-point 2-D boxcar filter works well. In the the
% final result the roll-off with index away from the perimeter is a
% smoothed  version of N-points wide linear sponge. Take maximum or
% perimeter values will fall due to zero padding in the conv2 operation.

for k=1:Nfilter
  tmp = conv2(P, ones([3 3])/9, 'same');
  P = max(tmp, P);
end

% Take one last pass of the filter to smooth the kink at the perimeter.

P = conv2(P, ones([3 3])/9, 'same');

% Extract the rho point grid subset

s = P(2:(end-1),2:(end-1));

% Normalize to be sure range is 0 to 1

s = s/max(s(:));

% Inflate the width at lower values of s

s = s.^(1/4);

% Convert s (0->1) to visc_factor and diff_factor which should range
% from 1 in the interior to a user defined value at the perimeter so that
% the interior value of viscosity is the value "uvnu2 set by ROMS
% ocean.in, i.e. spatially varying viscosity will be visc_factor*uvnu2

S.diff_factor = 1+(factor-1)*s;
S.visc_factor = S.diff_factor;

S

%--------------------------------------------------------------------------
% Plot sponge coefficients.
%--------------------------------------------------------------------------

if (Lplot),
  figure;
% contourf(G.lon_rho, G.lat_rho, nanland(S.diff_factor,G));
  pcolorjw(G.lon_rho, G.lat_rho, nanland(S.diff_factor,G));
  title('Sponge Diffusivity/Viscosity Factor')
  colorbar;
end

%--------------------------------------------------------------------------
% If appropriate, add/write sponge coefficient to GRID NetCDF file.
%--------------------------------------------------------------------------

if (Lwrite)
  add_sponge(Gname, S.visc_factor, S.diff_factor)
  sponge_info = [datestr(now,1) ' created with ' which(mfilename)       ...
                 ': factor = ', num2str(factor)                         ...             
                 ', Nfilter = ', num2str(Nfilter)];
  status = nc_attadd(Gname, 'sponge', sponge_info);
end

