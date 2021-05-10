function W = wrf_weights (Wname, Rname, Oname, varargin)

%
% WRF_WEIGHTS: Sets weight factors for merging DATA and ROMS components
%
% [W] = wrf_weights(Wname, Rname, Oname, Lsmooth, Nconv, Lplot)
%
% This function computes the WRF melding weights to combine fields from
% DATA and ROMS components. The merging between DATA and ROMS fields is
% done gradually.
%
% The DATA component supplies needed data to a particular ESM component.
% For example, it may export data to WRF at locations not covered by ROMS
% because of smaller grid coverage. If WRF and ROMS grids are incongruent,
% WRF needs to import sea surface temperature (SST) on those grid points
% not covered by ROMS. Thus, the weighting coefficients are used to merge
% the SST data:
%
%   SST_wrf(:,:) = Wroms(:,:) * SST_roms(;,;) + Wdata(:,:) * SST_data(:,:)
%
% where Wroms(:,:) + Wdata(:,:) = 1.
%
% On Input:
%
%    Wname       WRF initial/history NetCDF filename (string)
%
%    Rname       ROMS grid NetCDF filename (string)
%            or, an existing ROMS grid structure (struct)
%
%    Oname       Output weight factors NetCDF filename (string)
%
%    Lsmooth     Switch to compute gradually melding weights (Optional)
%                  (default = false)
%
%    Nconv       Number of smoothing convolutions iterations (Optional)
%                  (default = 50; higher values for wider melding zone)
%
%    Lplot       Switch to plot melding weights (Optional)
%                  (default = false)
%
% On Output:
%
%    W           Melding weights (struct)
%  

% svn $Id: wrf_weights.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license           Hernan G. Arango      %
%    See License_ROMS.txt                           John L. Wilkin        %
%=========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    Lsmooth = false;
    Nconv = 50;
    Lplot = false;
  case 1
    Lsmooth = varargin{1};
    Nconv = 50;
    Lplot = false;
  case 2
    Lsmooth = varargin{1};
    Nconv = varargin{2};
    Lplot = false;
  case 3
    Lsmooth = varargin{1};
    Nconv = varargin{2};
    Lplot = varargin{3};
end

if (Lsmooth)
  Lplot = true;
end

% Initialize out structure.

W = struct('lon', [],                                                   ...
           'lat', [],                                                   ...
           'mask', [],                                                  ...
           'XboxW', [],                                                 ...
           'YboxW', [],                                                 ...
           'XboxR', [],                                                 ...
           'YboxR', [],                                                 ...
           'data_weight_rigid', [],                                     ...
           'ocean_weight_rigid', [],                                    ...
           'data_weight_smooth', [],                                    ...
           'ocean_weight_smooth', []);

% Get WRF grid longitude, latitude, masks, and perimeter.

W.lon = nc_read(Wname, 'XLONG', 1);
W.lat = nc_read(Wname, 'XLAT', 1);

wland = nc_read(Wname, 'LANDMASK', 1);
wlake = nc_read(Wname, 'LAKEMASK', 1);

Imin=1;
Jmin=1;
[Imax,Jmax]=size(W.lon);

W.XboxW=[squeeze(W.lon(Imin:Imax,Jmin));                                ...
         squeeze(W.lon(Imax,Jmin+1:Jmax))';                             ...
         squeeze(flipud(W.lon(Imin:Imax-1,Jmax)));                      ...
         squeeze(fliplr(W.lon(Imin,Jmin:Jmax-1)))'];

W.YboxW=[squeeze(W.lat(Imin:Imax,Jmin));                                ...
         squeeze(W.lat(Imax,Jmin+1:Jmax))';                             ...
         squeeze(flipud(W.lat(Imin:Imax-1,Jmax)));                      ...
         squeeze(fliplr(W.lat(Imin,Jmin:Jmax-1)))'];

% Get ROMS grid structure and perimeter.

if (~isstruct(Rname)),
  G = get_roms_grid(Rname);
else
  G = Gname;
end
S = grid_perimeter(G);

W.XboxR=S.grid.perimeter.X_psi;
W.YboxR=S.grid.perimeter.Y_psi;

clear S

%--------------------------------------------------------------------------
% Compute land/sea mask (0:land, 1:sea) from the WRF masks:
%   landmask:  1: land, 0: elsewhere
%   lakemask:  1: lake, 0: elsewhere
%--------------------------------------------------------------------------

W.mask = ones(size(wland));
W.mask(wland == 1) = 0;           % zeroh out land
W.mask(wlake == 1) = 0;           % zeroh out lakes

% Find WRF grid cells inside of ROMS perimeter.

[IN,ON]=inpolygon(W.lon(:), W.lat(:), W.XboxR, W.YboxR);

%IN(ON) = true;                   % add cells on perimeter

% Compute DATA and ROMS uniform weights.

W.data_weight_rigid  = zeros(size(W.mask));
W.ocean_weight_rigid = zeros(size(W.mask));

W.data_weight_rigid(~IN) = 1;     % DATA component outside ROMS grid
W.ocean_weight_rigid(IN) = 1;     % ROMS contribution to field

% Zeroh out indices in WRF land grid cells.

W.data_weight_rigid (W.mask == 0) = 0;
W.ocean_weight_rigid(W.mask == 0) = 0;

%--------------------------------------------------------------------------
% Smooth gradually the transitio between DATA and OCEAN components.
%--------------------------------------------------------------------------

if (Lsmooth)
  [Dweight,Oweight] = smooth_weights(W, Nconv, Lplot);
  W.data_weight_smooth  = Dweight;
  W.ocean_weight_smooth = Oweight;
end

%--------------------------------------------------------------------------
% Create weight factors NetCDF file.
%--------------------------------------------------------------------------

[Im, Jm] = size(W.lon);
c_weights(Im, Jm, Oname);

% Set few global attributes.

Text = ['WRF coupling import field melding weights ',                   ...
        'between `DATA and ROMS components'];
ncwriteatt(Oname, '/', 'title', Text);

ncwriteatt(Oname, '/', 'wrf_grid', Wname);
ncwriteatt(Oname, '/', 'roms_grid', Rname);

if (Lsmooth)
  ncwriteatt(Oname, '/', 'convolutions', num2str(Nconv));
end

Text = ['Weights file created using Matlab script ',                    ...
        which(mfilename), blanks(1), datestr(now)];
ncwriteatt(Oname, '/', 'history', Text);

%--------------------------------------------------------------------------
% Write out 
%--------------------------------------------------------------------------

ncwrite(Oname, 'lon',  W.lon);
ncwrite(Oname, 'lat',  W.lat);
ncwrite(Oname, 'mask', W.mask);
ncwrite(Oname, 'data_weight_rigid',  W.data_weight_rigid);
ncwrite(Oname, 'ocean_weight_rigid', W.ocean_weight_rigid);

if (Lsmooth)
  ncwrite(Oname, 'data_weight_smooth',  W.data_weight_smooth);
  ncwrite(Oname, 'ocean_weight_smooth', W.ocean_weight_smooth);
end

return
