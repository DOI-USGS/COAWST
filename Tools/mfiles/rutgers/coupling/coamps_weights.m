function W = coamps_weights (Cname, Rname, Oname, varargin)

%
% COAMPS_WEIGHTS: Sets weight factors for merging DATA and ROMS components
%
% [W] = coamps_weights(Cname, Rname, Oname, Lsmooth, Nconv, Lplot)
%
% This function computes the COAMPS melding weights to combine fields from
% DATA and ROMS components. The merging between DATA and ROMS fields is
% done gradually.
%
% The DATA component supplies needed data to a particular ESM component.
% For example, it may export data to COAMPS at locations not covered by
% ROMS because of smaller grid coverage. If COAMPS and ROMS grids are
% incongruent, COAMPS needs to import sea surface temperature (SST) on
% those grid points not covered by ROMS. Thus, the weighting coefficients
% are used to merge the SST data:
%
% SST_coamps(:,:) = Wroms(:,:) * SST_roms(;,;) + Wdata(:,:) * SST_data(:,:)
%
% where Wroms(:,:) + Wdata(:,:) = 1.
%
% Both COAMPS and ROMS logitudes are wrapped to the interval [0 360]
% using Matlab intrinsic function wrapTo360 for easy use in applications
% with mixed conventions.
%
% On Input:
%
%    Cname       COAMPS history HDF5 filename (string)
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

% svn $Id: coamps_weights.m 996 2020-01-10 04:28:56Z arango $
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
           'XboxC', [],                                                 ...
           'YboxC', [],                                                 ...
           'XboxR', [],                                                 ...
           'YboxR', [],                                                 ...
           'data_weight_rigid', [],                                     ...
           'ocean_weight_rigid', [],                                    ...
           'data_weight_smooth', [],                                    ...
           'ocean_weight_smooth', []);

% Get COAMPS grid longitude, latitude, masks, and perimeter.
% If COAMPS longitudes is a module of 360 degrees, wrap angle in
% degrees to [-180 180] as in ROMS.

slon = 'longit_sfc_000000_000000_1';
slat = 'latitu_sfc_000000_000000_1';
smsk = 'lndsea_sfc_000000_000000_1';

V = nc_vnames(Cname);

Vindex = strncmp({V.Variables.Name}, slon, length(slon));
if (any(Vindex))
  Vname = V.Variables(Vindex).Name;
  W.lon = nc_read(Cname, Vname);
  W.lon = wrapTo360(W.lon);
else
  error(['Cannot find longitude variable ', slon, '...']);
end

Vindex = strncmp({V.Variables.Name}, slat, length(slat));
if (any(Vindex))
  Vname = V.Variables(Vindex).Name;
  W.lat = nc_read(Cname, Vname);
else
  error(['Cannot find latitute variable ', slat, '...']);
end

Vindex = strncmp({V.Variables.Name}, smsk, length(smsk));
if (any(Vindex))
  Vname = V.Variables(Vindex).Name;
  Cmask = nc_read(Cname, Vname);
else
  error(['Cannot find land/sea mask variable ', smsk, '...']);
end

Imin=1;
Jmin=1;
[Imax,Jmax]=size(W.lon);

W.XboxC=[squeeze(W.lon(Imin:Imax,Jmin));                                ...
         squeeze(W.lon(Imax,Jmin+1:Jmax))';                             ...
         squeeze(flipud(W.lon(Imin:Imax-1,Jmax)));                      ...
         squeeze(fliplr(W.lon(Imin,Jmin:Jmax-1)))'];

W.YboxC=[squeeze(W.lat(Imin:Imax,Jmin));                                ...
         squeeze(W.lat(Imax,Jmin+1:Jmax))';                             ...
         squeeze(flipud(W.lat(Imin:Imax-1,Jmax)));                      ...
         squeeze(fliplr(W.lat(Imin,Jmin:Jmax-1)))'];

% Get ROMS grid structure and perimeter.

if (~isstruct(Rname)),
  G = get_roms_grid(Rname);
else
  G = Gname;
end
G.lon_psi = wrapTo360(G.lon_psi);

S = grid_perimeter(G);

W.XboxR=S.grid.perimeter.X_psi;
W.YboxR=S.grid.perimeter.Y_psi;

clear S

%--------------------------------------------------------------------------
% Compute land/sea mask (0:land, 1:sea) from the COAMPS mask:
%   Cmask: -1: inland lake, 0: sea, 1: land, 2: ice, 3: seaice
%   lakemask:  1: lake, 0: elsewhere
%--------------------------------------------------------------------------

W.mask = zeros(size(Cmask));
W.mask(Cmask == 0) = 1;           % set 1=sea, 0:elsewhere

% Find COAMP grid cells inside of ROMS perimeter.

[IN,ON]=inpolygon(W.lon(:), W.lat(:), W.XboxR, W.YboxR);

%IN(ON) = true;                   % add cells on perimeter

% Compute DATA and ROMS uniform weights.

W.data_weight_rigid  = zeros(size(W.mask));
W.ocean_weight_rigid = zeros(size(W.mask));

W.data_weight_rigid(~IN) = 1;     % DATA component outside ROMS grid
W.ocean_weight_rigid(IN) = 1;     % ROMS contribution to field

% Zeroh out indices in COAMPS land grid cells.

W.data_weight_rigid (Cmask ~= 0) = 0;
W.ocean_weight_rigid(Cmask ~= 0) = 0;

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

Text = ['COAMPS coupling import field melding weights ',                ...
        'between DATA and ROMS components'];
ncwriteatt(Oname, '/', 'title', Text);

ncwriteatt(Oname, '/', 'coamps_grid', Cname);
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
