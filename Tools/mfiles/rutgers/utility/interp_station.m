function S = interp_station(Ginp, ncname, Vname, Tindex, Xsta, Ysta, varargin)

%
% INTERP_STATION: Interpolates a ROMS variable to requested location
%
% S = interp_station(G, ncname, Vname, Tindex, Xsta, Ysta, Zsta, method)
%
% Interpolate a ROMS variable to requested station location (X,Y,Z).
%
% On Input:
%
%    Ginp          ROMS Grid/History NetCDF file/URL name containing
%                    all grid variables (string)
%
%              or, an existing grid structure to which the vertical
%                    coordinates are added or updated (struct array)
%
%    ncname        ROMS NetCDF file data to process (string)
%
%    Vname         ROMS NetCDF variable name to process (string)
%
%    Tindex        Time record index to process (scalar)
%
%    Xsta          Station X-location(s) to interpolate (scalar or vector)
%                    If spherical  => longitude (degrees_east)
%                    If Cartesian  => X-coordinate (meters)
%
%    Ysta          Station Y-location(s) to interpolate (scalar or vector)
%                    If spherical  => latitude (degrees_north)
%                    If Cartesian  => Y-coordinate (meters)
%
%    Zsta          Station Z-location(s) depths (m) to interpolate
%                    (optional; scalar or vector; depths are negative)
%
%                    - IF Vname is a 2D field use Zsta = []  
%                    - IF Vname is a 3D field and Zeta = [], no vertical
%                      interpolation is carried out and all vertical
%                      level values are written in output structure
%
%    method        Horizontal interpolation method (optional; string)
%                    'natural'     natural neighbor interpolation
%                    'linear'      linear interpolation (default)
%                    'nearest'     nearest-neighbor interpolation
% On Output:
%
%    S             Interpolated data (struct)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Check input arguments.

if (~isstruct(Ginp))
  S.grid_name = Ginp;
  G = get_roms_grid(Ginp, ncname);
else
  G = Ginp;
end

vert_interp = false;

switch numel(varargin)
  case 0
    method = 'linear';
  case 1
    Zsta   = varargin{1};
    if (~isempty(Zsta))
      vert_interp = true;
    end
    method = 'linear';
  case 2
    Zsta   = varargin{1};
    if (~isempty(Zsta))
      vert_interp = true;
    end
    method = varargin{2};
end

Vmethod = 'linear';          % vertical interpolation method

% Inquire requested variable.

I = nc_inq(ncname);

if (~any(strcmp({I.Variables.Name}, Vname)))
  error(['INTERP_STATION - Cannot find variable: ', Vname,              ...
         ' in file ', ncname])
else
  V = nc_vinfo(ncname, Vname);
  nvdims = length(V.Dimensions);
end

% Deterrmine spherical switch.

spherical=false;

if (any(strcmp({I.Variables.Name}, 'spherical')))
  spherical = nc_read(ncname, 'spherical');
  if (ischar(spherical))
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
end

% Define output structure.


if (spherical)
  S = struct('ncfile'           , [],                                   ...
             'variable_name'    , [],                                   ...
             'number_levels'    , [],                                   ...
             'longitude'        , [],                                   ...
             'latitude'         , [],                                   ...
             'depth'            , [],                                   ...
             'time_record'      , [],                                   ...
             'date'             , [],                                   ...
             'values'           , [],                                   ...
             'min'              , [],                                   ...
             'max'              , [],                                   ...
             'mean'             , []);
 S.longitude = Xsta;
 S.latitude  = Ysta;
else
  S = struct('ncfile'           , [],                                   ...
             'variable_name'    , [],                                   ...
             'number_levels'    , [],                                   ...
             'x'                , [],                                   ...
             'y'                , [],                                   ...
             'depth'            , [],                                   ...
             'time_record'      , [],                                   ...
             'date'             , [],                                   ...
             'values'           , [],                                   ...
             'min'              , [],                                   ...
             'max'              , [],                                   ...
             'mean'             , []);
 S.x = Xsta;
 S.y = Ysta;
end

S.ncfile        = ncname;
S.variable_name = Vname;
S.time_record   = Tindex;

% Get ROMS time.

Tvalue = nc_read(ncname, 'ocean_time', Tindex);
Tunits = nc_getatt(ncname, 'units', 'ocean_time');

if (strfind(Tunits,'seconds since ')),
  basedate = Tunits(14:end);
  epoch    = datenum(basedate);
  S.date   = datestr(epoch + Tvalue/86400, 0);
else
  S = rmfield(S, 'date');  
end

%--------------------------------------------------------------------------
% Variable locations.
%--------------------------------------------------------------------------

% Get variable horizontal coordinates.

is3d  = false;
isr3d = false;
isw3d = false;

if (nvdims > 0)
  for n=1:nvdims,
    dimnam = char(V.Dimensions(n).Name);
    switch dimnam
      case 's_rho'
        is3d  = true;
        isr3d = true;
        N = G.N;     
      case 's_w'
        is3d  = true;
        isw3d = true;
        N = G.N + 1;     
      case {'xi_rho','eta_rho'}
        VM = G.mask_rho;
        if (spherical)
          VX = G.lon_rho;
          VY = G.lat_rho;
        else
          VX = G.x_rho;
          VY = G.y_rho;
        end
        VZ = G.z_r;     
      case {'xi_psi','eta_psi'}
        VM = G.mask_psi;
        if (spherical)
          VX = G.lon_psi;
          VY = G.lat_psi;
        else
          VX = G.x_psi;
          VY = G.y_psi;
        end
      case {'xi_u','eta_u'}
        VM = G.mask_u;
        if (spherical)
          VX = G.lon_u;
          VY = G.lat_u;
        else
          VX = G.x_u;
          VY = G.y_u;
        end
        VZ = G.z_u;     
      case {'xi_v','eta_v'}
        VM = G.mask_v;
        if (spherical)
          VX = G.lon_v;
          VY = G.lat_v;
        else
          VX = G.x_v;
          VY = G.y_v;
        end
        VZ = G.z_v;     
    end
  end
end
  
if (~is3d)
  vert_interp = false;
  S = rmfield(S, 'depth');  
end

if (isw3d && vert_interp)
  VZ = G.z_w;
end

%--------------------------------------------------------------------------
% Interpolate to requested locations.
%--------------------------------------------------------------------------

% Read ROMS variable.

ReplaceValue = NaN;
PreserveType = false;

Values = nc_read(ncname, Vname, Tindex, ReplaceValue, PreserveType);

% Interpolate with scatteredInterpolant (ROMS horizontal coordinates are
% not plaid).

%..........................................................................
if (is3d)                             % 3D field variable interpolation
%..........................................................................

  x = VX(:);                          % scatteredInterpolant wants 1-D
  y = VY(:);                          % vectors as input

  S.number_levels = N;
  
  Dind = find(VM < 0.5);
  if ~isempty(Dind)
    x(Dind) = [];                     % remove land points, if any
    y(Dind) = [];
  end  
  
  F = scatteredInterpolant(x, y, ones(size(x)), method);

  VK = NaN([length(Xsta), N]);        % initialize values
  if (vert_interp)
    S.values = ([length(Xsta), N]);
  end
  
  if (vert_interp)                    % initialize depths
    ZK = NaN([length(Xsta), N]);
    S.values = ([length(Xsta), 1]);
  end
  
  % Horizontal interpolation.

  for k = 1:N                         % level-by-level

    vk = squeeze(Values(:,:,k));
    vk = vk(:);

    if (~isempty(Dind))
      vk(Dind) = [];
    end

    F.Values = vk;
    VK(:,k) = F(Xsta, Ysta);          % interpolate at level k
    
    if (vert_interp)
      zk = squeeze(VZ(:,:,k));
      zk = zk(:); 
      if (~isempty(Dind))
        zk(Dind) = [];
      end

      F.Values = zk;
      ZK(:,k) = F(Xsta, Ysta);        % interpolate at level k
    end 
  
  end
  
  % If appropriate, perform vertical interpolation.

  if (vert_interp)
    S.depth = Zsta;
    for i = 1:length(Xsta)
      S.values(i) = interp1(ZK(i,:), VK(i,:), Zsta(i), Vmethod);    
    end
  else
    S.values = VK;
  end

%..........................................................................
else                                  % 2D field variable interpolation
%..........................................................................

  S = rmfield(S, 'number_levels');

  x = VX(:);                          % scatteredInterpolant wants 1-D
  y = VY(:);                          % vectors as input
  v = Values(:);

  Dind = find(VM < 0.5);
  if (~isempty(Dind))
    x(Dind) = [];                     % remove land points, if any
    y(Dind) = [];
    v(Dind) = [];
  end

  Dind = isnan(v);                    % remove NaN's
  if (any(Dind))
    x(Dind) = [];
    y(Dind) = [];
    v(Dind) = [];
  end  
  
  F = scatteredInterpolant(x, y, v, method);
  S.values = F(Xsta, Ysta);

end

% Statistics.

S.min  = min(S.values(:));
S.max  = max(S.values(:));
S.mean = mean(S.values(:), 'omitnan');

return
