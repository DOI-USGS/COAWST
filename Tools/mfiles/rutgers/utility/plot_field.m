function F=plot_field(Gname, Hname, Vname, Tindex, varargin)

%
% PLOT_FIELD:  Plot requested ROMS variable from input NetCDF file
%
% F=plot_field(Gname, Hname, Vname, Tindex, Level, Caxis, Mmap, wrtPNG)
%
% This function plots requested ROMS variable from input history
% NetCDF file. This function is very useful when debugging a ROMS
% application. The plotting is not that fancy but it provides enough
% information for browsing ROMS variables very quickly. The location
% of the minumum is marked with a filled magenta circle whereas the
% maximum is marked with a filled magenta square.
%
% It also draws the coastline if the data is available in the
% grid file or structure.
%
% On Input:
%
%    Gname         ROMS Grid NetCDF file/URL name (string)
%              or, an existing ROMS grid structure (struct array)
%
%    Hname         ROMS history NetCDF file/URL name (string)
%
%    Vname         ROMS NetCDF variable name to plot (string)
%
%    Tindex        Time record index used to process (scalar)
%                    (Use Inf or inf for last record)
%
%    Level         If 3D variable, vertical level to plot (scalar)
%                    (Optional, default: surface level)
%
%                     Level > 0,    terrain-following level
%                     Level < 0,    depth (field interpolation)
%
%    Caxis         Color axis (vector)
%                    (Optional, default: [Inf Inf], choosen internally)
%
%    Mmap          Switch to use m_map utility (true or false)
%                    (Optional, default: false)
%
%    wrtPNG        Switch to write out PNG file (true or false)
%                    (Optional, default: false)
%
% On Output:
%
%    F             Requested 2D or 3D variable (array)
%

% svn $Id: plot_field.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

got.coast = false;
got.Mname = false;
got.Xname = false;
got.Yname = false;
got.Zname = false;

isr3d = false;
isw3d = false;
isvec = false;

Tname = [];
Tsize = 0;
recordless = true;
MyProjection = 'mercator';

fill_land = true;                 % use patch (true) or draw coast (false)

%Land = [0.3412 0.2549 0.1843];   % dark brown
%Land = [0.4745 0.3765 0.2980];   % medium brown
%Land = [0.6706 0.5841 0.5176];   % light brown
 Land = [0.6 0.65 0.6];           % gray-green
 Lake = Land;
 
D = nc_dinfo(Hname);
N = D(strcmp({D.Name}, 's_rho')).Length;

%  Optional arguments.

switch numel(varargin)
  case 0
    Level  = N;
    Caxis  = [-Inf Inf];
    Mmap   = false;
    wrtPNG = false;
  case 1
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = [-Inf Inf];
    Mmap   = false;
    wrtPNG = false;
  case 2
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = false;
    wrtPNG = false;
  case 3
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = varargin{3};
    wrtPNG = false;
  case 4
    if (~isinf(varargin{1}))
      Level = varargin{1};
    else
      Level = N;
    end
    Caxis  = varargin{2};
    Mmap   = varargin{3};
    wrtPNG = varargin{4};
end

% Set ROMS Grid structure.
  
if (~isstruct(Gname)),
  G = get_roms_grid(Gname);
else
  G = Gname;
end

% Set colormap.

switch Vname
  case {'temp', 'temp-sur'}
    Cmap = cm_balance(512);
  case {'salt', 'salt_sur'}
    Cmap = cm_delta(512);
  case {'u_sur', 'u_sur_eastward', 'v_sur', 'v_sur_northward'}
    Cmap = cm_curl(512);
  case {'ubar', 'ubar_eastward', 'vbar', 'vbar_northward'}
    Cmap = cm_speed(512);
  case 'zeta'
    Cmap = cm_delta(512);
  case 'swrad'
    Cmap = cm_thermal(512);
  case 'lwrad'
    Cmap = flipud(cm_thermal(512));
  case {'latent', 'sensible', 'shflux'}
    Cmap = cm_balance(512);
  case {'EminusP', 'evaporation', 'ssflux'}
    Cmap = cm_delta(512);
  case 'rain'
    Cmap = flipud(cm_haline(512));
  case {'Uwind', 'Vwind'}
    Cmap = cm_curl(512);
  case {'sustr', 'svstr'}
    Cmap = cm_curl(512);
  case 'Pair'
    Cmap = cm_haline(512);
  case 'Qair'
    Cmap = cm_delta(512);
  case 'Hair'
    Cmap = cm_delta(512);
  otherwise
    Cmap = cm_balance(512);
end

%--------------------------------------------------------------------------
% Get Variable information.
%--------------------------------------------------------------------------

I = nc_vinfo(Hname, Vname);

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

nvdims = length(I.Dimensions);

if (nvdims > 0),
  for n=1:nvdims,
    dimnam = char(I.Dimensions(n).Name);
    switch dimnam
      case {'s_rho', 'level'}
        isr3d = true;
      case 's_w'
        isw3d = true;
      case {'xi_rho','lon_rho', 'lon'}
        Mname = 'mask_rho';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (G.spherical),
            Xname = 'lon_rho';
            Yname = 'lat_rho';
          else
            Xname = 'x_rho';
            Yname = 'y_rho';
          end
          Zname = 'z_r';
          got.Xname = true;
          got.Yname = true;
          got.Zname = true;
        end
      case {'xi_psi','lon_psi'}
        Mname = 'mask_psi';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (G.spherical),
            Xname = 'lon_psi';
            Yname = 'lat_psi';
          else
            Xname = 'x_psi';
            Yname = 'y_psi';
          end
          Zname = 'z_r';
          got.Xname = true;
          got.Yname = true;       
          got.Zname = true;
        end
      case {'xi_u','lon_u'}
        Mname = 'mask_u';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (G.spherical),
            Xname = 'lon_u';
            Yname = 'lat_u';
          else
            Xname = 'x_u';
            Yname = 'y_u';
          end 
          Zname = 'z_u';
          got.Xname = true;
          got.Yname = true;        
          got.Zname = true;
        end
        isvec = true;
      case {'xi_v','lon_v'}
        Mname = 'mask_v';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (G.spherical),
            Xname = 'lon_v';
            Yname = 'lat_v';
          else
            Xname = 'x_v';
            Yname = 'y_v';
          end
          Zname = 'z_v';
          got.Xname = true;
          got.Yname = true;        
          got.Zname = true;
        end
        isvec = true;
      case {'ocean_time', 'time', ~isempty(strfind(dimnam,'time'))}
        recordless = false;    
        Tsize = I.Dimensions(n).Length;
    end
  end
  if (isw3d),
    Zname = 'z_w';
  end  
end

is3d = isr3d || isw3d;

% Check if 'time' attribute is available.

itime = strcmp({I.Attributes.Name}, 'time');
if (any(itime)),
  Tname = I.Attributes(itime).Value;
end

%--------------------------------------------------------------------------
% Get coordinates.
%--------------------------------------------------------------------------

if (isfield(G,Xname)),
  if (~isempty(G.(Xname)))  
    X = G.(Xname);
  else
    error([' PLOT_FIELD - field '', Xname, ''',                          ...
           ' is empty in receiver grid structure: G']);
  end
else
  error([' PLOT_FIELD - unable to find field '', Xname, ''',             ...
         ' in receiver grid structure: G']);
end

if (isfield(G,Yname)),
  if (~isempty(G.(Yname)))
    Y = G.(Yname);
  else
    error([' PLOT_FIELD - field '', Yname, ''',                          ...
           ' is empty in receiver grid structure: G']);
  end
else
  error([' PLOT_FIELD - unable to find field '', Yname, ''',             ...
         ' in receiver grid structure: G']);
end

if (is3d),
  if (isfield(G,Zname)),
    if (~isempty(G.(Zname)))
      Z = G.(Zname);
    else
      error([' PLOT_FIELD - field '', Zname, ''',                        ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_FIELD - unable to find field '', Zname, ''',           ...
           ' in receiver grid structure: G']);
  end
end
  
if (isfield(G,Mname)),
  if (~isempty(G.(Mname)))
    mask = G.(Mname);
  else
    error([' PLOT_FIELD - field '', Mname, ''',                          ...
           ' is empty in receiver grid structure: G']);
  end
else
  error([' PLOT_FIELD - unable to find field '', Mname, ''',             ...
         ' in receiver grid structure: G']);
end

if (isfield(G,'lon_coast') && isfield(G,'lat_coast')),
  Clon = G.lon_coast;
  Clat = G.lat_coast;
  got.coast = true;
end

if (~G.spherical)
  X = 0.001 .* X;                 % km
  Y = 0.001 .* Y;                 % km
end

%--------------------------------------------------------------------------
% Read in requested variable from NetCDF file.
%--------------------------------------------------------------------------

if (~recordless && Tindex > Tsize),
  Tindex = Tsize;                     % process last time record available
end 

if (~isempty(Tname)),
  Tvalue = nc_read(Hname,Tname,Tindex);
  Tattr  = nc_getatt(Hname,'units',Tname);
  Tdays  = true;
  if (~isempty(strfind(Tattr, 'second'))),
    Tvalue = Tvalue/86400;                    % seconds to days
    Tdays  = false;
  end  
  iatt = strfind(Tattr, 'since');
  if (~isempty(iatt)),
    Torigin = Tattr(iatt+6:end);
    epoch   = datenum(Torigin,31);            % 'yyyy-mm-dd HH:MM:SS' 
    Tstring = datestr(epoch+Tvalue);
  else
    Tstring = num2str(Tvalue);    
  end
end

ReplaceValue = NaN;
PreserveType = false;

field = nc_read(Hname,Vname,Tindex,ReplaceValue,PreserveType);

if (is3d),
  if (Level > 0),
    F = squeeze(field(:,:,Level));
  else
    [~,~,F] = nc_slice(G,Hname,Vname,Level,Tindex);
  end
else
  F = field;
end

Fmin = min(F(:));  [Imin,Jmin] = find(F == min(F(:))); 
Fmax = max(F(:));  [Imax,Jmax] = find(F == max(F(:))); 

if (length(Imin) > 1 && length(Jmin > 1)),
  Imin = Imin(1);
  Jmin = Jmin(1);
end

if (length(Imax) > 1 && length(Jmax > 1)),
  Imax = Imax(1);
  Jmax = Jmax(1);
end

%--------------------------------------------------------------------------
% Plot requested field.
%--------------------------------------------------------------------------

figure;

if (Mmap)
  LonMin=min(X(:));   LonMax=max(X(:));
  LatMin=min(Y(:));   LatMax=max(Y(:));
  m_proj(MyProjection,'longitudes',[LonMin,LonMax],                     ...
                      'latitudes' ,[LatMin,LatMax]); 
  m_grid('tickdir','out','yaxisloc','left');
  hold on;
  m_pcolor(X, Y, nanland(F,G));
else
  pcolorjw(X, Y, nanland(F,G));
  hold on;
end

shading interp;
colorbar;
colormap(Cmap);
caxis(Caxis);

if (is3d),
  if (~isempty(Tname)),
    ht = title([Vname, ':', blanks(4),                                  ...
                'Level = ', num2str(Level), ',', blanks(4),             ...
                'Record = ', num2str(Tindex), ',', blanks(4),           ...
                'time = ', Tstring],                                    ...
               'FontSize', 14, 'FontWeight', 'bold' );
  else
    ht = title([Vname, ':', blanks(4),                                  ...
                'Level = ', num2str(Level), ',', blanks(4),             ...
                'Record = ', num2str(Tindex)],                          ...
               'FontSize', 14, 'FontWeight', 'bold' );
  end
else
  if (~isempty(Tname)),
    ht = title([Vname, ':', blanks(4),                                  ... 
                'Record = ', num2str(Tindex), ',', blanks(4),           ...
                'time = ', Tstring],                                    ...
               'FontSize', 14, 'FontWeight', 'bold' );
  else
    ht = title([Vname, ':', blanks(4),                                  ... 
                'Record = ', num2str(Tindex)],                          ...
               'FontSize', 14, 'FontWeight', 'bold' );
  end
end

hx = xlabel(['Min = ', num2str(Fmin), blanks(4),                        ...
             '(', num2str(Imin), ', ', num2str(Jmin), '),', blanks(8),  ...
             'Max = ', num2str(Fmax), blanks(4),                        ...
             '(', num2str(Imax), ', ', num2str(Jmax), ')'],             ...
            'FontSize', 14, 'FontWeight', 'bold' );

%  Mark minimum (down triangle) and maximum (up triangle) locations.

if (Mmap)
  if (fill_land)
    m_gshhs_i('patch', Land, 'edgecolor', Lake);
  else
    m_gshhs_i('color','k');
  end
  m_plot(X(Imin,Jmin), Y(Imin,Jmin), 'v',                               ...
         'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',12);
  m_plot(X(Imax,Jmax), Y(Imax,Jmax), '^',                               ...
         'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',12);
  m_plot(X(Imin,Jmin), Y(Imin,Jmin), 'kv',                              ...
         'MarkerSize',12);
  m_plot(X(Imax,Jmax), Y(Imax,Jmax), 'k^',                              ...
         'MarkerSize',12);
else
  plot(X(Imin,Jmin), Y(Imin,Jmin), 'v',                                 ...
       X(Imax,Jmax), Y(Imax,Jmax), '^',                                 ...
       'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',12);
  plot(X(Imin,Jmin), Y(Imin,Jmin), 'wv',                                ...
       X(Imax,Jmax), Y(Imax,Jmax), 'w^',                                ...
       'MarkerSize',12);
  if (got.coast),
    hc = plot(Clon,Clat,'k-');
  end
end

%  Write out PNG file.

if (wrtPNG)
  png_file=strcat(Vname,'_',num2str(Tindex, '%4.4i'),'.png');
  print(png_file, '-dpng', '-r300');
end

hold off;

return
