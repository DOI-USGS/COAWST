function F=plot_field(Gname, Hname, Vname, Tindex, Level)

%
% PLOT_FIELD:  Plot requested ROMS variable from input NetCDF file
%
% F=plot_field(Gname, Hname, Vname, Tindex, Level)
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
%
%                     Level > 0,    terrain-following level
%                     Level < 0,    depth (field interpolation)
%
% On Output:
%
%    F             Requested 2D or 3D variable (array)
%

% svn $Id: plot_field.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
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

% Set ROMS Grid structure.
  
if (~isstruct(Gname)),
  G = get_roms_grid(Gname);
else
  G = Gname;
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

%--------------------------------------------------------------------------
% Read in requested variable from donor NetCDF file.
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

%pcolor(X,Y,nanland(F,G)); shading interp; colorbar
pcolorjw(X,Y,nanland(F,G)); shading interp; colorbar

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

%  Mark minimum and maximum locations.

hold on;

hm = plot(X(Imin,Jmin), Y(Imin,Jmin), 'o',                              ...
          X(Imax,Jmax), Y(Imax,Jmax), 's');
set(hm,'MarkerEdgeColor','none','MarkerFaceColor','m','MarkerSize',12);

hb = plot(X(Imin,Jmin), Y(Imin,Jmin), 'ko',                             ...
          X(Imax,Jmax), Y(Imax,Jmax), 'ks');
set(hb,'MarkerSize',12);

if (got.coast),
  hc = plot(Clon,Clat,'k-');
end

hold off;

return
