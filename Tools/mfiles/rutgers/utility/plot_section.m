function F=plot_section(Gname, Hname, Vname, Tindex, orient, index,     ...
                        varargin)

%
% PLOT_SECTION:  Plots requested variable section from input NetCDF file
%
% F=plot_section(Gname, Hname, Vname, Tindex, orient, index,            ...
%                ptype, Caxis, wrtPNG)
%
% This function plots requested ROMS variable section from input
% history NetCDF file.
%
% On Input:
%
%    Gname         ROMS Grid NetCDF file/URL name (string)
%              or, an existing ROMS grid structure (struct array)
%
%    Hname         ROMS history NetCDF file/URL name (string)
%              or, a field structure to plot (struct array)
%
%    Vname         ROMS NetCDF variable name to plot (string)
%
%    Tindex        Time record index used to process (scalar)
%                    (Use Inf or inf for last record)
%
%    orient        Orientation for the extraction (string):
%                    orient='r'  row (west-east) extraction
%                    orient='c'  column (south-north) extraction
%
%    index         Row or column to extract (integer)
%                    if orient='r', then   1 <= index <= Mp
%                    if orient='c', then   1 <= index <= Lp
%
%    ptype         Plot type (integer)
%                     ptype < 0     use contourf with abs(ptype) colors
%                     ptype = 0     do not draw figure
%                     ptype = 1     use pcolor
%                     ptype = 2     use pcolorjw
%
%    Caxis         Color axis (vector)
%                    (Optional, default: [Inf Inf], choosen internally)
%
%    wrtPNG        Switch to write out PNG file (true or false)
%                    (Optional, default: false)
%
% On Output:
%
%    F             Requested 3D variable section (struc)
%
% To limit the depth axis to the upper 500m use:
%
%    axis([-Inf Inf -500 0)]

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

if (orient == 'c')
  F = struct('ncname'     , [], 'Vname'     , [],                       ...
             'Tindex'     , [], 'Tname'     , [], 'Tstring'   , [],     ...
             'column'     , [],                                         ...
             'value'      , [], 'min'       , [], 'max'       , [],     ...
             'X'          , [], 'Z'         , [],                       ...
             'Imin'       , [], 'Jmin'      , [],                       ...
             'Imax'       , [], 'Jmax'      , []);
  F.column = index;
else
  F = struct('ncname'     , [], 'Vname'     , [],                       ...
             'Tindex'     , [], 'Tname'     , [], 'Tstring'   , [],     ...
             'row'        , [],                                         ...
             'value'      , [], 'min'       , [], 'max'       , [],     ...
             'X'          , [], 'Z'         , [],                       ...
             'Imin'       , [], 'Jmin'      , [],                       ...
             'Imax'       , [], 'Jmax'      , []);
  F.row = index;
end

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

Land = [0.6 0.65 0.6];           % gray-green

%  Optional arguments.

switch numel(varargin)
  case 0
    ptype  = 0;
    Caxis  = [-Inf Inf];
    wrtPNG = false;
  case 1
    ptype  = varargin{1};
    Caxis  = [-Inf Inf];
    wrtPNG = false;
  case 2
    ptype  = varargin{1};
    Caxis  = varargin{2};
    wrtPNG = false;
  case 3
    ptype  = varargin{1};
    Caxis  = varargin{2};
    wrtPNG = varargin{3};
end

% Set ROMS Grid structure.

if (~isstruct(Gname))
  G = get_roms_grid(Gname);
else
  G = Gname;
end

if (~isstruct(Hname))
  getdata = true;
  D = nc_dinfo(Hname);
else
  getdata = false;

  H = Hname;
  X = H.X;
  Y = H.Y;
  Z = H.Z;
  mask = H.mask;

  field = H.field;
  [Im,Jm,Km] = size(field);
  Hname = [];
  Tstring = [];
end

% Set colormap.

switch Vname
  case {'AKs', 'AKt', 'AKv'}
    Cmap = cm_matter(255);
    scale = 1.0;
  case {'gls', 'tke'}
    Cmap = flipud(magma(255));
    scale = 1.0;
  case {'mud_01', 'mud_02', 'mud_03'}
    Cmap = cm_turbid(255);
    scale = 1.0;
  case {'salt'}
    Cmap = flipud(cividis(255));
    scale = 1.0;
  case {'temp'}
    Cmap = cm_balance(512);
    scale = 1.0;
  case {'u', 'v'}
    Cmap = cm_curl(512);
    scale = 1.0;
  case {'w','Ha', 'hz'}
    Cmap = cm_curl(512);
    scale = 86400.0;               % m/day
  case {'Hz', 'hz'}
    Cmap = cividis(512);
    scale = 1.0;
  otherwise
    Cmap = cm_balance(512);
    scale = 1.0;
end

%--------------------------------------------------------------------------
% Get Variable information.
%--------------------------------------------------------------------------

if (getdata)

  V = nc_vnames(Hname);

  switch (Vname)
    case {'Hz', 'hz'}
      if (any(strcmp({V.Variables(:).Name}, 'Hz')))
        I = nc_vinfo(Hname, 'Hz');
      else
        I = nc_vinfo(Hname, 'z_w');
      end
    otherwise
      I = nc_vinfo(Hname, Vname);
  end

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

  nvdims = length(I.Dimensions);

  if (nvdims > 0)
    for n=1:nvdims
      dimnam = char(I.Dimensions(n).Name);
      switch dimnam
        case {'s_rho', 'level'}
          isr3d = true;
        case 's_w'
          isw3d = true;
        case {'xi_rho','lon_rho', 'lon'}
          Mname = 'mask_rho';
          got.Mname = true;
          if (~(got.Xname || got.Yname))
            if (G.spherical)
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
          if (~(got.Xname || got.Yname))
            if (G.spherical)
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
          if (~(got.Xname || got.Yname))
            if (G.spherical)
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
          if (~(got.Xname || got.Yname))
            if (G.spherical)
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
    if (isw3d)
      Zname = 'z_w';
    end
  end

  is3d = isr3d || isw3d;

% Check if 'time' attribute is available.

  itime = strcmp({I.Attributes.Name}, 'time');
  if (any(itime))
    Tname = I.Attributes(itime).Value;
  end
end

%--------------------------------------------------------------------------
% Get coordinates.
%--------------------------------------------------------------------------

if (getdata)
  if (isfield(G,Xname))
    if (~isempty(G.(Xname)))
      X = G.(Xname);
    else
      error([' PLOT_SECTION - field '', Xname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_SECTION - unable to find field '', Xname, ''',        ...
           ' in receiver grid structure: G']);
  end

  if (isfield(G,Yname))
    if (~isempty(G.(Yname)))
      Y = G.(Yname);
    else
      error([' PLOT_SECTION - field '', Yname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_SECTION - unable to find field '', Yname, ''',        ...
           ' in receiver grid structure: G']);
  end

  if (is3d)
    if (isfield(G,Zname))
      if (~isempty(G.(Zname)))
        Z = G.(Zname);
      else
        error([' PLOT_SECTION - field '', Zname, ''',                   ...
               ' is empty in receiver grid structure: G']);
      end
    else
      error([' PLOT_SECTION - unable to find field '', Zname, ''',      ...
             ' in receiver grid structure: G']);
    end
    N = size(Z,3);
  else
    error([' PLOT_SECTION - cannot plot a section for a 2D field: ''',  ...
           Vname, ''])
  end

  if (isfield(G,Mname))
    if (~isempty(G.(Mname)))
      mask = G.(Mname);
    else
      error([' PLOT_SECTION - field '', Mname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_SECTION - unable to find field '', Mname, ''',        ...
         ' in receiver grid structure: G']);
  end

  if (~G.spherical)
    if (max(X(:)) > 5000 || max(Y(:)) > 5000)
      X = 0.001 .* X;                           % scale to kilometers
      Y = 0.001 .* Y;
    end
  end
end

%--------------------------------------------------------------------------
% Read in requested variable from NetCDF file.
%--------------------------------------------------------------------------

if (getdata)
  if (~recordless && Tindex > Tsize)
    Tindex = Tsize;                   % process last time record available
  end

  if (~isempty(Tname))
    Tvalue = nc_read(Hname,Tname,Tindex);
    Tattr  = nc_getatt(Hname,'units',Tname);
    Tdays  = true;
    if (~isempty(strfind(Tattr, 'second')))
      Tvalue = Tvalue/86400;                  % seconds to days
      Tdays  = false;
    end
    iatt = strfind(Tattr, 'since');
    if (~isempty(iatt))
      Torigin = Tattr(iatt+6:iatt+24);        % remove extra like GMT
      epoch   = datenum(Torigin,31);          % 'yyyy-mm-dd HH:MM:SS'
      Tstring = datestr(epoch+Tvalue);
    else
      Tstring = num2str(Tvalue);
    end
  end

  ReplaceValue = NaN;
  PreserveType = false;

  switch (Vname)
    case {'Hz', 'hz'}
      if (any(strcmp({V.Variables(:).Name}, 'Hz')))
        field = nc_read(Hname,Vname,Tindex,ReplaceValue,PreserveType);
      else
        Zw = nc_read(Hname,'z_w',Tindex,ReplaceValue,PreserveType);
        Np = size(Zw,3); N = Np-1;
        for k = 2:Np
          field(:,:,k-1) = Zw(:,:,k) - Zw(:,:,k-1);
        end
      end
      Z = G.z_r;
   otherwise
      field = nc_read(Hname,Vname,Tindex,ReplaceValue,PreserveType);
  end

end

% Extract data section to plot.

switch orient
  case 'c'
    V = squeeze(field(index,:,:)); [Im,Km]=size(V);
    m = squeeze(mask(index,:));
    M = repmat(m, [1 Km]);
    V = nanland(V, M);
    s = squeeze(Y(index,:));
    Z = squeeze(Z(index,:,:));
    if (G.spherical)                            % bathymetry
      x = squeeze(G.lat_rho(index,:));
    else
      x = squeeze(G.y_rho(index,:));
    end
    z = -squeeze(G.h(index,:));
 case 'r'
    V = squeeze(field(:,index,:)); [Im,Km]=size(V);
    m = squeeze(mask(:,index));
    M = repmat(m, [1 Km]);
    V = nanland(V, M);
    s = squeeze(X(:,index));
    S = repmat(s, [1 Km]);
    Z = squeeze(Z(:,index,:));
    if (G.spherical)                            % bathymetry
      x = squeeze(G.lon_rho(:,index));
    else
      x = squeeze(G.x_rho(:,index));
    end
    z = -squeeze(G.h(:,index));
end

V = V .* scale;

Fmin = min(V(:));  [Imin,Jmin] = find(V == min(V(:)));
Fmax = max(V(:));  [Imax,Jmax] = find(V == max(V(:)));

if (length(Imin) > 1 && length(Jmin > 1))
  Imin = Imin(1);
  Jmin = Jmin(1);
end

if (length(Imax) > 1 && length(Jmax > 1))
  Imax = Imax(1);
  Jmax = Jmax(1);
end

%  Load values to output structure.

F.ncname  = Hname;
F.Vname   = Vname;
F.Tindex  = Tindex;
F.Tname   = Tname;
F.Tstring = Tstring;
F.value   = V;
F.min     = Fmin;
F.max     = Fmax;
F.X       = S;
F.Z       = Z;
F.Imin    = Imin;
F.Jmin    = Jmin;
F.Imax    = Imax;
F.Jmax    = Jmax;

%--------------------------------------------------------------------------
% If applicable, plot requested field.
%--------------------------------------------------------------------------

if (ptype ~= 0)

  figure;

  if (ptype < 0)
   [c,hp] = contourf(S,Z,V,abs(ptype));
  elseif (ptype == 1)
   hp = pcolor(S,Z,V);
  elseif (ptype == 1)
   hp = pcolorjw(S,Z,V);
  end

  shading interp;
  colorbar; colormap(Cmap); caxis(Caxis);

  % plot bathymetry curve.

  hold on;
  area(x, z, min(z), 'FaceColor', Land, 'EdgeColor', Land);
  hold off;

  if (~isempty(Tname))
    ht = title([untexlabel(Vname), ':', blanks(4),                      ...
                'Record = ', num2str(Tindex), ',', blanks(4),           ...
                'time = ', Tstring],                                    ...
                'FontSize', 14, 'FontWeight', 'bold' );
  else
    ht = title([untexlabel(Vname), ':', blanks(4),                      ...
                'Record = ', num2str(Tindex)],                          ...
                'FontSize', 14, 'FontWeight', 'bold' );
  end

  hx = xlabel(['Min = ', num2str(Fmin), blanks(4),                      ...
               '(', num2str(Imin), ', ', num2str(Jmin), '),',           ...
               blanks(8), 'Max = ', num2str(Fmax), blanks(4),           ...
               '(', num2str(Imax), ', ', num2str(Jmax), ')'],           ...
               'FontSize', 14, 'FontWeight', 'bold' );

  hy = ylabel('Z (m)', 'FontSize', 14, 'FontWeight', 'bold');

%  Write out PNG file.

  if (wrtPNG)
    png_file=strcat(Vname,'_',num2str(Tindex, '%4.4i'),'.png');
    exportgraphics(gcf, png_file, 'resolution', 300);
  end

end

return
