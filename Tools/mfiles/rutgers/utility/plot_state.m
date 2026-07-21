function S = plot_state(Gname, Sname, rec, varargin)

% PLOT_STATE:  Plots data assimilation state vector fields
%
% F = plot_state(Gname, Sname, rec, level, ptype, Mmap, orient, index,
%                wrtPNG, PNGsuffix, R)
%
% This function plots horizontal or cross-section fields of the data
% assimilation state vector like increment, analysis, and trajectory.
%
% In cross-sections, we change the axis to the upper range, for example:
%
%  axis([-Inf Inf -500 0]);      % zoom of the uppr 500m
%
% On Input:
%
%    Gname         ROMS Grid NetCDF file/URL name (string)
%              or, an existing ROMS grid structure (struct array)
%
%    Iname         Increment NetCDF filename (string)
%
%    rec           Field/file time record to process (integer)
%                    If negative, it overlay parallel tile partition
%
%    level         horizontal or vertical section index (optional; scalar)
%                    level > 0   model level 1:N (default=N)
%                    level < 0   iso-depth slice (interpolation)
%
%    ptype         Plot type (optional; integer)
%                    ptype < 0   use contourf with abs(ptype) colors
%                    ptype = 1   use pcolor  (default)
%                    ptype = 2   use pcolorjw
%
%    Mmap          Switch to use map projection (optional; integer)
%                    Mmap = 0,  no map projection (default)
%                    Mmap = 1,  'm_map' utility
%                    Mmap = 2,  native Matlab toolbox
%
%    orient        Orientation of section extraction (optional; string):
%                    orient='r'  row (west-east) vertical section
%                    orient='c'  column (south-north) vertical section
%
%    index         Cross-section gris index (optional; integer)
%                    if orient='r', then   1 <= index <= Mp  west-east
%                    if orient='c', then   1 <= index <= Lp  south-north
%
%    wrtPNG        Switch to write out PNG files (optional; switch)
%                    (default: 0 or false)
%
%                    if abs(wrtPNG) > 1 and cross-sections, then
%                    the vertical axis is zoom at Zmin. Then,
%
%                      Zmin   = -abs(wrtPNG)
%                      wrtPNG = true
%                      axis([-Inf Inf Zmin 0])
%
%                    if wrtPNG < 1, ommit figure tile, doTitle = false
%
%    PNGsuffix     PNG filename suffix qualifier (string; OPTIONAL)
%
%    R             State fields color range values (struct)
%
%                    R.zeta  [min max]
%                    R.u     [min max]
%                    R.v     [min max]
%                    R.temp  [min max]
%                    R.salt  [min max]
%
%                    use [-Inf Inf] for full range
%
% On Output:
%
%    S             Processed state variable structure (array)
%

% Initialize.

F = struct('ncname'     , [], 'Vname'     , [],                       ...
           'Tindex'     , [], 'Tname'     , [], 'Tstring'     , [],   ...
           'Level'      , [], 'is3d'      , [], 'tiling'      , [],   ...
           'X'          , [], 'Y'         , [], 'value'       , [],   ...
           'min'        , [], 'max'       , [], 'checkval'    , [],   ...
           'Caxis'      , [], 'Cmap'      , [], 'ptype'       , [],   ...
           'doMap'      , [], 'projection', [],                       ...
           'orient'     , [], 'index'     , [],                       ...
           'gotCoast'   , [], 'lon_coast' , [], 'lat_coast'   , [],   ...
           'shading'    , [], 'pltHandle' , [], 'wrtPNG'      , []);

F.projection = 'mercator';

F.ncname   = Sname;
F.gotCoast = false;
F.shading  = 'interp';
F.Tname    = [];

Land = [0.6 0.65 0.6];           % gray-green

% State vector variables.

I        = nc_inq(Sname);
Vnames   = {I.Variables.Name};
N        = I.Dimensions(strcmp({I.Dimensions.Name}, 's_rho')).Length;
Tindex   = abs(rec);

Svarlist = {'zeta', 'u', 'v', 'u_eastward', 'v_northward', 'temp', 'salt'};

% Set tile partition, if 'tiling' global attribute exist.

tile = strcmp({I.Attributes.Name}, 'tiling');
if (any(tile))
  tiling = I.Attributes(tile).Value;
  tiling(strfind(tiling, 'x'))=' ';
  F.tiling = str2num(tiling);
end
draw_tiling = rec < 0;
doPNG = false;
doRange = false;
doTitle = true;
doZoom = false;

% Optional arguments.

switch numel(varargin)
  case 0
    level     = N;
    ptype     = 1;
    Mmap      = 0;
    orient    = [];
    index     = [];
    wrtPNG    = false;
    PNGsuffix = [];
    R         = [];
 case 1
    level     = varargin{1};
    ptype     = 1;
    Mmap      = 0;
    orient    = [];
    index     = [];
    wrtPNG    = false;
    PNGsuffix = [];
    R         = [];
 case 2
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = 0;
    orient    = [];
    index     = [];
    wrtPNG    = false;
    PNGsuffix = [];
    R         = [];
 case 3
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = varargin{3};
    orient    = [];
    index     = [];
    wrtPNG    = false;
    PNGsuffix = [];
    R         = [];
 case 4
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = varargin{3};
    orient    = varargin{4};
    index     = N;
    wrtPNG    = false;
    PNGsuffix = [];
    R         = [];
 case 5
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = varargin{3};
    orient    = varargin{4};
    index     = varargin{5};
    wrtPNG    = false;
    PNGsuffix = [];
    R         = [];
 case 6
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = varargin{3};
    orient    = varargin{4};
    index     = varargin{5};
    wrtPNG    = varargin{6};
    PNGsuffix = [];
    R         = [];
    doPNG     = true;
 case 7
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = varargin{3};
    orient    = varargin{4};
    index     = varargin{5};
    wrtPNG    = varargin{6};
    PNGsuffix = varargin{7};
    R         = [];
    doPNG     = true;
 case 8
    level     = varargin{1};
    ptype     = varargin{2};
    Mmap      = varargin{3};
    orient    = varargin{4};
    index     = varargin{5};
    wrtPNG    = varargin{6};
    PNGsuffix = varargin{7};
    R         = varargin{8};
    doPNG     = true;
    doRange   = true;
end

% Set parameters affecting the writing of PNG files.

if (doPNG)
  if (~islogical(wrtPNG))
    if (wrtPNG < 0)
      doTitle = false;
    end
    if (abs(wrtPNG) > 1)
      Zmin   = -abs(wrtPNG);
      doZoom = true;
      wrtPNG = true;
    end
  end
end

F.Caxis  = [-Inf Inf];
F.doMap  = Mmap;
F.index  = index;
F.orient = orient;
F.ptype  = ptype;
F.wrtPNG = false;

doSection = false;
if ~isempty(orient)
  if (orient == 'c' || orient == 'r')
    doSection = true;
  end
end

% Set ROMS Grid structure.

if (~isstruct(Gname))
  G = get_roms_grid(Gname);
else
  G = Gname;
end

%--------------------------------------------------------------------------
% Plot each avaliable state variable.
%--------------------------------------------------------------------------

nfield = 0;

for var = Svarlist

  field = char(var);
  ind   = strcmp(Vnames, field);
  isr3d = false;
  isw3d = false;

% Determine grid geometry variables.

  if (any(ind))

    nfield = nfield+1;
    Dnames = {I.Variables(ind).Dimensions.Name};
    nvdims = length(Dnames);

    for n = 1:nvdims
      dimnam = I.Variables(ind).Dimensions(n).Name;
      switch dimnam
        case {'s_rho'}
          isr3d = true;
        case 's_w'
          isw3d = true;
        case {'xi_rho', 'eta_rho'}
          Mname = 'mask_rho';
          if (G.spherical)
            Xname = 'lon_rho';
            Yname = 'lat_rho';
          else
            Xname = 'x_rho';
            Yname = 'y_rho';
          end
          Zname = 'z_r';
        case {'xi_psi', 'eta_psi'}
          Mname = 'mask_psi';
          if (G.spherical)
            Xname = 'lon_psi';
            Yname = 'lat_psi';
          else
            Xname = 'x_psi';
            Yname = 'y_psi';
          end
          Zname = 'z_r';
        case {'xi_u', 'eta_u'}
          Mname = 'mask_u';
          if (G.spherical)
            Xname = 'lon_u';
            Yname = 'lat_u';
          else
            Xname = 'x_u';
            Yname = 'y_u';
          end
          Zname = 'z_u';
        case {'xi_v', 'eta_v'}
          Mname = 'mask_v';
          if (G.spherical)
            Xname = 'lon_v';
            Yname = 'lat_v';
          else
            Xname = 'x_v';
            Yname = 'y_v';
          end
          Zname = 'z_v';
        case {'ocean_time', 'time', contains(dimnam,'time')}
          recordless = false;
          Tsize = I.Variables(ind).Dimensions(n).Length;
          Tname = dimnam;
      end
    end
    if (isw3d)
      Zname = 'z_w';
    end

    is3d = isr3d || isw3d;

    F.is3d = is3d;

% Check if 'time' attribute is available.

    itime = strcmp({I.Variables(ind).Attributes.Name}, 'time');
    if (any(itime))
      F.Tname = I.Variables(ind).Attributes(itime).Value;
    end

% Get grid coordinates.

    if (isfield(G,Xname))
      if (~isempty(G.(Xname)))
        X = G.(Xname);
      else
        error([' PLOT_STATE - field '', Xname, ''',                     ...
               ' is empty in receiver grid structure: G']);
      end
    else
      error([' PLOT_STATE - unable to find field '', Xname, ''',        ...
             ' in receiver grid structure: G']);
    end

    if (isfield(G,Yname))
      if (~isempty(G.(Yname)))
        Y = G.(Yname);
      else
        error([' PLOT_STATE - field '', Yname, ''',                     ...
               ' is empty in receiver grid structure: G']);
      end
    else
      error([' PLOT_STATE - unable to find field '', Yname, ''',        ...
             ' in receiver grid structure: G']);
    end

    if (is3d)
      if (isfield(G,Zname))
        if (~isempty(G.(Zname)))
          Z = G.(Zname);
        else
          error([' PLOT_STATE - field '', Zname, ''',                   ...
                 ' is empty in receiver grid structure: G']);
        end
      else
        error([' PLOT_STATE - unable to find field '', Zname, ''',      ...
             ' in receiver grid structure: G']);
      end
    end

    if (isfield(G,Mname))
      if (~isempty(G.(Mname)))
        mask = G.(Mname);
      else
        error([' PLOT_STATE - field '', Mname, ''',                     ...
               ' is empty in receiver grid structure: G']);
      end
    else
      error([' PLOT_STATE - unable to find field '', Mname, ''',        ...
             ' in receiver grid structure: G']);
    end

    if (isfield(G,'lon_coast') && isfield(G,'lat_coast'))
      F.lon_coast = G.lon_coast;
      F.lat_coast = G.lat_coast;
      F.gotCoast = true;
    end

    if (~G.spherical)
      X = 0.001 .* X;                 % km
      Y = 0.001 .* Y;                 % km
    end

    F.X = X;
    F.Y = Y;

% Read in and process time variable.

    if (~recordless && Tindex > Tsize)
      Tindex = Tsize;               % process last time record available
    end
    F.Tindex = Tindex;

    if (~isempty(F.Tname))
      Tvalue = nc_read(Sname, F.Tname, Tindex);
      Tattr  = nc_getatt(Sname, 'units', F.Tname);
      if (contains(Tattr, 'second'))
        Tvalue = Tvalue/86400;                  % seconds to days
      end
      iatt = strfind(Tattr, 'since');
      if (~isempty(iatt))
        Torigin = Tattr(iatt+6:end);
        epoch   = datenum(Torigin,31);          % 'yyyy-mm-dd HH:MM:SS'
        Tstring = datestr(epoch+Tvalue);
      else
        Tstring = num2str(Tvalue);
      end
    end
    F.Tstring = Tstring;

% Read in state variable.

    ReplaceValue = NaN;
    PreserveType = false;

    F.Vname = field;
    F.Level = level;

    f = nc_read(Sname, F.Vname, Tindex, ReplaceValue, PreserveType);

    if (is3d)
      Km = size(f, 3);

      if (doSection)

        Z0 = Z;                          % impose zero value at surface
        Z0(:,:,Km) = 0.0;                % for nicer cross-sections

        switch orient
          case 'c'
            V = squeeze(f(index,:,:)); [Im,Km]=size(V);
            m = squeeze(mask(index,:));
            M = repmat(m', [1 Km]);
            V = nanland(V, M);
            s = squeeze(F.Y(index,:));
            Q = repmat(s', [1 Km]);
            Z = squeeze(Z0(index,:,:));
            if (G.spherical)                          % bathymetry
              x = squeeze(G.lat_rho(index,:));
            else
              x = squeeze(G.y_rho(index,:));
            end
            z = -squeeze(G.h(index,:));
            sec_index = ['i=', num2str(index)];
          case 'r'
            V = squeeze(f(:,index,:)); [Im,Km]=size(V);
            m = squeeze(mask(:,index));
            M = repmat(m, [1 Km]);
            V = nanland(V, M);
            s = squeeze(F.X(:,index));
            Q = repmat(s, [1 Km]);
            Z = squeeze(Z0(:,index,:));
            if (G.spherical)                          % bathymetry
              x = squeeze(G.lon_rho(:,index));
            else
              x = squeeze(G.x_rho(:,index));
            end
            z = -squeeze(G.h(:,index));
            sec_index = ['j=', num2str(index)];
        end
        F.X     = Q;
        F.Y     = Z;
        F.value = V;
      else
        if (F.Level > 0)
          value = squeeze(f(:,:,F.Level));
        else
          [~,~,value] = nc_slice(G, Sname, field, F.Level, Tindex);
        end
        F.value = value;
      end
    else
      F.is3d  = false;
      F.Level = [];
      F.value = f;
    end

    F.min      = min(F.value(:));
    F.max      = max(F.value(:));
    F.checkval = bitcount(F.value(:));

    xlabel1    = ['Min = ', sprintf('%12.5e',F.min), blanks(3),         ...
                  'Max = ', sprintf('%12.5e',F.max)];

%--------------------------------------------------------------------------
% Plot state field
%--------------------------------------------------------------------------

% If plotting cross-sections, ignore 2D fields.

    if (~is3d && doSection)
      nfield = nfield -1;
      continue
    end

% Set color map.

    switch (field)
      case 'zeta'
        Cmap = mpl_Accent(256);
%       Cmap = cmap_odv('Odv');;
      case {'u', 'u_eastward'}
        Cmap = mpl_rainbow200;
%       Cmap = cmap_odv('BlueRed_473');
      case {'v', 'v_northward'}
        Cmap = flipud(mpl_rainbow200);
%       Cmap = cmap_odv('BlueRed_473');
      case 'temp'
        Cmap = flipud(mpl_Paired(256));
%       Cmap = cmap_odv('Odv_465');
      case 'salt'
        Cmap = mpl_Set3(256);
%       Cmap = flipud(mpl_gist_ncar(256));
      otherwise
        Cmap = mpl_Accent(256);
    end
    F.Cmap = Cmap;

% Plot field.


   if (is3d && doSection)
     figure;

     if (F.ptype < 0)
       NC = abs(F.ptype);
       [C,F.pltHandle] = contourf(F.X, F.Y, nanland(F.value,G), NC);
     elseif (F.ptype == 2)
       F.pltHandle = pcolorjw(F.X, F.Y, nanland(F.value,G));
     else
       F.pltHandle = pcolor(F.X, F.Y, nanland(F.value,G));
     end
               % bathymetry curve
     hold on;
     area(x, z, min(z), 'FaceColor', Land, 'EdgeColor', Land);
     plot(-128.292,-105.8667,'o','MarkerSize',8,'MarkerEdgeColor','k', ...
          'MarkerFaceColor',[0.8,0.8,0.80]);   % singe obs location
     hold off;

     colorbar;
     colormap(Cmap);
     shading interp;
     if (doZoom)
       axis([-Inf Inf Zmin 0]);
     end
     xlabel({xlabel1, F.Tstring});
     ylabel('Z (m)');
     if (doTitle)
       if (orient == 'c')
         title([untexlabel(field),                                     ...
                ',  section along i = ', num2str(F.index),             ...
                ',  Rec = ', num2st r(F.Tindex)]);
       else
         title([untexlabel(field),                                     ...
                ',  section along j = ', num2str(F.index),             ...
                ',  Rec = ', num2str(F.Tindex)]);
       end
     else
       title(blanks(2));
     end
   else

     F = hplot(G, F);

     if (draw_tiling)
       hold on;
       h = ptiles(F.tiling(1), F.tiling(2), F.ncname, false, 'r-', false);
       hold off
     end

     if (doTitle)
       if (is3d)
         title([untexlabel(field),                                     ...
                ',  Level = ', num2str(F.Level),                       ...
                ',  Rec = ', num2str(F.Tindex)]);
       else
         title([untexlabel(field),                                     ...
                ',  Rec = ', num2str(F.Tindex)]);
       end
     else
       title(blanks(2));
     end
     xlabel({xlabel1, F.Tstring});
   end

   if (doRange)
     caxis(R.(field));
   end

   if (wrtPNG)
     if (isempty(PNGsuffix))
       suffix=num2str(F.Tindex, '%3.3i');
     else
       suffix=PNGsuffix;
     end
     if (doSection)
       png_file=strcat(F.Vname,'_sec_',suffix,'.png');
     else
       png_file=strcat(F.Vname,'_',suffix,'.png');
     end
     exportgraphics(gcf, png_file, 'resolution', 300);
   end

   S(nfield) = F;

  end

end