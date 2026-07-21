function P = hplot(G, F)

%
% HPLOT:  Horizontal field pseudocolor or contour plots
%
% P = hplot(G, F)
%
% This function plots requested horizontal field. The input field structure
% is usually created in 'plot_field'.
%
% On Input:
%
%    G             Grid structure  (struct)
%
%    F             Field structure (struct)
%
%                    F.Vname        field variable name
%                    F.Tindex       field time record
%                    F.Tname        field time variable name
%                    F.Tstring      field time string
%                    F.Level        vertical level, if 3D field
%                    F.X            X-locations
%                    F.Y            Y-localitons
%                    F.value        field values
%                    F.Caxis        color axis
%                    F.doMap        0: no, 1: m_map, 2: mapping toolbox
%                    F.ptype        <0: 'contourf' with abs(ptype) contours
%                                   1:  'pcolor'
%                                   2:  'pcolorjw'
%                    F.projection   map projection
%                    F.gotCoast     coastline switch
%                    F.lon_coast    coastline longitude
%                    F.lat_coast    coastline latitude
%                    F.wrtPNG       create PNG file
%
% On Output:
%
%    P               Updated field F-structure (struct)
%
%                    P.min          field minimum value (land excluded)
%                    P.max          field maximum value (land excluded)
%                    P.Caxis        [F.min, F.max], if not provided
%                    P.geomap       Matlab's map toolbox structure, doMap=2
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Set output P structure, copy from input F structure.
%--------------------------------------------------------------------------

P = F;

%--------------------------------------------------------------------------
% Initialize.
%--------------------------------------------------------------------------

if (isfield(P, 'projection') & P.doMap > 0)
  if (~isempty(P.projection))
    MyProjection = P.projection;
  else
    MyProjection = 'mercator';
    P.projection = MyProjection;
  end
end

fill_land = true;                 % use patch (true) or draw coast (false)

Marks = false;                    % draw min/max marks

%Land = [0.3412 0.2549 0.1843];   % dark brown
%Land = [0.4745 0.3765 0.2980];   % medium brown
%Land = [0.6706 0.5841 0.5176];   % light brown
 Land = [0.6 0.65 0.6];           % gray-green
 Lake = Land;

%Cedge = 'none';                  % MarkerEdgeColor
 Cedge = 'y';

%--------------------------------------------------------------------------
% Set colormap.
%--------------------------------------------------------------------------

if (isfield(P, 'Cmap'))
  Cmap = P.Cmap;
else
  switch P.Vname
    case {'temp', 'temp_sur', 'temp_slice'}
      Cmap = flipud(mpl_Paired(256));
    case {'salt', 'salt_sur', 'salt_slice'}
      Cmap = mpl_Set3(256);
    case {'pvorticity', 'pvorticity_slice'}
      Cmap = mpl_Set1(256);
    case {'rvorticity', 'rvorticity_slice'}
      Cmap = cmap_odv('BlueRed_473');
%     Cmap = cmap_odv('Ferret_blue_orange',256);
    case {'u', 'u_sur', 'u_slice', 'u_eastward', 'u_eastward_slice'}
      Cmap = mpl_rainbow200;
    case {'v', 'v_sur', 'v_slice', 'v_northward', 'v_northward_slice'}
      Cmap = flipud(mpl_rainbow200);
    case {'ubar', 'ubar_eastward', 'vbar', 'vbar_northward'}
      Cmap = cm_speed(512);
    case 'zeta'
      Cmap = mpl_Accent(256);
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
      Cmap = cmap('R1');
  end
end

f = nanland(P.value, G);

Pmin = min(f(:));  [Imin,Jmin] = find(f == min(f(:)));
Pmax = max(f(:));  [Imax,Jmax] = find(f == max(f(:)));

if (length(Imin) > 1 && ~isempty(Jmin > 1))
  Imin = Imin(1);
  Jmin = Jmin(1);
end

if (length(Imax) > 1 && ~isempty(Jmax > 1))
  Imax = Imax(1);
  Jmax = Jmax(1);
end

%--------------------------------------------------------------------------
% Plot requested field.
%--------------------------------------------------------------------------

Xmin = min(P.X(:));
Xmax = max(P.X(:));
Ymin = min(P.Y(:));
Ymax = max(P.Y(:));

figure;

if (P.doMap == 1)                               % use m_map toolbox

  m_proj(MyProjection,'longitudes',[Xmin,Xmax],'latitudes',[Ymin,Ymax]);
  m_grid('tickdir','out','yaxisloc','left');
  hold on;
  if (P.ptype < 0)
    NC = abs(P.ptype);
    [C,H] = m_contourf(P.X, P.Y, f, NC);
  else
    H = m_pcolor(P.X, P.Y, f);
  end
elseif (P.doMap == 2)                           % Matlab mapping toolbox

  delta=max(fix((Xmax-Xmin)/5), fix((Ymax-Ymin)/5));
  x = wrapTo180(0:delta:360);
  indx = find(x >= Xmin & x <= Xmax);
  y = -90:delta:90;
  indy = find(y >= Ymin & y <= Ymax);

  geomap = defaultm(MyProjection);              % create mapping structure

  clf
  hax = axesm(geomap, 'MapLonLimit', [Xmin Xmax],                       ...
                      'MapLatLimit', [Ymin Ymax],                       ...
                      'Grid','on',                                      ...
                      'MLineLocation', delta,                           ...
                      'PLineLocation', delta,                           ...
                      'MeridianLabel', 'on',                            ...
                      'ParallelLabel', 'on',                            ...
                      'MlabelParallel', 'south',                        ...
                      'PLabelMeridian', 'west',                         ...
                      'MLabelLocation', x(indx),                        ...
                      'PLabelLocation', y(indy));

  geomap = defaultm(geomap);                    % adjust map structure
                                                % set empty properties
  if (P.ptype < 0)
    NC = abs(P.ptype);
    [C,H] = contourfm(P.Y, P.X, f, NC);
  else
    H = pcolorm(P.Y, P.X, f);
  end

  a = axis;                                     % coastline fill area
  geoshow('landareas.shp', 'FaceColor', Land)
  axis(a)

  tightmap

  P.geomap = geomap;                            % pass map structure

  hold on;

else                                            % no map

  if (P.ptype < 0)
    NC = abs(P.ptype);
    [C,H] = contourf(P.X, P.Y, f, NC);
  elseif (P.ptype == 2)
    H = pcolorjw(P.X, P.Y, f);
  else
    H = pcolor(P.X, P.Y, f);
  end
  hold on;

  if (P.gotCoast)
    hc = plot(P.lon_coast, P.lat_coast, 'k-');
    axis([Xmin Xmax Ymin Ymax]);
  end

end

if (isfield(P, 'shading'))
  switch (P.shading)
    case 'faceted'
      shading faceted;
    case 'flat'
      shading flat;
    case 'interp'
      shading interp;
  end
else
  shading flat
end
colorbar;
colormap(Cmap);

if (isfield(P, 'Caxis'))
  caxis(P.Caxis);
end

% Color fill land or draw coastlines.

if (P.doMap == 1)
 if (fill_land)
%   m_coast('patch', Land);
%   m_gshhs_i('patch', Land, 'edgecolor', Lake);
    m_gshhs_h('patch', Land, 'edgecolor', 'none');
 else
    m_gshhs_i('color','k');
  end
% [x,y]=m_ll2xy(-128.292,37.918);     % WC13 T,S observation
% plot(x, y, 'o','MarkerSize',8,                                       ...
%      'MarkerEdgeColor', 'r', 'MarkerFaceColor',[0.8,0.8,0.80]);
% x=[-134, -122.5];
% y=[37.666 37.666];
% [x,y]=m_ll2xy(x,y);                 % WC13 cross-section
% plot(x,y,'r:');
end

P.pltHandle = H;

%--------------------------------------------------------------------------
% Set figure title and xlabel.
%--------------------------------------------------------------------------

doTitle = true;
if (P.wrtPNG < 0)
  doTitle = false;
end

if (doTitle)
  if (P.is3d)
    if (~isempty(P.Tname))
      ht = title([untexlabel(P.Vname), ':', blanks(4),                  ...
                 'Level = ', num2str(P.Level), ',', blanks(4),          ...
                 'Record = ', num2str(P.Tindex), ',', blanks(4),        ...
                 'time = ', P.Tstring],                                 ...
                 'FontSize', 14, 'FontWeight', 'bold' );
    else
      ht = title([untexlabel(P.Vname), ':', blanks(4),                  ...
                 'Level = ', num2str(P.Level), ',', blanks(4),          ...
                 'Record = ', num2str(P.Tindex)],                       ...
                 'FontSize', 14, 'FontWeight', 'bold' );
    end
  else
    if (~isempty(P.Tname))
       ht = title([untexlabel(P.Vname), ':', blanks(4),                 ...
                  'Record = ', num2str(P.Tindex), ',', blanks(4),       ...
                  'time = ', P.Tstring],                                ...
                  'FontSize', 14, 'FontWeight', 'bold' );
    else
      ht = title([untexlabel(P.Vname), ':', blanks(4),                  ...
                 'Record = ', num2str(P.Tindex)],                       ...
                 'FontSize', 14, 'FontWeight', 'bold' );
    end
  end
end

hx = xlabel(['Min = ', num2str(Pmin), blanks(4),                        ...
             '(', num2str(Imin), ', ', num2str(Jmin), '),', blanks(8),  ...
             'Max = ', num2str(Pmax), blanks(4),                        ...
             '(', num2str(Imax), ', ', num2str(Jmax), ')'],             ...
            'FontSize', 14, 'FontWeight', 'bold' );

%--------------------------------------------------------------------------
%  Mark minimum (down triangle) and maximum (up triangle) locations.
%--------------------------------------------------------------------------

if (P.doMap == 1)                               % using m_map toolbox

  if (Marks)
    m_plot(P.X(Imin,Jmin), P.Y(Imin,Jmin), 'v',                         ...
           'MarkerEdgeColor',Cedge,                                     ...
           'MarkerFaceColor','k','MarkerSize',12);
    m_plot(P.X(Imax,Jmax), P.Y(Imax,Jmax), '^',                         ...
           'MarkerEdgeColor',Cedge,                                     ...
           'MarkerFaceColor','k','MarkerSize',12);
    m_plot(P.X(Imin,Jmin), P.Y(Imin,Jmin), 'kv', 'MarkerSize',12);
    m_plot(P.X(Imax,Jmax), P.Y(Imax,Jmax), 'k^', 'MarkerSize',12);
  end

elseif (P.doMap == 2)                           % using mapping  toolbox

  if (Marks)
    plotm(P.Y(Imin,Jmin), P.X(Imin,Jmin), 'v',                          ...
          'MarkerEdgeColor',Cedge,                                      ...
          'MarkerFaceColor','k','MarkerSize',12);
    plotm(P.Y(Imax,Jmax), P.X(Imax,Jmax), '^',                          ...
          'MarkerEdgeColor',Cedge,                                      ...
          'MarkerFaceColor','k','MarkerSize',12);
    plotm(P.Y(Imin,Jmin), P.X(Imin,Jmin), 'kv', 'MarkerSize',12);
    plotm(P.Y(Imax,Jmax), P.X(Imax,Jmax), 'k^', 'MarkerSize',12);
  end

else

  if (Marks)
    plot(P.X(Imin,Jmin), P.Y(Imin,Jmin), 'v',                           ...
         P.X(Imax,Jmax), P.Y(Imax,Jmax), '^',                           ...
         'MarkerEdgeColor',Cedge,'MarkerFaceColor','k',                 ...
         'MarkerSize',12);
    plot(P.X(Imin,Jmin), P.Y(Imin,Jmin), 'wv',                          ...
         P.X(Imax,Jmax), P.Y(Imax,Jmax), 'w^',                          ...
         'MarkerSize',12);
  end

end

%--------------------------------------------------------------------------
%  Write out PNG file.
%--------------------------------------------------------------------------

if (abs(P.wrtPNG))
  png_file=strcat(P.Vname,'_',num2str(P.Tindex, '%4.4i'),'.png');
  exportgraphics(gcf, png_file, 'resolution', 300);
end

hold off;
