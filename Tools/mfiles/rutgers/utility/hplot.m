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
%                    F.ptype        0:  'pcolor',
%                                   <0: 'contourf' with abs(ptype) contours
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

% svn $Id: hplot.m 1162 2023-03-27 21:08:44Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
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

switch P.Vname
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

P.min = Pmin;
P.max = Pmax;


if (any(isinf(P.Caxis)))
  P.Caxis = [Pmin Pmax];
end

%--------------------------------------------------------------------------
% Plot requested field.
%--------------------------------------------------------------------------

figure;

if (P.doMap == 1)                               % use m_map toolbox

  LonMin=min(P.X(:));   LonMax=max(P.X(:));
  LatMin=min(P.Y(:));   LatMax=max(P.Y(:));
  m_proj(MyProjection,'longitudes',[LonMin,LonMax],                     ...
                      'latitudes' ,[LatMin,LatMax]); 
  m_grid('tickdir','out','yaxisloc','left');
  hold on;
  if (P.ptype)
    NC = abs(P.ptype);
    [C,H] = m_contourf(P.X, P.Y, nanland(P.value,G), NC);
  else
    H = m_pcolor(P.X, P.Y, nanland(P.value,G));
  end

elseif (P.doMap == 2)                           % Matlab mapping toolbox

  LonLim=[min(P.X(:)) max(P.X(:))];
  LatLim=[min(P.Y(:)) max(P.Y(:))];

  delta=max(fix((LonLim(2)-LonLim(1))/5), fix((LatLim(2)-LatLim(1))/5));
  x = wrapTo180(0:delta:360);
  indx = find(x >= LonLim(1) & x <= LonLim(2)); 
  y = -90:delta:90;
  indy = find(y >= LatLim(1) & y <= LatLim(2)); 
  
  geomap = defaultm(MyProjection);              % create mapping structure

  clf
  hax = axesm(geomap, 'MapLonLimit', LonLim,                            ...
                      'MapLatLimit', LatLim,                            ...
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
  if (P.ptype)
    NC = abs(P.ptype);
    [C,H] = contourfm(P.Y, P.X, nanland(P.value,G), NC);
  else
    H = pcolorm(P.Y, P.X, nanland(P.value,G));
  end   
  
  a = axis;                                     % coastline fill area
  geoshow('landareas.shp', 'FaceColor', Land)
  axis(a)
  
  tightmap

  P.geomap = geomap;                            % pass map structure

  hold on;

else                                            % no map

  if (P.ptype)
    NC = abs(P.ptype);
    [C,H] = contourf(P.X, P.Y, nanland(P.value,G), NC);
  else
    H = pcolorjw(P.X, P.Y, nanland(P.value,G));
  end
  hold on;

  if (P.gotCoast)
    hc = plot(P.lon_coast, P.lat_coast, 'k-');
  end

end

%shading flat
shading interp;
colorbar;
colormap(Cmap);
if (P.Caxis(1) ~= P.Caxis(2)) 
  caxis(P.Caxis);
end

% Color fill land or draw coastlines.

if (P.doMap == 1)
 if (fill_land)
    m_gshhs_i('patch', Land, 'edgecolor', Lake);
  else
    m_gshhs_i('color','k');
  end
end

P.pltHandle = H;

%--------------------------------------------------------------------------
% Set figure title and xlabel.
%--------------------------------------------------------------------------

if (P.is3d)
  if (~isempty(P.Tname))
    ht = title([P.Vname, ':', blanks(4),                                ...
                'Level = ', num2str(P.Level), ',', blanks(4),           ...
                'Record = ', num2str(P.Tindex), ',', blanks(4),         ...
                'time = ', P.Tstring],                                  ...
               'FontSize', 14, 'FontWeight', 'bold' );
  else
    ht = title([P.Vname, ':', blanks(4),                                ...
                'Level = ', num2str(P.Level), ',', blanks(4),           ...
                'Record = ', num2str(P.Tindex)],                        ...
               'FontSize', 14, 'FontWeight', 'bold' );
  end
else
  if (~isempty(P.Tname))
    ht = title([P.Vname, ':', blanks(4),                                ... 
                'Record = ', num2str(P.Tindex), ',', blanks(4),         ...
                'time = ', P.Tstring],                                  ...
               'FontSize', 14, 'FontWeight', 'bold' );
  else
    ht = title([P.Vname, ':', blanks(4),                                ... 
                'Record = ', num2str(P.Tindex)],                        ...
               'FontSize', 14, 'FontWeight', 'bold' );
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

if (P.wrtPNG)
  png_file=strcat(P.Vname,'_',num2str(P.Tindex, '%4.4i'),'.png');
  print(png_file, '-dpng', '-r300');
end

hold off;
