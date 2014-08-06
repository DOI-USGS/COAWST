function F=plot_nesting(Gnames, Hnames, Vname, Tindex, Level, varargin)

%
% PLOT_NESTING:  Plots requested ROMS variable from nesting NetCDF files
%
% F=plot_nesting(Gnames, Hnames, Vname, Tindex, Level, iplot, Caxis)
%
% This function plots requested ROMS nesting variable from input
% history NetCDF files. This function is very useful when debugging
% a ROMS application. The plotting is not that fancy but it provides
% enough information for browsing ROMS variables very quickly. The
% location of the minumum is marked with a filled magenta circle
% whereas the maximum is marked with a filled magenta square.
%
% It also draws the coastline if the data is available in the
% grid file or structure.
%
% On Input:
%
%    Gnames        ROMS Grid NetCDF files/URL name (cell array)
%              or, an existing ROMS grid structure (struct array)
%
%    Hnames        ROMS history NetCDF files/URL name (cell array)
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
%    iplot         Plotting function to use (scalar, OPTIONAL)
%
%                     iplot = 0,    use PCOLOR
%                     iplot = 1,    use modified PCOLORJW (default)
%
%    Caxis         A two element vector [cmin cmax] to set manual
%                    scaling of pseudocolor for plotting objects
%                    (array, OPTIONAL)
%
% On Output:
%
%    F             Requested 2D or 3D plotting field (structure array)
%
%                    F(n).ncname       NetCDF file
%                    F(n).vname        Variable name
%                    F(n).record       Variable time record number
%                    F(n).time         Variable time
%                    F(n).value        Variable values
%                    F(n).min          Minimum plotted variable value
%                    F(n).max          Maximum plotted variable value
%                    F(n).x            Variable x-coordinates
%                    F(n).y            Variable y-coordinates
%                    F(n).hp           figure pcolor handler
%                    F(n).ht           figure title handler
%                    F(n).hc           figure coastlines handler
  
% svn $Id: plot_nesting.m 735 2014-04-28 23:15:37Z arango $
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
got.Caxis = false;

isr3d = false;
isw3d = false;
isvec = false;

iplot = 1;             % use PCOLORJW to locate correctly C-grid variables

Tname = [];
Tsize = 0;
recordless = true;

switch numel(varargin)
  case 1
    iplot = varargin{1};
  case 2
    iplot = varargin{1};
    Caxis = varargin{2};
    got.Caxis = true;
end

if (~isstruct(Gnames)),
  if (length(Gnames) ~= length(Hnames)),
    error([' PLOT_NESTING - need to provide the same number of files',  ...
           ' in ''Gnames'' and ''Hnames''']);
  end
end

Nfiles = length(Hnames); 

%--------------------------------------------------------------------------
% Get Variable information.
%--------------------------------------------------------------------------

if (~isstruct(Gnames)),
  parent = {'parent_grid',                                              ...
            'parent_Imin', 'parent_Imax',                               ...
            'parent_Jmin', 'parent_Jmax'};
  for n=1:Nfiles,
    g = get_roms_grid(char(Gnames(n)), char(Hnames(n)));
    if (isfield(g, 'parent_grid')),
      G(n) = rmfield(g, parent);
    else
      G(n) = g;
    end
  end
else
  G = Gnames;
end

% Set ROMS Grid structure. If the grid structure have the parent fields,
% remove them to have an array of similar structures.
  

% Get variable information.

I = nc_vinfo(char(Hnames(1)), Vname);

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
      case {'xi_rho','lon_rho'}
        Mname = 'mask_rho';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (G(1).spherical),
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
          if (G(1).spherical),
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
          if (G(1).spherical),
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
          if (G(1).spherical),
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
      case {'ocean_time', 'time'}
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
% Plot requested field.
%--------------------------------------------------------------------------

figure;

for n=1:Nfiles,

% Get coordinates.

  if (isfield(G(n),Xname)),
    if (~isempty(G(n).(Xname)))  
      X = G(n).(Xname);
    else
      error([' PLOT_NESTING - field '', Xname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_NESTING - unable to find field '', Xname, ''',        ...
           ' in receiver grid structure: G']);
  end

  if (isfield(G(n),Yname)),
    if (~isempty(G(n).(Yname)))
      Y = G(n).(Yname);
    else
      error([' PLOT_NESTING - field '', Yname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_NESTING - unable to find field '', Yname, ''',        ...
           ' in receiver grid structure: G']);
  end

  if (is3d),
    if (isfield(G(n),Zname)),
      if (~isempty(G(n).(Zname)))
        Z = G(n).(Zname);
      else
        error([' PLOT_NESTING - field '', Zname, ''',                   ...
               ' is empty in receiver grid structure: G']);
      end
    else
      error([' PLOT_NESTING - unable to find field '', Zname, ''',      ...
             ' in receiver grid structure: G']);
    end
  end
  
  if (isfield(G(n),Mname)),
    if (~isempty(G(n).(Mname)))
      mask = G(n).(Mname);
    else
      error([' PLOT_NESTING - field '', Mname, ''',                     ...
             ' is empty in receiver grid structure: G']);
    end
  else
    error([' PLOT_NESTING - unable to find field '', Mname, ''',        ...
           ' in receiver grid structure: G']);
  end

  if (isfield(G(n),'lon_coast') && isfield(G(n),'lat_coast')),
    Clon = G(n).lon_coast;
    Clat = G(n).lat_coast;
    got.coast = true;
  end

% Read in requested variable (V) from history NetCDF file.

  if (~recordless && Tindex > Tsize),
    Tindex = Tsize;                   % process last time record available
  end 

  if (~isempty(Tname)),
    Tvalue = nc_read(char(Hnames(n)),Tname,Tindex);
    Tvalue = Tvalue/86400;                  % seconds to days
  end

  ReplaceValue = NaN;
  PreserveType = false;

  field = nc_read(char(Hnames(n)),Vname,Tindex,ReplaceValue,PreserveType);

  if (is3d),
    if (Level > 0),
      V = squeeze(field(:,:,Level));
    else
      [~,~,V] = nc_slice(G,Hname,Vname,Level,Tindex);
    end
  else
    V = field;
  end

  Vmin = min(V(:));  [Imin,Jmin] = find(V == min(V(:))); 
  Vmax = max(V(:));  [Imax,Jmax] = find(V == max(V(:))); 

  if (~got.Caxis && n ==1),
    Caxis = [Vmin Vmax];
  end
  
  if (length(Imin) > 1 && length(Jmin > 1)),
    Imin = Imin(1);
    Jmin = Jmin(1);
  end

  if (length(Imax) > 1 && length(Jmax > 1)),
    Imax = Imax(1);
    Jmax = Jmax(1);
  end

  disp([blanks(5), 'Grid = ', num2str(n), blanks(4),                    ...
        'Min = ', num2str(Vmin), blanks(4),                             ...
        '(', num2str(Imin), ', ', num2str(Jmin), '),', blanks(8),       ...
        'Max = ', num2str(Vmax), blanks(4),                             ...
        '(', num2str(Imax), ', ', num2str(Jmax), ')']);
  
% Plot requested field.

  if (iplot == 0),
    hp = pcolor  (X,Y,nanland(V,G(n))); shading faceted; colorbar
  elseif (iplot == 1),
    hp = pcolorjw(X,Y,nanland(V,G(n))); shading faceted; colorbar
  end

  F(n).ncname = char(Hnames(n));
  F(n).vname  = Vname;
  F(n).record = Tindex;
  F(n).time   = Tvalue;
  F(n).values = V;
  F(n).min    = Vmin;
  F(n).max    = Vmax;
  F(n).x      = X;
  F(n).y      = Y;
  F(n).hp     = hp;                         % figure pcolor handler
  
  caxis(Caxis);
  
  if (n == 1),
    if (is3d),
      if (~isempty(Tname)),
        ht = title([Vname, ':', blanks(4),                              ...
                    'Level = ', num2str(Level), ',', blanks(4),         ...
                    'Record = ', num2str(Tindex), ',', blanks(4),       ...
                    'time = ', num2str(Tvalue)],                        ...
                    'FontSize', 14, 'FontWeight', 'bold' );
      else
        ht = title([Vname, ':', blanks(4),                              ...
                    'Level = ', num2str(Level), ',', blanks(4),         ...
                    'Record = ', num2str(Tindex)],                      ...
                    'FontSize', 14, 'FontWeight', 'bold' );
      end
    else
      if (~isempty(Tname)),
        ht = title([Vname, ':', blanks(4),                              ... 
                   'Record = ', num2str(Tindex), ',', blanks(4),        ...
                   'time = ', num2str(Tvalue)],                         ...
                   'FontSize', 14, 'FontWeight', 'bold' );
      else
        ht = title([Vname, ':', blanks(4),                              ... 
                   'Record = ', num2str(Tindex)],                       ...
                   'FontSize', 14, 'FontWeight', 'bold' );
      end
    end

    F(n).ht = ht;                           % figure title handdler

    set (gcf, 'renderer', 'OpenGL');
    alpha (hp, 1);
    
    hold on;

    if (got.coast),
      hc = plot(Clon,Clat,'k-');
      F(n).hc = hc;                          % figure coastline handler
    else
      F(n).hc = [];
    end
  end
end

hold off;

return
