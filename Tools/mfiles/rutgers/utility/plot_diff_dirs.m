function F = plot_diff_dirs (G, dir1, dir2, ncname, vname, rec1, varargin)
  
% PLOT_DIFF:  Plots field difference between two directories
%
% F = plot_diff (G, dir1, dir2, ncname, vname, rec1, rec2, index, orient,
%                Caxis, Mmap, ptype, wrtPNG);  
%
% This function plots the field difference between two directories having
% the same NetCDF filename at the specified time records.
%
% The plot can be a horizontal or vertical section.
%
% On Input:
%
%    G             A existing ROMS grid structure (struct array)
%
%    dir1          1st directory of ROMS NetCDF file (string)
%
%    dir2          2nd directory of ROMS NetCDF file (string)
%
%    ncname        ROMS NetCDF filename (string)
%  
%    Vname         ROMS NetCDF variable name to process (string)
%
%    rec1          1st field/file time record to process (scalar)
%
%    rec2          2nd field/file time record to process (optional; scalar)
%                    (default rec1 = rec2 for comparing files)
%
%    index         horizontal or vertical section index (optional; integer)
%                    if horizontal, then   1 <= index <= N  
%                    if orient='r', then   1 <= index <= Mp
%                    if orient='c', then   1 <= index <= Lp
%
%    orient        Orientation of section extraction (optional; string):
%                    orient='r'  row (west-east) extraction
%                    orient='c'  column (south-north) extraction
%
%    Caxis         Color axis range (optional; vector)
%                    (default: [-Inf Inf])
%
%    Mmap          Switch to use m_map utility (optional; integer)
%                    Mmap = 0,  no map projection (default)
%                    Mmap = 1,  'm_map' utility
%                    Mmap = 2,  native Matlab toolbox
%
%    ptype         Plot type (integer)
%                     ptype < 0     use contourf with abs(ptype) colors
%                     ptype = 1     use 'pcolor'
%                     ptype = 2     use 'pcolorjw'
%
%    wrtPNG        Switch to write out PNG file (true or false)
%                    (Optional, default: false)A
%
% On Output:
%
%    F             Requested variable structure (array)
%

% svn $Id$
%=======================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

% Initialize.

F = struct('dir1'       , [], 'dir2'      , [], 'ncname'      , [],   ...
           'Vname'      , [], 'rec1'      , [], 'rec2'        , [],   ...
           'Tname'      , [], 'Tstring'   , [],                       ...
           'Level'      , [], 'is3d'      , [],                       ...
           'X'          , [], 'Y'         , [],                       ...
           'value1'     , [], 'min1'      , [], 'max1'        , [],   ...
           'checkval1'  , [],                                         ...
           'value2'     , [], 'min2'      , [], 'max2'        , [],   ...
           'checkval2'  , [],                                         ...
           'diff'       , [], 'diffmin'   , [], 'diffmax'     , [],   ...
           'Caxis'      , [], 'doMap'     , [], 'projection'  , [],   ...
           'ptype'      , [], 'orient'    , [], 'index'       , [],   ...
           'gotCoast'   , [], 'lon_coast' , [], 'lat_coast'   , [],   ...
           'shading'    , [], 'pltHandle' , [], 'wrtPNG'      , []);

% Optional arguments.

switch numel(varargin)
  case 0
    rec2      = rec1;
    index     = [];
    orient    = [];
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false; 
  case 1
    rec2      = varargin{1};
    index     = [];
    orient    = [];
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false; 
  case 2
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = [];
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false; 
  case 3
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = varargin{3};
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false; 
  case 4
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = varargin{3};
    Caxis     = varargin{4};
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false; 
  case 5
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = varargin{3};
    Caxis     = varargin{4};
    Mmap      = varargin{5};
    ptype     = 1;
    wrtPNG    = false; 
  case 6
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = varargin{3};
    Caxis     = varargin{4};
    Mmap      = varargin{5};
    ptype     = varargin{6};
    wrtPNG    = false; 
  case 7
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = varargin{3};
    Caxis     = varargin{4};
    Mmap      = varargin{5};
    ptype     = varargin{6};
    wrtPNG    = varargin{7}; 
end

F.ncname1 = strcat(dir1,'/',ncname);
F.ncname2 = strcat(dir2,'/',ncname);
F.Vname   = vname;
F.Tname   = [];
F.Tindex  = rec1;
F.rec1    = rec1;
F.rec2    = rec2;

recordless = true;

% Set spatial coordinates.

I = nc_vinfo(F.ncname1, vname);
nvdims = length(I.Dimensions);

isr3d = false;
isw3d = false;

if (nvdims > 0)
  for n=1:nvdims
    dimnam = char(I.Dimensions(n).Name);
    switch dimnam
      case {'s_rho'}
        isr3d = true;
      case {'s_w'}
        isw3d = true;
      case {'xi_rho','lon_rho'}
        mask = G.mask_rho;
        if (G.spherical)
          F.X = G.lon_rho;
          F.Y = G.lat_rho;
        else
          F.X = G.x_rho./1000;
          F.Y = G.y_rho./1000;
        end
        if (isfield(G, 'z_r'))
          Z = G.z_r;
        end
     case {'xi_u','lon_u'}
        mask = G.mask_u;
        if (G.spherical)
          F.X = G.lon_u;
          F.Y = G.lat_u;
        else
          F.X = G.x_u./1000;
          F.Y = G.y_u./1000;
        end
        if (isfield(G, 'z_u'))
          Z = G.z_u;     
        end
     case {'xi_v','lon_v'}
        mask = G.mask_v;
        if (G.spherical)
          F.X = G.lon_v;
          F.Y = G.lat_v;
        else
          F.X = G.x_v./1000;
          F.Y = G.y_v./1000;
        end
        if (isfield(G, 'z_v'))
          Z = G.z_v;     
        end
     case {'ocean_time', 'time', ~contains(dimnam,'time')}
        recordless = false;    
        Tsize = I.Dimensions(n).Length;
    end
  end
end

is3d = isr3d || isw3d;
if (isw3d)
  Z = G.z_w;
end

%------------------------------------------------------------------------
% Read in requested variable from NetCDF file.
%------------------------------------------------------------------------

% Get time string.

if (~recordless)
  if (rec1 > Tsize || rec2 > Tsize)
    Tindex1 = Tsize;                 % process last time record available
    Tindex2 = Tsize;
  else
    Tindex1 = rec1;
    Tindex2 = rec2;
  end
else
  Tindex1 = [];
  Tindex2 = [];
end 
F.Tindex = Tindex1;

if (~isempty(F.Tname))
  Tvalue = nc_read(F.ncname1, F.Tname, F.Tindex);
  Tattr  = nc_getatt(F.ncname1, 'units', F.Tname);
  Tdays  = true;
  if (~contains(Tattr, 'second'))
    Tvalue = Tvalue/86400;                  % seconds to days
    Tdays  = false;
  end  
  iatt = strfind(Tattr, 'since');
  if (~isempty(iatt))
    iatt=iatt+6;
    Torigin = Tattr(iatt:iatt+18);
    epoch   = datenum(Torigin,31);          % 'yyyy-mm-dd HH:MM:SS' 
    Tstring = datestr(epoch+Tvalue);
  else
    Tstring = num2str(Tvalue);    
  end
  F.Tstring = Tstring;
end

% Process field difference.

F.value1    = nc_read(F.ncname1, vname, Tindex1);
F.min1      = min(F.value1(:));
F.max1      = max(F.value1(:));
F.checkval1 = bitcount(F.value1(:));
F.value2    = nc_read(F.ncname2, vname, Tindex2);
F.min2      = min(F.value2(:));
F.max2      = max(F.value2(:));
F.checkval2 = bitcount(F.value2(:));
F.diff      = F.value1 - F.value2;

F.Caxis     = [-Inf Inf];
F.Cmap      = red_blue(256);
F.doMap     = 0;
F.ptype     = ptype;
F.shading   = 'flat';
F.wrtPNG    = false;

doSection = false;
if ~isempty(orient)
  if (orient == 'c' || orient == 'r')
    doSection = true;
  end
end

if (isfield(G,'lon_coast') && isfield(G,'lat_coast'))
  F.lon_coast = G.lon_coast;
  F.lat_coast = G.lat_coast;
  F.gotCoast = true;
else
  F.gotCoast = false;
end

%------------------------------------------------------------------------
%  Extract horizontal or vertical section of the field difference.
%------------------------------------------------------------------------

if (is3d)
  Km = size(F.diff, 3);

  if (doSection)
    switch orient
      case 'c'
        V = squeeze(F.diff(index,:,:)); [Im,Km]=size(V);
        m = squeeze(mask(index,:));
        M = repmat(m', [1 Km]);
        V = nanland(V, M);
        s = squeeze(F.Y(index,:));
        S = repmat(s', [1 Km]);
        Z = squeeze(Z(index,:,:));
        sec_index = ['i=', num2str(index)];
      case 'r'
        V = squeeze(F.diff(:,index,:)); [Im,Km]=size(V);
        m = squeeze(mask(:,index));
        M = repmat(m, [1 Km]);
        V = nanland(V, M);
        s = squeeze(F.X(:,index));
        S = repmat(s, [1 Km]);
        Z = squeeze(Z(:,index,:));  
        sec_index = ['j=', num2str(index)];
    end    
    F.value = V;  
  else
    F.is3d  = false;
    if (isempty(index))
      Level = Km
    else
      Level = min(Km,index);
    end      
    F.Level = Level;
    for k=1:Km
      Dmin(k)=min(min(F.diff(:,:,k)));
      Dmax(k)=max(max(F.diff(:,:,k)));
    end
    F.diffmin = Dmin;
    F.diffmax = Dmax;
    F.value   = squeeze(F.diff(:,:,Level));
  end  
else
  F.is3d    = false;
  F.Level   = 1;
  F.value   = F.diff;
  F.diffmin = min(F.value(:));
  F.diffmax = max(F.value(:));
end
Dmin = min(F.value(:));
Dmax = max(F.value(:));

%------------------------------------------------------------------------
%  Plot horizontal or vertical section field difference.
%------------------------------------------------------------------------

if (doSection)
  figure;
  pcolorjw(S,Z,F.value); colorbar; shading flat;
  ylabel('Z (m)');
else
  P = hplot(G, F);
end

[~,name,~] = fileparts(ncname);

if (is3d)
  if (doSection)
    title(['File: ', untexlabel(name), blanks(4),                     ...
           'Var = ', vname,                                           ...
           ',  Along ', sec_index,                                    ...
           ',  Rec = ', num2str(rec1),'/', num2str(rec2) ]);
  else
    title(['File: ', untexlabel(name), blanks(4),                      ...
           'Var = ', vname,                                            ...
           ',  Level = ', num2str(Level),                              ...
           ',  Rec = ', num2str(rec1),'/', num2str(rec2) ]);
  end  
else
  title(['File: ', untexlabel(name), blanks(4),                         ...
         'Var = ', vname,                                               ...
         ',  Rec = ', num2str(rec1),'/', num2str(rec2) ]);
end
xlabel(['Min = ', num2str(Dmin), blanks(4), 'Max = ', num2str(Dmax)]);

return
