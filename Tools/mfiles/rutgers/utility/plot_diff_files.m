function F = plot_diff_files (G, ncname1, ncname2, vname, rec, varargin)

% PLOT_DIFF:  Plots field difference between two records or files
%
% F=plot_diff_files(G, ncname1, ncname2, vname, rec, index, orient,
%                   Caxis, Mmap, ptype, wrtPNG);
%
% This function plots the field difference between two files at the
% specified time records.
%
% The plot can be a horizontal or vertical section.
%
% On Input:
%
%    G             A existing ROMS grid structure (struct array)
%
%    ncname1       ROMS 1st NetCDF filename (string)
%
%    ncname2       ROMS 2nd NetCDF filename (string)
%
%    Vname         ROMS NetCDF variable name to process (string)
%
%    rec           Field/file time record to process (scalar or vector)
%                    if vector, rec(1) is used in ncname1
%                               rec(2) is used in ncname2
%
%    index         horizontal or vertical section index (optional; integer)
%                    if horizontal, then   1 <= index <= N
%                    if orient='r', then   1 <= index <= Mp
%                    if orient='c', then   1 <= index <= Lp
%
%    orient        Orientation of section extraction (optional; string):
%                    orient='h'  horizontal extraction
%                    orient='r'  row (west-east) vertical extraction
%                    orient='c'  column (south-north) vertical extraction
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
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

F = struct('ncname1'    , [], 'ncname2'   , [], 'Vname'       , [],   ...
           'Tindex'     , [], 'Tname'     , [], 'Tstring'     , [],   ...
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
    index     = [];
    orient    = 'h';
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false;
 case 1
    index     = varargin{1};
    orient    = 'h';
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false;
 case 2
    index     = varargin{1};
    orient    = varargin{2};
    Caxis     = [-Inf Inf];
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false;
 case 3
    index     = varargin{1};
    orient    = varargin{2};
    Caxis     = varargin{3};;
    Mmap      = 0;
    ptype     = 1;
    wrtPNG    = false;
 case 4
    index     = varargin{1};
    orient    = varargin{2};
    Caxis     = varargin{3};
    Mmap      = varargin{4};
    ptype     = 1;
    wrtPNG    = false;
 case 5
    index     = varargin{1};
    orient    = varargin{2};
    Caxis     = varargin{3};
    Mmap      = varargin{4};
    ptype     = varargin{5};
    wrtPNG    = false;
 case 6
    index     = varargin{1};
    orient    = varargin{2};
    Caxis     = varargin{3};
    Mmap      = varargin{4};
    ptype     = varargin{5};
    wrtPNG    = varargin{6};
end

if (orient == 'c' || orient == 'r')
  doSection = true;
else
  doSection = false;
end

F.projection = 'mercator';

F.ncname1 = ncname1;
F.ncname2 = ncname2;
F.Vname   = vname;
F.Tname   = 'ocean_time';
F.Tindex  = rec(1);
F.rec1    = rec(1);
if (length(rec) > 1)
  F.rec2  = rec(2);
else
  F.rec2  = rec(1);
end

% Get time string.

Tvalue = nc_read(F.ncname1, F.Tname, F.Tindex);
Tattr  = nc_getatt(F.ncname1, 'units', F.Tname);
Tdays  = true;
if (~isempty(strfind(Tattr, 'second')))
  Tvalue = Tvalue/86400;                    % seconds to days
  Tdays  = false;
end
iatt = strfind(Tattr, 'since');
if (~isempty(iatt))
  iatt=iatt+6;
  Torigin = Tattr(iatt:iatt+18);
  epoch   = datenum(Torigin,31);            % 'yyyy-mm-dd HH:MM:SS'
  Tstring = datestr(epoch+Tvalue);
else
  Tstring = num2str(Tvalue);
end
F.Tstring = Tstring;

% Get available variables.

V1 = nc_vnames(F.ncname1);
V2 = nc_vnames(F.ncname2);


% Process field difference.

switch (vname)
  case {'Hz', 'hz'}
    if (any(strcmp({V1.Variables(:).Name}, 'Hz')))
      F.value1 = nc_read(F.ncname1, 'Hz', F.rec1);
    else
      F1 = nc_read(F.ncname1, 'z_w', F.rec1);
      Np = size(F1,3); N = Np-1;
      for k = 2:Np
        F.value1(:,:,k-1) = F1(:,:,k) - F1(:,:,k-1);
      end
    end

    if (any(strcmp({V2.Variables(:).Name}, 'Hz')))
      F.value2 = nc_read(F.ncname2, 'Hz', F.rec2);
    else
      F2 = nc_read(F.ncname2, 'z_w', F.rec2);
      Np = size(F2,3); N = Np-1;
      for k = 2:Np
        F.value2(:,:,k-1) = F2(:,:,k) - F2(:,:,k-1);
      end
    end
  otherwise
    ivar=strcmp({V1.Variables(:).Name}, vname);
    if (any(contains({V1.Variables(ivar).Dimensions.Name},'time')))
      F.value1 = nc_read(F.ncname1, vname, F.rec1);
      F.value2 = nc_read(F.ncname2, vname, F.rec2);
    else
      F.value1 = nc_read(F.ncname1, vname);
      F.value2 = nc_read(F.ncname2, vname);
    end
end

F.min1      = min(F.value1(:));
F.max1      = max(F.value1(:));
F.checkval1 = bitcount(F.value1(:));

F.min2      = min(F.value2(:));
F.max2      = max(F.value2(:));
F.checkval2 = bitcount(F.value2(:));

F.diff      = F.value1 - F.value2;

F.Caxis     = Caxis;
F.Cmap      = mpl_sstanom(256);
F.doMap     = Mmap;
F.ptype     = ptype;
F.shading   = 'flat';
F.wrtPNG    = wrtPNG;

if (isfield(G,'lon_coast') && isfield(G,'lat_coast'))
  F.lon_coast = G.lon_coast;
  F.lat_coast = G.lat_coast;
  F.gotCoast = true;
else
  F.gotCoast = false;
end

% Set spatial coordinates.

switch (vname)
  case {'Hz', 'hz'}
    if (any(strcmp({V1.Variables(:).Name}, 'Hz')))
      I = nc_vinfo(F.ncname1, 'Hz');
    else
      I = nc_vinfo(F.ncname1, 'z_w');
    end
  otherwise
    I = nc_vinfo(F.ncname1, vname);
end
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
      case {'xi_psi','lon_psi'}
        mask = G.mask_psi;
        if (G.spherical)
          F.X = G.lon_psi;
          F.Y = G.lat_psi;
        else
          F.X = G.x_psi./1000;
          F.Y = G.y_psi./1000;
        end
        if (isfield(G, 'z_psi'))
          Z = G.z_psi;
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
    end
  end
end

is3d = isr3d || isw3d;
if (isw3d)
  Z = G.z_w;
end

if (~is3d)
  doSection = false;
end

%--------------------------------------------------------------------------
%  Extract horizontal or vertical section of the field difference.
%--------------------------------------------------------------------------

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
      Level = Km;
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

%--------------------------------------------------------------------------
%  Plot horizontal or vertical section field difference.
%--------------------------------------------------------------------------

if (doSection)
  figure;
  pcolorjw(S,Z,F.value); colorbar; shading flat;
  ylabel('Z (m)');
else
  P = hplot(G, F);
end

[~,name1,~] = fileparts(ncname1);
[~,name2,~] = fileparts(ncname2);

if (is3d)
  if (doSection)
    title(['File1: ', untexlabel(name1), blanks(4),                       ...
           'Var = ', vname,                                               ...
           ',  Along ', sec_index,                                        ...
           ',  Rec = ', num2str(F.rec1),'/', num2str(F.rec2) ]);
  else
    title(['File1: ', untexlabel(name1), blanks(4),                       ...
           'Var = ', vname,                                               ...
           ',  Level = ', num2str(Level),                                 ...
           ',  Rec = ', num2str(F.rec1),'/', num2str(F.rec2) ]);
  end
else
  title(['File1: ', untexlabel(name1), blanks(4),                         ...
         'Var = ', vname,                                                 ...
         ',  Rec = ', num2str(F.rec1),'/', num2str(F.rec2) ]);
end


return
