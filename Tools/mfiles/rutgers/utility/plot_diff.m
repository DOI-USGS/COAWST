function F = plot_diff (G, dir1, dir2, ncname, vname, rec1, varargin)
  
% PLOT_DIFF:  Plots field difference between two records or files
%
% F = plot_diff (G, dir1, dir2, ncname, vname, rec1, rec2, index, orient)
%
% This function plots the field difference between specified time records
% in the same NetCDF file or between specified files with the same record.
% The plot can be a horizontal or vertical section.
%
% On Input:
%
%    G             A existing ROMS grid structure (struct array)
%
%    dir1          1st directory of ROMS NetCDF file (string)
%
%    dir1          2nd directory of ROMS NetCDF file (string)
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
% On Output:
%
%    F             Requested variable structure (array)
%

% svn $Id: plot_diff.m 1162 2023-03-27 21:08:44Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.

doSection = false;

switch numel(varargin)
  case 0
    rec2      = rec1;
    index     = [];
    orient    = [];
  case 1
    rec2      = varargin{1};
    index     = [];
    orient    = [];
  case 2
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = [];
  case 3
    rec2      = varargin{1};
    index     = varargin{2};
    orient    = varargin{3};
    doSection = true;
end

F.ncname1 = strcat(dir1,'/',ncname);
F.ncname2 = strcat(dir2,'/',ncname);
F.Vname   = vname;
F.Tname   = 'ocean_time';
F.Tindex  = rec1;
F.rec1    = rec1;
F.rec2    = rec2;

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

% Process field difference.

F.value1    = nc_read(F.ncname1, vname, rec1);
F.min1      = min(F.value1(:));
F.max1      = max(F.value1(:));
F.checkval1 = bitcount(F.value1(:));
F.value2    = nc_read(F.ncname2, vname, rec2);
F.min2      = min(F.value2(:));
F.max2      = max(F.value2(:));
F.checkval2 = bitcount(F.value2(:));
F.diff      = F.value1 - F.value2;

F.Caxis     = [-Inf Inf];
F.doMap     = 0;
F.ptype     = 0;
F.gotCoast  = 0;
F.wrtPNG    = false;

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
    end
  end
end

is3d = isr3d || isw3d;
if (isw3d)
  Z = G.z_w;
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
  
if (is3d)
  if (doSection)
    title(['File: ', untexlabel(ncname), blanks(4),                       ...
           'Var = ', vname,                                               ...
           ',  Along ', sec_index,                                        ...
           ',  Rec = ', num2str(rec1),'/', num2str(rec2) ]);
  else
    title(['File: ', untexlabel(ncname), blanks(4),                       ...
           'Var = ', vname,                                               ...
           ',  Level = ', num2str(Level),                                 ...
           ',  Rec = ', num2str(rec1),'/', num2str(rec2) ]);
  end  
else
  title(['File: ', untexlabel(ncname), blanks(4),                         ...
         'Var = ', vname,                                                 ...
         ',  Rec = ', num2str(rec1),'/', num2str(rec2) ]);
end
xlabel(['Min = ', num2str(Dmin), blanks(4), 'Max = ', num2str(Dmax)]);

return


