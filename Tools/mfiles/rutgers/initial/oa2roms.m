function F = oa2roms(G,OAname,Vname,varargin)

% OA2ROMS:  Vertically interpolates OA variable to ROMS grid
%
% F = oa2roms(Gname,OAname,Vname,Tindex,method)
%
% This function vertically interpolates requested OA variable to
% ROMS terrain-following vertical grid. The OA package generates
% fields on a constant depth: standard depth levels.
%
% On Input:
%
%    G           ROMS Grid structure containing depth fields
%                  (struct array, see 'get_roms_grid.m')
%
%    OAname      Objective Analysis NetCDF file name (string)
%
%    Vname       ROMS NetCDF variable name to plot (string)
%
%    Tindex      Time record index used to process (Optional, default=1)
%
%    method      Interpolation method (Optional, string):
%
%                  'nearest'       Nearest neighbor
%                  'linear'        Linear (default)
%                  'spline'        Piecewise cubic spline
%                  'cubic'         Cubic
%
% On Output:
%
%    F             Vertically interpolated field (3D array)
%

% svn $Id: oa2roms.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set optional arguments.

Tindex = 1;
method = 'linear';

switch numel(varargin)
  case 1
    Tindex = varargin{1};
  case 2
    Tindex = varargin{1};
    method = varargin{2};
end

% Check grid structure.
  
if (~isstruct(G)),
  error('OA2ROMS: Input ''G'' argument is not a ROMS grid structure.');
else
  if (~isfield(G, 'z_r')),
    error('OA2ROMS: Input ''G'' does not contain depth arrays.');
  end
end
  
% Inquire about the OA NetCDF file.

I = nc_inq(OAname);

% Check requested variable.

Vind = strcmp({I.Variables.Name}, Vname);

if (any(Vind)),
  Vsize = [I.Variables(Vind).Dimensions.Length];
  igrid = I.Variables(Vind).Cgridtype.Value;
else
  error(['OA2ROMS: unable to find ', Vname, ' in:', OAname]);
end

%--------------------------------------------------------------------------
% Get ROMS Grid depth array.
%--------------------------------------------------------------------------

switch (igrid)
  case 1
    Zname = 'z_r'; 
  case 3
    Zname = 'z_u'; 
  case 4
    Zname = 'z_v'; 
  case 5
    Zname = 'z_w'; 
end

if (isfield(G, Zname)),
  Z = G.(Zname);
  Zsize = size(Z);
else
  error(['OA2ROMS: unable to find '',Zname,'' in ROMS Grid structure G']);
end

for n=1:2,
  if (Vsize(n) ~= Zsize(n)),
    error(['OA2ROMS: Dimensions mismatch, ',                            ...
           I.Variables(Vind).Dimensions(n).Name, ' = ',                 ...
           numstr(Vsize(n)), '  ', num2str(Zsize(n))]);
  end
end  

%--------------------------------------------------------------------------
% Vertically interpolate OA field to ROMS terrain-following coordinates.
%--------------------------------------------------------------------------

ReplaceValue = NaN;
PreserveType = false;

Zoa = nc_read(OAname,'zout');

V = nc_read(OAname,Vname,Tindex,ReplaceValue,PreserveType);

[Im,Jm,Km] = size(V);

for j=1:Jm,
  for i=1:Im,
    Voa  = squeeze(V(i,j,:));
    Zwrk = squeeze(Z(i,j,:));
    F(i,j,:) = interp1(Zoa,Voa,Zwrk,method);
  end
end

return
