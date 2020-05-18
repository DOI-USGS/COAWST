function V = roms2roms(ncname,D,R,Vname,Tindex,Rvector,varargin)

%
% ROMS2ROMS: Interpolates requested ROMS variable to specified ROMS grid
%
% V = roms2roms(ncname,D,R,Vname,Tindex,Rvector,Hmethod,offset,RemoveNaN);
%
% This function interpolates requested 2D/3D variable between two ROMS
% application grids. The receiver grid must be inside of the donor grid.
%
% This function is intended for down-scaling or nesting applications.
% The horizontal/vertical coordinates for the donor and the receiver
% grids are specified with array structures 'D' and 'R', which are
% builded elsewhere using the script 'get_roms_grid.m' for efficiency
% and functionality.  It uses 'TriScatteredInterp' for interpolating
% a 2D or 3D variable.
% 
% On Input:
%
%    ncname        Donor NetCDF file/URL name (string) containing
%                    variable to process
%
%    D             Donor grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    R             Receiver grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    Vname         Field variable name to process (string)
%
%    Tindex        Time record index to process (scalar).
%
%    Rvector       Switch to Interpolate U- and V-points variables
%                    to RHO-points to facilitate rotation in
%                    in curvilinear grid applications somewhere
%                    else (logical)
%
%    Hmethod       Interpolation method for 'TriScatteredInterp'
%                    (string):
%
%                    'natural'     natural neighbor interpolation
%                    'linear'      linear interpolation (default)
%                    'nearest'     nearest-neighbor interpolation
%
%    offset        Number of extra points to use when sampling the
%                    donor grid so it is large enough to contain
%                    the receiver grid  (default 5)
%
%    RemoveNaN     Switch to remove NaN values from the interpolated 
%                    variable with a second interpolation step
%                    using the nearest-neighbor method
%                    (default false)
%
% On Output:
%
%    V             Interpolated requested 2D/3D variable
%
 
% svn $Id: roms2roms.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set optional arguments.

Hmethod = 'linear';
offset = 5;
RemoveNaN = false;

switch numel(varargin)
  case 1
    Hmethod = varargin{1};
  case 2
    Hmethod = varargin{1};
    offset = varargin{2};
  case 3
    Hmethod = varargin{1};
    offset = varargin{2};
    RemoveNaN = varargin{3};
end

%  Set vertical interpolation method.

Vmethod = 'linear';               

%  Initialize.

got.Mname = false;
got.Xname = false;
got.Yname = false;
got.Zname = false;

isr3d = false;
isw3d = false;
isvec = false;
recordless = true;

V = [];

%  Get information about requested variable.
 
Info = nc_vinfo(ncname,Vname);

nvdims = length(Info.Dimensions);

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

if (nvdims > 0),
  for n=1:nvdims,
    dimnam = char(Info.Dimensions(n).Name);
    switch dimnam
      case 's_rho'
        isr3d = true;
      case 's_w'
        isw3d = true;
      case {'xi_rho','eta_rho'}
        Mname = 'mask_rho';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (D.spherical),
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
      case {'xi_psi','eta_psi'}
        Mname = 'mask_psi';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (D.spherical),
            Xname = 'lon_psi';
            Yname = 'lat_psi';
          else
            Xname = 'x_psi';
            Yname = 'y_psi';
          end
          got.Xname = true;
          got.Yname = true;       
        end
      case {'xi_u','eta_u'}
        Mname = 'mask_u';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (D.spherical),
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
     case {'xi_v','eta_v'}
        Mname = 'mask_v';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (D.spherical),
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
      case 'ocean_time'
        recordless = false;
    end
  end
  if (isw3d),
    Zname = 'z_w';
  end  
end

is3d = isr3d || isw3d;

%--------------------------------------------------------------------------
%  Get horizontal and vertical coordinates from donor and receiver grids.
%--------------------------------------------------------------------------

%  Donor grid.

if (isfield(D,Xname)),
  if (~isempty(D.(Xname)))
    XD = D.(Xname);
  else
    error([' ROMS2ROMS - field '', Xname, ''',                          ...
           ' is empty in donor grid structure: D']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Xname, ''',             ...
         ' in donor grid structure: D']);
end

if (isfield(D,Yname)),
  if (~isempty(D.(Yname)))
    YD = D.(Yname);
  else
    error([' ROMS2ROMS - field '', Yname, ''',                          ...
           ' is empty in donor grid structure: D']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Yname, ''',             ...
         ' in donor grid structure: D']);
end

if (is3d),
  if (isfield(D,Zname)),
    if (~isempty(D.(Zname)))
      ZD = D.(Zname);
    else
      error([' ROMS2ROMS - field '', Zname, ''',                        ...
             ' is empty in donor grid structure: D']);
    end
  else
    error([' ROMS2ROMS - unable to find field '', Zname, ''',           ...
           ' in donor grid structure: D']);
  end
end

if (isfield(D,Mname)),
  if (~isempty(D.(Mname)))
    Dmask = D.(Mname);
  else
    error([' ROMS2ROMS - field '', Mname, ''',                          ...
           ' is empty in donor grid structure: D']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Mname, ''',             ...
         ' in donor grid structure: D']);
end

%  Receiver grid.

if (isvec && Rvector),
  if (D.spherical),
    Xname = 'lon_rho';           % If requested, interpolate
    Yname = 'lat_rho';           % U- and V-points variables
  else                           % to RHO-points instead to
    Xname = 'x_rho';             % facilitate the curvilinear
    Yname = 'y_rho';             % rotation to true East and
  end                            % North.  This rotation is
  Mname = 'mask_rho';            % done somewhere else.
  if (is3d),
    Zname = 'z_r';
  end
end

if (isfield(R,Xname)),
  if (~isempty(R.(Xname)))  
    XR = R.(Xname);
  else
    error([' ROMS2ROMS - field '', Xname, ''',                          ...
           ' is empty in receiver grid structure: R']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Xname, ''',             ...
         ' in receiver grid structure: R']);
end

if (isfield(R,Yname)),
  if (~isempty(R.(Yname)))
    YR = R.(Yname);
  else
    error([' ROMS2ROMS - field '', Yname, ''',                          ...
           ' is empty in receiver grid structure: R']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Yname, ''',             ...
         ' in receiver grid structure: R']);
end

if (is3d),
  if (isfield(R,Zname)),
    if (~isempty(R.(Zname)))
      ZR = R.(Zname);
    else
      error([' ROMS2ROMS - field '', Zname, ''',                        ...
             ' is empty in receiver grid structure: R']);
    end
  else
    error([' ROMS2ROMS - unable to find field '', Zname, ''',           ...
           ' in receiver grid structure: R']);
  end
end
  
if (isfield(R,Mname)),
  if (~isempty(R.(Mname)))
    Rmask = R.(Mname);
  else
    error([' ROMS2ROMS - field '', Mname, ''',                          ...
           ' is empty in receiver grid structure: R']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Mname, ''',             ...
         ' in receiver grid structure: R']);
end

%--------------------------------------------------------------------------
%  Read in requested variable from donor NetCDF file.
%--------------------------------------------------------------------------

ReplaceValue = NaN;
PreserveType = false;

VD = nc_read(ncname,Vname,Tindex,ReplaceValue,PreserveType);

%--------------------------------------------------------------------------
%  Set donor grid sampling indices to accelerate the interpolation.
%  The horizontally sampled donor grid is large enough to contain
%  the receiver grid. The parameter 'offset' is used to add extra
%  points when computing the sampling indices (Istr:Iend,Jstr:Jend).
%  That is, the sampled grid is 'offset' points larger in all sides.
%  This is done to resolve well the interpolation near the boundaries
%  of the receiver grid.
%--------------------------------------------------------------------------

[Istr,Iend,Jstr,Jend] = sample_grid(XD,YD,XR,YR,offset);

%--------------------------------------------------------------------------
%  Interpolate requested variable to receiver grid: Build interpolation
%  data structure, I.
%--------------------------------------------------------------------------

I.Vname = Vname;   

%  Determine if interpolating 2D or 3D variables.  Notice that it is
%  possible to interpolate recordless 2D or 3D variables.  That is,
%  variables with no time dimension.

if (recordless),
  I.nvdims = nvdims;
else
  I.nvdims = nvdims-1;
end

%  Set interpolation data.

switch (I.nvdims),

 case 2

   [ImD,JmD]=size(XD);
   [ImR,JmR]=size(XR);

   disp(' ');
   disp(['Interpolating 2D variable: ', Vname,                          ...
          ' (', num2str(ImR), 'x', num2str(JmR),') from donor ',        ...
           '(', num2str(ImD), 'x', num2str(JmD),') ...']);
   disp(' ');
   
   I.VD    = VD(Istr:1:Iend,Jstr:1:Jend);  

   I.Dmask = Dmask(Istr:1:Iend,Jstr:1:Jend);
   I.XD    = XD(Istr:1:Iend,Jstr:1:Jend);
   I.YD    = YD(Istr:1:Iend,Jstr:1:Jend);
   I.ZD    = [];
   
   I.Rmask = Rmask;
   I.XR    = XR;
   I.YR    = YR;
   I.ZR    = [];

   I.Zsur  = [];
   I.Zbot  = [];

 case 3

   [ImD,JmD,KmD]=size(ZD);
   [ImR,JmR,KmR]=size(ZR);

   disp(' ');
   disp(['Interpolating 3D variable: ', Vname,                          ...
          ' (', num2str(ImR), 'x', num2str(JmR), 'x',                   ...
	        num2str(KmR), ') from donor ',                          ...
           '(', num2str(ImD), 'x', num2str(JmD), 'x',                   ...
	        num2str(KmD), ') ...']);
   disp(' ');
   
   I.VD    = VD(Istr:1:Iend,Jstr:1:Jend,:);

   I.Dmask = Dmask(Istr:1:Iend,Jstr:1:Jend);
   I.XD    = XD(Istr:1:Iend,Jstr:1:Jend);
   I.YD    = YD(Istr:1:Iend,Jstr:1:Jend);
   I.ZD    = ZD(Istr:1:Iend,Jstr:1:Jend,:);

   I.Rmask = Rmask;
   I.XR    = XR;
   I.YR    = YR;
   I.ZR    = ZR;

   I.Zsur  = max(R.z_w(:))+eps;
   I.Zbot  = min(R.z_w(:))-eps;
end

%  Interplate requested variable.

V = interp_field(I,Hmethod,Vmethod,RemoveNaN);

return
