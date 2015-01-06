function F = coarse2fine(Ginp,Gout,Gfactor,varargin)

%
% COARSE2FINE:  Creates a finer resolution ROMS Grid NetCDF file
%
% F = coarse2fine(Ginp,Gout,Gfactor,Imin,Imax,Jmin,Jmax)
%
% Given a coarse resolution Grid NetCDF file (Ginp), this function creates
% a finer resolution grid in the region specified by the coarser grid
% coordinates (Imin,Jmin) and (Imax,Jmax). Notice that (Imin,Jmin), and
% (Imax,Jmax) indices are in terms of the PSI-points because it actually
% define the physical boundaries of the refined grid. The grid refinement
% coefficient is specified with Gfactor.
%
% On Input:
%
%    Ginp       Input  coaser Grid NetCDF file name (character string)
%    Gout       Output finer  Grid NetCDF file name (character string)
%    Gfactor    Grid refinement factor (3,5,7,9,11,13,15,...)
%    Imin       Coarse grid lower-left  I-coordinate (PSI-point)
%    Imax       Coarse grid upper-right I-coordinate (PSI-point)
%    Jmin       Coarse grid lower-left  J-coordinate (PSI-point)
%    Jmax       Coarse grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    F          Fine resolution Grid structure
%

% svn $Id: coarse2fine.m 738 2014-10-14 21:49:14Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Check refinement factor. The values below (3:2:27) are somewhat
% rediculous for an application in ROMS. However, these values are
% possible if the strategy is extract the final grids to be used
% in ROMS from an intermediate very fine resolution grid.

legal = abs([3:2:27] - Gfactor) < 4*eps;
if (~any(legal)),
  error([' COARSE2FINE: illegal refinement factor, Gfactor = ',         ...
         num2str(Gfactor)]);
end  
  
% Set TriScatteredInterp method.  Use 'natural' for Natural Neighbor
% Interpolation Method.  This does NOT implies Nearest Neighbor
% Interpolation.  Natural Neighbor Interpolation is a triangulantion
% based method that has an area-of-influence weighting associated
% with each sample data point.  It is C1 continous except at the
% sample locations.
  
method = 'natural';

% Get coarse grid structure.

C = get_roms_grid(Ginp);

% Set coarse grid refinement PSI-indices.

[Lp,Mp] = size(C.h);

L = Lp-1;
M = Mp-1;

% Set offset "half" factor which is used to extract boundary points
% outside of the PSI-points perimeter defined by Imin, Imax, Jmin,
% and Jmax.  This also takes into account the exterior contact points.

half = floor(Gfactor-1)/2;

Imin = half;                   % Default PSI-points range to apply
Imax = Lp-half;                % refinement to the full input grid
Jmin = half;
Jmax = Mp-half;

switch numel(varargin)
  case 1
    Imin = varargin{1};
  case 2
    Imin = varargin{1};
    Imax = varargin{2};
  case 3
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
  case 4
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
    Jmax = varargin{4};
end

% Set spherical switch.

spherical = C.spherical;

if (spherical),
  F.spherical = 1;
else
  F.spherical = 0;
end

% Set uniform grid switch.

if isfield(C, 'uniform'),
  F.uniform = C.uniform;
else
  F.uniform = 0;
end

% Set curvilinear switch.

curvilinear = C.curvilinear;

if (spherical),
  F.curvilinear = 1;
else
  F.curvilinear = 0;
end

% Get grid lengths.

F.xl = 0;
F.el = 0;

% Check available grid variables. Get input file information
% structure, I.

I = nc_inq(Ginp);
vnames = {I.Variables.Name};

got_list = {'lon_rho'  , 'lat_rho'  , 'lon_psi'  , 'lat_psi'  ,         ...
            'lon_u'    , 'lat_u'    , 'lon_v'    , 'lat_v'    ,         ...
            'x_rho'    , 'y_rho'    , 'x_psi'    , 'y_psi'    ,         ...
            'x_u'      , 'y_u'      , 'x_v'      , 'y_v'      ,         ...
            'hraw'     , 'angle'    ,                                   ...
            'mask_rho' , 'mask_psi' , 'mask_u'   , 'mask_v'   ,         ...
            'lon_coast', 'lat_coast'};

for value = got_list,
  field = char(value);
  got.(field) =  any(strcmp(vnames, field));
end
	    
% Set fields to process.

field_list = {'f', 'h', 'pm', 'pn'};

if (curvilinear),
  field_list = [field_list, 'dmde', 'dndx'];
end

if (got.angle),
  field_list = [field_list, 'angle'];
end

if (got.x_rho && got.y_rho),
  field_list = [field_list, 'x_rho', 'y_rho'];
end

if (got.x_psi && got.y_psi),
  field_list = [field_list, 'x_psi', 'y_psi'];
end

if (got.x_u && got.y_u),
  field_list = [field_list, 'x_u', 'y_u'];
end

if (got.x_v && got.y_v),
  field_list = [field_list, 'x_v', 'y_v'];
end

if (spherical),
  field_list = [field_list, 'lon_rho', 'lat_rho',                       ...
                            'lon_psi', 'lat_psi',                       ...
                            'lon_u',   'lat_u',                         ...
                            'lon_v',   'lat_v'];
end

if (got.mask_rho || got.mask_psi || got.mask_u || got.mask_v),
  field_list = [field_list, 'mask_rho', 'mask_psi',                     ...
                            'mask_u',   'mask_v'];
end

% Set coarser grid fractional coordinates.

[Lp,Mp] = size(C.h);          % RHO-points

L = Lp-1;
M = Mp-1;

[YrC,XrC] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
[YpC,XpC] = meshgrid(1.0:1:M     , 1.0:1:L     );
[YuC,XuC] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );
[YvC,XvC] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);

C.Lr = Lp;
C.Mr = Mp;

% Check refinement region.

if (Imin >= Imax),
  error([' COARSE2FINE: Imin >= Imax,    ',                             ...
         ' Imin = ', num2str(Imin), ',',                                ...
         ' Imax = ', num2str(Imax)]);
end

if (Jmin >= Jmax),
  error([' COARSE2FINE: Jmin >= Jmax,    ',                             ...
         ' Jmin = ', num2str(Jmin), ','                                 ...
         ' Jmax = ', num2str(Jmax)]);
end

if (Imin < 1),
  error([' COARSE2FINE: Imin < 1,   ',                                  ...
         ' Imin = ', num2str(Imax)]);
end

if (Imax > L),
  error([' COARSE2FINE: Imax > L,    ',                                 ...
         ' Imax = ', num2str(Imax), ',',                                ...
         ' L = ', num2str(Lp-half)]);
end

if (Jmin < 1),
  error([' COARSE2FINE: Jmin < 1,   ',                                  ...
         ' Jmin = ', num2str(Imax)]);
end

if (Jmax > M),
  error([' COARSE2FINE: Jmax > M,    ',                                 ...
         ' Jmax = ', num2str(Jmax), ',',                                ...
         ' M = ', num2str(Mp-half)]);
end

% Set finer grid fractional coordinates.

delta = 1.0/Gfactor;
halfdelta  = 0.5*delta;

IpF = (Imin:delta:Imax);                           % PSI-points
JpF = (Jmin:delta:Jmax);                           % PSI-points
IrF = [IpF(1)-halfdelta IpF+halfdelta];            % RHO-points
JrF = [JpF(1)-halfdelta JpF+halfdelta];            % RHO-points

[YrF,XrF] = meshgrid(JrF,IrF);                     % RHO-points
[YpF,XpF] = meshgrid(JpF,IpF);                     % PSI-points
[YuF,XuF] = meshgrid(JrF,IpF);                     % U-points
[YvF,XvF] = meshgrid(JpF,IrF);                     % V-points

F.Lr = 0;
F.Mr = 0;

%--------------------------------------------------------------------------
% Interpolate grid variables.
%--------------------------------------------------------------------------
%
% The default strategy is to interpolate the Cartesian coordinates
% (x_rho, y_rho, ...) in meters using the nondimensional fractional
% coordinates (XrC, YrC, ...).  Then, the other variables are
% interpolated in term of the Cartesian coordinates between
% coarse and finer grids.
%
% However, some application grid NetCDF files may not have the
% Cartesian coordinates (x_rho, y_rho, ...) variables. Then, the
% coarse to fine interpolation is done in terms of spherical
% coordinates (lon_rho, lat_rho, ...).
%

disp(' ');
disp('Interpolating from coarse to fine ...');

Rutgers=0; COAWST=1;
if (Rutgers)

% Grid locations at RHO-points.

if (got.x_rho && got.y_rho),

  RCr = TriScatteredInterp(XrC(:),YrC(:),C.x_rho(:),method);

                        F.x_rho = RCr(XrF, YrF);
  RCr.V = C.y_rho(:);   F.y_rho = RCr(XrF, YrF);

  if (got.lon_rho && got.lat_rho),
    RSr = TriScatteredInterp(C.x_rho(:),C.y_rho(:),C.lon_rho(:),method);

                            F.lon_rho = RSr(F.x_rho, F.y_rho);
    RSr.V = C.lat_rho(:);   F.lat_rho = RSr(F.x_rho, F.y_rho);
  end
  
elseif (got.lon_rho && got.lat_rho),

  RSr = TriScatteredInterp(XrC(:),YrC(:),C.lon_rho(:),method);

                          F.lon_rho = RSr(XrF, YrF);
  RSr.V = C.lat_rho(:);   F.lat_rho = RSr(XrF, YrF);
  
  if (got.x_rho && got.x_rho),
    RCr = TriScatteredInterp(C.lon_rho(:),C.lat_rho(:),C.x_rho(:),method);

                          F.x_rho = RCr(F.lon_rho, F.lat_rho);
    RCr.V = C.y_rho(:);   F.y_rho = RCr(F.lon_rho, F.lat_rho);
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at RHO-points')
  
end

% Grid locations at PSI-points.

if (got.x_psi && got.y_psi),

  RCp = TriScatteredInterp(XpC(:),YpC(:),C.x_psi(:),method);

                        F.x_psi = RCp(XpF, YpF);
  RCp.V = C.y_psi(:);   F.y_psi = RCp(XpF, YpF);

  if (got.lon_psi && got.lat_psi),
    RSp = TriScatteredInterp(C.x_psi(:),C.y_psi(:),C.lon_psi(:),method);

                            F.lon_psi = RSp(F.x_psi, F.y_psi);
    RSp.V = C.lat_psi(:);   F.lat_psi = RSp(F.x_psi, F.y_psi);
  end

elseif (got.lon_psi && got.lat_psi),

  RSp = TriScatteredInterp(XpC(:),YpC(:),C.lon_psi(:),method);

                          F.lon_psi = RSp(XpF, YpF);
  RSp.V = C.lat_psi(:);   F.lat_psi = RSp(XpF, YpF);

  if (got.x_psi && got.x_psi),
    RCp = TriScatteredInterp(C.lon_psi(:),C.lat_psi(:),C.x_psi(:),method);

                          F.x_psi = RCp(F.lon_psi, F.lat_psi);
    RCp.V = C.y_psi(:);   F.y_psi = RCp(F.lon_psi, F.lat_psi);
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at PSI-points')

end

% Grid locations at U-points.

if (got.x_u && got.y_u),

  RCu = TriScatteredInterp(XuC(:),YuC(:),C.x_u(:),method);

                      F.x_u = RCu(XuF, YuF);
  RCu.V = C.y_u(:);   F.y_u = RCu(XuF, YuF);

  if (got.lon_u && got.lat_u),
    RSu = TriScatteredInterp(C.x_u(:),C.y_u(:),C.lon_u(:),method);

                          F.lon_u = RSu(F.x_u, F.y_u);
    RSu.V = C.lat_u(:);   F.lat_u = RSu(F.x_u, F.y_u);
  end
  
elseif (got.lon_u && got.lat_u),

  RSu = TriScatteredInterp(XuC(:),YuC(:),C.lon_u(:),method);

                        F.lon_u = RSu(XuF, YuF);
  RSu.V = C.lat_u(:);   F.lat_u = RSu(XuF, YuF);

  if (got.x_u && got.x_u),
    RCu = TriScatteredInterp(C.lon_u(:),C.lat_u(:),C.x_u(:),method);

                        F.x_u = RCu(F.lon_u, F.lat_u);
    RCu.V = C.y_u(:);   F.y_u = RCu(F.lon_u, F.lat_u);
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at U-points')

end

% Grid locations at V-points.

if (got.x_v && got.y_v),

  RCv = TriScatteredInterp(XvC(:),YvC(:),C.x_v(:),method);

                      F.x_v = RCv(XvF, YvF);
  RCv.V = C.y_v(:);   F.y_v = RCv(XvF, YvF);

  if (got.lon_v && got.lat_v),
    RSv = TriScatteredInterp(C.x_v(:),C.y_v(:),C.lon_v(:),method);

                          F.lon_v = RSv(F.x_v, F.y_v);
    RSv.V = C.lat_v(:);   F.lat_v = RSv(F.x_v, F.y_v);
  end
  
elseif (got.lon_v && got.lat_v),

  RSv = TriScatteredInterp(XvC(:),YvC(:),C.lon_v(:),method);

                        F.lon_v = RSv(XvF, YvF);
  RSv.V = C.lat_v(:);   F.lat_v = RSv(XvF, YvF);

  if (got.x_v && got.x_v),
    RCv = TriScatteredInterp(C.lon_v(:),C.lat_v(:),C.x_v(:),method);

                        F.x_v = RCv(F.lon_v, F.lat_v);
    RCv.V = C.y_v(:);   F.y_v = RCv(F.lon_v, F.lat_v);
  end

else
  
  error(' COARSE2FINE: unable to find coordinates at V-points')
  
end


else %COAWST

% Grid locations at RHO-points.

if (got.x_rho && got.y_rho),
  F.x_rho=interp2(XrC',YrC',C.x_rho',XrF',YrF','spline')';
  F.y_rho=interp2(XrC',YrC',C.y_rho',XrF',YrF','spline')';
  if (got.lon_rho && got.lat_rho),
    F.lon_rho=interp2(XrC',YrC',C.lon_rho',XrF',YrF','spline')';
    F.lat_rho=interp2(XrC',YrC',C.lat_rho',XrF',YrF','spline')';
  end
elseif (got.lon_rho && got.lat_rho),
  F.lon_rho=interp2(XrC',YrC',C.lon_rho',XrF',YrF','spline')';
  F.lat_rho=interp2(XrC',YrC',C.lat_rho',XrF',YrF','spline')';
  if (got.x_rho && got.y_rho),
    F.x_rho=interp2(XrC',YrC',C.x_rho',XrF',YrF','spline')';
    F.y_rho=interp2(XrC',YrC',C.y_rho',XrF',YrF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at RHO-points')
end

% Grid locations at PSI-points.

if (got.x_psi && got.y_psi),
  F.x_psi=interp2(XpC',YpC',C.x_psi',XpF',YpF','spline')';
  F.y_psi=interp2(XpC',YpC',C.y_psi',XpF',YpF','spline')';
  if (got.lon_psi && got.lat_psi),
    F.lon_psi=interp2(XpC',YpC',C.lon_psi',XpF',YpF','spline')';
    F.lat_psi=interp2(XpC',YpC',C.lat_psi',XpF',YpF','spline')';
  end
elseif (got.lon_psi && got.lat_psi),
  F.lon_psi=interp2(XpC',YpC',C.lon_psi',XpF',YpF','spline')';
  F.lat_psi=interp2(XpC',YpC',C.lat_psi',XpF',YpF','spline')';
  if (got.x_psi && got.y_psi),
    F.x_psi=interp2(XpC',YpC',C.x_psi',XpF',YpF','spline')';
    F.y_psi=interp2(XpC',YpC',C.y_psi',XpF',YpF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at PSI-points')
end

% Grid locations at U-points.

if (got.x_u && got.y_u),
  F.x_u=interp2(XuC',YuC',C.x_u',XuF',YuF','spline')';
  F.y_u=interp2(XuC',YuC',C.y_u',XuF',YuF','spline')';
  if (got.lon_u && got.lat_u),
    F.lon_u=interp2(XuC',YuC',C.lon_u',XuF',YuF','spline')';
    F.lat_u=interp2(XuC',YuC',C.lat_u',XuF',YuF','spline')';
  end
elseif (got.lon_u && got.lat_u),
  F.lon_u=interp2(XuC',YuC',C.lon_u',XuF',YuF','spline')';
  F.lat_u=interp2(XuC',YuC',C.lat_u',XuF',YuF','spline')';
  if (got.x_u && got.y_u),
    F.x_u=interp2(XuC',YuC',C.x_u',XuF',YuF','spline')';
    F.y_u=interp2(XuC',YuC',C.y_u',XuF',YuF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at U-points')
end

% Grid locations at V-points.

if (got.x_v && got.y_v),
  F.x_v=interp2(XvC',YvC',C.x_v',XvF',YvF','spline')';
  F.y_v=interp2(XvC',YvC',C.y_v',XvF',YvF','spline')';
  if (got.lon_v && got.lat_v),
    F.lon_v=interp2(XvC',YvC',C.lon_v',XvF',YvF','spline')';
    F.lat_v=interp2(XvC',YvC',C.lat_v',XvF',YvF','spline')';
  end
elseif (got.lon_v && got.lat_v),
  F.lon_v=interp2(XvC',YvC',C.lon_v',XvF',YvF','spline')';
  F.lat_v=interp2(XvC',YvC',C.lat_v',XvF',YvF','spline')';
  if (got.x_v && got.y_v),
    F.x_v=interp2(XvC',YvC',C.x_v',XvF',YvF','spline')';
    F.y_v=interp2(XvC',YvC',C.y_v',XvF',YvF','spline')';
  end
else
  error(' COARSE2FINE: unable to find coordinates at V-points')
end

end


% Get grid lengths.

if (got.x_psi && got.y_psi),
  F.xl = max(F.x_psi(:)) - min(F.x_psi(:));
  F.el = max(F.y_psi(:)) - min(F.y_psi(:));
else
  F.xl = 0;
  F.el = 0;
end

% Other grid variables. The inverse metrics pm and pn cannot be
% interpolated. They need to be recomputed.

Rr  = TriScatteredInterp(C.x_rho(:),C.y_rho(:),C.f(:),method);
F.f = Rr(F.x_rho, F.y_rho);

if (got.angle),
  Rr.V = C.angle(:);   F.angle = Rr(F.x_rho, F.y_rho);
end

% Bathymetry.  Make sure that the interpolation method is linear
% to insure that the global integral of bathymetry is conserved.
% Perhaps, we need to have here a conservation interpolation.
% See Clark and Farley (1984) equations 30-36.
%
% Clark, T.L. and R.D. Farley, 1984:  Severe Downslope Windstorm
%   Calculations in Two and Three Spatial Dimensions Using Anelastic
%   Interative Grid Nesting: A Possible Mechanism for Gustiness,
%   J. Atmos. Sci., 329-350.

Rr.V = C.h(:);    Rr.Method = 'linear';

F.h  = Rr(F.x_rho, F.y_rho);

if (got.hraw),
  C.hraw = nc_read(Ginp,'hraw',1);
  
  Rr.V = C.hraw(:);   F.hraw = Rr(F.x_rho, F.y_rho);
end

% Land/sea masking: use nondimensional fractional coordinates.

if (got.mask_rho || got.mask_psi || got.mask_u || got.mask_v),
  F.mask_rho = interp2(XrC', YrC', C.mask_rho', XrF, YrF, 'nearest');
  
  [F.mask_u, F.mask_v, F.mask_psi]=uvp_masks(F.mask_rho);
end

% Recompute some variables.

deg2rad = pi / 180.0;
omega   = 2.0 * pi * 366.25 / (24.0 * 3600.0 * 365.25);

if (spherical),
  F.f_new = 2.0 * omega * sin(deg2rad * F.lat_rho);
end

% Set grid metrics.

if (C.uniform),

% Coarse has a uniform grid distribution.

  dx = 1/unique(C.pm(:));
  dy = 1/unique(C.pn(:));

  F.pm = ones(size(F.h)) .* (Gfactor/dx);
  F.pn = ones(size(F.h)) .* (Gfactor/dy);
  
  F.dndx = zeros(size(F.h));
  F.dmde = zeros(size(F.h));

else

% Recompute metrics at refined resolution.  We cannot interpolate.
% The computation of metrics from discrete point is subject to
% roundoff.  There is not much that we can do here.  The roundoff
% is small and of the order 1.0E-16 (eps value).

  disp(' ');
  if (spherical),
    GreatCircle = true;
    disp('Computing grid spacing: great circle distances');
  else
    GreatCircle = false;
    disp('Computing grid spacing: Cartesian distances');
  end

  [F.pm, F.pn, F.dndx, F.dmde]=grid_metrics(F, GreatCircle);
end

%--------------------------------------------------------------------------
% Create finer resolution Grid NetCDF file and write out data.
%--------------------------------------------------------------------------

% Set number of grid points.

[LpF,MpF] = size(F.x_rho);              % RHO-points

disp(' ');
disp(['Number of points:',                                              ...
      ' Coarse = ', num2str(Lp) ,' x ', num2str(Mp),                    ...
      ',  fine = ', num2str(LpF),' x ', num2str(MpF)]);

% Check information structure for compliance. If applicable, change
% spherical variable to integer. Modify the Land/Sea masking attributes
% for compliance. If necessary, add the "coordinates" attribute.

I = nc_check(I);

% Create ROMS Grid NetCDF file. Rewrite dimension in input file
% information structure, I.

I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length = LpF;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length = MpF;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_psi' )).Length = LpF-1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_psi')).Length = MpF-1;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_u'   )).Length = LpF-1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_u'  )).Length = MpF;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_v'   )).Length = LpF;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_v'  )).Length = MpF-1;

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

nc_create(Gout, mode, I);

% Set global attributes.

status = nc_attadd(Gout, 'parent_grid', Ginp);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Imin', int32(Imin));
if (status ~= 0), return, end
  
status = nc_attadd(Gout, 'parent_Imax', int32(Imax));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Jmin', int32(Jmin));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Jmax', int32(Jmax));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'refine_factor', int32(Gfactor));
if (status ~= 0), return, end

status = nc_attdel(Gout, 'history');
if (status ~= 0), return, end

history = ['GRID file created using Matlab script: coarse2fine, ',      ...
           date_stamp];
status = nc_attadd(Gout, 'history', history);
if (status ~= 0), return, end

% Write out fine resolution grid variables.

disp(['Writing finer grid variables into: ', Gout]);
disp(' ');

% Write out fine resolution grid variables.

status = nc_write (Gout, 'spherical', F.spherical);
if (status ~= 0), return, end

status = nc_write (Gout, 'xl', F.xl);
if (status ~= 0), return, end

status = nc_write (Gout, 'el', F.el);
if (status ~= 0), return, end

if (got.hraw),
  status = nc_write (Gout, 'hraw', F.hraw, 1);
  if (status ~= 0), return, end
end

for value = field_list,
  field = char(value);
  status = nc_write (Gout, field, F.(field));
  if (status ~= 0), return, end
end

% If appropriate, write coastline data.

if (spherical && got.lon_coast && got.lat_coast),
  add_coastline (Gout, C.lon_coast, C.lat_coast);
  if (status ~= 0), return, end
end

%--------------------------------------------------------------------------
% Get full refinement grid structure.
%--------------------------------------------------------------------------

F = get_roms_grid(Gout);

return


