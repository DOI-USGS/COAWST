function R = grid_rotated(Ginp,Gout,Xcenter,Ycenter,theta,dx,dy,Im,Jm,Lplt)

%
% GRID_ROTATED:  Extracts a rotated subdomain Grid NetCDF file
%
% R = grid_rotated(Ginp,Gout,Xcenter,Ycenter,theta,dx,dy,Im,Jm,Lplt)
%
% Given a larger Grid NetCDF file (Ginp), this function extracts and
% creates a rotated subdomain Grid NetCDF file (Gout) centered at
% (Xmid,Ymid) and rotated by angle theta. 
%
% In Cartesian coordinates, the rotation is defined as:
%
%    Xrot = X * cosd(theta) + Y * sind(theta);
%    Yrot = Y * cosd(theta) - X * sind(theta);
%   
% On Input:
%
%    Ginp       Input  larger  Grid NetCDF file name (string)
%
%    Gout       Output rotated Grid NetCDF file name (string)
%
%    Xcenter    Center of Cartesian (meter) or longitude (degrees_east)
%                 of the rotated grid at PSI-points
%
%    Ycenter    Center of Cartesian (meter) or latitude (degrees_north)
%                 of the rotated grid at PSI-points
%
%    theta      Rotation angle between positive XI-axis relative to EAST
%                (degrees)
%
%    dx         Rotated grid spacing in the XI-direction (meter)
%
%    dy         Rotated grid spacing in the ETA-direction (meter)
%
%    Im         Number of points in the XI-direction
%
%    Jm         Number of points in the ETA-direction
%
%    Lplt       Switch to plot rotated grid (OPTIONAL, default=false)
%
% On Output:
%
%    R          Rotated Grid structure
%
%
% For Example, the rotated grid in the LAKE_JERSEY application is
% extracted using the following command:
%
%    R = grid_rotated('lake_jersey_grd_a.nc','lake_jersey_grd_f.nc', ...
%                     25000, 22000, -24.5, 100, 100, 164, 86);
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2019 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set switch to plot rotated grid.
  
Lplot = false;
if (nargin == 10)
  Lplot = Lplt;
end

% Get larger donor grid structure and inquire about its variables.

D = get_roms_grid(Ginp);
I = nc_inq(Ginp);
vnames = {I.Variables.Name};

% Initialize parameters: We need to use scatteredInterpolant for spherical
% because the grid are not extrictly monotonic and plaid.

if D.spherical
  Smethod = 'linear';         % scatteredInterpolant: C0 continuity
  Gmethod = 'spline';         % griddedInterpolant:   C2 continuity
else
  Gmethod = 'linear';         % griddedInterpolant:   C0 continuity
end

% Set grid variables to process.

grd_vars = {'h', 'f', 'angle', 'pm', 'pn', 'dndx', 'dmde',              ...
            'x_rho', 'y_rho', 'x_psi', 'y_psi',                         ...
            'x_u', 'y_u', 'x_v', 'y_v'};

new_vars = {'xi_rho', 'eta_rho', 'xi_psi', 'eta_psi',                   ...
            'xi_u',   'eta_u',   'xi_v',   'eta_v'};

if (D.spherical)
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u', 'lat_u', 'lon_v', 'lat_v',             ...
                        'lon_coast', 'lat_coast'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

grd_vars = [grd_vars, new_vars];

for value = grd_vars
  field = char(value);
  got.(field) =  any(strcmp(vnames, field));
end
got.hraw = any(strcmp(vnames, 'hraw'));        % special record variable

%--------------------------------------------------------------------------
% Extract rotated subdomain grid.
%--------------------------------------------------------------------------

% Set larger donor grid fractional (XI,ETA) coordinates. The fractional
% coordinated will be added to the rotated grid NetCDF file to facilitate
% computing nesting connectivity elsewhere.

[Lp,Mp] = size(D.h);               % RHO-points

L = Lp-1;
M = Mp-1;

[D.eta_rho, D.xi_rho] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
[D.eta_psi, D.xi_psi] = meshgrid(1.0:1:M     , 1.0:1:L     );
[D.eta_u  , D.xi_u  ] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );
[D.eta_v  , D.xi_v  ] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);

% Compute rotated coordinates.

if (D.spherical)
  nx = Im - 1;                     % XI-number  of interior RHO-points (Lm)
  ny = Jm - 1;                     % ETA-number of interior RHO-points (Mm)
  Clat = Ycenter;   % ?? jcw

  R = spherical_grid (D, nx, ny, dx, dy, Clat, Xcenter, Ycenter, theta);

else

  R = cartesian_grid (D, Im, Jm, dx, dy, Xcenter, Ycenter, theta);

end

% Set rotated grid fractional (XI,ETA) coordinates.

if (D.spherical)

  Drho = scatteredInterpolant(D.lon_rho(:), D.lat_rho(:), D.xi_rho(:),  ...
                              Smethod);
                              R.xi_rho  = Drho(R.lon_rho, R.lat_rho);
  Drho.Values = D.eta_rho(:); R.eta_rho = Drho(R.lon_rho, R.lat_rho);

  Dpsi = scatteredInterpolant(D.lon_psi(:), D.lat_psi(:), D.xi_psi(:),  ...
                              Smethod);
                              R.xi_psi  = Dpsi(R.lon_psi, R.lat_psi);
  Dpsi.Values = D.eta_psi(:); R.eta_psi = Dpsi(R.lon_psi, R.lat_psi);

  Du   = scatteredInterpolant(D.lon_u(:), D.lat_u(:), D.xi_u(:),        ...
                              Smethod);
                              R.xi_u  = Du(R.lon_u, R.lat_u);
  Du.Values = D.eta_u(:);     R.eta_u = Du(R.lon_u, R.lat_u);

  Dv   = scatteredInterpolant(D.lon_v(:), D.lat_v(:), D.xi_v(:),        ...
                              Smethod);
                              R.xi_v  = Dv(R.lon_v, R.lat_v);
  Dv.Values = D.eta_v(:);     R.eta_v = Dv(R.lon_v, R.lat_v);
 
else
 
  Drho = griddedInterpolant(D.x_rho, D.y_rho, D.xi_rho, Gmethod);
                            R.xi_rho  = Drho(R.x_rho, R.y_rho);
  Drho.Values = D.eta_rho;  R.eta_rho = Drho(R.x_rho, R.y_rho);

  Dpsi = griddedInterpolant(D.x_psi, D.y_psi, D.xi_psi, Gmethod);
                            R.xi_psi  = Dpsi(R.x_psi, R.y_psi);
  Dpsi.Values = D.eta_psi;  R.eta_psi = Dpsi(R.x_psi, R.y_psi);
  
  Du   = griddedInterpolant(D.x_u, D.y_u, D.xi_u, Gmethod);
                            R.xi_u  = Du(R.x_u, R.y_u);
  Du.Values = D.eta_u;      R.eta_u = Du(R.x_u, R.y_u);

  Dv   = griddedInterpolant(D.x_v, D.y_v, D.xi_v, Gmethod);
                            R.xi_v  = Dv(R.x_v, R.y_v);
  Dv.Values = D.eta_v;      R.eta_v = Dv(R.x_v, R.y_v);

end

% If spherical grid, interpolate Cartesian coordinates it available
% in donor grid.

if (D.spherical)

  if (got.x_rho && got.y_rho)
    Drho = scatteredInterpolant(D.xi_rho(:), D.eta_rho(:), D.x_rho(:),  ...
                                Smethod);
                                R.x_rho = Drho(R.xi_rho, R.eta_rho);
    Drho.Values = D.y_rho(:);   R.y_rho = Drho(R.xi_rho, R.eta_rho);
  end

  if (got.x_psi && got.y_psi)
    Dpsi = scatteredInterpolant(D.xi_psi(:), D.eta_psi(:), D.x_psi(:),  ...
                                Smethod);
                                R.x_psi = Dpsi(R.xi_psi, R.eta_psi);
    Dpsi.Values = D.y_psi(:);   R.y_psi = Dpsi(R.xi_psi, R.eta_psi);
  end

  if (got.x_u && got.y_u)
    Du   = scatteredInterpolant(D.xi_u(:), D.eta_u(:), D.x_u(:),        ...
                                Smethod);
                                R.x_u = Du(R.xi_u, R.eta_u);
    Du.Values = D.y_u(:);       R.y_u = Du(R.xi_u, R.eta_u);
  end

  if (got.x_v && got.y_v)
    Dv   = scatteredInterpolant(D.xi_v(:), D.eta_v(:), D.x_v(:),        ...
                                Smethod);
                                R.x_v = Dv(R.xi_v, R.eta_v);
    Dv.Values = D.y_v(:);       R.y_v = Dv(R.xi_v, R.eta_v);
  end
  
end

%--------------------------------------------------------------------------
% Set other grid variables. 
%--------------------------------------------------------------------------

% Set application bathymetry (XI-ETA coordinates interpolation).

if (D.spherical)
  Rrho = scatteredInterpolant(D.xi_rho(:), D.eta_rho(:), D.h(:), Smethod);
else
  Rrho = griddedInterpolant(D.xi_rho, D.eta_rho, D.h, Gmethod);
end
R.h  = Rrho(R.xi_rho, R.eta_rho);

% Raw bathymetry (XI-ETA coordinates interpolation).

if (got.hraw)
  try
    D.hraw = nc_read(Ginp, 'hraw', 1);
    Rrho.Values = D.hraw;
    R.hraw  = Rrho(R.xi_rho, R.eta_rho);
  catch
    got.hraw = false;
  end
end
 
% Coriolis parameter (XI-ETA coordinates interpolation).

Rrho.Values = D.f(:);
R.f = Rrho(R.xi_rho, R.eta_rho);

% Land/sea masking: use nondimensional fractional coordinates. Use
% nearest point (XI-ETA) interpolation.

if (got.mask_rho || got.mask_psi || got.mask_u || got.mask_v)
  Rrho.Values = D.mask_rho(:);
  Rrho.Method = 'nearest';
  R.mask_rho  = Rrho(R.xi_rho, R.eta_rho);

  [R.mask_u, R.mask_v, R.mask_psi]=uvp_masks(R.mask_rho);
end

%--------------------------------------------------------------------------
% If requested, plot extracted rotated grid.
%--------------------------------------------------------------------------

if (Lplot)
  figure
  if (D.spherical)
    plot(D.lon_psi(:,1:end),  D.lat_psi(:,1:end),  'k-',                ...
         D.lon_psi(1:end,:)', D.lat_psi(1:end,:)', 'k-',                ...
         R.lon_psi(:), R.lat_psi(:), 'r+');
    axis tight;
    xlabel ('\bf\fontsize{18}Lon_{\psi}')
    ylabel ('\bf\fontsize{18}Lat_{\psi}')
    title (['Xcenter = ', num2str(Xcenter), blanks(4),                  ...
            'Ycenter = ', num2str(Ycenter), blanks(4),                  ...
            'dx = ', num2str(dx/1000), ' km', blanks(4),                ...
            'dy = ', num2str(dy/1000), ' km', blanks(4),                ...
            'theta = ', num2str(theta), blanks(4),                      ...
            'Im = ', num2str(Im), blanks(4),                            ...
            'Jm = ', num2str(Jm)]);
  else
    plot(D.x_psi(:,1:end)/1000,  D.y_psi(:,1:end)/1000,  'k-',          ...
         D.x_psi(1:end,:)'/1000, D.y_psi(1:end,:)'/1000, 'k-',          ...
         R.x_psi(:)/1000, R.y_psi(:)/1000, 'r+');
    xlabel ('\bf\fontsize{18}X_{\psi} (km)')
    ylabel ('\bf\fontsize{18}Y_{\psi} (km)')
    title (['Xcenter = ', num2str(Xcenter/1000), ' km', blanks(4),      ...
            'Ycenter = ', num2str(Ycenter/1000), ' km', blanks(4),      ...
            'dx = ', num2str(dx/1000), ' km', blanks(4),                ...
            'dy = ', num2str(dy/1000), ' km', blanks(4),                ...
            'theta = ', num2str(theta), blanks(4),                      ...
            'Im = ', num2str(Im), blanks(4),                            ...
            'Jm = ', num2str(Jm)]);
  end

  figure
  plot(D.xi_psi(:,1:end),  D.eta_psi(:,1:end),  'k-',                   ...
       D.xi_psi(1:end,:)', D.eta_psi(1:end,:)', 'k-',                   ...
       R.xi_psi(:), R.eta_psi(:), 'r+');
  axis tight;
  xlabel ('\bf\fontsize{18}\xi_{i,j}')
  ylabel ('\bf\fontsize{18}\eta_{i,j}')
  if (D.spherical)
    title (['Xcenter = ', num2str(Xcenter), blanks(4),                    ...
            'Ycenter = ', num2str(Ycenter), blanks(4),                    ...
            'dx = ', num2str(dx/1000),' km', blanks(4),                   ...
            'dy = ', num2str(dy/1000),' km', blanks(4),                   ...
            'theta = ', num2str(theta), blanks(4),                        ...
            'Im = ', num2str(Im), blanks(4),                              ...
            'Jm = ', num2str(Jm)]);
  else
    title (['Xcenter = ', num2str(Xcenter/1000), ' km', blanks(4),        ...
            'Ycenter = ', num2str(Ycenter/1000), ' km', blanks(4),        ...
            'dx = ', num2str(dx/1000),' km', blanks(4),                   ...
            'dy = ', num2str(dy/1000),' km', blanks(4),                   ...
            'theta = ', num2str(theta), blanks(4),                        ...
            'Im = ', num2str(Im), blanks(4),                              ...
            'Jm = ', num2str(Jm)]);
  end    
    
  figure;
  if (D.spherical)
    plot(R.lon_psi(:), R.lat_psi(:), 'r+');

    hold on;
    pcolorjw(D.lon_psi, D.lat_psi, D.mask_psi);
    shading flat; colorbar;
    alpha (0.5)   
    xlabel ('\bf\fontsize{18}Lon_{\psi}')
    ylabel ('\bf\fontsize{18}Lat_{\psi}')
    title (['Xcenter = ', num2str(Xcenter), blanks(4),                  ...
            'Ycenter = ', num2str(Ycenter), blanks(4),                  ...
            'dx = ', num2str(dx/1000), ' km', blanks(4),                ...
            'dy = ', num2str(dy/1000), ' km', blanks(4),                ...
            'theta = ', num2str(theta), blanks(4),                      ...
            'Im = ', num2str(Im), blanks(4),                            ...
            'Jm = ', num2str(Jm)]);
    hold off;
  else
    plot(R.x_psi(:)/1000, R.y_psi(:)/1000, 'r+');

    hold on;
    pcolorjw(D.x_psi/1000, D.y_psi/1000, D.mask_psi);
    shading flat; colorbar;
    alpha (0.5)   
    xlabel ('\bf\fontsize{18}X_{\psi} (km)')
    ylabel ('\bf\fontsize{18}Y_{\psi} (km)')
    title (['Xcenter = ', num2str(Xcenter/1000), ' km', blanks(4),      ...
            'Ycenter = ', num2str(Ycenter/1000), ' km', blanks(4),      ...
            'dx = ', num2str(dx/1000), ' km', blanks(4),                ...
            'dy = ', num2str(dy/1000), ' km', blanks(4),                ...
            'theta = ', num2str(theta), blanks(4),                      ...
            'Im = ', num2str(Im), blanks(4),                            ...
            'Jm = ', num2str(Jm)]);
    hold off;
  end
end

%--------------------------------------------------------------------------
% Create rotated grid NetCDF file.
%--------------------------------------------------------------------------

% Set number of grid points.

[LpR,MpR] = size(R.h);                  % RHO-points

disp(' ');
disp(['Number of points:',                                              ...
      ' Donor = ',     num2str(Lp) ,' x ', num2str(Mp),                 ...
      ',  Rotated = ', num2str(LpR),' x ', num2str(MpR)]);

% Create ROMS Grid NetCDF file. Rewrite dimension in input file
% information structure, I.

I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length = LpR;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length = MpR;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_psi' )).Length = LpR-1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_psi')).Length = MpR-1;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_u'   )).Length = LpR-1;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_u'  )).Length = MpR;

I.Dimensions(strcmp({I.Dimensions.Name},'xi_v'   )).Length = LpR;
I.Dimensions(strcmp({I.Dimensions.Name},'eta_v'  )).Length = MpR-1;

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% Add new grid fractional (XI,ETA) coordinates variables.  In rotated
% grids, the fractionam (XI,ETA) coordinates are not trivial.  They are
% needed whe processing rotated nested contact points.

if (~got.angle)
  got.angle = true;
  new_vars = [new_vars, 'angle'];        a  % add angle if not in donor grid
end

ic = length(I.Variables);

for value = new_vars
  field = char(value);
  got.(field) =  true;
  ic = ic + 1;
  I.Variables(ic) = roms_metadata(field,  D.spherical);
end

% Check metadata structure for compliance. If applicable, change
% spherical variable to integer. Modify the Land/Sea masking attributes
% for compliance. If necessary, add the "coordinates" attribute.

I = nc_check(I);
I = check_metadata(I);

% Create rotated grid NetCDF file.

nc_create(Gout, mode, I);

% Set global attributes.  The extracting factors global attributes names
% depend on "Gfactor" since refinement applications (Gfactor > 1) they
% are used to process the contact points.  The Gfactor=1 extracted grid
% may be used as the coaser grid in a nested application or in a non-
% nested application. It is always good idea to keep this information
% in the extrated file global attributes.

status = nc_attadd(Gout, 'donor_grid', Ginp);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'donor_Xcenter', Xcenter);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'donor_Ycenter', Ycenter);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'rotated_angle', theta);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'rotated_dx', dx);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'rotated_dy', dy);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'rotated_Im', int32(Im));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'rotated_Jm', int32(Jm));
if (status ~= 0), return, end

status = nc_attdel(Gout, 'history');
if (status ~= 0), return, end

history = ['GRID file created using Matlab script ',                    ...
           which(mfilename), blanks(1), date_stamp];
status = nc_attadd(Gout, 'history', history);
if (status ~= 0), return, end

%--------------------------------------------------------------------------
% Write out rotated grid variables.
%--------------------------------------------------------------------------

disp(['Writing finer grid variables into: ', Gout]);
disp(' ');

status = nc_write (Gout, 'spherical', double(R.spherical));
if (status ~= 0), return, end

status = nc_write (Gout, 'xl', R.xl);
if (status ~= 0), return, end

status = nc_write (Gout, 'el', R.el);
if (status ~= 0), return, end

if got.hraw
  status = nc_write (Gout, 'hraw', R.hraw, 1);
  if (status ~= 0), return, end
end

for value = grd_vars
  field = char(value);
  if (isfield(R, field))
    status = nc_write (Gout, field, R.(field));
    if (status ~= 0), return, end
  end
end

% If appropriate, write coastline data.

if D.spherical && got.lon_coast && got.lat_coast
  add_coastline (Gout, D.lon_coast, D.lat_coast);
  if (status ~= 0), return, end
end

%--------------------------------------------------------------------------
% Get full rotated grid structure.
%--------------------------------------------------------------------------

R = get_roms_grid(Gout);

return
