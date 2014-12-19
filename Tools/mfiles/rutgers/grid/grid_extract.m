function S = grid_extract(Ginp,Gout,Imin,Imax,Jmin,Jmax)

%
% GRID_EXTRACT:  Extracts a subdomain Grid NetCDF file
%
% S = grid_extract(Ginp,Gout,Imin,Imax,Jmin,Jmax)
%
% Given a larger Grid NetCDF file (Ginp), this function extracts and
% creates a subdomain Grid NetCDF file (Gout) for the sampled region
% in terms of larger grid coordinates (Imin,Jmin) and (Imax,Jmax).
% Notice that the (Imin,Jmin) and (Imax,Jmax) indices are located at
% PSI-points. They actually define the physical boundaries of the
% smaller subdomain grid.
%
% On Input:
%
%    Ginp       Input  larger  Grid NetCDF file name (string)
%           or, an existing grid structure (struct array)
%
%    Gout       Output smaller Grid NetCDF file name (string)
%
%    Imin       Larger grid lower-left  I-coordinate (PSI-point)
%
%    Imax       Larger grid upper-right I-coordinate (PSI-point)
%
%    Jmin       Larger grid lower-left  J-coordinate (PSI-point)
%
%    Jmax       Larger grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    S          Smaller Grid structure
%

% svn $Id: grid_extract.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% If applicable, get larger grid structure.

if (~isstruct(Ginp)),
  L = get_roms_grid(Ginp);
else
  L = Ginp;
end

% Check extraction region.

if (Imin >= Imax),
  error([' GRID_EXTRACT: Imin >= Imax,    ',                            ...
         ' Imin = ', num2str(Imin),                                     ...
         ' Imax = ', num2str(Imax)]);
end

if (Jmin >= Jmax),
  error([' GRID_EXTRACT: Jmin >= Jmax,    ',                            ...
         ' Jmin = ', num2str(Jmin),                                     ...
         ' Jmax = ', num2str(Jmax)]);
end

if (Imax > L.Lm+1),
  error([' GRID_EXTRACT: Imax > L,    ',                                ...
         ' Imax = ', num2str(Imax),                                     ...
         ' L = ', num2str(L.Lm+1)]);
end

if (Jmax > L.Mm+1),
  error([' GRID_EXTRACT: Jmax > M,    ',                                ...
         ' Jmax = ', num2str(Jmax),                                     ...
         ' M = ', num2str(L.Mm+1)]);
end

% Set grid variables to process.

grd_vars = {'h', 'f', 'angle', 'pm', 'pn', 'dndx', 'dmde',              ...
            'x_rho', 'y_rho', 'x_psi', 'y_psi',                         ...
            'x_u', 'y_u', 'x_v', 'y_v'};

if (L.spherical),
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u', 'lat_u', 'lon_v', 'lat_v'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

grd_vars = [grd_vars, 'hraw'];               % needs to be last

%--------------------------------------------------------------------------
% Extract subdomain grid.  Only valid for C-grid...
%--------------------------------------------------------------------------

% Set indices to extract.  We need to account for the zero shift of
% indices in Matlab. It is not that intuitive but the logic below is
% correct. It is tricky... For example in Matlab, the Imin corresponds
% to Cartesian coordinates to x_psi(Ip)=x_rho(Ir) since Ip=Ir. However,
% in Fortran Ip=Ir+1.  Therefore, the zero shift is done at the end.

delta = 1;

Ip = (Imin:delta:Imax  );
Ir = (Imin:delta:Imax+1);
Iu = (Imin:delta:Imax  );
Iv = (Imin:delta:Imax+1);

Jp = (Jmin:delta:Jmax  );
Jr = (Jmin:delta:Jmax+1);
Ju = (Jmin:delta:Jmax+1);
Jv = (Jmin:delta:Jmax  );

% Initialize several subsomain structure paameters.

S.grid_name = Gout;
S.Lm = Imax-Imin+1;
S.Mm = Jmax-Jmin+1;
S.spherical = L.spherical;

% Extract grid variables.

for value = grd_vars,
  field = char(value);
  got.(field) = false;
  if (isfield(L,field)),
    switch (field)
      case {'lon_psi', 'lat_psi', 'mask_psi', 'x_psi', 'y_psi'}
        S.(field) = L.(field)(Ip,Jp);
        got.(field) = true;     
      case {'lon_u', 'lat_u', 'mask_u', 'x_u', 'y_u'}
        S.(field) = L.(field)(Iu,Ju);
        got.(field) = true;     
      case {'lon_v', 'lat_v', 'mask_v', 'x_v', 'y_v'}
        S.(field) = L.(field)(Iv,Jv);
          got.(field) = true;     
      case {'hraw'}     
        S.(field) = L.(field)(Ir,Jr,:);
        got.(field) = true;     
      otherwise
        S.(field) = L.(field)(Ir,Jr);
        got.(field) = true;     
    end
  end  
end

% Get grid lengths.

if (got.x_psi && got.y_psi),
  S.xl = max(S.x_psi(:)) - min(S.x_psi(:));
  S.el = max(S.y_psi(:)) - min(S.y_psi(:));
else
  S.xl = 0;
  S.el = 0;
end

%--------------------------------------------------------------------------
% Create subdomain Grid NetCDF file and write out data.
%--------------------------------------------------------------------------

% Set number of grid points.

[Lp,Mp] = size(S.h);                 % RHO-points

disp(' ');
disp(['Number of points:',                                              ...
      ' Large Grid = ', num2str(L.Lm+2) ,' x ', num2str(L.Mm+2),        ...
      ',  Small Grid = ', num2str(S.Lm+2),' x ', num2str(S.Mm+2)]);

% Create ROMS Grid NetCDF file.

NewFile = true;

status = c_grid(Lp, Mp, Gout, NewFile, S.spherical);
if (status ~= 0), return, end

% Set global attributes.

status = nc_attdel(Gout, 'history');
if (status ~= 0), return, end

history = ['GRID file created using Matlab script: grid_extract, ',     ...
           date_stamp];
status = nc_attadd(Gout, 'history', history);
if (status ~= 0), return, end

% Write out fine resolution grid variables.

disp(['Writing subdomain grid variables into: ', Gout]);
disp(' ');

status = nc_write (Gout, 'spherical', S.spherical);
if (status ~= 0), return, end

status = nc_write (Gout, 'xl', S.xl);
if (status ~= 0), return, end

status = nc_write (Gout, 'el', S.el);
if (status ~= 0), return, end

if (got.hraw),
  bath = size(S.hraw,3);
  for rec=1:bath,
    status = nc_write (Gout, 'hraw', S.hraw, rec);
    if (status ~= 0), return, end
  end
else
  status = nc_write (Gout, 'hraw', S.h, 1);
  if (status ~= 0), return, end
end

for value = 1:length(grd_vars)-1,
  field = char(grd_vars(value));
  if (got.(field)),
    status = nc_write (Gout, field, S.(field));
    if (status ~= 0), return, end
  end  
end

%--------------------------------------------------------------------------
% Add coastline data if available.
%--------------------------------------------------------------------------

if (isfield(L, 'lon_coast') && isfield(L, 'lat_coast')),
  add_coastline (Gout, L.lon_coast, L.lat_coast);
end

%--------------------------------------------------------------------------
% Get full extracted grid structure.
%--------------------------------------------------------------------------

S = get_roms_grid(Gout);

return
