function add_obslonlat (obs_ncfile, Ginp)

%
% ADD_OBSLONLAT:  Adds (lon,lat) location of observations
%  
% add_obslonlat (obs_ncfile, Ginp)
%
% Adds the observations locations (lon,lat) to an existing NetCDF file by
% interoplating from the fractional coordinates (obs_Xgrid, obs_Ygrid). All
% the observations are assumed to be at RHO-points in ROMS.
%
% Although the observations (lon,lat) coordinates are not needed in ROMS,
% this information is useful for other post-processing applications.
%
% On Input:
%
%    obs_ncfile        Observations NetCDF file name (string)
%
%    Ginp              Application GRID NetCDF file name (string)
%                  or, an existing grid structure (struct array)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire grid NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(obs_ncfile);

got.lon = any(strcmp({I.Variables.Name}, 'obs_lon'));
got.lat = any(strcmp({I.Variables.Name}, 'obs_lat'));

%--------------------------------------------------------------------------
% If appropriate, define (lon,lat) variables.
%--------------------------------------------------------------------------

% The strategy here is to copy the same metadata structure as in the
% other variables in the observation file and edit the appropiate
% values using the intrinsic interface.

S.Dimensions = I.Dimensions;

ic = 0;

append_vars = false;

if (~got.lon),
  index = strcmp({I.Variables.Name}, 'obs_Xgrid');
  ic = ic + 1;
  S.Variables(ic) = I.Variables(index); 
  S.Variables(ic).Name = 'obs_lon';
  S.Variables(ic).Attributes(1).Value = 'observation longitude';
  S.Variables(ic).Attributes(2).Name  = 'units';
  S.Variables(ic).Attributes(2).Value = 'degrees_east';
  append_vars = true;
end

if (~got.lat),
  index = strcmp({I.Variables.Name}, 'obs_Ygrid');
  ic = ic + 1;
  S.Variables(ic) = I.Variables(index); 
  S.Variables(ic).Name = 'obs_lat';
  S.Variables(ic).Attributes(1).Value = 'observation latitude';
  S.Variables(ic).Attributes(2).Name  = 'units';
  S.Variables(ic).Attributes(2).Value = 'degrees_north';
  append_vars = true;
end

if (append_vars),
  nc_append(obs_ncfile, S);
end
  
%--------------------------------------------------------------------------
% Write out observation (lon,lat).
%--------------------------------------------------------------------------

% Get observations fractional coordinates.

Xobs = nc_read(obs_ncfile, 'obs_Xgrid');      % fractional I-grid
Yobs = nc_read(obs_ncfile, 'obs_Ygrid');      % fractional J-grid

% Get application spherical coordinates at RHO-points

if (~isstruct(Ginp))
  rlon = nc_read(Ginp, 'lon_rho');
  rlat = nc_read(Ginp, 'lat_rho');
else
  rlon = Ginp.lon_rho;
  rlat = Ginp.lat_rho;
end
[Lp,Mp] = size(rlon);

% Interpolate (lon,lat) from fractional grid coordinates.  Notice
% that in 4D-Var, the arbitrary origin of the fractional coordinates
% is (0.0, 0.0).

[Jgrid, Igrid]=meshgrid(0:Mp-1, 0:Lp-1);

method = 'linear';

F = griddedInterpolant(Igrid, Jgrid, rlon, method);

                    obs_lon = F(Xobs, Yobs);
F.Values = rlat;    obs_lat = F(Xobs, Yobs);

% Write out observations (lon,lat) locations.

status = nc_write(obs_ncfile, 'obs_lon', obs_lon);
status = nc_write(obs_ncfile, 'obs_lat', obs_lat);

return
