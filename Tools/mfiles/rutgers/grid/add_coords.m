function add_coords(ncfile, gname)

%
% ADD_COORDS:  Adds horizontal coordinates to a NetCDF file
%
% add_coords(ncfile, gname)
%
% Adds horizontal spherical or cartesian coordinates to specified
% NetCDF file.  The coordinates are read from given grid NetCDF file.
% In some cases, we need CF compliant NetCDF files for third party
% software.
%
% On Input:
%
%    ncfile       NetCDF filename to process (string)
%    gname        Source NetCDF containing grid variables (string)
%

% svn $Id: add_coords.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire NetCDF files.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);              % NetCDF file to process
G = nc_inq(gname);               % Source NetCDF conatining grid variables

if (any(strcmp({G.Variables.Name}, 'spherical')))
  spherical = nc_read(gname, 'spherical');
  if (ischar(spherical))
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
else
  spherical = true;
end

if (spherical)
  Coords = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v'};
  for var = Coords
    field = char(var);
    [got.(field)] = any(strcmp({I.Variables.Name}, field));
  end
else
  Coords = {'x_rho', 'y_rho', 'x_u', 'y_u', 'x_v', 'y_v'};
  for var = Coords
    field = char(var);
    [got.(field)] = any(strcmp({I.Variables.Name}, field));
  end
end

%--------------------------------------------------------------------------
%  If appropriate, define sponge variables.
%--------------------------------------------------------------------------

ic = 0;

append_vars = false;

for var = Coords
  field = char(var);
  if (~got.(field))
    ic = ic + 1;
    S.Dimensions = I.Dimensions;
    S.Variables(ic) = roms_metadata(field, spherical);
    append_vars = true;
  else
    disp(['Variable ''',field,''' already exists. Updating its value']);
  end
end

if (append_vars),
  check_metadata(S);
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write out coordinate variables into NetCDF file.
%--------------------------------------------------------------------------

for var = Coords
  field = char(var);
  f = nc_read(gname, field);
  status = nc_write(ncfile, field, f)
end

return
