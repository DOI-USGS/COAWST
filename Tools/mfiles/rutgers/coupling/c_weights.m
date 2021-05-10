function c_weights(Im, Jm, ncname)

%
% C_WEIGHTS:  Creates weights factors NetCDF for coupling
%
% c_weights(Im, Jm, ncname)
%
% This function creates a new weight factor NetCDF file for melding
% DATA and OCEAN components during ESMF coupling. They are usually
% by the atmosphere component when its grid is incongruent with that
% of the OCEAN. The DATA model is needed to provide values in grid
% cell not covered by the OCEAN model.
%
% On Input:
%
%    Im          Atmosphere model number of points in the west-east
%                  direction (scalar)
%
%    Jm          Atmosphere model number of points in the south-north
%                  direction (scalar)
%  
%    ncname      Ouptut weights factors NetCDF filename (string)
%
%

% svn $Id: c_weights.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Create weight factors NetCDF.
%--------------------------------------------------------------------------

% Create file.

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

ncid = netcdf.create(ncname, mode);

% Define dimensions.

Did.lon = netcdf.defDim(ncid, 'lon', Im);
Did.lat = netcdf.defDim(ncid, 'lat', Jm);

% Define global attributes.

varid = netcdf.getConstant('nc_global');

Text = 'CF-1.0';
netcdf.putAtt(ncid, varid, 'Conventions', Text);

Text = 'Coupling field melding weights between data and ocean components';
netcdf.putAtt(ncid, varid, 'title', Text);

%--------------------------------------------------------------------------
% Define variables.
%--------------------------------------------------------------------------

Dims = [Did.lon Did.lat];

% Longitude.

varid = netcdf.defVar(ncid, 'lon', netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name', 'longitude');
netcdf.putAtt(ncid, varid, 'standard_name', 'longitude');
netcdf.putAtt(ncid, varid, 'units', 'degrees_east');

% Latitude.

varid = netcdf.defVar(ncid, 'lat', netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name', 'latitude');
netcdf.putAtt(ncid, varid, 'standard_name', 'latitude');
netcdf.putAtt(ncid, varid, 'units', 'degrees_north');

% Land/Sea mask.

varid = netcdf.defVar(ncid, 'mask', netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name', 'land-sea mask');
netcdf.putAtt(ncid, varid, 'flag_values', [0, 1]);
netcdf.putAtt(ncid, varid, 'flag_meanings', 'land sea');
netcdf.putAtt(ncid, varid, 'coordinates', 'lon lat');

% DATA component rigid (uniform) weights.

varid = netcdf.defVar(ncid, 'data_weight_rigid',                        ...
                      netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'DATA component rigid melding weights');
netcdf.putAtt(ncid, varid, 'valid_min', 0.0d0);
netcdf.putAtt(ncid, varid, 'valid_max', 1.0d0);
netcdf.putAtt(ncid, varid, 'coordinates', 'lon lat');

% OCEAN component rigid (uniform) weights.

varid = netcdf.defVar(ncid, 'ocean_weight_rigid',                       ...
                      netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'OCEAN component rigid melding weights');
netcdf.putAtt(ncid, varid, 'valid_min', 0.0d0);
netcdf.putAtt(ncid, varid, 'valid_max', 1.0d0);
netcdf.putAtt(ncid, varid, 'coordinates', 'lon lat');

% DATA component smooth (gradually melding) weights.

varid = netcdf.defVar(ncid, 'data_weight_smooth',                       ...
                      netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'DATA component smooth melding weights');
netcdf.putAtt(ncid, varid, 'valid_min', 0.0d0);
netcdf.putAtt(ncid, varid, 'valid_max', 1.0d0);
netcdf.putAtt(ncid, varid, 'coordinates', 'lon lat');

% OCEAN component smooth (gradually melding) weights.

varid = netcdf.defVar(ncid, 'ocean_weight_smooth',                      ...
                      netcdf.getConstant('nc_double'), Dims);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'OCEAN component smooth melding weights');
netcdf.putAtt(ncid, varid, 'valid_min', 0.0d0);
netcdf.putAtt(ncid, varid, 'valid_max', 1.0d0);
netcdf.putAtt(ncid, varid, 'coordinates', 'lon lat');

%--------------------------------------------------------------------------
% Close file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);

netcdf.close(ncid);

return
