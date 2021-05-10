function [F]=roms_field(ncname, rec, FieldName)

% ROMS_field:  Computes requested secondary ROMS field.
%
% [F]=roms_field(ncname, rec, FieldName)
%
% This function computes requested secondary 2D or 3D field from ROMS
% state output history/average data.
%
% On Input:
%
%    ncname      ROMS history/average NetCDF filename (string)
%
%    rec         Time record to process from history/average file (scalar)
%
%    FieldName   Secondary field to compute (string):
%
%                  'rvorticity2d'       2D relative vorticity
%                  'rvorticity3d'       3D relative vorticity
% On Output:
%
%    F           Requested field structure (struct):
%
%                  F.ncname            input NetCDF
%                  F.record            time record
%                  F.name              secondary field name
%                  F.units             field units
%                  F.value             field values
%                  F.mask              land/sea mask
%                  F.X                 field X-locations
%                  F.Y                 field Y-locations
%                  F.Z                 field Z-locations
%
% Currently, relative vorticity is supported.  Other fields will be added
% in the future.
%

% svn $Id: roms_field.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Get ROMS grid structure.

G=get_roms_grid(ncname, ncname, rec);

% Compute requested field.

F.ncname = ncname;
F.record = rec;

switch FieldName
  case 'rvorticity2d'
    F.name  = '2D relative vorticity at PSI-points';
    F.units = 's-1';
    ubar = nc_read(ncname, 'ubar', rec);
    vbar = nc_read(ncname, 'vbar', rec);
    F.value = rvorticity(ubar, vbar, G.pm, G.pn, G.mask_psi);
    F.mask = G.mask_psi;
    F.X = G.lon_psi;
    F.Y = G.lat_psi;  
  case 'rvorticity3d'
    F.name  = '3D relative vorticity at PSI-points';
    F.units = 's-1';
    u = nc_read(ncname, 'u', rec);
    v = nc_read(ncname, 'v', rec);
    F.value = rvorticity(u, v, G.pm, G.pn, G.mask_psi);
    F.mask = G.mask_psi;
    F.X = G.lon_psi;
    F.Y = G.lat_psi;
    F.Z = depths(ncname, ncname, 2, 0, rec);
  otherwise
    disp(['Secondary field:',FieldName, ' is not available...'])
end

return
