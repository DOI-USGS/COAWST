function add_mask(ncfile, spherical, mask)

%
% ADD_MASK:  Adds land/sea mask arrays to a NetCDF file
%
% add_mask(ncfile, spherical, mask)
%
% Adds land/sea masking array to specified NetCDF file that it is not a
% ROMS grid NetCDF file.  It primarily used to add mask field to datasets
% to aid on spatial interpolations.
%
% On Input:
%
%    ncfile       NetCDF filename to process (string)
%    spherical    Spherical switch: true or false (logical)
%    mask         Land/Sea mask (array)
%

% svn $Id: add_mask.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire NetCDF files.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);              % NetCDF file to process

Masks = {'mask'};
for var = Masks
  field = char(var);
  [got.(field)] = any(strcmp({I.Variables.Name}, field));
end

%--------------------------------------------------------------------------
%  If appropriate, define mask variable.
%--------------------------------------------------------------------------

ic = 0;

append_vars = false;

for var = Masks
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

if (append_vars)
  check_metadata(S);
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write out coordinate variables into NetCDF file.
%--------------------------------------------------------------------------

for var = Masks
  field = char(var);
  nc_write(ncfile, field, mask)
end

return
