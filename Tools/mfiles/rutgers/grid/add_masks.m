function add_masks(ncfile, gname)

%
% ADD_MASKS:  Adds land/sea mask arrays to a NetCDF file
%
% add_masks(ncfile, gname)
%
% Adds land/sea masking array to specified NetCDF file.  The mask are
% read from given grid NetCDF file. In some cases, we need CF compliant
% NetCDF files for third party software.
%
% On Input:
%
%    ncfile       NetCDF filename to process (string)
%    gname        Source NetCDF containing land/sea mask variables (string)
%

% svn $Id: add_masks.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire NetCDF files.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);              % NetCDF file to process
G = nc_inq(gname);               % Source NetCDF conatining masks variables

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

Masks = {'mask_rho', 'mask_u', 'mask_v'};
for var = Masks
  field = char(var);
  [got.(field)] = any(strcmp({I.Variables.Name}, field));
end

%--------------------------------------------------------------------------
%  If appropriate, define sponge variables.
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

if (append_vars),
  check_metadata(S);
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write out coordinate variables into NetCDF file.
%--------------------------------------------------------------------------

for var = Masks
  field = char(var);
  f = nc_read(gname, field);
  status = nc_write(ncfile, field, f)
end

return
