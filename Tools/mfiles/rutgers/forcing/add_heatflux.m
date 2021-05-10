function add_heatflux (ncfile)

%
% ADD_HEATFLUX:  Adds heat flux componets to existing NetCDF file
%  
% add_heatflux (ncfile)
%
% Adds the heat flux components (shortwave, longwave, latent, and sensible)
% to an existent NetCdF file if not present.  The file must contains the
% net surface heat flux variable so it can be used as metadata for the
% heat flux components.
%
% On Input:
%
%    ncfile            Existing ROMS NetCDF filename (string)
%

% svn $Id: add_heatflux.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire input NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);

got.shflux   = any(strcmp({I.Variables.Name}, 'shflux'));
got.swrad    = any(strcmp({I.Variables.Name}, 'swrad'));
got.lwrad    = any(strcmp({I.Variables.Name}, 'lwrad'));
got.latent   = any(strcmp({I.Variables.Name}, 'latent'));
got.sensible = any(strcmp({I.Variables.Name}, 'sensible'));

%--------------------------------------------------------------------------
% If appropriate, define heat flux component variables.
%--------------------------------------------------------------------------

% The strategy here is to copy the same metadata structure as in the
% surface net heat flux 'shflux' and edit the appropiate values using
% the intrinsic interface.

S.Dimensions = I.Dimensions;

ic = 0;

append_vars = false;

if (got.shflux)

  if (~got.swrad),
    index = strcmp({I.Variables.Name}, 'shflux');
    ic = ic + 1;
    S.Variables(ic) = I.Variables(index); 
    S.Variables(ic).Name = 'swrad';
    S.Variables(ic).Attributes(1).Value = 'solar shortwave radiation flux';
    disp(['Adding ', S.Variables(ic).Name, ': ',                        ...
          S.Variables(ic).Attributes(1).Value]);
    append_vars = true;
  end

  if (~got.lwrad),
    index = strcmp({I.Variables.Name}, 'shflux');
    ic = ic + 1;
    S.Variables(ic) = I.Variables(index); 
    S.Variables(ic).Name = 'lwrad';
    S.Variables(ic).Attributes(1).Value = 'longwave radiation flux';
    disp(['Adding ', S.Variables(ic).Name, ': ',                        ...
          S.Variables(ic).Attributes(1).Value]);
    append_vars = true;
  end

  if (~got.latent),
    index = strcmp({I.Variables.Name}, 'shflux');
    ic = ic + 1;
    S.Variables(ic) = I.Variables(index); 
    S.Variables(ic).Name = 'latent';
    S.Variables(ic).Attributes(1).Value = 'surface latent heat flux';
    disp(['Adding ', S.Variables(ic).Name, ': ',                        ...
          S.Variables(ic).Attributes(1).Value]);
    append_vars = true;
  end

  if (~got.sensible),
    index = strcmp({I.Variables.Name}, 'shflux');
    ic = ic + 1;
    S.Variables(ic) = I.Variables(index); 
    S.Variables(ic).Name = 'sensible';
    S.Variables(ic).Attributes(1).Value = 'surface sensible heat flux';
    disp(['Adding ', S.Variables(ic).Name, ': ',                        ...
          S.Variables(ic).Attributes(1).Value]);
    append_vars = true;
  end
  
else
  
  error(['Cannot find variable: shflux']);

end

if (append_vars),
  nc_append(ncfile, S);
end

return
