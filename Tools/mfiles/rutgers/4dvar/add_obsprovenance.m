function add_obsprovenance (obs_ncfile, Ginp)

%
% ADD_OBSPROVENANCE:  Adds obs_provenance to observation file
%  
% add_obsprovenance (obs_ncfile)
%XJ
% Adds the observation origin (obs_prevenance) to an existing NetCDF file.
%
% On Input:
%
%    obs_ncfile        Observations NetCDF file name (string)
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

got.provenance = any(strcmp({I.Variables.Name}, 'obs_provenance'));

%--------------------------------------------------------------------------
% If appropriate, define (lon,lat) variables.
%--------------------------------------------------------------------------

% The strategy here is to copy the same metadata structure as in the
% other variables in the observation file and edit the appropiate
% values using the intrinsic interface.

S.Dimensions = I.Dimensions;

ic = 0;

append_vars = false;

if (~got.provenance),
  index = strcmp({I.Variables.Name}, 'obs_type');
  ic = ic + 1;
  S.Variables(ic) = I.Variables(index); 
  S.Variables(ic).Name = 'obs_provenance';
  S.Variables(ic).Attributes(1).Value = 'observation origin';
  S.Variables(ic).Attributes(2).Name  = 'flag_values';
  S.Variables(ic).Attributes(2).Value = 0;
  S.Variables(ic).Attributes(2).Name  = 'flag_meaning';
  S.Variables(ic).Attributes(2).Value = 'model generated';
  append_vars = true;
end

if (append_vars),
  nc_append(obs_ncfile, S);
end
  

return
