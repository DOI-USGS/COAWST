function add_tide_date (ncfile, tide_ref)

%
% ADD_TIDE_DATE:  Adds zero phase date to existing tidal NetCDF file
%  
% add_tide_date (ncfile, tide_ref)
%
% Adds the "zero_start_date" variable and the "base_date" global attribute
% to and existing tidal forcing NetCDF. The "zero_phase_date" reference is
% specified when processing the tidal data and corresponds to the date
% of zero tidal phase. Usually, it is the same as ROMS reference time but
% it can be different.
%
% The NetCDF metadata for "zero_phase_date" variable is:
%
%   double zero_phase_date ;
%     zero_phase_date:long_name = "tide reference date for zero phase" ;
%     zero_phase_date:units = "days as %Y%m%d.%f" ;
%     zero_phase_date:C_format = "%13.4f" ;
%     zero_phase_date:FORTRAN_format = "(f13.4)" ;  
%
% On Input:
%
%    ncfile            Existing ROMS tidal forcing NetCDF filename (string)
%    tide_ref          Tide reference date number value
%
% For example:
%
%    add_tide_date('frc_tides.nc', datenum(2006,1,1))
%  

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire input NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);

got.zero_phase_date = any(strcmp({I.Variables.Name}, 'zero_phase_date'));

%--------------------------------------------------------------------------
% If appropriate, define heat flux component variables.
%--------------------------------------------------------------------------

% The strategy here is to copy the same metadata structure as in the
% surface net heat flux 'shflux' and edit the appropiate values using
% the intrinsic interface.

S.Dimensions = I.Dimensions;

append_vars = false;

if (~got.zero_phase_date)
  ic = 1;
  long_name = 'tidal reference date for zero phase';
  S.Variables(ic).Name = 'zero_phase_date';
  S.Variables(ic).Dimensions = [];
  S.Variables(ic).Size = 1;
  S.Variables(ic).Datatype = 'double';
  S.Variables(ic).Attributes(1).Name  = 'long_name';
  S.Variables(ic).Attributes(1).Value = long_name;
  S.Variables(ic).Attributes(2).Name  = 'units';
  S.Variables(ic).Attributes(2).Value = 'days as %Y%m%d.%f';
  S.Variables(ic).Attributes(3).Name  = 'C_format';
  S.Variables(ic).Attributes(3).Value = '%13.4f';  
  S.Variables(ic).Attributes(4).Name  = 'FORTRAN_format';
  S.Variables(ic).Attributes(4).Value = '(f13.4)';
  
  disp(['Adding variable: ', S.Variables(ic).Name, ' = ',               ...
          S.Variables(ic).Attributes(1).Value]);
  append_vars = true;
end

if (append_vars),
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write 
%--------------------------------------------------------------------------

if (~got.zero_phase_date)
  dayfrac=tide_ref-fix(tide_ref);
  tidedate=round(str2double(datestr(fix(tide_ref),'yyyymmdd'))+dayfrac,4);
  ncwrite(ncfile, 'zero_phase_date', tidedate);
end

return
