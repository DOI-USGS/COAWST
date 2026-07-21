function c_observations(S, varargin)

%
% C_OBSERVATIONS:  Creates 4D-Var observations NetCDF file
%
% [status]=c_observations(S, file)
%
% This function creates ROMS 4D-Var observation NetCDF file using specified
% in structure array, S.
%
% On Input:
%
%    S   Observations file creation parameters (structure array):
%
%          S.ncfile                NetCDF file name (string)
%          S.spherical             spherical grid switch
%          S.Nstate                number of state variables
%          S.Nsurvey               number of data surveys 
%          S.Ndatum                number of observations
%          S.state_flag_values     obs_type 'flag_values' attribute
%          S.state_flag_meanings   obs_type 'flag_meanings attribute
%          S.origin_flag_values    obs_provenance 'flag_values' attribute
%          S.origin_flag_meanings  obs_provenance 'flag_meaning' attribute
%
%        Optional fields:
%
%          S.title                 application title
%          S.grd_file              'grd_file' global attribute
%          S.global_variables      'state_variables' global attribute
%          S.global_provenance     'obs_provenance' global attribute
%          S.global_sources        'obs_sources' global attribute
%
%    file  Output observation file name (string, OPTIONAL). If specified,
%            it creates this file instead the one in S.ncfile.
%
%
% Examples showing how to assign some of the structure fields:
%
%    newline=sprintf('\n');
%
%    S.state_flag_values=[1:1:7];
%
%    S.state_flag_meanings={'zeta', ...
%                           'ubar', ...
%                           'vbar', ...
%                           'u', ...
%                           'v', ...
%                           'temperature', ...
%                           'salinity'};
%
%    S.origin_flag_values=[1:1:11];
%
%    S.origin_flag_meanings={'gridded_AVISO_SLA', ...
%                            'blended_SST', ...
%                            'XBT_Met_Office', ...
%                            'CTD_temperature_Met_Office', ...
%                            'CTD_salinity_Met_Office', ...
%                            'ARGO_temperature_Met_Office', ...
%                            'ARGO_salinity_Met_Office', ...
%                            'CTD_temperature_CalCOFI', ...
%                            'CTD_salinity_CalCOFI', ...
%                            'CTD_temperature_GLOBEC', ...
%                            'CTD_salinity_GLOBEC'};
%
%    S.state_variables=[newline, ...
%      '1: free-surface (m) ', newline, ...
%      '2: vertically integrated u-momentum component (m/s) ', newline, ...
%      '3: vertically integrated v-momentum component (m/s) ', newline, ...
%      '4: u-momentum component (m/s) ', newline, ...
%      '5: v-momentum component (m/s) ', newline, ...
%      '6: potential temperature (Celsius) ', newline, ...
%      '7: salinity (nondimensional)'];
%
%    S.obs_provenance=[newline, ...
%       ' 1: gridded AVISO sea level anomaly ', newline, ...
%       ' 2: blended satellite SST ', newline, ...
%       ' 3: XBT temperature from Met Office ', newline, ...
%       ' 4: CTD temperature from Met Office ', newline, ...
%       ' 5: CTD salinity from Met Office ', newline, ...
%       ' 6: ARGO floats temperature from Met Office ', newline, ...
%       ' 7: ARGO floats salinity from Met Office ', newline, ...
%       ' 8: CTD temperature from CalCOFI ', newline, ...
%       ' 9: CTD salinity from CalCOFI ', newline, ...
%       '10: CTD temperature from GLOBEC ', newline, ...
%       '11: CTD salinity from GLOBEC'];
%
  
% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Get error covariance standard deviation creation parameters.
%--------------------------------------------------------------------------

switch numel(varargin)
  case 0
    if (isfield(S,'ncfile'))
      ncfile=S.ncfile;
    else
      error(['C_OBSERVATIONS - Cannot find file name field: ncfile, ',  ...
             'in structure array S']);
    end
  case 1
    ncfile=varargin{1};    
end

%--------------------------------------------------------------------------
% Set dimensions.
%--------------------------------------------------------------------------

Dname.survey ='survey';          Dsize.survey = S.Nsurvey;
Dname.state  ='state_variable';  Dsize.state  = S.Nstate;
Dname.datum  ='datum';           Dsize.datum  = S.Ndatum;

%--------------------------------------------------------------------------
% Create 4D-Var observation NetCDF file.
%--------------------------------------------------------------------------

disp(' ');
disp(['*** Creating observations file:  ', ncfile]);

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode, netcdf.getConstant('64BIT_OFFSET'));

ncid = netcdf.create(ncfile,mode);

%--------------------------------------------------------------------------
% Define dimensions.
%--------------------------------------------------------------------------

Did.survey = netcdf.defDim(ncid, 'survey',         S.Nsurvey);
Did.state  = netcdf.defDim(ncid, 'state_variable', S.Nstate);
Did.datum  = netcdf.defDim(ncid, 'datum',          S.Ndatum);

%--------------------------------------------------------------------------
% Define global attributes.
%--------------------------------------------------------------------------

varid = netcdf.getConstant('nc_global');

netcdf.putAtt(ncid, varid, 'type', 'ROMS Observations');

if (isfield(S,'title'))
  netcdf.putAtt(ncid, varid, 'title', strip(S.title));
end

netcdf.putAtt(ncid, varid, 'Conventions', 'CF-1.4');

if (isfield(S,'grd_file'))
  netcdf.putAtt(ncid, varid, 'grd_file', strip(S.grd_file));
end
  
if (isfield(S,'grid_Lm_Mm_N')),
  netcdf.putAtt(ncid, varid, 'grid_Lm_Mm_N', int32(S.grid_Lm_Mm_N));
end
  
if (isfield(S,'global_variables'))
  netcdf.putAtt(ncid, varid, 'global_variables', S.global_variables);
end
  
if (isfield(S,'global_provenance'))
  netcdf.putAtt(ncid, varid, 'global_provenance', S.global_provenance);
end

if (isfield(S,'global_sources'))
  netcdf.putAtt(ncid, varid, 'global_sources', S.global_sources);
end  
  
netcdf.putAtt(ncid, varid, 'history',                                   ...
              ['4D-Var observations, ', date_stamp]);

%--------------------------------------------------------------------------
% Define variables.
%--------------------------------------------------------------------------

% Spherical switch.

varid = netcdf.defVar(ncid, 'spherical',                                ...
                      netcdf.getConstant('nc_int'), []);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'grid type logical switch') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'Cartesian spherical');

% Define number of observation per survey.

varid = netcdf.defVar(ncid, 'Nobs',                                     ...
                      netcdf.getConstant('nc_int'), Did.survey);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'number of observations with the same survey time');

% Define observation survey time.

varid = netcdf.defVar(ncid, 'survey_time',                              ...
                      netcdf.getConstant('nc_double'), Did.survey);
netcdf.putAtt(ncid, varid, 'long_name', 'survey time');
netcdf.putAtt(ncid, varid, 'units', 'days');

% Define observation global variance (obsolete).

varid = netcdf.defVar(ncid, 'obs_variance',                             ...
                      netcdf.getConstant('nc_double'), Did.state);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'global temporal and spatial observation variance');
netcdf.putAtt(ncid, varid, 'units', 'squared state variable units');

% Define observation type.

varid = netcdf.defVar(ncid, 'obs_type',                                 ...
                      netcdf.getConstant('nc_int'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'model state variable associated with observations');
if (isfield(S,'state_flag_values'))
  netcdf.putAtt(ncid, varid, 'flag_values', int32(S.state_flag_values));
end
if (isfield(S,'state_flag_meanings'))
  netcdf.putAtt(ncid, varid, 'flag_meanings', S.state_flag_meanings);
end

% Define observation provenance.

varid = netcdf.defVar(ncid, 'obs_provenance',                           ...
                      netcdf.getConstant('nc_int'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name', 'observation origin');
if (isfield(S,'origin_flag_values'))
  netcdf.putAtt(ncid, varid, 'flag_values', int32(S.origin_flag_values));
end
if (isfield(S,'origin_flag_meanings'))
  netcdf.putAtt(ncid, varid, 'flag_meanings', S.origin_flag_meanings);
end

% Define observation time.

varid = netcdf.defVar(ncid, 'obs_time',                                 ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name', 'time of observation');
netcdf.putAtt(ncid, varid, 'units', 'days');

% Define observation longitude.

if (isfield(S,'lon'))
  varid = netcdf.defVar(ncid, 'obs_lon',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
  netcdf.putAtt(ncid, varid, 'long_name', 'observation longitude');
  netcdf.putAtt(ncid, varid, 'units', 'degrees_east');
  netcdf.putAtt(ncid, varid, 'standard_name', 'longitude');
end

% Define observation latitude.

if (isfield(S,'lat'))
  varid = netcdf.defVar(ncid, 'obs_lat',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
  netcdf.putAtt(ncid, varid, 'long_name', 'observation latitude');
  netcdf.putAtt(ncid, varid, 'units', 'degrees_north');
  netcdf.putAtt(ncid, varid, 'standard_name', 'latitude');
end

% Define observation depth.

varid = netcdf.defVar(ncid, 'obs_depth',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name', 'depth of observation');
netcdf.putAtt(ncid, varid, 'units', 'meters');
netcdf.putAtt(ncid, varid, 'negative', 'downwards');

% Define observation fractional grid X-location.

varid = netcdf.defVar(ncid, 'obs_Xgrid',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'observation fractional x-grid location');

% Define observation fractional grid Y-location.

varid = netcdf.defVar(ncid, 'obs_Ygrid',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'observation fractional y-grid location');

% Define observation fractional grid Z-location.

varid = netcdf.defVar(ncid, 'obs_Zgrid',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'observation fractional z-grid location');

% Define observation error covariance.

varid = netcdf.defVar(ncid, 'obs_error',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name', 'observation error covariance');
netcdf.putAtt(ncid, varid, 'units', 'squared state variable units');

% Define observation value.

varid = netcdf.defVar(ncid, 'obs_value',                                ...
                      netcdf.getConstant('nc_double'), Did.datum);
netcdf.putAtt(ncid, varid, 'long_name', 'observation value');
netcdf.putAtt(ncid, varid, 'units', 'associated state variable units');

% Define observation meta value: additional data qualifiers for the
% extra-observation value specified in obs_value, which can be used
% in complex operators like HF radar radials, etc.

if (isfield(S,'meta'))
  varid = netcdf.defVar(ncid, 'obs_meta',                               ...
                        netcdf.getConstant('nc_double'), Did.datum);
  netcdf.putAtt(ncid, varid, 'long_name', 'observation meta value');
  netcdf.putAtt(ncid, varid, 'units', 'associated state variable units');
end

%--------------------------------------------------------------------------
% Close NetCDF file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);

netcdf.close(ncid);

return

