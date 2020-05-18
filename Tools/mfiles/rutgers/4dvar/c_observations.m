function [status]=c_observations(S,file)

%
% C_OBSERVATIONS:  Creates 4D-Var observations NetCDF file
%
% [status]=c_observations(S,file)
%
% This function creates ROMS 4D-Var observation NetCDF file using specified
% in structure array, S.
%
% On Input:
%
%    S     Observations file creation parameters (structure array):
%
%            S.ncfile                NetCDF file name (string)
%            S.spherical             spherical grid switch
%            S.Nstate                number of state variables
%            S.Nsurvey               number of data surveys 
%            S.Ndatum                number of observations
%            S.state_flag_values     obs_type 'flag_values' attribute
%            S.state_flag_meanings   obs_type 'flag_meanings attribute
%            S.origin_flag_values    obs_provenance 'flag_values' attribute
%            S.origin_flag_meanings  obs_provenance 'flag_meaning' attribute
%
%          Optional fields:
%
%            S.title                 application title
%            S.grd_file              'grd_file' global attribute
%            S.global_variables      'state_variables' global attribute
%            S.global_provenance     'obs_provenance' global attribute
%            S.global_sources        'obs_sources' global attribute
%
%            S.do_longitude          switch to define longitude locations
%            S.do_latitude           switch to define latitude  locations
%
%    file  Output observation file name (string, OPTIONAL). If specified,
%            it creates this file instead the one in S.ncfile.
%
% On Output:
%
%    status  Error flag.
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
  
% svn $Id: c_observations.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  Get error covariance standard deviation creation parameters.
%--------------------------------------------------------------------------

if (nargin < 2),
  if (isfield(S,'ncfile')),
    ncfile=S.ncfile;
  else
    error(['C_OBSERVATIONS - Cannot find file name field: ncfile, ',    ...
          'in structure array S']);
  end
else
  ncfile=file;    
end

if (isfield(S,'spherical')),
  spherical=S.spherical;
else
  spherical=0;
end

%--------------------------------------------------------------------------
%  Set dimensions.
%--------------------------------------------------------------------------

Dname.survey ='survey';          Dsize.survey = S.Nsurvey;
Dname.state  ='state_variable';  Dsize.state  = S.Nstate;
Dname.datum  ='datum';           Dsize.datum  = S.Ndatum;

%--------------------------------------------------------------------------
%  Set all posible variables names.
%--------------------------------------------------------------------------

Vname.spherical  = 'spherical';
Vname.Nobs       = 'Nobs';
Vname.survey     = 'survey_time';
Vname.variance   = 'obs_variance';
Vname.type       = 'obs_type';
Vname.provenance = 'obs_provenance';
Vname.time       = 'obs_time';
Vname.depth      = 'obs_depth';
Vname.lon        = 'obs_lon';
Vname.lat        = 'obs_lat';
Vname.Xgrid      = 'obs_Xgrid';
Vname.Ygrid      = 'obs_Ygrid';
Vname.Zgrid      = 'obs_Zgrid';
Vname.error      = 'obs_error';
Vname.value      = 'obs_value';

%--------------------------------------------------------------------------
%  Create 4D-Var observation NetCDF file.
%--------------------------------------------------------------------------

disp(' ');
disp(['*** Creating observations file:  ', ncfile]);

[ncid,status]=mexnc('create',ncfile,'nc_clobber');
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_OBSERVATIONS: CREATE - unable to create file: ', ncfile,    ...
	 sprintf('\n'), blanks(16), 'In current directory: ', pwd]);
end

%--------------------------------------------------------------------------
%  Define dimensions.
%--------------------------------------------------------------------------

[did.survey,status]=mexnc('def_dim',ncid,Dname.survey,Dsize.survey); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: DEF_DIM - unable to define dimension: ',     ...
	  Dname.survey]);
end

[did.state,status]=mexnc('def_dim',ncid,Dname.state,Dsize.state); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: DEF_DIM - unable to define dimension: ',     ...
	  Dname.state]);
end

[did.datum,status]=mexnc('def_dim',ncid,Dname.datum,Dsize.datum); 
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: DEF_DIM - unable to define dimension: ',     ...
	  Dname.datum]);
end

%--------------------------------------------------------------------------
%  Create global attributes.
%--------------------------------------------------------------------------

str='ROMS observations';
lstr=length(str);
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'type',nc_constant('nc_char'),lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',  ...
        ' type.']);
end

if (isfield(S,'title')),
  lstr=length(S.title);
  [status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),          ...
                 'title',nc_constant('nc_char'),lstr,S.title);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',...
          ' title.']);
  end
end

str='CF-1.4';
lstr=length(str);
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'Conventions',nc_constant('nc_char'),lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:', ...
	  ' Conventions.']);
end

if (isfield(S,'grd_file')),
  lstr=length(S.grd_file);
  [status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),          ...
                 'grd_file',nc_constant('nc_char'),lstr,S.grd_file);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',...
          ' grd_file.']);
  end
end

if (isfield(S,'grid_Lm_Mm_N')),
  nval=length(S.grid_Lm_Mm_N);
  [status]=mexnc('put_att_int',ncid,nc_constant('nc_global'),           ...
                 'grid_Lm_Mm_N',nc_constant('nc_int'),nval,             ...
		 int32(S.grid_Lm_Mm_N));
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',...
           ' grd_Lm_Mm_N.']);
  end
end

if (isfield(S,'global_variables')),
  lstr=length(S.global_variables);
  [status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),          ...
                 'state_variables',nc_constant('nc_char'),              ...
                 lstr,S.global_variables);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',...
           ' state_variables.']);
  end
end

if (isfield(S,'global_provenance')),
  lstr=length(S.global_provenance);
  [status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),          ...
                 'obs_provenance',nc_constant('nc_char'),               ...
                 lstr,S.global_provenance);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',...
           ' obs_provenance.']);
  end
end

str='squared state variable units';
lstr=length(str);
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'variance_units',nc_constant('nc_char'),lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribure:',  ...
        ' variance_units.']);
end

if (isfield(S,'global_sources')),
  lstr=length(S.global_sources);
  [status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),          ...
                 'obs_sources',nc_constant('nc_char'),                  ...
                 lstr,S.global_sources);
  if (status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',...
	  ' obs_sources.']);
  end
end

str=['4D-Var observations, ',date_stamp];
lstr=length(str);
[status]=mexnc('put_att_text',ncid,nc_constant('nc_global'),            ...
               'history',nc_constant('nc_char'),lstr,str);
if (status ~= 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_OBSERVATIONS: PUT_ATT_TEXT - unable to global attribute:',  ...
        ' history.']);
end

%--------------------------------------------------------------------------
%  Define configuration variables.
%--------------------------------------------------------------------------

% Define spherical switch.

Var.name          = Vname.spherical;
Var.type          = nc_constant('nc_int');
Var.dimid         = [];
Var.long_name     = 'grid type logical switch';
Var.flag_values   = [0 1];
Var.flag_meanings = ['Cartesian', blanks(1), ...
                     'spherical'];
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

% Define observation variables.

Var.name          = Vname.Nobs;
Var.type          = nc_constant('nc_int');
Var.dimid         = [did.survey];
Var.long_name     = 'number of observations with the same survey time';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.survey;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.survey];
Var.long_name     = 'survey time';
Var.units         = 'days';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.variance;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.state];
Var.long_name     = 'global temporal and spatial observation variance';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.type;
Var.type          = nc_constant('nc_int');
Var.dimid         = [did.datum];
Var.long_name     = 'model state variable associated with observations';
Var.flag_values   = S.state_flag_values;
Var.flag_meanings = S.state_flag_meanings;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.provenance;
Var.type          = nc_constant('nc_int');
Var.dimid         = [did.datum];
Var.long_name     = 'observation origin';
Var.flag_values   = S.origin_flag_values;
Var.flag_meanings = S.origin_flag_meanings;
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.time;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'time of observation';
Var.units         = 'days';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

if (isfield(S,'do_longitude')),
  if (S.do_longitude),
    Var.name          = Vname.lon;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.datum];
    Var.long_name     = 'observation longitude';
    Var.units         = 'degrees_east';
    Var.standard_name = 'longitude';
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_latitude')),
  if (S.do_latitude),
    Var.name          = Vname.lat;
    Var.type          = nc_constant('nc_double');
    Var.dimid         = [did.datum];
    Var.long_name     = 'observation latitude';
    Var.units         = 'degrees_north';
    Var.standard_name = 'latitude';
    [~,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

Var.name          = Vname.depth;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'depth of observation';
Var.units         = 'meters';
Var.negative      = 'downwards';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.Xgrid;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'observation fractional x-grid location';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.Ygrid;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'observation fractional y-grid location';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.Zgrid;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'observation fractional z-grid location';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.error;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'observation error covariance';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name          = Vname.value;
Var.type          = nc_constant('nc_double');
Var.dimid         = [did.datum];
Var.long_name     = 'observation value';
[~,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%--------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%--------------------------------------------------------------------------

[status]=mexnc('enddef',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error('C_OBSERVATIONS: ENDDEF - unable to leave definition mode.');
end

[status]=mexnc('close',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error([ 'C_OBSERVATIONS: CLOSE - unable to close NetCDF file: ',      ...
         ncname]);
end

return

