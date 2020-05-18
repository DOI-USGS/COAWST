function [S]=obs_read(ncfile);

%
% OBS_READ:  Reads ROMS 4D-Var observation NetCDF file
%
% [S]=obs_read(ncfile)
%
% This function reads ROMS 4D-Var observation NetCDF file and stores all
% the variables in structure array, S.
%
% On Input:
%
%    ncfile  4D-Var Observations NetCDF file name (string)
%
% On Output:
%
%    S       Observations data (structure array):
%
%              S.ncfile         NetCDF file name (string)
%              S.Nsurvey        number of observations surveys
%              S.Nstate         number of state variables
%              S.Ndatum         total number of observations
%              S.spherical      spherical grid switch
%              S.Nobs           number of observations per survey
%              S.survey_time    time for each survey time
%              S.variance       global variance per state variable
%              S.type           state variable associated with observation
%              S.time           time for each observation
%              S.depth          depth of observation
%              S.Xgrid          observation fractional x-grid location
%              S.Ygrid          observation fractional y-grid location
%              S.Zgrid          observation fractional z-grid location
%              S.error          observation error
%              S.value          observation value
%
%            The following optional variables will be read if available:
%
%              S.provenance     observation origin
%              S.lon            observation longitude
%              S.lat            observation latitude
% 
%            The following variable attributes will be read if available:
%
%              S.state_flag_values     obs_type 'flag_values' attribute
%              S.state_flag_meanings   obs_type 'flag_meanings attribute
%              S.origin_flag_values    obs_provenance 'flag_values' attribute
%              S.origin_flag_meanings  obs_provenance 'flag_meaning' attribute
%
%            The following global attributes will be read if available:
%
%              S.global_variables      'state_variables' global attribute
%              S.global_provenance     'obs_provenance' global attribute
%              S.global_sources        'obs_sources' global attribute
%

% svn $Id: obs_read.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Initialize output structure.

S=[];

%----------------------------------------------------------------------------
%  Inquire input observations NetCDF file.
%----------------------------------------------------------------------------

I=nc_inq(ncfile);

%----------------------------------------------------------------------------
%  Read in all available variables.
%----------------------------------------------------------------------------

V=nc_vnames(ncfile);
nvars=length(V.Variables);

S.ncfile=ncfile;
S.Nsurvey=0;
S.Nstate=0;
S.Ndatum=0;

got.provenance=0;
got.label=0;

%  Insure that read vector variables have the singleton in the first
%  dimension to allow vector concatenation when merging several datasets.

for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch name
    case 'spherical'
      S.spherical=nc_read(ncfile,'spherical');
    case 'Nobs'
      S.Nobs=nc_read(ncfile,'Nobs');
      if (size(S.Nobs,1) > 1),
        S.Nobs=transpose(S.Nobs);
      end,
    case 'survey_time'
      S.survey_time=nc_read(ncfile,'survey_time');
      if (size(S.survey_time,1) > 1),
        S.survey_time=transpose(S.survey_time);
      end,
      S.Nsurvey=length(S.survey_time);
    case 'obs_variance'
      S.variance=nc_read(ncfile,'obs_variance');
      if (size(S.variance,1) > 1),
        S.variance=transpose(S.variance);
      end,
      S.Nstate=length(S.variance);
    case 'obs_type'
      S.type=nc_read(ncfile,'obs_type');
      if (size(S.type,1) > 1),
        S.type=transpose(S.type);
      end,
    case 'obs_provenance'
      S.provenance=nc_read(ncfile,'obs_provenance');
      if (size(S.provenance,1) > 1),
        S.provenance=transpose(S.provenance);
      end,
      got.provenance=1;
    case 'obs_label'
      S.label=nc_read(ncfile,'obs_label');
      if (size(S.label,1) > 1),
        S.label=transpose(S.label);
      end,
      got.label=1;
    case 'obs_time'
      S.time=nc_read(ncfile,'obs_time');
      if (size(S.time,1) > 1),
        S.time=transpose(S.time);
      end,
    case 'obs_lon'
      S.lon=nc_read(ncfile,'obs_lon');
      if (size(S.lon,1) > 1),
        S.lon=transpose(S.lon);
      end,
    case 'obs_lat'
      S.lat=nc_read(ncfile,'obs_lat');
      if (size(S.lat,1) > 1),
        S.lat=transpose(S.lat);
      end,
    case 'obs_depth'
      S.depth=nc_read(ncfile,'obs_depth');
      if (size(S.depth,1) > 1),
        S.depth=transpose(S.depth);
      end,
    case 'obs_Xgrid'
      S.Xgrid=nc_read(ncfile,'obs_Xgrid');
      if (size(S.Xgrid,1) > 1),
        S.Xgrid=transpose(S.Xgrid);
      end,
    case 'obs_Ygrid'
      S.Ygrid=nc_read(ncfile,'obs_Ygrid');
      if (size(S.Ygrid,1) > 1),
        S.Ygrid=transpose(S.Ygrid);
      end,
    case 'obs_Zgrid'
      S.Zgrid=nc_read(ncfile,'obs_Zgrid');
      if (size(S.Zgrid,1) > 1),
        S.Zgrid=transpose(S.Zgrid);
      end,
    case 'obs_error'
      S.error=nc_read(ncfile,'obs_error');
      if (size(S.error,1) > 1),
        S.error=transpose(S.error);
      end,
    case 'obs_value'
      S.value=nc_read(ncfile,'obs_value');
      if (size(S.value,1) > 1),
        S.value=transpose(S.value);
      end,
      S.Ndatum=length(S.value);
  end,
end,

%----------------------------------------------------------------------------
%  Read in relevant variable and global attributes.
%----------------------------------------------------------------------------

%  Read in 'flag_values' and 'flag_meanings attributes for variable
%  'obs_type'.

ivar=strcmp({I.Variables.Name},'obs_type');
if (strcmp({I.Variables(ivar).Attributes.Name},'flag_values')),
  Avalue=nc_getatt(ncfile,'flag_values','obs_type');
  if (~isempty(Avalue)),
    S.state_flag_values=Avalue;
  end
end

ivar=strcmp({I.Variables.Name},'obs_type');
if (strcmp({I.Variables(ivar).Attributes.Name},'flag_meanings')),
  Avalue=nc_getatt(ncfile,'flag_meanings','obs_type');
  if (~isempty(Avalue)),
    S.state_flag_meanings=Avalue;
  end
end

%  Read in 'flag_values' and 'flag_meanings attributes for variable
%  'obs_provenace'.

if (got.provenance), 
  ivar=strcmp({I.Variables.Name},'obs_provenance');
  if (strcmp({I.Variables(ivar).Attributes.Name},'flag_values')),
    Avalue=nc_getatt(ncfile,'flag_values','obs_provenance');
    if (~isempty(Avalue)),
      S.origin_flag_values=Avalue;
    end
  end

  ivar=strcmp({I.Variables.Name},'obs_provenance');
  if (strcmp({I.Variables(ivar).Attributes.Name},'flag_meanings')),
    Avalue=nc_getatt(ncfile,'flag_meanings','obs_provenance');
    if (~isempty(Avalue)),
      S.origin_flag_meanings=Avalue;
    end
  end
end
  
% Read in 'title' global attribute.

if (strcmp({I.Attributes.Name},'title')),
  Avalue=nc_getatt(ncfile,'title');
  if (~isempty(Avalue)),
    S.title=Avalue;
  end
end

% Read in 'grd_file' global attribute.

if (strcmp({I.Attributes.Name},'grd_file')),
  Avalue=nc_getatt(ncfile,'grd_file');
  if (~isempty(Avalue)),
    S.grd_file=Avalue;
  end
end

% Read in 'grid_Lm_Mm_N' global attribute.

if (strcmp({I.Attributes.Name},'grd_Lm_Mm_N')),
  Avalue=nc_getatt(ncfile,'grid_Lm_Mm_N');
  if (~isempty(Avalue)),
    S.grid_Lm_Mm_N=int32(Avalue);
  end
end

% Read in 'state_variables' global attribute.

if (strcmp({I.Attributes.Name},'state_variables')),
  Avalue=nc_getatt(ncfile,'state_variables');
  if (~isempty(Avalue)),
    S.global_variables=Avalue;
  end
end

% Read in 'obs_provenance' global attribute.

if (strcmp({I.Attributes.Name},'obs_provenance')),
  Avalue=nc_getatt(ncfile,'obs_provenance');
  if (~isempty(Avalue)),
    S.global_provenance=Avalue;
  end
end

% Read in 'obs_sources' global attribute.

if (strcmp({I.Attributes.Name},'obs_sources')),
  Avalue=nc_getatt(ncfile,'obs_sources');
  if (~isempty(Avalue)),
    S.global_sources=Avalue;
  end
end

return
