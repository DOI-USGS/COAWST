function [S]=obs_depth(ncfile, S, zflag);

%
% OBS_DEPTH:  Computes the observations fractional z-grid locations
%
% [S]=obs_depth(ncfile, S, zflag)
%
% This function computes the observation fractional z-grid locations. It
% removes vertical outliers. The application vertical grid parameters are
% extracted from provided ROMS history file.
%
% On Input:
%
%    ncfile  ROMS history file name (string)
%
%    S         Observations structure data or NetCDF file name:
%
%    zflag   Flag to compute ROMS depths (integer):
%
%              zflag = 0        use zero free-surface
%              zflag > 0        read zflag history record for free-surface
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

% svn $Id: obs_depth.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Read in observation if S is a 4D-Var observation NetCDF file.

if (ischar(S)),
  S=obs_read(S);
end,

%  Check if 'provenace', 'lon', and 'lat' fields are available.

has.lonlat=false;
if (isfield(S,'lon') && isfield(S,'lat')),
  has.lonlat=true;
end,

has.provenance=false;
if (isfield(S,'provenance')),
  has.provenance=true;
end,

%  Get number of data surveys.

Nsurvey=length(S.survey_time);

%  Set observations dynamical fields (cell array) in structure, S.

field_list = {'type', 'time', 'Xgrid', 'Ygrid', 'Zgrid', 'depth', ...
              'error', 'value'};

if (has.lonlat),
  field_list = [field_list, 'lon', 'lat'];
end,

if (has.provenance),
  field_list = [field_list, 'provenance'];
end,

%----------------------------------------------------------------------------
%  Compute depths of ROMS application.
%----------------------------------------------------------------------------

igrid = 5;                                  % depth of W-points
tindex = abs(zflag);                        % use zflag record for zeta

switch ( zflag ),
  case 0,
    Zw=depths(ncfile, ncfile, igrid, 0, 0);
  otherwise,
    Zw=depths(ncfile, ncfile, igrid, 0, tindex);
end,

[Lr Mr Nw]=size(Zw);
Nr = Nw -1;

%----------------------------------------------------------------------------
%  Compute observations fractional z-grid location. Use W-points depths
%  to include observations in the bottom half and top half of the grid
%  cells.
%----------------------------------------------------------------------------

%  Initialize.

S.Zgrid = S.depth;

%  Check zero values. It is assume that zero value maybe assigned to
%  surface observations like satellite data.

ind = find(S.depth == 0);

if (~isempty(ind)),
  S.Zgrid(ind) = Nr;
end,

%  Interpolate negative depth values to model fractional z-grid location.

ind = find(S.depth < 0);

if (~isempty(ind));
  for n=1:length(ind),
    iobs = ind(n);
    I = 1.0 + floor(S.Xgrid(iobs));
    J = 1.0 + floor(S.Ygrid(iobs));
    z = reshape(Zw(I,J,:),1,Nw);
    S.Zgrid(iobs) = interp1(z, [0:1:Nr], S.depth(iobs));
  end,
end,

%  Since the state variables are at vertical RHO-points, assign bottom half
%  locations to k=1.

ind = find((0.0 < S.Zgrid) & (S.Zgrid <= 0.5));
if (~isempty(ind));
  S.Zgrid(ind) = 1.0;
end,

%  Since the state variables are at vertical RHO-points, assign top half
%  locations to k=Km-1.

ind = find((Nr < S.Zgrid) & (S.Zgrid <= Nr+0.5));
if (~isempty(ind));
  S.Zgrid(ind) = Nr;
end,

%  Remove NaNs from the vertical interpolation in all dynamical fields
%  in the observation structure.

ind = find(isnan(S.Zgrid));

if (~isempty(ind)),
  for value = field_list,
    field = char(value);             % convert from cell to string
    S.(field)(ind) = [];             % remove outliers
  end,

  disp(' ');
  disp([' Number of vertical outlier observations removed = ', ...
        num2str(length(ind))]);
  disp(' ');
  
  for i=1:S.Nsurvey,
    ind=find(S.time == S.survey_time(i));
    S.Nobs(i)=length(ind);
  end,
end,

S.Ndatum = length(S.time);

return
