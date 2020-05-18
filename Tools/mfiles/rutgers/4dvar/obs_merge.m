function [S]=obs_merge(files, dt);

%
% OBS_MERGE:  Merges data from several 4D-Var observation NetCDF files
%
% [S]=obs_merge(files, dt)
%
% Given a cell array of 4D-Var NetCDF observation files, this function
% combines them into a single observation data structure.
%
% On Input:
%
%    files   Observations NetCDF file names (string cell array)
%
%    dt      Minimum amount of time (days) for merging surveys from input
%              files. It must be greater than the ROMS application time
%              step. Surveys that are closer than "dt" together will be
%              combined.
%
% On Output:
%
%    S       Observations data (structure array):
%
%              S.ncfile         NetCDF file name (string)
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
% NOTE: On exit, a new structure, S, is returned which can be written to
%       disk using 'c_observations' and 'obs_write'.
%

% svn $Id: obs_merge.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license           Brian Powell            %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Initialize internal parameters.

check_time=true;
if (nargin < 2),
  check_time=false;
end,

Nfiles=length(files);

%----------------------------------------------------------------------------
%  Process each observation file and build observation structure.
%----------------------------------------------------------------------------

%  Set observations dynamical fields (cell array) in structure, S.

field_list = {'type', 'time', 'Xgrid', 'Ygrid', 'Zgrid', 'depth', ...
              'error', 'value'};

%  Initialize structure with first file.

[S]=obs_read(char(files(1)));

%  If appropriate, add other dynamical fields to structure list.

if (isfield(S,'lon') && isfield(S,'lat')),
  field_list = [field_list, 'lon', 'lat'];
end,

if (isfield(S,'provenance')),
  field_list = [field_list, 'provenance'];
end,

%  Proccess remaining files.

for m=2:Nfiles,

  [N]=obs_read(char(files(m)));

%  Add new data to S structure.

  for value = field_list,
    field = char(value);
    S.(field) = [S.(field), N.(field)];
  end,
 
  S.Nobs        = [S.Nobs,         N.Nobs       ];
  S.survey_time = [S.survey_time,  N.survey_time];
  S.variance    = (S.variance    + N.variance   ) / 2.0;

%  Order according to ascending 'survey_time'.

  [l,ind]=sort(S.survey_time);

  S.Nobs        = S.Nobs(ind);
  S.survey_time = S.survey_time(ind);

  S.Nsurvey     = length(S.Nobs);

%  Order according to ascending 'obs_time'.

  [l,ind]=sort(S.time);

  for value = field_list,
    field = char(value);
    S.(field) = S.(field)(ind);
  end,

end,

%  Loop over 'survey_times' until they are well-spaced.

if (check_time),

  while (true),
    ntime = S.survey_time;
    sdt = diff(S.survey_time);
    l = find(sdt ~= 0 & sdt < dt);
    if (isempty( l )),
      break;
    end,
    otime = S.time;
    for i=1:length(l),
      lo = find(S.time == S.survey_time(l(i)   ) | ...
                S.time == S.survey_time(l(i)+1));
      ntime(l(i):l(i)+1) = nanmedian(S.time(lo));
      otime(lo) = ntime(l(i));
    end,
    S.time = otime;
    S.survey_time = ntime;
  end,

end,

%  Combine all surveys with similar times.

[l,i] = unique(S.survey_time);

if (length(l) ~= length(S.survey_time)),
  survey_time = S.survey_time;
  S.survey_time = l;
  Nobs = S.Nobs;
  S.Nobs = [];
  for i=1:length(l)
    m = find(survey_time == l(i));
    S.Nobs(i,1) = sum(Nobs(m));
  end,
end,

S.Nsurvey = length(S.survey_time);

%  Sort the observations by the 'obs_type' for each survey.

for i=1:S.Nsurvey,

  l = find(S.time == S.survey_time(i));
  [m,i]=sort(S.type(l));

  for value = field_list,
    field = char(value);
    S.(field)(l) = S.(field)(l(i));
  end,

end,

% Lastly, remove any NaN's that may have slipped through.

l=find(isnan(S.value) | isnan(S.time)   | ...
       isnan(S.Xgrid) | isnan(S.Ygrid)  | ...
       isnan(S.Zgrid) | isnan(S.depth));

if (~isempty(l)),
  for value = field_list,
    field = char(value);
    S.(field)(l) = [];
  end,

  for i=1:S.Nsurvey,
    ind=find(S.time == S.survey_time(i));
    S.Nobs(i)=length(ind);
  end,
end,

S.Ndatum = length(S.time);
S.files  = files;

return
