function [S]=plot_super(GRDfile, OBSfile, state_var, survey)
  
%
% PLOT_SUPER:  Computes and plot super observations
%
% [S]=plot_super(GRDfile, OBSfile, state_var, survey)
%
% This script can be used to check the effect of binning 4D-Var
% observations after processing with 'super_obs'.
%
% On Input:
%
%    GRDfile      Application NetCDF grid file name (string)
%
%    OBSfile      Observation NetCDF file name (string)
%
%    state_var    Associate state variable:
%
%                   state_var = 1     free-surface
%                   state_var = 2     vertically integrated u-momentum
%                   state_var = 3     vertically integrated v-momentum
%                   state_var = 4     u-momentum
%                   state_var = 5     v-momentum
%                   state_var = 6     temperature
%                   state_var = 7     salinity
%
%    survey       Survey record in observation file to plot (OPTIONAL).
%                   If not provided, it will plot the first survey.
%
% On Output:
%
%    S            Binned observations data (structure array):
%
%                   S.ncfile       NetCDF file name (string)
%                   S.Ndatum       total number of observations
%                   S.spherical    spherical grid switch
%                   S.Nobs         number of observations per survey
%                   S.survey_time  time for each survey time
%                   S.variance     global variance per state variable
%                   S.type         state variable associated with observation
%                   S.time         time for each observation
%                   S.depth        depth of observation
%                   S.Xgrid        observation fractional x-grid location
%                   S.Ygrid        observation fractional y-grid location
%                   S.Zgrid        observation fractional z-grid location
%                   S.error        observation error
%                   S.value        observation value
%
%                   S.std          binning observation standard deviation
%
%                 The following optional variables will be processed if
%                 available:
%
%                   S.provenance   observation origin
%                   S.lon          observation longitude
%                   S.lat          observation latitude
%

% svn $Id: plot_super.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

if (nargin < 3),
  error([' PLOT_SUPER: You must specify a grid and observation file' ...
         ' along with starting and ending times']);
end,

if (nargin < 4),
  survey = 1;                % plot first survey
end,
  
%  Read in application grid.

V = nc_vnames(GRDfile);
nvars = length(V.Variables);

got_coast = false;

for n=1:nvars,
  name = char(V.Variables(n).Name);
  switch name
    case 'lon_rho'
      rlon = nc_read(GRDfile, 'lon_rho');
    case 'lat_rho'
      rlat = nc_read(GRDfile, 'lat_rho');
    case 'lon_coast'
      clon = nc_read(GRDfile, 'lon_coast');
      clat = nc_read(GRDfile, 'lat_coast');
      got_coast = true;
  end,
end,

%  Read in provided observations into structure.

[O] = obs_read(OBSfile);

%  Create super observations, if necessary.

[S] = super_obs(O);

%  Plot original (O) and new (S) observations.

Nsurvey = length(O.survey_time);        % Number of surveys

survey_time = O.survey_time(survey);    % select requested survey

ind_o = find(O.type == state_var & O.time == survey_time);

if (~isempty(ind_o)),
  Olon = O.lon(ind_o);
  Olat = O.lat(ind_o);
  Oval = O.value(ind_o);
else
  disp(' ');
  disp(['   PlOT_SUPER: No observations found for state variable = ', ...
        num2str(state_var)]);
  disp(['               at survey ', num2str(survey), ' = ', ...
        num2str(survey_time)]);
  return
end,

ind_s = find(S.type == state_var & S.time == survey_time);

if (~isempty(ind_s)),
  Slon = S.lon(ind_s);
  Slat = S.lat(ind_s);
  Sval = S.value(ind_s);
end,

pcolor(rlon,rlat,ones(size(rlon)).*NaN);
hold on;
if (got_coast),
  plot(clon,clat,'k-');
end,
scatter(Olon,Olat,30,Oval,'ko');
scatter(Slon,Slat,30,Sval,'d','filled');
colorbar('fontsize',10,'fontweight','bold');
hold off;
title('Original (black circles), Super Obs (filled diamonds)', ...
      'fontsize',10,'fontweight','bold');

return
