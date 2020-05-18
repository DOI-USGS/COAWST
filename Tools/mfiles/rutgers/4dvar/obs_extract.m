function [out_files]=obs_extract(time_interval, obs_inp, obs_out, report)

% OBS_EXTRACT:  Extacts and creates observation NetCDF files
%
% [obs_files]=obs_extract(time_interval, obs_inp, obs_out, report)
%
% Extracts data from input observations NetCDF file at the
% requested time interval (days) and creates new observation
% NetCDF files at such intervals
%
% On Input:
%
%    time_interval    Sampling time interval for extraction (days)
%                       For esample, use time_interval=1.0 for daily
%
%    obs_inp          Input observation NetCDF file name to process
%                       (string)
%
%    obs_out          Output observation NetCDF base name (string)
%                       It is concatenated with day to create new files
%
%    report           Report extracting information (logical; OPTIONAL)
%
% On Output:
%
%    out_files        Created output file names (cell array)
%
% Usage:
%
%    out_files = obs_extract(1, 'obs_37623.nc', 'obs_new', true);
%
  
% svn $Id: obs_extract.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2017 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if (nargin < 4)
  report = false;
end
  
%--------------------------------------------------------------------------
% Inquire input observations NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(obs_inp);

%--------------------------------------------------------------------------
% Get input observation structure.
%--------------------------------------------------------------------------

S = obs_read(obs_inp);

%--------------------------------------------------------------------------
% Determine number of sample interval, number of surveys per sample
% interval, and number of observations per sample interval.  Load values
% to extract structure, E.
%--------------------------------------------------------------------------

survey_min = min(S.survey_time);
survey_max = max(S.survey_time);

survey_samples = floor(survey_min):time_interval:ceil(survey_max);
if (survey_samples(end) < survey_max)
  survey_samples = [survey_samples survey_samples(end)+time_interval];
end
Nsamples = length(survey_samples) - 1;

for n=1:Nsamples
  ind = find(S.survey_time >  survey_samples(n)     &                   ...
             S.survey_time <= survey_samples(n+1));

  ncfile = strcat(obs_out,'_',num2str(survey_samples(n)),'.nc');

  E.sample(n).ncfile = ncfile;
  E.sample(n).survey_time = S.survey_time(ind);
  E.sample(n).Nobs = S.Nobs(ind);
  E.sample(n).datum = sum(S.Nobs(ind));
end

out_files = {E.sample.ncfile};

%--------------------------------------------------------------------------
% Report extraction parameters.
%--------------------------------------------------------------------------

if (report)
  for n=1:Nsamples
    disp(' ');
    disp(['Time Interval: ', num2str(n),                                ...
          '   Time Range: ', num2str(survey_samples(n)), ' - '          ...
                             num2str(survey_samples(n+1))]);
    disp(['  Output File: ', E.sample(n).ncfile]); 
    disp([' Survey Times: ', num2str(E.sample(n).survey_time)]); 
    disp(['         Nobs: ', num2str(E.sample(n).Nobs)]);
    disp(['    New datum: ', num2str(E.sample(n).datum)]);
  end
end

%--------------------------------------------------------------------------
% Create new NetCDF files and write out extracted data.
%--------------------------------------------------------------------------

for n=1:Nsamples

  survey=length(E.sample(n).survey_time);
  datum=E.sample(n).datum;

  if (datum > 0)

    ind=find(S.time >  survey_samples(n)     &                          ...
             S.time <= survey_samples(n+1));

% Edit NetCDF information structure and change the value of the "datum"
% dimension.

    I.Dimensions(strcmp({I.Dimensions.Name},'survey' )).Length=survey;
    I.Dimensions(strcmp({I.Dimensions.Name},'datum' )).Length=datum;

% Create sampled observations NetCDF file.

    mode=netcdf.getConstant('CLOBBER');
    mode=bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

    ncfile=char(E.sample(n).ncfile);
    nc_create(ncfile, mode, I);

% Write out observation variables.

    nvars=length(I.Variables);

    notwritten=[];
  
    for m=1:nvars
      vname=char(I.Variables(m).Name);
    
      switch vname
        case 'spherical'
          if (isfield(S,'spherical'))
            nc_write(ncfile,'spherical',S.spherical);
          else
            disp('   field ''spherical''   not found in input structure');
            notwritten=[notwritten 'spherical '];
          end
        case 'Nobs'
          var=E.sample(n).Nobs;
          nc_write(ncfile,vname,var);
        case 'survey_time'
          var=E.sample(n).survey_time;
          nc_write(ncfile,vname,var);
        case 'obs_variance'
          if (isfield(S,'variance'))
            nc_write(ncfile,vname,S.variance);
          else
            disp('   field ''variance''    not found in input structure');
            notwritten=[notwritten 'obs_variance '];
          end
        case 'obs_type'
          if (isfield(S,'type'))
            var=S.type(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''type''        not found in input structure');
            notwritten=[notwritten 'obs_type '];
          end
        case 'obs_provenance'
          if (isfield(S,'provenance'))
            var=S.provenance(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''provenance''  not found in input structure');
            notwritten=[notwritten 'obs_provenance '];
          end
        case 'obs_label'
          if (isfield(S,'label'))
            var=S.label(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''label''  not found in input structure');
            notwritten=[notwritten 'obs_label '];
          end
        case 'obs_time'
          if (isfield(S,'time'))
            var=S.time(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''time''        not found in input structure');
            notwritten=[notwritten 'obs_time '];
          end
        case 'obs_lon'
          if (isfield(S,'lon'))
            var=S.lon(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''lon''         not found in input structure');
            notwritten=[notwritten 'obs_lon '];
          end
        case 'obs_lat'
          if (isfield(S,'lat'))
            var=S.lat(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''lat''         not found in input structure');
            notwritten=[notwritten 'obs_lat '];
          end
        case 'obs_depth'
          if (isfield(S,'depth'))
            var=S.depth(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''depth''       not found in input structure');
            notwritten=[notwritten 'obs_depth '];
          end
        case 'obs_Xgrid'
          if (isfield(S,'Xgrid'))
            var=S.Xgrid(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''Xgrid''       not found in input structure');
            notwritten=[notwritten 'obs_Xgrid '];
          end
        case 'obs_Ygrid'
          if (isfield(S,'Ygrid'))
            var=S.Ygrid(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''Ygrid''       not found in input structure');
            notwritten=[notwritten 'obs_Ygrid '];
          end
        case 'obs_Zgrid'
          if (isfield(S,'Zgrid'))
            var=S.Zgrid(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''Zgrid''       not found in input structure');
            notwritten=[notwritten 'obs_Zgrid '];
          end
        case 'obs_error'
          if (isfield(S,'error'))
            var=S.error(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''error''       not found in input structure');
            notwritten=[notwritten 'obs_error '];
          end
        case 'obs_value'
          if (isfield(S,'value'))
            var=S.value(ind);
            nc_write(ncfile,vname,var);
          else
            disp('   field ''value''       not found in input structure');
            notwritten=[notwritten 'obs_value '];
          end
      end
    end

    if (~isempty(notwritten))
      disp(' ');
      disp(['   Following variables were not written:  ', notwritten]);
    end

% Add extract attribute.

    nc_attadd(ncfile, 'obs_source', obs_inp);

    history=['Extracted from Matlab script:  ',mfilename,               ...
             '  on ',date_stamp];
    nc_attadd(ncfile, 'history', history);
  else
    disp(' ');
    disp(['WARNING: No observations found for the period: ',            ...
          num2str(survey_samples(n)), ' to ',                           ...
          num2str(survey_samples(n+1)),                                 ...
          '  (', E.sample(n).ncfile, ')']);
  end

end

return
