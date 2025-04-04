%
%  D_STD_frc:  Driver script to compute and write 4D-Var surface forcing
%              standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  surface forcing standard deviations which are used to convert modeled
%  error correlations to error covariances.
%
%  The first step is to run the model application with air/sea bulk
%  fluxes (BULK_FLUXES option) for a period that is long enough to
%  compute meaningful circulation statistics, like mean and standard
%  deviations for surface forcing state variables: sustr, svstr,
%  shflux, and ssflux.
%
%  As an example, we compute the 4D-Var standard deviations for the WC13
%  application from daily history files from 1/1/2000 to 12/25/2004.
%

% svn $Id: d_std_frc.m 1156 2023-02-18 01:44:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
%  User tunable parameters.
%--------------------------------------------------------------------------

%  Set standard deviation NetCDF file. The file name is edited and the
%  month will be appended as *i_jan.nc:

STDfile = 'wc13_std_f.nc';

mstr = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};

%  Set input grif file.

my_root = '/home/arango/ocean/toms/repository/test';

GRDfile = strcat(my_root, '/WC13/Data/wc13_grd.nc');

%  Set input history files (string cell structure).

my_root = '/home/arango/ocean/toms/repository/test';

HISdir  = strcat(my_root, '/WC13/STD/Data');

HISfile = dir(fullfile(HISdir, 'wc*.nc'));

nfiles = length(HISfile);

%  Set surface forcing state variables dynamical fields (cell array)
%  to process.

field_list = {'sustr', 'svstr', 'shflux', 'ssflux'};

%  Set grid variables dynamical fields (cell array) to write in
%  output file(s).

grid_list  = {'theta_s', 'theta_b', 'Tcline' , 'hc'     , 's_rho'  , ...
              's_w'    , 'Cs_r'   , 'Cs_w'   , 'h'      , 'lon_rho', ...
              'lat_rho', 'lon_u'  , 'lat_u'  , 'lon_v'  , 'lat_v'};

%  Initialize working structure.

S.title = 'California Current System, 1/3 degree resolution (WC13)';

S.grd_file  = GRDfile;              % application grid

S.do_sustr  = true;                 % surface u-momentum stress
S.do_svstr  = true;                 % surface v-momentum stress
S.do_shflux = true;                 % surface net heat flux
S.do_ssflux = true;                 % surface salt flux (E-P)*SALT

%--------------------------------------------------------------------------
%  Compute monthly averages and standard deviations.
%--------------------------------------------------------------------------

%  Extract first history file name.

HisFile1 = fullfile(HISdir, HISfile(1).name);

%  Set NetCDF file creation parameters.

V = nc_vnames(HisFile1);
nvars = length(V.Variables);
 
S.curvilinear = false;
S.masking     = false;
S.Vtransform  = 1;
S.Vstretching = 1;

for n=1:nvars,
  name = char(V.Variables(n).Name);
  switch name
    case 'spherical'
      S.spherical = nc_read(HisFile1, 'spherical');
      if (ischar(S.spherical)),
        if (S.spherical == 'T' || S.spherical == 't');
          S.spherical = 1;
        else
          S.spherical = 0;
        end
      end
    case 'Vtransform'
      S.Vtransform  = nc_read(HisFile1, 'Vtransform');
    case 'Vstretching'
      S.Vstretching = nc_read(HisFile1, 'Vstretching');
    case 'angle'
      S.angle = nc_read(HisFile1, 'angle');
      S.curvilinear = true;
    case 'mask_rho'
      S.rmask = nc_read(HisFile1, 'mask_rho');
      S.umask = nc_read(HisFile1, 'mask_u');
      S.vmask = nc_read(HisFile1, 'mask_v');
      S.masking = true;
  end
end

if (S.curvilinear),
  grid_list = [grid_list, 'angle'];
end,
  
if (S.masking),
  grid_list = [grid_list, 'mask_rho', 'mask_u', 'mask_v'];
end

%  Get grid size.

[Lr,Mr,Nr] = size(nc_read(HisFile1, 'temp', 1));

S.Lm = Lr - 2;                      % number of interior RHO x-points
S.Mm = Mr - 2;                      % number of interior RHO y-points
S.N  = Nr;                          % number of vertical RHO levels

%  Read in grid.

S.rlon = nc_read(HisFile1, 'lon_rho');
S.rlat = nc_read(HisFile1, 'lat_rho');

S.ulon = nc_read(HisFile1, 'lon_u');
S.ulat = nc_read(HisFile1, 'lat_u');

S.vlon = nc_read(HisFile1, 'lon_v');
S.vlat = nc_read(HisFile1, 'lat_v');

%  Process month-by-month.

for m=1:12,

%  Set standard deviation file name.

  lstr = length(STDfile)-3;
  S.ncname = char(strcat(STDfile(1:lstr), '_', mstr(m), '.nc'));
  
%  Initialize mean and variance arrays.

  Rcount = 0;                         % record counter
  rec    = 1;                         % initialization record

  S.month = m;

  for fval = field_list,
    field     = char(fval);
    field_avg = [field, '_avg'];
    field_std = [field, '_std'];
  
    try
      S.(field_avg) = zeros(size(nc_read(HisFile1, field, rec))); 
      S.(field_std) = S.(field_avg);
    catch
      disp([' D_STD_FRC: error while processing, rec = ', num2str(rec)]);
      disp(['            for variable : ', field]);
      didp(['            in file: ', HisFile1]);
      return
    end,
  end,

  disp(' ');
  disp([ 'Computing mean and standard deviation fields, month = ',      ...
        num2str(m), ' ...']);
  disp(' ');

%  Accumulate montly fields.

  for n=1:nfiles,

    ncfile = fullfile(HISdir, HISfile(n).name);

    time = nc_read(ncfile,'ocean_time');
    Nrec = length(time);

    for rec=1:Nrec,

      [year,month,day]=datevec(datenum(1968,5,23) + time(rec)/86400);
      
      if (month == m),

        mydate=datestr(datenum(1968,5,23) + time(rec)/86400);
        disp([ '*** Processing Fields: ', mydate]);

        Rcount = Rcount + 1;

        for fval = field_list,
          field     = char(fval);           % convert cell to string
          field_avg = [field, '_avg'];      % average field 
          field_std = [field, '_std'];      % standard deviation field

          try
            F = nc_read(ncfile, field, rec);

          catch
            disp([' D_STD_FRC: error while processing, rec = ',         ...
                  num2str(rec)]);
            disp(['            for variable : ', field]);
            didp(['            in file: ', HisFile1]);
            return
          end

          S.(field_avg) = S.(field_avg) + F;
          S.(field_std) = S.(field_std) + F.^2;
        end
      
      end
    
    end
  
  end

%  Compute monthly mean and standard deviation fields. Use an
%  unbiased estimate for variance:
%
%    var = [(sum(Xi ^2), i=1:N) / (N-1)] - N * Xmean / (N-1)

  fac1 = 1 / max(1,Rcount-1);
  fac2 = fac1 * Rcount;

  for fval = field_list,
    field     = char(fval);                 % convert cell to string
    field_avg = [field, '_avg'];            % average field 
    field_std = [field, '_std'];            % standard deviation field

    S.(field_avg) = S.(field_avg) ./ Rcount;
    S.(field_std) = sqrt(fac1 * S.(field_std) - fac2 * S.(field_avg) .^ 2);
  end

%--------------------------------------------------------------------------
%  Write out surface forcing standard deviation fields.
%--------------------------------------------------------------------------

  disp(' ');
  disp([ 'Writting out standard deviation, file = ', S.ncname]);
  disp(' ');

%  Create surface standard deviations NetCDF file.

  s = c_std_frc(S);

%  Write out grid data.

  v = 'spherical';    s = nc_write(S.ncname, v, S.spherical);
  v = 'Vtransform';   s = nc_write(S.ncname, v, S.Vtransform);
  v = 'Vstretching';  s = nc_write(S.ncname, v, S.Vstretching);

  for fval = grid_list,
    field = char(fval);                     % convert cell to string

    f = nc_read(HisFile1, field);  s = nc_write(S.ncname, field, f);
  end,
    
% Write out surface forcing standard deviation data.

  rec = 1;
  
  s = nc_write(S.ncname, 'ocean_time', 0, rec);

  for fval = field_list,
    field     = char(fval);                 % convert cell to string
    field_std = [field, '_std'];            % standard deviation field

    s = nc_write(S.ncname, field, S.(field_std), rec);
  end,

% Process next month.

end,

disp(' ');
disp('Done.');
disp(' ');
