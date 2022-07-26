% get_coawst_data
%
% to get data from the COAWST archive
% US East coast + Gulf of Mexico, 200
% http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd
% use to create init and clm files for another grid
%
% this will rotate and interp in z as needed.
%
% jcw 05Jul2022
%

% Step 1: Enter start and end times
  tstr=datenum(2014,02,01);
  tend=datenum(2014,05,01);

% Step 2: Enter name of grid that you want init and clm files for
  grd_name='E:\data3\Projects\NOPP_hurricanes_modeling\grids\NYBight_grd5.nc';

% Step 3: Enter names of new init and clm files.
  clm_file='NYBight_grd5_clm.nc';
  ini_file='NYBight_grd5_ini.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   end of user input     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% path for coawst thredds
%
  url='http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd'
%
% get some grid info
  lonr=ncread(url,'lon_rho');
  latr=ncread(url,'lat_rho');
  lonu=ncread(url,'lon_u');
  latu=ncread(url,'lat_u');
  lonv=ncread(url,'lon_v');
  latv=ncread(url,'lat_v');
%
%get some time info
  coawst_hours=ncread(url,'time');
  time_att=ncreadatt(url,'time','units');
  time_start=datenum(str2num(time_att(13:16)),str2num(time_att(18:19)),str2num(time_att(21:22)),str2num(time_att(24:25)),str2num(time_att(27:28)),str2num(time_att(30:31)));
  coawst_time=time_start+coawst_hours/24;
%
%get index for start time
  zz=(tstr-coawst_time);
  sindxs=find(zz>=0);
  sindxs=sindxs(end);  %take last index that was positive
%
%get index for end time
  zz=(tend-coawst_time);
  sindxe=find(zz>=0);
  sindxe=sindxe(end);  %take last index that was positive
%
%what times are we doing?
  numsteps=sindxe-sindxs+1;
  coawst_str=datestr(coawst_time(sindxs),'yyyymmdd.HHMMSS');
  coawst_end=datestr(coawst_time(sindxe),'yyyymmdd.HHMMSS');
  disp(['getting data from ', coawst_str,' to ', coawst_end])
%
%create ini file




