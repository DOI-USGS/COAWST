delete('D:\Temporary\COAWST\SWAN\filesexist')

% grid to model on
gridname='USeast_grd3';
modelgrid='D:\Temporary\COAWST\SWAN\forcing\USeast_grd3.nc';%used to get info, size, lat, lon, etc.

%time constraints
timestart=str2num(datestr(now,'yyyymmdd'));
timeend2=str2num(datestr(now+3,'yyyymmdd'));

% load ww3 data using OpenDap
 loaddap_ww3_2TPAR(modelgrid)

% [lat,lon,time,dp,hs,tp,info,indicator]=loaddap_ww3;
% %create tpar files in D:\Temporary\COAWST\SWAN\forcing
% if indicator==1
% [timestr]=timefix(time);
% ww3totparnew(dp.dirpwsfc,hs.htsgwsfc,tp.perpwsfc,timestr,lat.lat,lon.lon,modelgrid,gridname,time);
% end

cd D:\Temporary\COAWST\SWAN\forcing

% load wind data and create wind files in  D:\Temporary\COAWST\SWAN\forcing
[windstart_gfs,windend_gfs,indicator_gfs]=loaddap_gfs(modelgrid);
[windstart,windend,indicator2]=loaddap_wind(modelgrid);
mash_dat(modelgrid,windstart,windstart_gfs);

%info needed for INPUT files
load D:\Temporary\COAWST\SWAN\forcing\specpts
[r c]=size(specpts);

%create INPUT files
SWANInput(timestart,timeend2,gridname,r,specpts,windstart,windend,modelgrid);%uses hot file
SWANInstat(timestart,timeend2,gridname,r,specpts,windstart,windend,modelgrid);%creates hot file
ncclose

%Make sure filesexist is removed
fname = tempname; fid = fopen(fname, 'wt'); fprintf(fid, '%s', 'cd /raid1/barmstrong/Projects/COAWST/swan/test; rm filesexist; exit');
fclose(fid); system(['ssh -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu < ' fname]);

%copies swan input and tpar files over
cd D:/Temporary/COAWST/SWAN/run
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" INPUT barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/nonstat
cd ../runstat
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" INPUT barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/stat
cd ../forcing
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" *.txt barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/forcings
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" comb_wind_update.dat barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/forcings

cd D:\Temporary\COAWST
% start batch file on nemo
fname = tempname; fid = fopen(fname, 'wt'); fprintf(fid, '%s', 'cd /raid1/barmstrong/Projects/COAWST/swan/test; qsub run_swan; exit');
fclose(fid); system(['ssh -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu < ' fname]);

% log run information
fid=fopen('D:\Temporary\COAWST\coawstlog.txt','a');
fprintf(fid,datestr(now-1));
fprintf(fid,' start ');
fprintf(fid,' %8.6f ',timestart);
fprintf(fid,' finish ');
fprintf(fid,'%8.6f ',timeend2);
fclose(fid);

% fname = tempname; fid = fopen(fname, 'wt'); fprintf(fid, '%s', 'cd /raid1/barmstrong/Projects/COAWST/swan/test; qsub run_swan1; exit');
% fclose(fid); system(['ssh -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu < ' fname]);


% !ssh -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu -c 'cd /raid1/barmstrong/Projects/COAWST/swan/test; qsub run_swan'
% t = timer('TimerFcn','disp(''time to wait...'')','StopFcn','output_phase','StartDelay',3600,'TasksToExecute',2);
% start(t);
output_phase;
%exit