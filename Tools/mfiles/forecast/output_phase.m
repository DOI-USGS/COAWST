%Set up a timer 
t2=timer('TimerFcn','disp(''time to wait...'')','StartDelay',(3600/4),...
    'StopFcn','!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/filesexist .',...
    'TasksToExecute',1); 
%check for special output file
cd D:/Temporary/COAWST/SWAN/
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/filesexist .
X=exist('filesexist');
while X==0 
   start(t2)
   wait(t2)
   !scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/filesexist .
   X=exist('filesexist')
end

cd run
%run .sh script to ssh files back to local machine
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" barmstrong@nemo.whoi.edu:/raid1/barmstrong/Projects/COAWST/swan/test/nonstat/*.mat .

%converts SWAN output files to netcdf
cd 'D:\Temporary\COAWST\SWAN\run'
%uses INPUT from second, nonstat run
swanmat2roms_bwave('INPUT','hsig.nc',0)
swanmat2roms_bwind('INPUT','wind.nc',0)
%converts netcdf files to png movie files
googlemap_movie_swan

%run .sh script to ssh png files to capecodder 
cd 'D:\Temporary\COAWST\SWAN\run\pngfiles'
!scp -i "C:/Documents and Settings/barmstrong/.ssh/COAWSTid_rsa" *.gif barmstrong@capecodder.er.usgs.gov:/Volumes/web/user/cccp/public/COAWST

%write to log file
fid=fopen('D:\Temporary\COAWST\coawstlog.txt','a');
fprintf(fid,' finished ');
fprintf(fid,datestr(now-1));
fprintf(fid,' \n');
fclose(fid);

exit

