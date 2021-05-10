% create_delilah_surface_forc
%
% 21 Dec 2020
%

cd E:\data\models\InWave\readswan\Projects\Delilah\coawst
%cd C:\work\Projects\Delilah
%
%  load InWave grid
%
netcdf_load('InWave_delilah_grd2.nc');
%
%  load the observed met data
%
fid=fopen('FRF_Oct1990.met');
line=fgetl(fid);
for mm=1:100000
  line=fgetl(fid);
  if (line==-1); break; end
  data(mm,:)=str2num(line);
end
fclose(fid);
%
%  here is the header
%  1    2    3      4        5        6          7        8       9         10        11 
% Month Day Time   DecDay    M616     M624     M932_s    M933    M932_x    M932_y    M604
%   M616   = Barometric pressure in mbar
%   M624   = Air temperature in degrees C
%   M932_s = Wind speed in m/s
%   M933   = direction the wind comes from, degrees clockwise from true north 
%   M932_x = Cross-shore wind component in m/s, +x is onshore
%   M932_y = Longshore wind component in m/s, +y is southerly
%   M604   = mm of rain fall in the 34-minute period 

%time_start=datenum(1990,10,01);
%time=datenum(1990,data(:,1),data(:,2),floor(data(:,3)/100),data(:,2)*0,data(:,2)*0);
%time=time-time_start;
time=data(:,4)-274.;

Uwind=data(:,7).*cos(-data(:,8)*pi/180+3*pi/2);
Vwind=data(:,7).*sin(-data(:,8)*pi/180+3*pi/2);

su=1.22*0.0013*Uwind.^2;
sv=1.22*0.0013*Vwind.^2;

% try to make a roms forc file with just 1 point of data, since it is the
%  same spatially.  Also, 20 degrees rotation.
angle_r=angle(1,1);
%sustr=ones(size(x_u,1),size(x_u,2),length(time(1:1:end)));
%svstr=ones(size(x_v,1),size(x_v,2),length(time(1:1:end)));
sustr=ones(2,2,length(time(1:1:end)));
svstr=ones(2,2,length(time(1:1:end)));

ct=0;
for mm=1:length(time)
  ct=ct+1;
  sustr(:,:,ct)=su(ct).*cos(angle_r)+sv(ct).*sin(angle_r);
  svstr(:,:,ct)=sv(ct).*cos(angle_r)-su(ct).*sin(angle_r);
end

lon=[x_rho(1,1) x_rho(end,1); x_rho(1,end) x_rho(end,end)];
lat=[y_rho(1,1) y_rho(end,1); y_rho(1,end) y_rho(end,end)];
create_roms_forcings(lon, lat, time, 'frc_delilah2.nc', 'sustr', 'svstr');

%now create swan wind file, if needed
%first interpolate to a constant dt
swan_time=[time(1):1/24:21-1/24];
swan_uwind=interp1(time,Uwind,swan_time);
swan_vwind=interp1(time,Vwind,swan_time);
fid = fopen('wind_swan_delilah.dat','w');
for i=1:length(swan_time)
% disp(['Writing winds for SWAN at ',datestr(Time(i)+datenum(1858,11,17,0,0,0))])
  uswan=squeeze(swan_uwind(i)); uswan=repmat(uswan,1,4);
  vswan=squeeze(swan_vwind(i)); vswan=repmat(vswan,1,4);
  fprintf(fid,'%10.2f\n',uswan');
  fprintf(fid,'%10.2f\n',vswan');
end
fclose(fid);
%then you need to use these wind commands
%&& KEYWORD TO CREATE WIND GRID &&
%INPGRID WIND REGULAR 0 0 0 1 1 1000 1500 &
%       NONSTATIONARY 19901001.010000 1 HR 19901021.230000
%READINP WIND 1 'Projects/Delilah/wind_swan_delilah.dat' 4 0 FREE




