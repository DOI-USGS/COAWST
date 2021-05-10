% create_delilah_surface_forc
%
% 21 Dec 2020
%

cd E:\data\models\InWave\readswan\Projects\Delilah\coawst
%
%  load InWave grid
%
netcdf_load('InWave_delilah_grd2.nc');
%
%  load the observed water levels
%
fid=fopen('FRF_Oct1990.tid');
line=fgetl(fid);
for mm=1:100000
  line=fgetl(fid);
  if (line==-1); break; end
  data(mm,:)=str2num(line);
end
fclose(fid);
%
%  here is the header
%  1    2     3    4    5      6      7        8
% Year Month Day Time DecDay Time_m Tide_p  TideDiff
%
%   Tide_m  = Measured tide in meters  
%   Tide_p  = Predicted tide in meters  
%   TideDif = Measured - Predicted tide in meters  
time=data(:,5)-data(1,5);
zeta=data(:,6);
fn='bndry_delilah.nc';

%
% now create a file for zeta and fill it.
%

% get some grid info
[LP,MP]=size(x_rho);
L=LP-1;
Lm=L-1;
M=MP-1;
Mm=M-1;
L  = Lm+1;
M  = Mm+1;
xpsi  = L;
xrho  = LP;
xu    = L;
xv    = LP;
epsi = M;
erho = MP;
eu   = MP;
ev   = M;
s   = 8;   %need to get this from the user

%% create bndry file
nc_bndry=netcdf.create(fn,'clobber');
if isempty(nc_bndry), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by updatclim on ' datestr(now)]);
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'type', 'climate forcing file from http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis');
%% Dimensions:

disp(' ## Defining Dimensions...')
 
psidimID = netcdf.defDim(nc_bndry,'xpsi',L);
xrhodimID = netcdf.defDim(nc_bndry,'xrho',LP);
xudimID = netcdf.defDim(nc_bndry,'xu',L);
xvdimID = netcdf.defDim(nc_bndry,'xv',LP);

epsidimID = netcdf.defDim(nc_bndry,'epsi',M);
erhodimID = netcdf.defDim(nc_bndry,'erho',MP);
eudimID = netcdf.defDim(nc_bndry,'eu',MP);
evdimID = netcdf.defDim(nc_bndry,'ev',M);
s_rhodimID = netcdf.defDim(nc_bndry,'s_rho',s);

zttdimID = netcdf.defDim(nc_bndry,'zeta_time',length(time));
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
ztID = netcdf.defVar(nc_bndry,'zeta_time','double',zttdimID);
netcdf.putAtt(nc_bndry,ztID,'long_name','zeta_time');
%netcdf.putAtt(nc_bndry,ztID,'units','days');
netcdf.putAtt(nc_bndry,ztID,'units','days since 1990-10-01 00:00:00 UTC');
netcdf.putAtt(nc_bndry,ztID,'field','zeta_time, scalar, series');

zeteID = netcdf.defVar(nc_bndry,'zeta_east','double',[erhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zeteID,'long_name','free-surface eastern boundary condition');
netcdf.putAtt(nc_bndry,zeteID,'units','meter');
netcdf.putAtt(nc_bndry,zeteID,'field','zeta_east, scalar, series');

zetnID = netcdf.defVar(nc_bndry,'zeta_north','double',[xrhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zetnID,'long_name','free-surface northern boundary condition');
netcdf.putAtt(nc_bndry,zetnID,'units','meter');
netcdf.putAtt(nc_bndry,zetnID,'field','zeta_north, scalar, series');

zetsID = netcdf.defVar(nc_bndry,'zeta_south','double',[xrhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zetsID,'long_name','free-surface southern boundary condition');
netcdf.putAtt(nc_bndry,zetsID,'units','meter');
netcdf.putAtt(nc_bndry,zetsID,'field','zeta_south, scalar, series');

%close file
netcdf.close(nc_bndry)

%%% now fill it with data

nc_bndry=netcdf.open(fn,'NC_WRITE');

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

%% Time
zettimeid = netcdf.inqVarID(nc_bndry,'zeta_time');
netcdf.putVar(nc_bndry,zettimeid,time);

%% zeta
zetbdyide = netcdf.inqVarID(nc_bndry,'zeta_east');%get id
zeta_east=repmat(zeta,1,MP)';
netcdf.putVar(nc_bndry,zetbdyide,zeta_east);%set variable
clear zeta_east

zetbdyidn = netcdf.inqVarID(nc_bndry,'zeta_north');%get id
zeta_north=repmat(zeta,1,LP)';
netcdf.putVar(nc_bndry,zetbdyidn,zeta_north);%set variable
clear zeta_north

zetbdyids = netcdf.inqVarID(nc_bndry,'zeta_south');%get id
zeta_south=repmat(zeta,1,LP)';
netcdf.putVar(nc_bndry,zetbdyids,zeta_south);%set variable
clear zeta_south

netcdf.close(nc_bndry);

