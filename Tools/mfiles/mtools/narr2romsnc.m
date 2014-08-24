% narr2romsnc_mw.m
%
% This routine :
%   input-  reads two 1-year-long 3-hour-averaged NARR Wind data netcdf files
%           ('ncu' and 'ncv' for eastward and northward winds, respectively)
%           from ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/
%
%   user-   selects time intercal [narr_wind_start:narr_wind_end] and
%                   spatial intervel [roms grid or generic grid]
%
%   output- uses native matlab netcdf to create a wind forcing file for ROMS
%           and creates ascii wind forcing files for SWAN
%
% 2012 - 04: I. Safak converted from the non-native netcdf setup of N.Kumar
% 23May2012 - jcwarner, modify slightly to be user defined.
%

%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) enter path and names of NARR wind netcdf files:
ncu = netcdf.open('uwnd.10m.2012.nc','NC_NOWRITE');
ncv = netcdf.open('vwnd.10m.2012.nc','NC_NOWRITE');

%(2) Enter name of output ROMS netcdf wind forcing file
ROMS_NARR_name='roms_narr_Oct2012.nc';
%ROMS_NARR_name='roms_narr_ref3_Oct2012.nc';

%(3) Enter name of output SWAN ASCII wind forcing file
SWAN_NARR_name='swan_narr_Oct2012.dat';
%SWAN_NARR_name='swan_narr_ref3_Oct2012.dat';

%(4) Enter start and end dates
wind_start = datenum('27-Oct-2012');
wind_end   = datenum('31-Oct-2012');

%(5) Need to interpolate winds to a roms grid or a user defined grid.
% Set one of these to a 1, the other to a 0.
interpto_roms_grid = 1;
interpto_user_grid = 0;
if (interpto_roms_grid)
  model_grid='Sandy_roms_grid.nc';
% model_grid='Sandy_roms_grid_ref3.nc';
else
  lon_rho=[255:0.25:310]-360;
  lat_rho=[ 10:0.25:50 ];  % Create a 1/4 degree lat-lon grid
end

%%%%%%%%%%%%%%%%%%%%%   END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

varid = netcdf.inqVarID(ncu,'time'); 
NARR_time = netcdf.getVar(ncu,varid);
NARR_time = NARR_time./24+julian(1800,1,1,0);
NARR_time = datenum(gregorian(NARR_time));

%compute start and end wind indices
narr_wind_start = floor(interp1(NARR_time,[1:1:length(NARR_time)],wind_start));
narr_wind_end   = floor(interp1(NARR_time,[1:1:length(NARR_time)],wind_end));
ntimes=narr_wind_end-narr_wind_start+1;

%get some data from the NARR files
disp('going to get narr data');

varid = netcdf.inqVarID(ncu,'lon'); 
nlon = netcdf.getVar(ncu,varid);
%nlon2 = nlon;nlon2(nlon<0) = nlon2(nlon<0)+360;  % Make all longitudes have positive values
%nlon=nlon2;
numx=size(nlon,1);
numy=size(nlon,2);

varid = netcdf.inqVarID(ncu,'lat'); 
nlat = netcdf.getVar(ncu,varid);

varid = netcdf.inqVarID(ncu,'uwnd');
uwnd = netcdf.getVar(ncu,varid,[0 0 narr_wind_start-1],[numx numy ntimes]);
uwnd_factor= netcdf.getAtt(ncu,varid,'scale_factor');

varid = netcdf.inqVarID(ncv,'vwnd');
vwnd = netcdf.getVar(ncv,varid,[0 0 narr_wind_start-1],[numx numy ntimes]);
vwnd_factor= netcdf.getAtt(ncv,varid,'scale_factor');

%open the roms grid, if needed
if (interpto_roms_grid)
  nc_model_grid= netcdf.open(model_grid);
  varid = netcdf.inqVarID(nc_model_grid,'lon_rho'); 
  lon_rho = netcdf.getVar(nc_model_grid,varid);
  varid = netcdf.inqVarID(nc_model_grid,'lat_rho'); 
  lat_rho = netcdf.getVar(nc_model_grid,varid);
else
  lon_rho=repmat(lon_rho',1,length(lat_rho));
  lat_rho=repmat(lat_rho,size(lon_rho,1),1);
end

% Pre-allocate variables that will contain wind data in the ROMS grid
[Lp,Mp]=size(lon_rho);
L=Lp-1;
M=Mp-1;
uroms=zeros(Lp,Mp,ntimes);
vroms=zeros(Lp,Mp,ntimes);

netcdf.close(ncu);
netcdf.close(ncv);

count = 0;
% Interpolate wind information from NARR onto the ROMS grid
for i = 1:ntimes
    count=count+1;
    disp(['Winds for ROMS at ',datestr(NARR_time(narr_wind_start-1+i))])
    un = double(uwnd(:,:,i)).*uwnd_factor;
        
    A=un>999;
    un(A)=0;
    
    F  = TriScatteredInterp(double(nlon(:)),double(nlat(:)),double(un(:)));
    uroms(:,:,count) = F(double(lon_rho),double(lat_rho));
    clear F A
    
    vn = double(vwnd(:,:,i)).*vwnd_factor;
    
    B=vn>999;
    vn(B)=0;
    
    F  = TriScatteredInterp(double(nlon(:)),double(nlat(:)),double(vn(:)));
    vroms(:,:,count) = F(double(lon_rho),double(lat_rho));
    
    clear F B
end

Time = NARR_time(narr_wind_start:narr_wind_end)-datenum(1858,11,17,0,0,0);

% Creation of NetCDF file for ROMS Winds

ncid = netcdf.create(ROMS_NARR_name,'nc_clobber');

% Global variables
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'type','Gridpak file');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'gridid','combined grid');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',['Created by "' mfilename '" on ' datestr(now)]);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'CPP_options','DCOMPLEX, DBLEPREC, NCARG_32, PLOTS,');

% Dimensions
dimid_time=netcdf.defDim(ncid,'wind_time',length(Time));
dimid_xr=netcdf.defDim(ncid,'xr',Lp);
dimid_er=netcdf.defDim(ncid,'er',Mp);

% Variables
varid = netcdf.defVar(ncid,'wind_time','double',dimid_time);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,Time);
netcdf.reDef(ncid);
netcdf.putAtt(ncid,varid,'long_name','wind time');
netcdf.putAtt(ncid,varid,'units','days since 1858-11-17 00:00:00 UTC');
netcdf.putAtt(ncid,varid,'field','scalar');

varid = netcdf.defVar(ncid,'lat','double',[dimid_xr dimid_er]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,lat_rho);
netcdf.reDef(ncid);
netcdf.putAtt(ncid,varid,'long_name','latitude of rho-points');
netcdf.putAtt(ncid,varid,'units','degree_north');

varid = netcdf.defVar(ncid,'lon','double',[dimid_xr dimid_er]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,lon_rho);
netcdf.reDef(ncid);                                                  
netcdf.putAtt(ncid,varid,'long_name','longitude of rho-points');                   
netcdf.putAtt(ncid,varid,'units','degree_east'); 

varid = netcdf.defVar(ncid,'Uwind','double',[dimid_xr dimid_er dimid_time]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,uroms);
netcdf.reDef(ncid);
netcdf.putAtt(ncid,varid,'long_name','Uwind');
netcdf.putAtt(ncid,varid,'units','meter second-1');
netcdf.putAtt(ncid,varid,'coordinates','lon lat');

varid = netcdf.defVar(ncid,'Vwind','double',[dimid_xr dimid_er dimid_time]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,vroms);
netcdf.reDef(ncid);
netcdf.putAtt(ncid,varid,'long_name','Vwind');
netcdf.putAtt(ncid,varid,'units','meter second-1');
netcdf.putAtt(ncid,varid,'coordinates','lon lat');

netcdf.close(ncid);
disp(['------------ wrote ',ROMS_NARR_name,' ------------']);

% Creation of ASCII files for SWAN Winds
fid = fopen(SWAN_NARR_name,'w');
for i=1:length(Time)
    disp(['Winds for SWAN at ',datestr(Time(i)+datenum(1858,11,17,0,0,0))])
    uswan=squeeze(uroms(:,:,i)');
    vswan=squeeze(vroms(:,:,i)');
    fprintf(fid,'%10.2f\n',uswan');
    fprintf(fid,'%10.2f\n',vswan');
end
fclose(fid);


