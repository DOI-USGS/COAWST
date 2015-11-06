function create_roms_netcdf_clm_mwUL(fn,gn,t_clim)

[xi eta]=size(gn.lon_rho);
%Write NetCDF file using netcdf builtins for 2010a
nc=netcdf.create(fn,'clobber');
if isempty(nc), return, end

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by updatclim on ' datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'climate forcing file from http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis');

% Dimensions:
disp(' ## Defining Dimensions...')
%dimid = netcdf.defDim(ncid,dimname,dimlen)
LP=xi;
MP=eta;
L=LP-1;
M=MP-1;
N=0;

psidimID = netcdf.defDim(nc,'xpsi',L);
xrhodimID = netcdf.defDim(nc,'xrho',LP);
xudimID = netcdf.defDim(nc,'xu',L);
xvdimID = netcdf.defDim(nc,'xv',LP);

epsidimID = netcdf.defDim(nc,'epsi',M);
erhodimID = netcdf.defDim(nc,'erho',MP);
eudimID = netcdf.defDim(nc,'eu',MP);
evdimID = netcdf.defDim(nc,'ev',M);
s_rhodimID = netcdf.defDim(nc,'s_rho',gn.N);

%octdimID = netcdf.defDim(nc,'ocean_time',netcdf.getConstant('NC_UNLIMITED'));
octdimID = netcdf.defDim(nc,'ocean_time',t_clim);
zttdimID = netcdf.defDim(nc,'zeta_time',t_clim);
v2tdimID = netcdf.defDim(nc,'v2d_time',t_clim);
v3tdimID = netcdf.defDim(nc,'v3d_time',t_clim);
sltdimID = netcdf.defDim(nc,'salt_time',t_clim);
tptdimID = netcdf.defDim(nc,'temp_time',t_clim);
onedimID = netcdf.defDim(nc,'one',1);

% Variables and attributes:
disp(' ## Defining Variables, and Attributes...')
%varid = netcdf.defVar(ncid,varname,xtype,dimids)
%netcdf.putAtt(ncid,varid,attrname,attrvalue)

ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
netcdf.putAtt(nc,ocID,'long_name','wind field time');
netcdf.putAtt(nc,ocID,'units','days');
netcdf.putAtt(nc,ocID,'field','wave_time, scalar, series');

ztID = netcdf.defVar(nc,'zeta_time','double',zttdimID);
netcdf.putAtt(nc,ztID,'long_name','zeta_time');
netcdf.putAtt(nc,ztID,'units','days');
netcdf.putAtt(nc,ztID,'field','zeta_time, scalar, series');

v2ID = netcdf.defVar(nc,'v2d_time','double',v2tdimID);
netcdf.putAtt(nc,v2ID,'long_name','v2d_time');
netcdf.putAtt(nc,v2ID,'units','days');
netcdf.putAtt(nc,v2ID,'field','v2d_time, scalar, series');

v3ID = netcdf.defVar(nc,'v3d_time','double',v3tdimID);
netcdf.putAtt(nc,v3ID,'long_name','v3d_time');
netcdf.putAtt(nc,v3ID,'units','days');
netcdf.putAtt(nc,v3ID,'field','v3d_time, scalar, series');

slID = netcdf.defVar(nc,'salt_time','double',sltdimID);
netcdf.putAtt(nc,slID,'long_name','salt_time');
netcdf.putAtt(nc,slID,'units','days');
netcdf.putAtt(nc,slID,'field','salt_time, scalar, series');

tpID = netcdf.defVar(nc,'temp_time','double',tptdimID);
netcdf.putAtt(nc,tpID,'long_name','temp_time');
netcdf.putAtt(nc,tpID,'units','days');
netcdf.putAtt(nc,tpID,'field','temp_time, scalar, series');

lonID = netcdf.defVar(nc,'lon_rho','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,lonID,'long_name','lon_rho');
netcdf.putAtt(nc,lonID,'units','degrees');
netcdf.putAtt(nc,lonID,'FillValue_',100000.);
netcdf.putAtt(nc,lonID,'missing_value',100000.);
netcdf.putAtt(nc,lonID,'field','xp, scalar, series');

latID = netcdf.defVar(nc,'lat_rho','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,latID,'long_name','lon_rho');
netcdf.putAtt(nc,latID,'units','degrees');
netcdf.putAtt(nc,latID,'FillValue_',100000.);
netcdf.putAtt(nc,latID,'missing_value',100000.);
netcdf.putAtt(nc,latID,'field','yp, scalar, series');

zetID = netcdf.defVar(nc,'zeta','double',[xrhodimID erhodimID zttdimID]);
netcdf.putAtt(nc,zetID,'long_name','zeta');
netcdf.putAtt(nc,zetID,'units','meter');
netcdf.putAtt(nc,zetID,'field','zeta, scalar, series');

salID = netcdf.defVar(nc,'salt','float',[xrhodimID erhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc,salID,'long_name','salt');
netcdf.putAtt(nc,salID,'units','psu');
netcdf.putAtt(nc,salID,'field','salt, scalar, series');

tmpID = netcdf.defVar(nc,'temp','float',[xrhodimID erhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc,tmpID,'long_name','temp');
netcdf.putAtt(nc,tmpID,'units','C');
netcdf.putAtt(nc,tmpID,'field','temp, scalar, series');

uID = netcdf.defVar(nc,'u','float',[xudimID eudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc,uID,'long_name','velx');
netcdf.putAtt(nc,uID,'units','meter second-1');
netcdf.putAtt(nc,uID,'field','velx, scalar, series');

vID = netcdf.defVar(nc,'v','float',[xvdimID evdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc,vID,'long_name','vely');
netcdf.putAtt(nc,vID,'units','meter second-1');
netcdf.putAtt(nc,vID,'field','vely, scalar, series');

ubID = netcdf.defVar(nc,'ubar','float',[xudimID eudimID v2tdimID]);
netcdf.putAtt(nc,ubID,'long_name','mean velx');
netcdf.putAtt(nc,ubID,'units','meter second-1');
netcdf.putAtt(nc,ubID,'field','mean velx, scalar, series');

vbID = netcdf.defVar(nc,'vbar','float',[xvdimID evdimID v2tdimID]);
netcdf.putAtt(nc,vbID,'long_name','mean vely');
netcdf.putAtt(nc,vbID,'units','meter second-1');
netcdf.putAtt(nc,vbID,'field','mean vely, scalar, series');

netcdf.close(nc)

