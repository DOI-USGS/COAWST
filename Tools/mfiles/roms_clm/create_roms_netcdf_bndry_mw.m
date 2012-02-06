function create_roms_netcdf_bndry_mw(fn,gn,t_clim)
%
% create the boundary file for US East
% so this is special and only has certain arrays.
%
% jcw April 18, 2009
%
%

% Get some grid info. 
  [LP,MP]=size(gn.lon_rho);
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
  s       = gn.N;

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

zttdimID = netcdf.defDim(nc_bndry,'zeta_time',length(t_clim));
v2tdimID = netcdf.defDim(nc_bndry,'v2d_time',length(t_clim));
v3tdimID = netcdf.defDim(nc_bndry,'v3d_time',length(t_clim));
sltdimID = netcdf.defDim(nc_bndry,'salt_time',length(t_clim));
tptdimID = netcdf.defDim(nc_bndry,'temp_time',length(t_clim));
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
ztID = netcdf.defVar(nc_bndry,'zeta_time','double',zttdimID);
netcdf.putAtt(nc_bndry,ztID,'long_name','zeta_time');
netcdf.putAtt(nc_bndry,ztID,'units','days');
netcdf.putAtt(nc_bndry,ztID,'field','zeta_time, scalar, series');

v2ID = netcdf.defVar(nc_bndry,'v2d_time','double',v2tdimID);
netcdf.putAtt(nc_bndry,v2ID,'long_name','v2d_time');
netcdf.putAtt(nc_bndry,v2ID,'units','days');
netcdf.putAtt(nc_bndry,v2ID,'field','v2d_time, scalar, series');

v3ID = netcdf.defVar(nc_bndry,'v3d_time','double',v3tdimID);
netcdf.putAtt(nc_bndry,v3ID,'long_name','v3d_time');
netcdf.putAtt(nc_bndry,v3ID,'units','days');
netcdf.putAtt(nc_bndry,v3ID,'field','v3d_time, scalar, series');

slID = netcdf.defVar(nc_bndry,'salt_time','double',sltdimID);
netcdf.putAtt(nc_bndry,slID,'long_name','salt_time');
netcdf.putAtt(nc_bndry,slID,'units','days');
netcdf.putAtt(nc_bndry,slID,'field','salt_time, scalar, series');

tpID = netcdf.defVar(nc_bndry,'temp_time','double',tptdimID);
netcdf.putAtt(nc_bndry,tpID,'long_name','temp_time');
netcdf.putAtt(nc_bndry,tpID,'units','days');
netcdf.putAtt(nc_bndry,tpID,'field','temp_time, scalar, series');

zetsID = netcdf.defVar(nc_bndry,'zeta_south','double',[xrhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zetsID,'long_name','free-surface southern boundary condition');
netcdf.putAtt(nc_bndry,zetsID,'units','meter');
netcdf.putAtt(nc_bndry,zetsID,'field','zeta_south, scalar, series');

zeteID = netcdf.defVar(nc_bndry,'zeta_east','double',[erhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zeteID,'long_name','free-surface eastern boundary condition');
netcdf.putAtt(nc_bndry,zeteID,'units','meter');
netcdf.putAtt(nc_bndry,zeteID,'field','zeta_east, scalar, series');

ubsID = netcdf.defVar(nc_bndry,'ubar_south','float',[xudimID v2tdimID]);
netcdf.putAtt(nc_bndry,ubsID,'long_name','2D u-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,ubsID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ubsID,'field','ubar_south, scalar, series');

ubeID = netcdf.defVar(nc_bndry,'ubar_east','float',[eudimID v2tdimID]);
netcdf.putAtt(nc_bndry,ubeID,'long_name','2D u-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,ubeID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ubeID,'field','ubar_east, scalar, series');

vbsID = netcdf.defVar(nc_bndry,'vbar_south','float',[xvdimID v2tdimID]);
netcdf.putAtt(nc_bndry,vbsID,'long_name','2D v-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,vbsID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vbsID,'field','vbar_south, scalar, series');

vbeID = netcdf.defVar(nc_bndry,'vbar_east','float',[evdimID v2tdimID]);
netcdf.putAtt(nc_bndry,vbeID,'long_name','2D v-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,vbeID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vbeID,'field','vbar_east, scalar, series');

usID = netcdf.defVar(nc_bndry,'u_south','float',[xudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,usID,'long_name','3D u-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,usID,'units','meter second-1');
netcdf.putAtt(nc_bndry,usID,'field','u_south, scalar, series');

ueID = netcdf.defVar(nc_bndry,'u_east','float',[eudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,ueID,'long_name','3D u-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,ueID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ueID,'field','u_east, scalar, series');

vsID = netcdf.defVar(nc_bndry,'v_south','float',[xvdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,vsID,'long_name','3D v-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,vsID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vsID,'field','v_south, scalar, series');

veID = netcdf.defVar(nc_bndry,'v_east','float',[evdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,veID,'long_name','3D v-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,veID,'units','meter second-1');
netcdf.putAtt(nc_bndry,veID,'field','v_east, scalar, series');

tmpsID = netcdf.defVar(nc_bndry,'temp_south','float',[xrhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc_bndry,tmpsID,'long_name','3D temperature southern boundary condition');
netcdf.putAtt(nc_bndry,tmpsID,'units','C');
netcdf.putAtt(nc_bndry,tmpsID,'field','temp_south, scalar, series');

tmpeID = netcdf.defVar(nc_bndry,'temp_east','float',[erhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc_bndry,tmpeID,'long_name','3D temperature eastern boundary condition');
netcdf.putAtt(nc_bndry,tmpeID,'units','C');
netcdf.putAtt(nc_bndry,tmpeID,'field','temp_east, scalar, series');

salsID = netcdf.defVar(nc_bndry,'salt_south','float',[xrhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc_bndry,salsID,'long_name','3D salinity southern boundary condition');
netcdf.putAtt(nc_bndry,salsID,'units','psu');
netcdf.putAtt(nc_bndry,salsID,'field','salt_south, scalar, series');

saleID = netcdf.defVar(nc_bndry,'salt_east','float',[erhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc_bndry,saleID,'long_name','3D salinity eastern boundary condition');
netcdf.putAtt(nc_bndry,saleID,'units','psu');
netcdf.putAtt(nc_bndry,saleID,'field','salt_east, scalar, series');

%close file
netcdf.close(nc_bndry)

