function create_inwave_ini(LP,MP,Nbins,Bindir,Bindir_c,pd,Ac,Cx,Cy,Ct,TA,inifile)

disp(' ')
disp(['## Creating the file : ',inifile])

initime=0;

%
%  Create the initial file
%

type = 'INITIAL file' ; 
history = 'InWave' ;
nc_bndry=netcdf.create(inifile,'clobber');
if isempty(nc_bndry), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by create_InWave_ini on ' datestr(now)]);

%% Dimensions:
disp(' ## Defining Dimensions...')
L=LP-1;
M=MP-1;

psidimID = netcdf.defDim(nc_bndry,'xpsi',L);
xrhodimID = netcdf.defDim(nc_bndry,'xrho',LP);
xudimID = netcdf.defDim(nc_bndry,'xu',L);
xvdimID = netcdf.defDim(nc_bndry,'xv',LP);

epsidimID = netcdf.defDim(nc_bndry,'epsi',M);
erhodimID = netcdf.defDim(nc_bndry,'erho',MP);
eudimID = netcdf.defDim(nc_bndry,'eu',MP);
evdimID = netcdf.defDim(nc_bndry,'ev',M);

etime_dimID = netcdf.defDim(nc_bndry,'energy_time',length(initime));
eangle_dimID = netcdf.defDim(nc_bndry,'energy_angle',Nbins+1);
eanglec_dimID = netcdf.defDim(nc_bndry,'energy_angle_c',Nbins);
TA_dimID = netcdf.defDim(nc_bndry,'TA_dim',1);


%% Variables and attributes:
disp(' ## Defining Variables and Attributes...')

etID = netcdf.defVar(nc_bndry,'energy_time','double',etime_dimID);
netcdf.putAtt(nc_bndry,etID,'long_name','time for energy envelope');
netcdf.putAtt(nc_bndry,etID,'units','seconds');
netcdf.putAtt(nc_bndry,etID,'field','energy_time, scalar, series');

eaID = netcdf.defVar(nc_bndry,'energy_angle','double',eanglec_dimID);
netcdf.putAtt(nc_bndry,eaID,'long_name','direction respect to the north of the bin');
netcdf.putAtt(nc_bndry,eaID,'units','degrees');
netcdf.putAtt(nc_bndry,eaID,'field','energy_angle, scalar, series');

ecID = netcdf.defVar(nc_bndry,'energy_angle_c','double',eangle_dimID);
netcdf.putAtt(nc_bndry,ecID,'long_name','direction respect to the north of the bin');
netcdf.putAtt(nc_bndry,ecID,'units','degrees');
netcdf.putAtt(nc_bndry,ecID,'field','energy_angle_c, scalar, series');

TAID = netcdf.defVar(nc_bndry,'TA_dim','double',TA_dimID);
netcdf.putAtt(nc_bndry,TAID,'long_name','representative absolute peak period');
netcdf.putAtt(nc_bndry,TAID,'units','seconds');
netcdf.putAtt(nc_bndry,TAID,'field','TA_dim, scalar, series');

ACID = netcdf.defVar(nc_bndry,'AC','double',[xrhodimID erhodimID eanglec_dimID]);
netcdf.putAtt(nc_bndry,ACID,'long_name','wave energy envelope');
netcdf.putAtt(nc_bndry,ACID,'units','Joules');
netcdf.putAtt(nc_bndry,ACID,'field','AC, scalar, series');

cxID = netcdf.defVar(nc_bndry,'Cx','double',[xudimID erhodimID eanglec_dimID]);
netcdf.putAtt(nc_bndry,cxID,'long_name','x component of the group celerity');
netcdf.putAtt(nc_bndry,cxID,'units','m s-1');
netcdf.putAtt(nc_bndry,cxID,'field','cx, scalar, series');

cyID = netcdf.defVar(nc_bndry,'Cy','double',[xrhodimID evdimID eanglec_dimID]);
netcdf.putAtt(nc_bndry,cyID,'long_name','y component of the group celerity');
netcdf.putAtt(nc_bndry,cyID,'units','m s-1');
netcdf.putAtt(nc_bndry,cyID,'field','cy, scalar, series');

ctID = netcdf.defVar(nc_bndry,'Ct','double',[xrhodimID erhodimID eangle_dimID]);
netcdf.putAtt(nc_bndry,ctID,'long_name','directional component of the group celerity');
netcdf.putAtt(nc_bndry,ctID,'units','rad s-1');
netcdf.putAtt(nc_bndry,ctID,'field','ct, scalar, series');

TaID = netcdf.defVar(nc_bndry,'Ta','double',TA_dimID);
netcdf.putAtt(nc_bndry,TaID,'long_name','representative absolute peak period');
netcdf.putAtt(nc_bndry,TaID,'units','seconds');
netcdf.putAtt(nc_bndry,TaID,'field','Ta, scalar, series');

netcdf.close(nc_bndry)

%
% Write variables
%
ncwrite(inifile,'energy_time',initime); 
ncwrite(inifile,'energy_angle',Bindir_c); 
ncwrite(inifile,'energy_angle_c',Bindir); 
ncwrite(inifile,'TA_dim',1); 

ncwrite(inifile,'AC',Ac);
ncwrite(inifile,'Cx',Cx); 
ncwrite(inifile,'Cy',Cy); 
ncwrite(inifile,'Ct',Ct); 
ncwrite(inifile,'Ta',TA);

disp(['Created initial file -->   ',inifile])

end




