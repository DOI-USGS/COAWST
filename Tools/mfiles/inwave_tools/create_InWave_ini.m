function create_inwave_ini(LP,MP,Nbins,Bindir,Bindir_c,pd,Ac,Cx,Cy,Ct,TA,inifile)

disp(' ')
disp(['## Creating the file : ',inifile])

initime=0;

%
%  Create the initial file
%

type = 'INITIAL file for the InWave model' ; 
history = 'InWave' ;
nc_ini=netcdf.create(inifile,'clobber');
if isempty(nc_ini), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_ini,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by create_InWave_ini on ' datestr(now)]);

%% Dimensions:
disp(' ## Defining Dimensions...')
L=LP-1;
M=MP-1;

psidimID = netcdf.defDim(nc_ini,'xpsi',L);
xrhodimID = netcdf.defDim(nc_ini,'xrho',LP);
xudimID = netcdf.defDim(nc_ini,'xu',L);
xvdimID = netcdf.defDim(nc_ini,'xv',LP);

epsidimID = netcdf.defDim(nc_ini,'epsi',M);
erhodimID = netcdf.defDim(nc_ini,'erho',MP);
eudimID = netcdf.defDim(nc_ini,'eu',MP);
evdimID = netcdf.defDim(nc_ini,'ev',M);

etime_dimID = netcdf.defDim(nc_ini,'energy_time',length(initime));
eangle_dimID = netcdf.defDim(nc_ini,'energy_angle',Nbins+1);
eanglec_dimID = netcdf.defDim(nc_ini,'energy_angle_c',Nbins);
TA_dimID = netcdf.defDim(nc_ini,'TA_dim',1);


%% Variables and attributes:
disp(' ## Defining Variables and Attributes...')

etID = netcdf.defVar(nc_ini,'energy_time','double',etime_dimID);
netcdf.putAtt(nc_ini,etID,'long_name','time for energy envelope');
netcdf.putAtt(nc_ini,etID,'units','seconds');
netcdf.putAtt(nc_ini,etID,'field','energy_time, scalar, series');

eaID = netcdf.defVar(nc_ini,'energy_angle','double',eanglec_dimID);
netcdf.putAtt(nc_ini,eaID,'long_name','direction respect to the north of the bin');
netcdf.putAtt(nc_ini,eaID,'units','degrees');
netcdf.putAtt(nc_ini,eaID,'field','energy_angle, scalar, series');

ecID = netcdf.defVar(nc_ini,'energy_angle_c','double',eangle_dimID);
netcdf.putAtt(nc_ini,ecID,'long_name','direction respect to the north of the bin');
netcdf.putAtt(nc_ini,ecID,'units','degrees');
netcdf.putAtt(nc_ini,ecID,'field','energy_angle_c, scalar, series');

TAID = netcdf.defVar(nc_ini,'TA_dim','double',TA_dimID);
netcdf.putAtt(nc_ini,TAID,'long_name','representative absolute peak period');
netcdf.putAtt(nc_ini,TAID,'units','seconds');
netcdf.putAtt(nc_ini,TAID,'field','TA_dim, scalar, series');

ACID = netcdf.defVar(nc_ini,'AC','double',[xrhodimID erhodimID eanglec_dimID]);
netcdf.putAtt(nc_ini,ACID,'long_name','wave action envelope');
netcdf.putAtt(nc_ini,ACID,'units','Joules s m-2 rad-1');
netcdf.putAtt(nc_ini,ACID,'field','AC, scalar, series');

cxID = netcdf.defVar(nc_ini,'Cx','double',[xudimID erhodimID eanglec_dimID]);
netcdf.putAtt(nc_ini,cxID,'long_name','x component of the group celerity');
netcdf.putAtt(nc_ini,cxID,'units','m s-1');
netcdf.putAtt(nc_ini,cxID,'field','cx, scalar, series');

cyID = netcdf.defVar(nc_ini,'Cy','double',[xrhodimID evdimID eanglec_dimID]);
netcdf.putAtt(nc_ini,cyID,'long_name','y component of the group celerity');
netcdf.putAtt(nc_ini,cyID,'units','m s-1');
netcdf.putAtt(nc_ini,cyID,'field','cy, scalar, series');

ctID = netcdf.defVar(nc_ini,'Ct','double',[xrhodimID erhodimID eangle_dimID]);
netcdf.putAtt(nc_ini,ctID,'long_name','directional component of the group celerity');
netcdf.putAtt(nc_ini,ctID,'units','rad s-1');
netcdf.putAtt(nc_ini,ctID,'field','ct, scalar, series');

TaID = netcdf.defVar(nc_ini,'Ta','double',TA_dimID);
netcdf.putAtt(nc_ini,TaID,'long_name','representative absolute peak period');
netcdf.putAtt(nc_ini,TaID,'units','seconds');
netcdf.putAtt(nc_ini,TaID,'field','Ta, scalar, series');

netcdf.close(nc_ini)

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




