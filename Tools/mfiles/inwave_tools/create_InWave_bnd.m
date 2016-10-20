function create_inwave_bnd(LP, MP, Nangle_bnd, Dir_bnd, obc, ...
                           AcN, AcE, AcS, AcW, TA, time, bndfile)

disp(' ')
disp(['## Creating the file : ',bndfile])

L=LP-1;
M=MP-1;

%  Create the boundary file
%

type = 'BOUNDARY file for the InWave model' ; 
nc_bndry=netcdf.create(bndfile,'clobber');
if isempty(nc_bndry), return, end

nbin=[1:1:Nangle_bnd];
dir=Dir_bnd;

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by create_InWave_ini on ' datestr(now)]);

%
%  Create dimensions
%

psidimID = netcdf.defDim(nc_bndry,'xpsi',L);
xrhodimID = netcdf.defDim(nc_bndry,'xrho',LP);
xudimID = netcdf.defDim(nc_bndry,'xu',L);
xvdimID = netcdf.defDim(nc_bndry,'xv',LP);

epsidimID = netcdf.defDim(nc_bndry,'epsi',M);
erhodimID = netcdf.defDim(nc_bndry,'erho',MP);
eudimID = netcdf.defDim(nc_bndry,'eu',MP);
evdimID = netcdf.defDim(nc_bndry,'ev',M);

etime_dimID = netcdf.defDim(nc_bndry,'energy_time',length(time));
%eangle_dimID = netcdf.defDim(nc_bndry,'energy_angle',Nbins+1);
eanglec_dimID = netcdf.defDim(nc_bndry,'energy_angle_c',length(nbin));
TA_dimID = netcdf.defDim(nc_bndry,'TA_dim',1);

NT=length(time);

%
%  Create variables and attributes
%
%% Variables and attributes:
disp(' ## Defining Variables and Attributes...')

etID = netcdf.defVar(nc_bndry,'energy_time','double',etime_dimID);
netcdf.putAtt(nc_bndry,etID,'long_name','time for energy envelope');
netcdf.putAtt(nc_bndry,etID,'units','seconds');
netcdf.putAtt(nc_bndry,etID,'field','energy_time, scalar, series');

ecID = netcdf.defVar(nc_bndry,'energy_angle_c','double',eanglec_dimID);
netcdf.putAtt(nc_bndry,ecID,'long_name','direction respect to the north of the bin');
netcdf.putAtt(nc_bndry,ecID,'units','degrees');
netcdf.putAtt(nc_bndry,ecID,'field','energy_angle_c, scalar, series');

TAID = netcdf.defVar(nc_bndry,'TA_dim','double',TA_dimID);
netcdf.putAtt(nc_bndry,TAID,'long_name','representative absolute peak period');
netcdf.putAtt(nc_bndry,TAID,'units','seconds');
netcdf.putAtt(nc_bndry,TAID,'field','TA_dim, scalar, series');

%

if obc(3)==1
%
%   Southern boundary
%
  ACID = netcdf.defVar(nc_bndry,'AC_south','double',[xrhodimID eanglec_dimID etime_dimID]);
  netcdf.putAtt(nc_bndry,ACID,'long_name','southern boundary wave action envelope');
  netcdf.putAtt(nc_bndry,ACID,'units','Joules s m-2 rad-1');
  netcdf.putAtt(nc_bndry,ACID,'field','AC, scalar, series');

  TaID = netcdf.defVar(nc_bndry,'Ta_south','double',TA_dimID);
  netcdf.putAtt(nc_bndry,TaID,'long_name','southern boundary representative absolute peak period');
  netcdf.putAtt(nc_bndry,TaID,'units','seconds');
  netcdf.putAtt(nc_bndry,TaID,'field','Ta, scalar, series');
end

if obc(2)==1
%
%   Eastern boundary
%
  ACID = netcdf.defVar(nc_bndry,'AC_east','double',[erhodimID eanglec_dimID etime_dimID]);
  netcdf.putAtt(nc_bndry,ACID,'long_name','eastern boundary wave action envelope');
  netcdf.putAtt(nc_bndry,ACID,'units','Joules s m-2 rad-1');
  netcdf.putAtt(nc_bndry,ACID,'field','AC, scalar, series');

  TaID = netcdf.defVar(nc_bndry,'Ta_east','double',TA_dimID);
  netcdf.putAtt(nc_bndry,TaID,'long_name','eastern boundary representative absolute peak period');
  netcdf.putAtt(nc_bndry,TaID,'units','seconds');
  netcdf.putAtt(nc_bndry,TaID,'field','Ta, scalar, series');
end

if obc(1)==1
%
%   Northern boundary
%
  ACID = netcdf.defVar(nc_bndry,'AC_north','double',[xrhodimID eanglec_dimID etime_dimID]);
  netcdf.putAtt(nc_bndry,ACID,'long_name','northern boundary wave action envelope');
  netcdf.putAtt(nc_bndry,ACID,'units','Joules s m-2 rad-1');
  netcdf.putAtt(nc_bndry,ACID,'field','AC, scalar, series');

  TaID = netcdf.defVar(nc_bndry,'Ta_north','double',TA_dimID);
  netcdf.putAtt(nc_bndry,TaID,'long_name','northern boundary representative absolute peak period');
  netcdf.putAtt(nc_bndry,TaID,'units','seconds');
  netcdf.putAtt(nc_bndry,TaID,'field','Ta, scalar, series');
end
%
if obc(4)==1
%
%   Western boundary
%
  ACID = netcdf.defVar(nc_bndry,'AC_west','double',[erhodimID eanglec_dimID etime_dimID]);
  netcdf.putAtt(nc_bndry,ACID,'long_name','western boundary wave action envelope');
  netcdf.putAtt(nc_bndry,ACID,'units','Joules s m-2 rad-1');
  netcdf.putAtt(nc_bndry,ACID,'field','AC, scalar, series');

  TaID = netcdf.defVar(nc_bndry,'Ta_west','double',TA_dimID);
  netcdf.putAtt(nc_bndry,TaID,'long_name','western boundary representative absolute peak period');
  netcdf.putAtt(nc_bndry,TaID,'units','seconds');
  netcdf.putAtt(nc_bndry,TaID,'field','Ta, scalar, series');
end

netcdf.close(nc_bndry)

%
% Create global attributes
%

%
% Write variables
%

ncwrite(bndfile,'energy_time',time(1,1:NT));
ncwrite(bndfile,'energy_angle_c',dir);
ncwrite(bndfile,'TA_dim',1);

if obc(3)==1
  ncwrite(bndfile,'AC_south',AcS);
  ncwrite(bndfile,'Ta_south',TA);
end

if obc(2)==1
  ncwrite(bndfile,'AC_east',AcE);
  ncwrite(bndfile,'Ta_east',TA);
end

if obc(1)==1
  ncwrite(bndfile,'AC_north',AcN);
  ncwrite(bndfile,'Ta_north',TA);
end

if obc(4)==1
  ncwrite(bndfile,'AC_west',AcW);
  ncwrite(bndfile,'Ta_west',TA);
end

disp(['Created boundary file -->   ',bndfile])
disp(' ')

return

