function roms_tide_forc_file(Fname,Lp,Mp,periods)

%create init file
nc_init=netcdf.create(Fname,'clobber');
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by create_roms_tide_forc on ' datestr(now)]);
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'type', 'tide forcing file');


%% Dimensions:

disp(' ## Defining Dimensions...')
 
%get some grid info
xrho = Lp;
erho = Mp;
ntide = length(periods);

xrhodimID = netcdf.defDim(nc_init,'xrho',Lp);
erhodimID = netcdf.defDim(nc_init,'erho',Mp);
ntidedimID = netcdf.defDim(nc_init,'tide_period',ntide);

periodsID = netcdf.defVar(nc_init,'tide_period','double',[ntidedimID]);
netcdf.putAtt(nc_init,periodsID,'long_name','tide angular period');
netcdf.putAtt(nc_init,periodsID,'units','hours');
netcdf.putAtt(nc_init,periodsID,'field','tide_period, scalar, series');

ephaseID = netcdf.defVar(nc_init,'tide_Ephase','double',[xrhodimID erhodimID ntidedimID]);
netcdf.putAtt(nc_init,ephaseID,'long_name','tidal elevation phase angle');
netcdf.putAtt(nc_init,ephaseID,'units','degrees, time of maximum elevation with respect chosen time origin');
netcdf.putAtt(nc_init,ephaseID,'field','tide_Ephase, scalar, series');

eampID = netcdf.defVar(nc_init,'tide_Eamp','double',[xrhodimID erhodimID ntidedimID]);
netcdf.putAtt(nc_init,eampID,'long_name','tidal elevation amplitude');
netcdf.putAtt(nc_init,eampID,'units','meter');
netcdf.putAtt(nc_init,eampID,'field','tide_Eamp, scalar, series');

cphaseID = netcdf.defVar(nc_init,'tide_Cphase','double',[xrhodimID erhodimID ntidedimID]);
netcdf.putAtt(nc_init,cphaseID,'long_name','tidal current phase angle');
netcdf.putAtt(nc_init,cphaseID,'units','degrees, time of maximum velocity with respect chosen time origin');
netcdf.putAtt(nc_init,cphaseID,'field','tide_Cphase, scalar, series');

cangleID = netcdf.defVar(nc_init,'tide_Cangle','double',[xrhodimID erhodimID ntidedimID]);
netcdf.putAtt(nc_init,cangleID,'long_name','tidal current inclination angle');
netcdf.putAtt(nc_init,cangleID,'units','degrees between semi-major axis and East');
netcdf.putAtt(nc_init,cangleID,'field','tide_Cangle, scalar, series');

cminID = netcdf.defVar(nc_init,'tide_Cmin','double',[xrhodimID erhodimID ntidedimID]);
netcdf.putAtt(nc_init,cminID,'long_name','minimum tidal current, ellipse semi-minor axis');
netcdf.putAtt(nc_init,cminID,'units','meter second-1');
netcdf.putAtt(nc_init,cminID,'field','tide_Cmin, scalar, series');

cmaxID = netcdf.defVar(nc_init,'tide_Cmax','double',[xrhodimID erhodimID ntidedimID]);
netcdf.putAtt(nc_init,cmaxID,'long_name','maximum tidal current, ellipse semi-minor axis');
netcdf.putAtt(nc_init,cmaxID,'units','meter second-1');
netcdf.putAtt(nc_init,cmaxID,'field','tide_Cmax, scalar, series');

netcdf.close(nc_init)



