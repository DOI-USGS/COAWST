function updatbdry_coawst_mw(fn,gn,bdy,wdr)

% written by Mingkui Li, May 2008
% further adaptions by jcw, April 18, 2009
% fn is the clm file
% gn is the working grid
% bdy is the bdy file name
% wdr is the working directory

%Use clm file to get the data
nc_clm=netcdf.open(fn,'NC_NOWRITE');
%get number of time steps in clm file
timeid = netcdf.inqVarID(nc_clm,'ocean_time');
ocean_time=netcdf.getVar(nc_clm,timeid);

%file name to be created
bndry_file=[bdy];

%call to create this bndry file.
create_roms_netcdf_bndry_mwUL(bndry_file,gn,length(ocean_time))
nc_bndry=netcdf.open(bndry_file,'NC_WRITE');

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

%% Time
zettimeid = netcdf.inqVarID(nc_bndry,'zeta_time');
netcdf.putVar(nc_bndry,zettimeid,ocean_time);

v2dtimeid = netcdf.inqVarID(nc_bndry,'v2d_time');
netcdf.putVar(nc_bndry,v2dtimeid,ocean_time);

v3dtimeid = netcdf.inqVarID(nc_bndry,'v3d_time');
netcdf.putVar(nc_bndry,v3dtimeid,ocean_time);

saltimeid = netcdf.inqVarID(nc_bndry,'salt_time');
netcdf.putVar(nc_bndry,saltimeid,ocean_time);

tmptimeid = netcdf.inqVarID(nc_bndry,'temp_time');
netcdf.putVar(nc_bndry,tmptimeid,ocean_time);

%% zeta
zetclmid = netcdf.inqVarID(nc_clm,'zeta');%get id
zeta=netcdf.getVar(nc_clm,zetclmid);%get variable
zeta_south=zeta(:,1);
zetbdyids = netcdf.inqVarID(nc_bndry,'zeta_south');%get id
netcdf.putVar(nc_bndry,zetbdyids,zeta_south);%set variable
clear zeta_south

zeta_east=zeta(end,:);
zetbdyide = netcdf.inqVarID(nc_bndry,'zeta_east');%get id
netcdf.putVar(nc_bndry,zetbdyide,zeta_east);%set variable
clear zeta_east

%% ubar
ubarclmid = netcdf.inqVarID(nc_clm,'ubar');
ubarclm=netcdf.getVar(nc_clm,ubarclmid);

ubar_south=ubarclm(:,1);
ubarbdyids = netcdf.inqVarID(nc_bndry,'ubar_south');
netcdf.putVar(nc_bndry,ubarbdyids,ubar_south);
clear ubar_south

ubar_east=ubarclm(end,:);
ubarbdyide = netcdf.inqVarID(nc_bndry,'ubar_east');
netcdf.putVar(nc_bndry,ubarbdyide,ubar_east);
clear ubar_east

%% vbar
vbarclmid = netcdf.inqVarID(nc_clm,'vbar');
vbarclm=netcdf.getVar(nc_clm,vbarclmid);

vbar_south=vbarclm(:,1);
vbarbdyids = netcdf.inqVarID(nc_bndry,'vbar_south');
netcdf.putVar(nc_bndry,vbarbdyids,vbar_south);
clear vbar_south

vbar_east=vbarclm(end,:);
vbarbdyide = netcdf.inqVarID(nc_bndry,'vbar_east');
netcdf.putVar(nc_bndry,vbarbdyide,vbar_east);
clear vbar_east


%% u
uclmid = netcdf.inqVarID(nc_clm,'u');
uclm=netcdf.getVar(nc_clm,uclmid);

u_south=uclm(:,1,:);
ubdyids = netcdf.inqVarID(nc_bndry,'u_south');
netcdf.putVar(nc_bndry,ubdyids,u_south);
clear u_south

u_east=uclm(end,:,:);
ubdyide = netcdf.inqVarID(nc_bndry,'u_east');
netcdf.putVar(nc_bndry,ubdyide,u_east);
clear u_east

%% v
vclmid = netcdf.inqVarID(nc_clm,'v');
vclm=netcdf.getVar(nc_clm,vclmid);

v_south=vclm(:,1,:);
vbdyids = netcdf.inqVarID(nc_bndry,'v_south');
netcdf.putVar(nc_bndry,vbdyids,v_south);
clear v_south

v_east=vclm(end,:,:);
vbdyide = netcdf.inqVarID(nc_bndry,'v_east');
netcdf.putVar(nc_bndry,vbdyide,v_east);
clear v_east

%% temp
tempclmid = netcdf.inqVarID(nc_clm,'temp');
tempclm=netcdf.getVar(nc_clm,tempclmid);

temp_south=tempclm(:,1,:);
tempbdyids = netcdf.inqVarID(nc_bndry,'temp_south');
netcdf.putVar(nc_bndry,tempbdyids,temp_south);
clear temp_south

temp_east=tempclm(end,:,:);
tempbdyide = netcdf.inqVarID(nc_bndry,'temp_east');
netcdf.putVar(nc_bndry,tempbdyide,temp_east);
clear temp_east

%% salt
saltclmid = netcdf.inqVarID(nc_clm,'salt');
saltclm=netcdf.getVar(nc_clm,saltclmid);

salt_south=saltclm(:,1,:);
saltbdyids = netcdf.inqVarID(nc_bndry,'salt_south');
netcdf.putVar(nc_bndry,saltbdyids,salt_south);
clear salt_south

salt_east=saltclm(end,:,:);
saltbdyide = netcdf.inqVarID(nc_bndry,'salt_east');
netcdf.putVar(nc_bndry,saltbdyide,salt_east);
clear salt_east

%%
netcdf.close(nc_bndry);
netcdf.close(nc_clm);
