function create_roms_netcdf_grid_file(theROMSFile,LP,MP)

%
% creates a netcdf file with roms grid variables.
% jcw 07Aug2012
%

% Open the ROMS File.
nc = netcdf.create(theROMSFile, 'clobber');
if isempty(nc)
	disp([' ## Unable to open ROMS NetCDF output file.'])
	return
end

% Populate the ROMS File.

%% Global attributes:

disp(' ## Defining Global Attributes...')

netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','ROMS GRID file');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'gridid','theGridTitle');
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by ', mfilename ', on ', datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','ROMS Application')

% Dimensions:

% The xi direction (left-right):
L = LP-1; % The psi dimension.

% The eta direction (up-down):
M = MP-1; % The psi dimension.

disp(' ## Defining Dimensions...')

xi_psi = netcdf.defDim(nc,'xi_psi',L);
xi_rho = netcdf.defDim(nc,'xi_rho',LP);
xi_u   = netcdf.defDim(nc,'xi_u',L);
xi_v   = netcdf.defDim(nc,'xi_v',LP);

eta_psi = netcdf.defDim(nc,'eta_psi',M);
eta_rho = netcdf.defDim(nc,'eta_rho',MP);
eta_u   = netcdf.defDim(nc,'eta_u',MP);
eta_v   = netcdf.defDim(nc,'eta_v',M);

one = netcdf.defDim(nc,'one',1);
two = netcdf.defDim(nc,'two',2);
%bath = netcdf.defDim(nc,'bath',netcdf.getConstant('NC_UNLIMITED'));
bath = netcdf.defDim(nc,'bath',1);

%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')

v1 = netcdf.defVar(nc,'xl','double',one);
netcdf.putAtt(nc,v1,'long_name','domain length in the XI-direction')
netcdf.putAtt(nc,v1,'units','meter')

v2 = netcdf.defVar(nc,'el','double',one);
netcdf.putAtt(nc,v2,'long_name','domain length in the ETA-direction')
netcdf.putAtt(nc,v2,'units','meter')

v3 = netcdf.defVar(nc,'JPRJ','char',two);
netcdf.putAtt(nc,v3,'long_name','Map projection type')
netcdf.putAtt(nc,v3,'option_ME','Mercator')
netcdf.putAtt(nc,v3,'option_ST','Stereographic')
netcdf.putAtt(nc,v3,'option_LC','Lambert conformal conic')

v4 = netcdf.defVar(nc,'spherical','char',one);
netcdf.putAtt(nc,v4,'long_name','Grid type logical switch')
netcdf.putAtt(nc,v4,'option_T','spherical')
netcdf.putAtt(nc,v4,'option_F','Cartesian')

v5 = netcdf.defVar(nc,'depthmin','short',one);
netcdf.putAtt(nc,v5,'long_name','domain length in the XI-direction')
netcdf.putAtt(nc,v5,'units','meter')

v6 = netcdf.defVar(nc,'depthmax','short',one);
netcdf.putAtt(nc,v6,'long_name','Deep bathymetry clipping depth')
netcdf.putAtt(nc,v6,'units','meter')

v7 = netcdf.defVar(nc,'hraw','double',[xi_rho eta_rho bath]);
netcdf.putAtt(nc,v7,'long_name','Working bathymetry at RHO-points')
netcdf.putAtt(nc,v7,'units','meter')
netcdf.putAtt(nc,v7,'field','bath, scalar')

v8 = netcdf.defVar(nc,'h','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v8,'long_name','Final bathymetry at RHO-points')
netcdf.putAtt(nc,v8,'units','meter')
netcdf.putAtt(nc,v8,'field','bath, scalar')

v9 = netcdf.defVar(nc,'f','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v9,'long_name','Coriolis parameter at RHO-points')
netcdf.putAtt(nc,v9,'units','second-1')
netcdf.putAtt(nc,v9,'field','Corilis, scalar')

v10 = netcdf.defVar(nc,'pm','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v10,'long_name','curvilinear coordinate metric in XI')
netcdf.putAtt(nc,v10,'units','meter-1')
netcdf.putAtt(nc,v10,'field','pm, scalar')

v11 = netcdf.defVar(nc,'pn','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v11,'long_name','curvilinear coordinate metric in ETA')
netcdf.putAtt(nc,v11,'units','meter-1')
netcdf.putAtt(nc,v11,'field','pn, scalar')

v12 = netcdf.defVar(nc,'dndx','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v12,'long_name','xi derivative of inverse metric factor pn')
netcdf.putAtt(nc,v12,'units','meter')
netcdf.putAtt(nc,v12,'field','dndx, scalar')

v13 = netcdf.defVar(nc,'dmde','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v13,'long_name','eta derivative of inverse metric factor pm')
netcdf.putAtt(nc,v13,'units','meter')
netcdf.putAtt(nc,v13,'field','dmde, scalar')

v14 = netcdf.defVar(nc,'x_rho','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v14,'long_name','x location of RHO-points')
netcdf.putAtt(nc,v14,'units','meter')

v15 = netcdf.defVar(nc,'y_rho','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v15,'long_name','y location of RHO-points')
netcdf.putAtt(nc,v15,'units','meter')

v16 = netcdf.defVar(nc,'x_psi','double',[xi_psi eta_psi]);
netcdf.putAtt(nc,v16,'long_name','x location of PSI-points')
netcdf.putAtt(nc,v16,'units','meter')

v17 = netcdf.defVar(nc,'y_psi','double',[xi_psi eta_psi]);
netcdf.putAtt(nc,v17,'long_name','y location of PSI-points')
netcdf.putAtt(nc,v17,'units','meter')

v18 = netcdf.defVar(nc,'x_u','double',[xi_u eta_u]);
netcdf.putAtt(nc,v18,'long_name','x location of U-points')
netcdf.putAtt(nc,v18,'units','meter')

v19 = netcdf.defVar(nc,'y_u','double',[xi_u eta_u]);
netcdf.putAtt(nc,v19,'long_name','y location of U-points')
netcdf.putAtt(nc,v19,'units','meter')
 
v20 = netcdf.defVar(nc,'x_v','double',[xi_v eta_v]);
netcdf.putAtt(nc,v20,'long_name','x location of V-points')
netcdf.putAtt(nc,v20,'units','meter')

v21 = netcdf.defVar(nc,'y_v','double',[xi_v eta_v]);
netcdf.putAtt(nc,v21,'long_name','y location of V-points')
netcdf.putAtt(nc,v21,'units','meter')
 
v22 = netcdf.defVar(nc,'lat_rho','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v22,'long_name','latitude of RHO-points')
netcdf.putAtt(nc,v22,'units','degree_north')

v23 = netcdf.defVar(nc,'lon_rho','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v23,'long_name','longitude of RHO-points')
netcdf.putAtt(nc,v23,'units','degree_east')

v24 = netcdf.defVar(nc,'lat_psi','double',[xi_psi eta_psi]);
netcdf.putAtt(nc,v24,'long_name','latitude of PSI-points')
netcdf.putAtt(nc,v24,'units','degree_north')

v25 = netcdf.defVar(nc,'lon_psi','double',[xi_psi eta_psi]);
netcdf.putAtt(nc,v25,'long_name','longitude of PSI-points')
netcdf.putAtt(nc,v25,'units','degree_east')

v26 = netcdf.defVar(nc,'lat_u','double',[xi_u eta_u]);
netcdf.putAtt(nc,v26,'long_name','latitude of U-points')
netcdf.putAtt(nc,v26,'units','degree_north')

v27 = netcdf.defVar(nc,'lon_u','double',[xi_u eta_u]);
netcdf.putAtt(nc,v27,'long_name','longitude of U-points')
netcdf.putAtt(nc,v27,'units','degree_east')

v28= netcdf.defVar(nc,'lat_v','double',[xi_v eta_v]);
netcdf.putAtt(nc,v28,'long_name','latitude of V-points')
netcdf.putAtt(nc,v28,'units','degree_north')

v29 = netcdf.defVar(nc,'lon_v','double',[xi_v eta_v]);
netcdf.putAtt(nc,v29,'long_name','longitude of V-points')
netcdf.putAtt(nc,v29,'units','degree_east')

v30 = netcdf.defVar(nc,'mask_rho','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v30,'long_name','mask on RHO-points')
netcdf.putAtt(nc,v30,'option_0','land')
netcdf.putAtt(nc,v30,'option_1','water')

v31 = netcdf.defVar(nc,'mask_u','double',[xi_u eta_u]);
netcdf.putAtt(nc,v31,'long_name','mask on U-points')
netcdf.putAtt(nc,v31,'option_0','land')
netcdf.putAtt(nc,v31,'option_1','water')

v32 = netcdf.defVar(nc,'mask_v','double',[xi_v eta_v]);
netcdf.putAtt(nc,v32,'long_name','mask on U-points')
netcdf.putAtt(nc,v32,'option_0','land')
netcdf.putAtt(nc,v32,'option_1','water')

v33 = netcdf.defVar(nc,'mask_psi','double',[xi_psi eta_psi]);
netcdf.putAtt(nc,v33,'long_name','mask on PSI-points')
netcdf.putAtt(nc,v33,'option_0','land')
netcdf.putAtt(nc,v33,'option_1','water')

v34 = netcdf.defVar(nc,'angle','double',[xi_rho eta_rho]);
netcdf.putAtt(nc,v34,'long_name','angle between xi axis and east')
netcdf.putAtt(nc,v34,'units','degree')

netcdf.close(nc)
 
