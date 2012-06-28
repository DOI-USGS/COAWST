function mat2roms_mw(theMatFile, theROMSFile)

% mat2roms_mw('theMatFile', 'theROMSFile')
% mat2roms_mw -- Convert grid data from a Mat-file to a ROMS file
% using the native mathworks (mw) netcdf interface.
 
% Version 28-Dec-2011 Stephen Nicholls
% Updated 28-Jan-2012 John C. Warner

if nargin < 1, theMatFile = '*.mat'; end
if nargin < 2, theROMSFile = 'roms_model_grid.nc'; end

% Get the file names.
if any(theMatFile == '*')
	help(mfilename)
	theFilterSpec = theMatFile;
	thePrompt = 'Select a Mat-File';
	[theFile, thePath] = uigetfile(theFilterSpec, thePrompt);
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath(end+1) = filesep; end
	theMatFile = [thePath theFile];
end

if any(theROMSFile == '*')
	theFilterSpec = theROMSFile;
	thePrompt = 'Save As ROMS File';
	[theFile, thePath] = uiputfile(theFilterSpec, thePrompt);
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath(end+1) = filesep; end
	theROMSFile = [thePath theFile];
end

if isequal(theMatFile, theROMSFile)
	disp([' ## Must not select same file for input and output.'])
	return
end

% Load the Mat-File.
s = load(theMatFile);
if isempty(s)
	disp([' ## Mat-File is empty.'])
	return
end

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

[m, n] = size(s.rho.depth);

% The xi direction (left-right):

LP = m;   % The rho dimension.
L = LP-1; % The psi dimension.

% The eta direction (up-down):

MP = n;   % The rho dimension.
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

netcdf.endDef(nc)
 
% Fill the variables with data.

disp(' ## Filling Variables...')

projection = s.projection;
switch lower(projection)
case 'mercator'
	theProjection = 'ME';
case 'stereographic'
	theProjection = 'ST';
case 'lambert conformal conic'
	theProjection = 'LC';
otherwise
	theProjection = '??';
end
netcdf.putVar(nc,v3,theProjection)
netcdf.putVar(nc,v4,s.spherical)

% Determine coordinates.
m_proj(s.projection);
if ((s.spherical=='T') || (s.spherical=='t'))
  geogrid_lon = interp2(s.rho.lon, 'bilinear');
  geogrid_lat = interp2(s.rho.lat, 'bilinear');
% compute x and y
  [x, y] = m_ll2xy(s.rho.lon, s.rho.lat);
  x=x*6371000;
  y=y*6371000;
else
  [s.rho.lon, s.rho.lat] = m_xy2ll(s.rho.x/6371000, s.rho.y/6371000);
  geogrid_lon = interp2(s.rho.lon, 'bilinear');
  geogrid_lat = interp2(s.rho.lat, 'bilinear');
  x = s.rho.x;
  y = s.rho.y;
end
grid_x = interp2(x, 'bilinear');
grid_y = interp2(y, 'bilinear');

xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));

netcdf.putVar(nc,v1,xl)
netcdf.putVar(nc,v2,el)

% Depths at RHO points.
bathymetry = s.rho.depth;
if ~isempty(bathymetry)
  netcdf.putVar(nc,v5,min(min(bathymetry)));
  netcdf.putVar(nc,v6,max(max(bathymetry)));
  netcdf.putVar(nc,v7,[0 0 0],[LP MP 1],bathymetry);
  netcdf.putVar(nc,v8,bathymetry);
end

% Coriolis
f = sw_f(s.rho.lat);
netcdf.putVar(nc,v9,f)

% Handy indices.
[m, n] = size(grid_x);

i_rho = 1:2:m;
j_rho = 1:2:n;
i_psi = 2:2:m;
j_psi = 2:2:n;
i_u = 2:2:m;
j_u = 1:2:n;
i_v = 1:2:m;
j_v = 2:2:n;

% Locations.
netcdf.putVar(nc,v14,x);
netcdf.putVar(nc,v15,y);

netcdf.putVar(nc,v16,grid_x(i_psi, j_psi));
netcdf.putVar(nc,v17,grid_y(i_psi, j_psi));

netcdf.putVar(nc,v18,grid_x(i_u, j_u));
netcdf.putVar(nc,v19,grid_y(i_u, j_u));

netcdf.putVar(nc,v20,grid_x(i_v, j_v));
netcdf.putVar(nc,v21,grid_y(i_v, j_v));

netcdf.putVar(nc,v22,geogrid_lat(i_rho, j_rho));
netcdf.putVar(nc,v23,geogrid_lon(i_rho, j_rho));

netcdf.putVar(nc,v24,geogrid_lat(i_psi, j_psi));
netcdf.putVar(nc,v25,geogrid_lon(i_psi, j_psi));

netcdf.putVar(nc,v26,geogrid_lat(i_u, j_u));
netcdf.putVar(nc,v27,geogrid_lon(i_u, j_u));

netcdf.putVar(nc,v28,geogrid_lat(i_v, j_v));
netcdf.putVar(nc,v29,geogrid_lon(i_v, j_v));


% Compute dx, dy, angle.
if ((s.spherical=='T') || (s.spherical=='t'))
  lat_u=geogrid_lat(i_u, j_u);
  lon_u=geogrid_lon(i_u, j_u);
  [dx,ang]=sw_dist(lat_u(:),lon_u(:),'km');

  dx=[dx(:); dx(end)]*1000;  % km==> m
  dx=reshape(dx,LP-1,MP);
  dx=[dx(1,:); dx(1:end-1,:); dx(end-1,:)];

  ang=[ang(:); ang(end)];
  ang=reshape(ang,LP-1,MP);
  ang=[ang(1,:); ang(1:end-1,:); ang(end-1,:)];
  ang=ang*pi/180;

  lat_v=geogrid_lat(i_v, j_v).';
  lon_v=geogrid_lon(i_v, j_v).';
  dy=sw_dist(lat_v(:),lon_v(:),'km');

  dy=[dy(:); dy(end)]*1000;   % km ==> m 
  dy=reshape(dy,MP-1,LP);
  dy=[dy(1,:); dy(1:end-1,:); dy(end-1,:)];
  dy=dy.';
else
  x_u=grid_x(i_u, j_u);
  dx=diff(x_u);
  dx=[dx(1,:); dx; dx(end,:)];

  y_v=grid_y(i_v, j_v).';
  dy=diff(y_v);
  dy=dy.';
  dy=[dy(:,1) dy dy(:,end)];

  y_v=grid_y(i_psi, j_psi);
  y_v=diff(y_v); 
  x_v=dx(2:end-1,2:end);
  ang=angle(x_v+y_v*sqrt(-1));
  ang=[ang(end,:); ang; ang(1,:)];
  ang=[ang ang(:,end)];
end

pm=1./dx;
pn=1./dy;
netcdf.putVar(nc,v10,pm);
netcdf.putVar(nc,v11,pn);

dmde = zeros(LP, MP);
dndx = zeros(LP, MP);
dmde(:,2:end-1) = 0.5*(dx(:,3:end) - dx(:,1:end-2));
dmde(:,1)=dmde(:,2);
dmde(:,end)=dmde(:,end-1);
dndx(2:end-1,:) = 0.5*(dy(3:end,:) - dy(1:end-2,:));
dmde(1,:)=dmde(2,:);
dmde(end,:)=dmde(end-1,:);
netcdf.putVar(nc,v12,dndx);
netcdf.putVar(nc,v13,dmde);

% Masking.
mask = s.rho.mask;

water = double(mask);
netcdf.putVar(nc,v30,water);

u_mask = water(1:end-1,:) & water(2:end,:);
netcdf.putVar(nc,v31,double(u_mask));

v_mask= water(:,1:end-1) & water(:,2:end);
netcdf.putVar(nc,v32,double(v_mask));

psi_mask= water(1:end-1,1:end-1) & water(1:end-1,2:end) & water(2:end,1:end-1) & water(2:end,2:end);
netcdf.putVar(nc,v33,double(psi_mask));

% Angle.
netcdf.putVar(nc,v34,ang);  % Degrees.

netcdf.close(nc)
